import hts/bam
import osproc
import strformat
import random
import algorithm
import hts/vcf
import os
import sets
import tables


type read = object
  #qname:string
  sequence:string
  base_qualities: seq[uint8]

type pair* = object
  qname*: string
  read1*: read
  read2*: read
#[
create a bam with only reads from the regions that overlap the bases used by somalier
]#

proc complement(s:char): char {.inline.} =
    if s == 'C':
        return 'G'
    elif s == 'G':
        return 'C'
    elif s == 'A':
        return 'T'
    elif s == 'T':
        return 'A'
    else:
        return s

proc reverse_complement(xs: string): string =
  result = newString(xs.len)
  for i, x in xs:
    result[xs.high-i] = complement(x)

proc simulate_variant(aln:Record, site:Variant, freq: float): read =
  var r = read()
  discard aln.sequence(r.sequence)
  discard aln.base_qualities(r.base_qualities)
  if site == nil or random(1'f) > freq or aln.chrom != $site.CHROM or aln.start > site.start or aln.stop < site.stop:
    if aln.flag.reverse:
      reverse(r.base_qualities)
      r.sequence = reverse_complement(r.sequence)
    return r

  var
    off = aln.start
    qoff = 0
    r_off_only = 0

  for event in aln.cigar:
    var cons = event.consumes
    if cons.query:
      qoff += event.len
    if cons.reference:
      off += event.len
      if not cons.query:
        roff_only += event.len
    if off <= site.start: continue

    if not (cons.query and cons.reference): continue

    var over = off - site.start - roff_only
    if over > qoff: break
    if over < 0: continue
    doAssert qoff - over >= 0

    if site.start mod 2 == 0:
      r.sequence[qoff - over] = site.REF[0]
    else:
      r.sequence[qoff - over] = site.ALT[0][0]
    break

  if aln.flag.reverse:
    reverse(r.base_qualities)
    r.sequence = reverse_complement(r.sequence)
  return r

proc extract(bam_path: string, vcf_path: string, fasta:string): seq[pair] =
  var s = initSet[string]()

  var ivcf:VCF
  var ibam:Bam
  if not open(ivcf, vcf_path):
    quit "couldn't open vcf"

  if not open(ibam, bam_path, index=true, fai=fasta):
    quit "couldn't open bam"

  doAssert 0 == ibam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, 8191 - SAM_SEQ.int - SAM_QUAL.int - SAM_AUX.int - SAM_RGAUX.int)

  for v in ivcf:
    for aln in ibam.query($v.CHROM, v.start, v.stop):
      s.incl(aln.qname)

  ibam.close()
  ibam = nil
  GC_fullCollect()

  if not open(ibam, bam_path, fai=fasta, threads=2):
    quit "couldn't open bam"
  doAssert 0 == ibam.set_option(FormatOption.CRAM_OPT_REQUIRED_FIELDS, 8191 - SAM_AUX.int - SAM_RGAUX.int)

  var pairs = newTable[string, pair](256)
  var n = 0

  for aln in ibam:
    n += 1
    if aln.qname notin s: continue
    if aln.flag.secondary or aln.flag.supplementary: continue

    var p = pairs.mgetOrPut(aln.qname, pair(qname: aln.qname))
    # TODO: hwo to get simulated variant?)
    # move adding to pairs above? yes.
    var r = simulate_variant(aln, nil, 1)

    if aln.flag.read1:
      p.read1 = r
    else:
      p.read2 = r

    pairs[aln.qname] = p

  result = newSeqOfCap[pair](pairs.len)
  for qname, p in pairs:
    if p.read1.sequence != "" and p.read2.sequence != "":
      result.add(p)

  echo "got ", result.len, " pairs from ", n, " total reads"

proc encode(bqs: seq[uint8]): string =
  result = newString(bqs.len)
  for i, b in bqs:
    result[i] = char(b + 33)

proc write(reads: seq[pair], base:string) =

  var fh1:File
  var fh2:File

  if not open(fh1, &"{base}_R1.fastq", mode=fmWrite):
    quit "couldn't open fastq"
  if not open(fh2, &"{base}_R2.fastq", mode=fmWrite):
    quit "couldn't open fastq"

  for r in reads:
    doAssert r.read2.sequence != ""
    doAssert r.read1.sequence != ""
    fh1.write(&"@{r.qname}\n{r.read1.sequence}\n+\n{encode(r.read1.base_qualities)}\n")
    fh2.write(&"@{r.qname}\n{r.read2.sequence}\n+\n{encode(r.read2.base_qualities)}\n")

  fh1.close()
  fh2.close()

proc contaminate(a: var seq[pair], cont: var seq[pair], out_fastq_base: string, level: float, fasta:string) =
  shuffle(a)
  shuffle(cont)

  var n_a = int(0.5 + (1 - level) * a.len.float)
  var n_c = int(0.5 + level * a.len.float)

  var reads = a[0..<n_a]
  reads.add(cont[0..<n_c])
  shuffle(reads)

  echo "writing ", reads.len, " pairs"
  reads.write(out_fastq_base)

  var name = splitFile(out_fastq_base).name

  var cmd = &"""bwa mem -R "@RG\tID:{name}\tSM:{name}{level:.1g}" -t 5 {fasta} {out_fastq_base}_R1.fastq {out_fastq_base}_R2.fastq"""
  cmd &= &"""| samblaster | samtools sort -@ 5 -o {out_fastq_base}_c{level:.3f}.bam && samtools index {out_fastq_base}_c{level:.3f}.bam"""
  doAssert execCmd(cmd) == 0
  echo &"wrote {out_fastq_base}_c{level:.3f}.bam and .bai"

  removeFile(&"{out_fastq_base}_R1.fastq")
  removeFile(&"{out_fastq_base}_R2.fastq")


when isMainModule:
  randomize(42)

  #  ../sites.chr.hg38.vcf.gz NA12891.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram NA12892.alt_bwamem_GRCh38DH.20150826.CEU.exome.cram /data/human/GRCh38_full_analysis_set_plus_decoy_hla.fa NA12891

  var vcf_path = paramStr(1) # path to sites vcf used by somalier
  var base_bam_path = paramStr(2)
  var contam_bam_path = paramStr(3)
  var fasta = paramStr(4)
  var out_fastq_base = paramStr(5)

  var base_read_pairs = extract(base_bam_path, vcf_path, fasta)
  var contam_read_pairs = extract(contam_bam_path, vcf_path, fasta)

  var contam_levels = @[0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.2, 0.3]

  for l in contam_levels:
    echo "contaminating at:", l
    contaminate(base_read_pairs, contam_read_pairs, out_fastq_base, l, fasta)

  echo "writing subset of contamination file"
  contam_read_pairs.write("contamination")
  var cmd = &"""bwa mem -R "@RG\tID:contam\tSM:contam" -t 5 {fasta} contamination_R1.fastq contamination_R2.fastq"""
  cmd &= &"""| samblaster | samtools sort -@3 -m4G -o contamination.bam && samtools index contamination.bam"""
  doAssert execCmd(cmd) == 0
  echo &"wrote contamination.bam and .bai"
  removeFile("contamination_R1.fastq")
  removeFile("contamination_R2.fastq")
