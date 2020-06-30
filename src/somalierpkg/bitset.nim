{.passC:"-mpopcnt".}
import bitops

type bitset* = seq[uint64]

proc create_bitset*(size:SomeInteger): bitset =
  ## create a new bitset able to hold at least `size`
  let L = (size.int + sizeof(uint64) * 8 - 1) div (sizeof(uint64) * 8)
  doAssert L * 8 * sizeof(uint64) >= size, $(size, L * sizeof(uint64), L)
  result = newSeq[uint64](L)

proc set*(b:var bitset, i:SomeInteger) {.inline.} =
  ## set the ith bit.
  let idx = i shr 6
  b[idx] = b[idx] or (1'u64 shl (i mod 64))

type genotypes* = tuple[hom_ref:bitset, het:bitset, hom_alt:bitset]
type IBSResult* = tuple[IBS0: int32, IBS2:int32, N:int32, shared_hets:int32, shared_hom_alts:int32, het_ab:int32]

proc IBS*(a: genotypes, b: genotypes): IBSResult =
  doAssert a.hom_ref.len == b.hom_ref.len
  var sh_het:uint64 # shared het
  var sh_ha:uint64 # shared hom-alt
  var n:uint64
  for k in 0..a.hom_ref.high:
    result.IBS0 += countSetBits((a.hom_ref[k] and b.hom_alt[k]) or (a.hom_alt[k] and b.hom_ref[k])).int32
    sh_het = a.het[k] and b.het[k]
    sh_ha = a.hom_alt[k] and b.hom_alt[k]
    result.IBS2 += countSetBits(sh_ha or sh_het or (a.hom_ref[k] and b.hom_ref[k])).int32
    result.shared_hets += countSetBits(sh_het).int32
    result.shared_hom_alts += countSetBits(sh_ha).int32
    n = (a.hom_ref[k] or a.het[k] or a.hom_alt[k]) and (b.hom_ref[k] or b.het[k] or b.hom_alt[k])
    # het_ab is count of sites that are het in either sample and not unknown in
    # the other sample. see equation 9 and 7 from king paper (PMID: 20926424) where they
    # state:
    # where NAa,Aa, NAa(i) and NAa(j) are the total numbers of SNPs at which both individuals of the pair are heterozygous, and the total number of heterozygotes for the i-th and j-th individual, respectively, excluding those SNPs with missing genotypes in either individual of the pair.
    # Here, NAa(i) and NAa(j) are a.het[k] and b.het[k].
    result.het_ab += (countSetBits(b.het[k] and n) + countSetBits(a.het[k] and n)).int32
    result.N += countSetBits(n).int32

proc XIBS*(a: genotypes, b: genotypes): IBSResult =
  doAssert a.hom_ref.len == b.hom_ref.len
  var sh_het:uint64 # shared het
  var sh_ha:uint64 # shared hom-alt
  for k in 0..a.hom_ref.high:
    result.IBS0 += countSetBits((a.hom_ref[k] and b.hom_alt[k]) or (a.hom_alt[k] and b.hom_ref[k])).int32
    sh_het = a.het[k] and b.het[k]
    sh_ha = a.hom_alt[k] and b.hom_alt[k]
    result.IBS2 += countSetBits(sh_ha or sh_het or (a.hom_ref[k] and b.hom_ref[k])).int32
