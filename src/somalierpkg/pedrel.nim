import strutils
import strformat
import argparse
import pedfile
import sequtils
import algorithm

proc pedrel_main*() =
  var argv = commandLineParams()
  if argv[0] == "pedrel":
    argv = argv[1..argv.high]

  var p = newParser("somalier pedrel"):
    help("report pairwise relationships from pedigree file")
    option("-o", "--output", help="output file path", default="stdout")
    option("-m", "--min-relatedness", help="minimum relatedness to report", default="0.01")
    arg("pedfile", help="pedigry (fam) file path")

  let opts = p.parse(argv)
  if opts.help:
    quit 0

  if opts.pedfile == "":
    echo p.help
    quit "[somalier] pedfile argument required"

  # Parse the pedigree file
  var samples = parse_ped(opts.pedfile)

  let min_rel = parseFloat(opts.min_relatedness)

  # Calculate pairwise relatedness for all sample pairs
  var output_file: File
  if opts.output == "stdout":
    output_file = stdout
  else:
    if not open(output_file, opts.output, fmWrite):
      quit "[somalier] couldn't open output file: " & opts.output

  # Write header
  output_file.write_line("sample_a\tsample_b\trelatedness")

  # Iterate through all sample pairs and report those with min_rel threshold
  for i in 0..<samples.len:
    let sampleA = samples[i]
    for j in (i+1)..<samples.len:
      let sampleB = samples[j]
      let rel = sampleA.relatedness(sampleB)

      if rel > min_rel:
        # Ensure consistent ordering (alphabetical by sample ID)
        let (sample1, sample2) = if sampleA.id < sampleB.id:
                                  (sampleA.id, sampleB.id)
                                else:
                                  (sampleB.id, sampleA.id)
        output_file.write_line(&"{sample1}\t{sample2}\t{rel:.6f}")

  if opts.output != "stdout":
    output_file.close()
    stderr.write_line(&"[somalier] wrote {samples.len} sample relationships to: {opts.output}")