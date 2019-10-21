include somalierpkg/relate
import slivarpkg/pedfile

import unittest

suite "somalier groups":
  test "that pairs added in groups update relatedness for those added by pedigree file":
    var fh:File
    doAssert open(fh, "_grps.txt", fmWrite)
    fh.write_line("normal0,tumor0")
    fh.close

    doAssert open(fh, "_ped.txt", fmWrite)
    fh.write("""FAM001	normal0	normal1	normal2	2	-9
FAM001	normal1	0	0	1	-9
FAM001	normal2	0	0	2	-9
""")
    fh.close

    var samples = parse_ped("_ped.txt")
    var sample_names: seq[string]
    for s in samples: sample_names.add(s.id)

    var groups: seq[pair]
    groups.add_ped_samples(samples, sample_names)
    groups.add(readGroups("_grps.txt", groups))

    check (a: "normal1", b: "tumor0", rel: 0.5) in groups
    check (a: "normal2", b: "tumor0", rel: 0.5) in groups




