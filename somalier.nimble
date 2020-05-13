#import somalierpkg/version as _

version       = "0.2.11" #somalierVersion
author        = "Brent Pedersen"
description   = "sample-swap checking directly on BAMs/CRAMs for cancer data"
license       = "academic only"



# Dependencies

requires "https://github.com/brentp/zip#dev"
requires "nim >= 0.20.0", "hts >= 0.3.4", "https://github.com/brentp/pedfile", "https://github.com/brentp/hileup", "argparse", "lapper", "arraymancer#head"
requires "https://github.com/brentp/slivar#head"
srcDir = "src"

#bin = @["./somalier"]
#bin = @["somalier"]

task test, "run the tests":
  exec "nim c  -d:useSysAssert -d:useGcAssert --lineDir:on --debuginfo -r tests/test_groups"
  exec "bash tests/functional-tests.sh"
