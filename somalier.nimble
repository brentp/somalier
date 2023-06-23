#import somalierpkg/version as _

version       = "0.2.17" #somalierVersion
author        = "Brent Pedersen"
description   = "sample-swap checking directly on BAMs/CRAMs for cancer data"
license       = "MIT"



# Dependencies

requires "https://github.com/brentp/zip#dev"
requires "nim >= 1.2.0", "hts >= 0.3.20", "https://github.com/brentp/pedfile >= 0.0.3", "https://github.com/brentp/hileup", "argparse >=0.10.1 & < 0.10.2", "lapper", "arraymancer#head"
requires "https://github.com/brentp/slivar#head"
srcDir = "src"

bin = @["somalier"]

task test, "run the tests":
  exec "nim c  -d:useSysAssert -d:useGcAssert --lineDir:on --debuginfo -r tests/test_groups"
  exec "bash tests/functional-tests.sh"
