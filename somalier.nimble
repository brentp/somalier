# Package
import ospaths
template thisModuleFile: string = instantiationInfo(fullPaths = true).filename

when fileExists(thisModuleFile.parentDir / "src/somalier.nim"):
  # In the git repository the Nimble sources are in a ``src`` directory.
  import src/somalierpkg/version as _
else:
  # When the package is installed, the ``src`` directory disappears.
  import somalierpkg/version as _

version       = somalierVersion
author        = "Brent Pedersen"
description   = "sample-swap checking directly on BAMs/CRAMs for cancer data"
license       = "academic only"


# Dependencies

requires "nim >= 0.19.0", "hts >= 0.2.13", "https://github.com/brentp/bpbio", "https://github.com/brentp/hileup", "argparse"
srcDir = "src"

#bin = @["./somalier.nim"]

