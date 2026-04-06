nim c -d:danger src/somalier.nim
DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib:/opt/homebrew/opt/openblas/lib ./src/somalier relate ../somalier_out/*.somalier -o isabl-test -g ../grps.txt
