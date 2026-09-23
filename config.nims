when defined(somalier_usearch):
  # USearch includes a C++ adapter. Nim's default C++ compiler is g++, but some
  # build environments (including Clang-only containers) do not provide it.
  if findExe("g++").len == 0 and findExe("clang++").len > 0:
    switch("cc", "clang")
