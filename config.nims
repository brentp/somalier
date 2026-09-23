when defined(somalier_usearch):
  # The query phase creates one Nim worker per reserved USearch search slot.
  # Configure this here so nimble and direct compiler builds behave alike.
  switch("threads", "on")

  # USearch includes a C++ adapter. Nim's default C++ compiler is g++, but some
  # build environments (including Clang-only containers) do not provide it.
  if findExe("g++").len == 0 and findExe("clang++").len > 0:
    switch("cc", "clang")
