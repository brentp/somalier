import bitops

type bitset* = seq[uint64]

proc create_bitset*(capacity:SomeInteger): bitset =
  ## create a new bitset able to hold at least `capacity`
  let size = capacity * sizeof(uint64)
  let L = (size + sizeof(uint64) * 8 - 1) div (sizeof(uint64) * 8)
  result = newSeq[uint64](L)

proc set*(b:var bitset, i:SomeInteger) {.inline.} =
  ## set the ith bit.
  let idx = i shr 6
  b[idx] = b[idx] or (1'u64 shl (i mod 64))

proc intersection_count(b:bitset, o:bitset): int {.inline.} =
  doAssert b.len == o.len

  var j = 0
  for k in countup(0, b.len - 7, 8):
    result += countSetBits(b[k] and o[k])
    result += countSetBits(b[k+1] and o[k+1])
    result += countSetBits(b[k+2] and o[k+2])
    result += countSetBits(b[k+3] and o[k+3])
    result += countSetBits(b[k+4] and o[k+4])
    result += countSetBits(b[k+5] and o[k+5])
    result += countSetBits(b[k+6] and o[k+6])
    result += countSetBits(b[k+7] and o[k+7])
    result += countSetBits(b[k+8] and o[k+8])
    j = k + 7

  for k in countup(j + 1, b.len):
    result += countSetBits(b[k] and o[k])

type genotypes* = tuple[hom_ref:bitset, het:bitset, hom_alt:bitset]
type IBSResult* = tuple[IBS0: int, IBS2:int, N:int, shared_hets:int, shared_hom_alts:int]

proc IBS*(a: genotypes, b: genotypes): IBSResult =
  doAssert a.hom_ref.len == b.hom_ref.len
  #var j = 0
  #for k in countup(0, ahr.len - 7, 8):
  var sh_het:uint64 # shared het
  var sh_ha:uint64 # shared hom-alt
  for k in 0..a.hom_ref.high:
    # TODO: simd this
    result.IBS0 += countSetBits((a.hom_ref[k] and b.hom_alt[k]) or (a.hom_alt[k] and b.hom_ref[k]))
    sh_het = a.het[k] and b.het[k]
    sh_ha = a.hom_alt[k] and b.hom_alt[k]
    result.IBS2 += countSetBits(sh_ha or sh_het or (a.hom_ref[k] and b.hom_ref[k]))
    result.shared_hets += countSetBits(sh_het)
    result.shared_hom_alts += countSetBits(sh_ha)
    result.N += countSetBits((a.hom_ref[k] or a.het[k] or a.hom_alt[k]) and (b.hom_ref[k] or b.het[k] or b.hom_alt[k]))

#[
proc krelated*(alts: var seq[int8], ibs: var seq[uint16], n: var seq[uint16], hets: var seq[uint16], homs: var seq[uint16], shared_hom_alts: var seq[uint16], n_samples: int): int {.inline.} =

  if alts[n_samples - 1] == 1:
    hets[n_samples-1] += 1
  elif alts[n_samples - 1] == 2:
    homs[n_samples-1] += 1

  var is_het: bool
  var aj, ak: int8
  var nused = 0

  for j in 0..<(n_samples-1):
    aj = alts[j]
    if aj == -1: continue
    is_het = (aj == 1)

    if is_het:
      hets[j] += 1
    elif aj == 2:
      homs[j] += 1

    nused += 1

    for k in j+1..<n_samples:
      ak = alts[k]
      if ak == -1: continue
      n[j * n_samples + k] += 1
      if is_het:
        # shared hets
        if ak == 1:
          ibs[k * n_samples + j] += 1
      else:
        # ibs0
        if aj != ak and aj + ak == 2:
          ibs[j * n_samples + k] += 1
      # ibs2
      if aj == ak: #and not is_het:
        n[k * n_samples + j] += 1
        if aj == 2:
          shared_hom_alts[j * n_samples + k] += 1
  return nused

]#
