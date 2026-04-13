import tables
import math
import ./common

const
  maxOtherFraction* = 0.04
  defaultCharrHomRate* = 0.20
  defaultCharrHomAlpha* = 0.001
  minContamAlleleAf = 1e-6'f64

type CharrStats* = object
  contamination*: float64
  n_sites_usable*: int
  n_hom_ref_usable*: int
  n_hom_alt_usable*: int

proc max_hom_minor_reads*(dp: int, hom_minor_rate: float, alpha: float): int =
  ## Return the largest minor-read count still plausible under a
  ## homozygous-like binomial null. This gets stricter as depth increases.
  let q = 1.0 - hom_minor_rate
  var
    pmf = exp(dp.float64 * ln(q))
    tail = 1.0
  result = -1
  for k in 0..dp:
    if tail >= alpha:
      result = k
    else:
      break
    if k == dp:
      break
    tail -= pmf
    pmf = pmf * ((dp - k).float64 / (k + 1).float64) * (hom_minor_rate / q)

proc call_homozygous_like_genotype*(ac: allele_count, min_depth: int,
    hom_minor_rate: float, alpha: float,
    max_minor_by_depth: var Table[int, int]): int8 {.inline.} =
  let dp = int(ac.nref + ac.nalt)
  if dp < min_depth:
    return -1
  if ac.nother.float / (dp.float + ac.nother.float) > maxOtherFraction:
    return -1
  let max_minor = max_minor_by_depth.mgetOrPut(dp,
      max_hom_minor_reads(dp, hom_minor_rate, alpha))
  let minor = min(ac.nref, ac.nalt).int
  if minor > max_minor:
    return -1
  if ac.nref > ac.nalt:
    return 0
  if ac.nalt > ac.nref:
    return 2
  return -1

proc estimate_charr*(
    allele_counts: openArray[allele_count],
    pop_afs: openArray[float32],
    min_depth: int,
    hom_minor_rate: float,
    hom_alpha: float
  ): CharrStats =
  doAssert allele_counts.len == pop_afs.len

  var
    total = 0'f64
    max_minor_by_depth = initTable[int, int]()
  for i, ac in allele_counts:
    let gt = call_homozygous_like_genotype(ac, min_depth, hom_minor_rate,
        hom_alpha, max_minor_by_depth)
    let dp = (ac.nref + ac.nalt).float64
    let pB = pop_afs[i].float64
    var
      infiltrating = 0'f64
      contam_af = 0'f64

    case gt
    of 0'i8:
      result.n_hom_ref_usable.inc
      infiltrating = ac.nalt.float64
      contam_af = pB
    of 2'i8:
      result.n_hom_alt_usable.inc
      infiltrating = ac.nref.float64
      contam_af = 1'f64 - pB
    else:
      continue

    if contam_af <= minContamAlleleAf:
      if gt == 0'i8:
        result.n_hom_ref_usable.dec
      else:
        result.n_hom_alt_usable.dec
      continue

    total += infiltrating / (contam_af * dp)
    result.n_sites_usable.inc

  if result.n_sites_usable == 0:
    result.contamination = NaN
  else:
    result.contamination = total / result.n_sites_usable.float64
