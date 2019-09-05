import math

const arraySize = 1024

type counter = object
  counts: array[arraySize, int]

proc push(c:var counter, x:float) {.inline.} =
  doAssert 0 <= x and x <= 0.5
  var idx = min(c.counts.high, int(x * (arraySize * 2.0) + 0.5))
  c.counts[idx].inc

proc median(c:counter): float {.inline.} =
  var n = sum(c.counts)
  if n < 1: return 0
  var mid = int(n.float / 2.0 + 0.5)
  var fsum = 0
  for i, v in c.counts:
    fsum += v
    if fsum >= mid:
      return i.float / (arraySize * 2.0)

type AD = object
  nref: int16
  nalt: int16

template depth(a:AD): int =
  a.nref.int + a.nalt.int

template ab(a:AD): float32 =
  # allele balance
  if a.depth > 0: a.nalt.int / a.depth else: -1

const error_rate = 2e-3

proc estimate_contamination*(self_abs: seq[float32], other_abs: seq[float32]): (float, int) =
  ## estimate contamination of self, by other.
  var sites_used = 0
  var c: counter
  for i, a in self_abs:
    if a == -1: continue
    var b = other_abs[i]
    if b == -1: continue
    if abs(a - b) < 0.015: continue

    sites_used += 1
    # aa: 0.5, bb: 0 -> 2
    # aa: 1,   bb: 0 -> 1
    # aa: 0.01, bb: 0.03 -> 2
    if a > 0.25 and a < 0.75: continue

    var scaler = min(2'f32, 1'f32 / abs(a - b).float32)
    # a: 0.01, b: 0.0001 then b can't contribute to a
    if a < 0.5 and b < a:
      scaler = 0

    # a: 0.99, b: 0.999 then reads can't come from b
    elif a > 0.5 and b > a:
      scaler = 0

    var ax = if a > 0.5: 1-a else: a
    if ax > 0.25:
      # remove cases where direction of support does not match
      # e.g.:
      # self:0.4744310677051544 other:0.9432255029678345 ax:0.02556893229484558 scaler:2.0 evidence for contamination of:0.05113786458969116
      if a < 0.5 and b > 0.5:
        scaler = 0.0 # no way that change in AF can come from b
      ax = 0.5 - ax

    var contam = scaler * ax

    #echo "self:", a, " other:", b, " ax:", ax, " scaler:", scaler, " evidence for contamination of:", contam

    c.push(contam) #scaler * (if a < 0.5: a else: 1 - a))
  #echo sum(c.counts), " ", c.counts[0..<min(arraySize, 100)]
  return (c.median, sum(c.counts))


proc estimate_contamination(a:AD, b:AD, min_depth:int=7): float32 {.inline.} =
  if a.depth < min_depth: return -1
  if b.depth < min_depth: return -1

  var aab = a.ab
  # TODO: make this depend on depth
  if aab > 0.25 and aab < 0.75: return -1

  var bab = b.ab
  # TODO: make this cutoff depend on depth.
  if abs(aab - bab) < 0.015: return -1

  var scaler = min(2'f32, 1'f32 / abs(aab - bab).float32)

  # a: 0.01, b: 0.0001 then b can't contribute to a
  if aab < 0.5 and bab < aab:
    return 0
  # a: 0.99, b: 0.999 then reads can't come from b
  elif aab > 0.5 and bab > aab:
    return 0

  var ax = if aab > 0.5: 1-aab else: aab
  if ax > 0.25:
      # remove cases where direction of support does not match
      # e.g.:
      # self:0.4744310677051544 other:0.9432255029678345 ax:0.02556893229484558 scaler:2.0 evidence for contamination of:0.05113786458969116
      if aab < 0.5 and bab > 0.5:
        return 0
      ax = 0.5 - ax

  return scaler * ax


proc estimate_contamination(aads: seq[AD], bads: seq[AD]): (float, int) =
  ## estimate contamination of self, by other.
  var sites_used = 0
  var c: counter
  doAssert aads.len == bads.len
  for i, a in aads:

    var contam = estimate_contamination(a, bads[i])
    if contam < 0: continue
    c.push(contam)
  echo sum(c.counts), " ", c.counts[0..<min(arraySize, 1000)]
  return (c.median, sum(c.counts))


when isMainModule:
  import unittest
  import random
  import alea

  const n_sites = 4000
  # number of sites mutated in the tumor
  # make it so when this is 0, there is no differnece found
  const n_mutated_sites = 0
  const mutated_af = 0.15
  const mean_depth = 600
  randomize()

  var bern = bernoulli(0.49)
  var rnd : Random
  rnd.random = proc(): float = random(1'f)

  proc contaminate_with(bb:var seq[AD], aa:var seq[AD], p:float=0.05) =
    # contaminate b with a
    for i, a in aa:
      bb[i].nref = int16(0.5 + (1 - p) * a.nref.float + p * bb[i].nref.float)
      bb[i].nalt = int16(0.5 + p * a.nalt.float + p * bb[i].nalt.float)

  proc getad(gt:int, dp:int, af_adjust:float32): AD =
    var adj = bernoulli(af_adjust)
    if gt == 0:
      result = AD(nref: dp.int16, nalt: 0)
      if af_adjust > 0:
        result.nalt += int16(0.5 + dp.float * rnd.mean(adj, samples=dp))
    elif gt == 1:
      result = AD(nref: 0.int16, nalt: dp.int16)
      if af_adjust > 0:
        result.nref += int16(0.5 + dp.float * rnd.mean(adj, samples=dp))
    elif gt == -1:
      result = AD(nref: (dp/2).int16, nalt: (dp/2).int16)
    else:
      # het
      result = AD()
      result.nref = int16(0.5 + dp.float * rnd.mean(bern, samples=dp))
      result.nalt = dp.int16 - result.nref
      if af_adjust > 0:
        result.nref += int16(0.5 + dp.float * rnd.mean(adj, samples=dp))

    # simulate errors
    if random(1'f) * dp.float > error_rate / 2:
      if random(1'f) < 0.5:
        if result.nref > 0:
          result.nref.dec
        result.nalt.inc
      else:
        if result.nalt > 0:
          result.nalt.dec
        result.nref.inc


  var mut_prob = n_mutated_sites.float / n_sites.float
  var normal = newSeq[AD](n_sites)
  var tumor = newSeq[AD](n_sites)
  # 0, 1, 2, or -1
  var normal_genotype = newSeq[int8](n_sites)
  var gt_choices = @[0, 1, 1, 2]
  for i in 0..2:
    gt_choices.add(gt_choices)
  gt_choices.add(-1)
  var mutated_idx = newSeq[int]()

  for i in 0..<n_sites:
    var dp = int(0.5 + random(0.5*mean_depth..1.6*mean_depth))
    var gt = random(gt_choices)
    normal_genotype[i] = gt.int8
    normal[i] = getad(gt, dp, 0.0'f32)
    var adj = 0'f32
    if random(1'f) < mut_prob:
      adj = mutated_af
      mutated_idx.add(i)

    tumor[i] = getad(gt, int(0.5 + random(0.5*mean_depth..1.6*mean_depth)), adj)

  tumor.contaminate_with(normal, p=0.01)
  for i in mutated_idx:
    echo normal[i], " ", tumor[i]
  echo mutated_idx.len
  echo estimate_contamination(tumor, normal)

