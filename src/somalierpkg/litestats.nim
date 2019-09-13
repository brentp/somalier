import math

type
  RunningStat* = object             ## an accumulator for statistical data
    n*: int                         ## number of pushed data
    mean*:float
    mom2*:float

proc push*(s: var RunningStat, x: float) {.inline.} =
  ## pushes a value `x` for processing
  inc(s.n)
  # See Knuth TAOCP vol 2, 3rd edition, page 232
  let delta = x - s.mean
  let delta_n = delta / toFloat(s.n)
  s.mean += delta_n

  let term1 = delta * delta_n * toFloat(s.n - 1)
  s.mom2 += term1

proc excl*(s: var RunningStat, x: float) {.inline.} =
  # remove a value from the runnint total
  doAssert s.n > 0
  dec(s.n)
  let delta = s.mean - x
  let delta_n = delta / toFloat(s.n)
  s.mean += delta_n
  let term1 = delta * delta_n * toFloat(s.n + 1)
  s.mom2 -= term1


proc push*(s: var RunningStat, x: int) {.inline.} =
  s.push(x.toFloat)

proc standardDeviation*(s: RunningStat): float =
  result = sqrt(s.mom2 / toFloat(s.n))

when isMainModule:

  var s: RunningStat
  s.push(2);
  s.push(4);
  s.push(7);
  var std247 = s.standardDeviation

  s.push(6);

  s.excl(6);
  s.excl(7)
  doAssert abs(s.mean - 3'f32) < 1e-5
  s.push(7)
  doAssert abs(std247 - s.standardDeviation) < 1e-5
