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

proc push*(s: var RunningStat, x: int) {.inline.} =
  s.push(x.toFloat)

proc standardDeviation*(s: RunningStat): float =
  result = sqrt(s.mom2 / toFloat(s.n))
