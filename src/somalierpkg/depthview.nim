import os
import math
import times
import tables
import stats
import strformat
import arraymancer
import ./common
import json
import ./relate
import argparse

type Trace = object
  y: seq[float32]
  x: seq[int]
  name: string

const tmpl_html = staticRead("depthview.html")


proc sample_normalize(by_chrom: TableRef[string, seq[Trace]], cnt: var counts, sites: seq[Site]) =
  var S: float64
  var n: float64

  for s in cnt.sites:
    if s.depth > 0'u32:
      S += s.depth.float64
      n += 1
  var m = (S / n).float32

  var sample_by_chrom = newTable[string, Trace]()
  for i, c in cnt.sites:
    var s = sites[i]
    if s.chrom notin sample_by_chrom:
      sample_by_chrom[s.chrom] = Trace(name: cnt.sample_name)

    if c.depth == 0:
      sample_by_chrom[s.chrom].y.add(NaN)
    else:
      sample_by_chrom[s.chrom].y.add(c.depth.float32 / m)

  for chrom, trace in sample_by_chrom:
    if chrom notin by_chrom:
      by_chrom[chrom] = newSeq[Trace]()
    by_chrom[chrom].add(trace)


proc colNormalize(traces: var seq[Trace]) =
  for i in 0..<traces[0].y.len:
    var rs: RunningStat
    var rsc: RunningStat
    for t in traces:
      if t.y[i].classify != fcNaN:
        rs.push(t.y[i])
    for t in traces.mitems:
      if t.y[i].classify != fcNaN:
        t.y[i] /= rs.mean.float32
        rsc.push(t.y[i])

    #if rsc.n < 2: return
    #var lo = rsc.mean - 3.5 * rsc.standardDeviation
    #var hi = rsc.mean + 2.0 * rsc.standardDeviation
    #for t in traces.mitems:
    #  if t.y[i].classify != fcNaN and t.y[i] > lo and t.y[i] < hi:
    #    t.y[i] = NaN

proc colNormalize(by_chrom: var TableRef[string, seq[Trace]]) =
  for traces in by_chrom.mvalues:
    traces.colNormalize



proc smooth(trace: var Trace) =

  for i in 1..<trace.y.high:
    if trace.y[i - 1].classify == fcNaN: continue
    if trace.y[i].classify == fcNaN: continue
    if trace.y[i + 1].classify == fcNaN: continue

    if i == 1 or i == trace.y.high - 1 or trace.y[i - 2].classify == fcNaN or trace.y[i + 2].classify == fcNaN:
      trace.y[i] = (trace.y[i-1] + trace.y[i] + trace.y[i + 1])/3'f32
    else:
      trace.y[i] = (trace.y[i-2] + trace.y[i-1] + trace.y[i] + trace.y[i + 1] + trace.y[i + 2])/5'f32


proc smooth(by_chrom: var TableRef[string, seq[Trace]]) =
  for traces in by_chrom.mvalues:
    for t in traces.mitems:
      t.smooth()

proc drop(trace: var Trace, chrom: string, lo: float32, hi:float32) =
  var ys = newSeq[float32]()
  var xs = newSeq[int]()

  var ok = (proc(y: float32): bool =
    return y > lo and y < hi
  )
  var lastNa = false
  var lastHi = false
  var lastLo = false
  for i, y in trace.y:
    if y.classify == fcNaN or ok(y):
      if not lastNa:
        xs.add(i)
        ys.add(NaN)
        lastNa = true
      continue

    if y > 1:
      if not lastHi:
        lastHi = true
        lastLo = false
        xs.add(i)
        ys.add(NaN)
    else:
      if not lastLo:
        lastHi = false
        lastLo = true
        xs.add(i)
        ys.add(NaN)

    ys.add(y)
    xs.add(i)
    lastNa = false
  trace.x = xs
  trace.y = ys

proc drop(by_chrom: var TableRef[string, seq[Trace]], lo: float32, hi: float32) =
  for chrom, traces in by_chrom.mpairs:
    for t in traces.mitems:
      t.drop(chrom, lo, hi)

proc depth_main*() =

  var argv = commandLineParams()
  if argv.len == 0: argv = @["-h"]
  if argv[0] == "depthview": argv = argv[1..argv.high]

  var p = newParser("somalier depth"):
    help("depth plot on somalier-extracted data")
    option("-v", "--vcf", help="sites vcf file to indicate genomic positions")
    option("--lo", default="0.8", help="don't plot values above this threshold and below --hi (helps reduce size of html output)")
    option("--hi", default="1.3", help="don't plot values below this threshold and above --lo (helps reduce size of html output)")
    option("-o", default="somalier-depthview.html", help="html output file")
    arg("extracted", nargs= -1, help="$sample.somalier files for each sample.")

  var opts = p.parse(argv)
  if opts.help:
    quit 0

  if opts.extracted.len == 0:
    echo p.help
    quit "send argument for extracted files"

  var
    lo = parseFloat(opts.lo).float32
    hi = parseFloat(opts.hi).float32

  var sites = readSites(opts.vcf)
  var text_by_chrom = newTable[string, seq[string]]()
  for s in sites:
    if s.chrom notin text_by_chrom: text_by_chrom[s.chrom] = newSeq[string]()
    text_by_chrom[s.chrom].add(&"{s.chrom}:{s.position}")

  var cnt : counts
  var by_chrom = newTable[string, seq[Trace]]()
  for i, f in opts.extracted:
    read_extracted(f, cnt)
    by_chrom.sample_normalize(cnt, sites)

  by_chrom.colNormalize()
  by_chrom.smooth()
  by_chrom.drop(lo, hi)

  var s = %* by_chrom
  var t = tmpl_html.replace("<INPUT_TEXT>", $(%text_by_chrom))
  t = t.replace("<INPUT_JSON>", $s)

  var fh: File
  if not open(fh, opts.o, mode=fmWrite):
    quit "couldn't open output file"
  fh.write(t)
  fh.close()


when isMainModule:
  depth_main()





