import os
import math
import tables
import stats
import strformat
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
    for t in traces:
      if t.y[i].classify != fcNaN:
        rs.push(t.y[i])
    for t in traces.mitems:
      if t.y[i].classify != fcNaN:
        t.y[i] /= rs.mean.float32

proc colNormalize(by_chrom: var TableRef[string, seq[Trace]]) =
  for traces in by_chrom.mvalues:
    traces.colNormalize


proc smooth(trace: var Trace, rsc: var RunningStat, w:int=5) =
  var w = int(w.float / 2)

  for i in 1..<trace.y.high:
    #if trace.y[i - 1].classify == fcNaN: continue
    if trace.y[i].classify == fcNaN: continue
    #if trace.y[i + 1].classify == fcNaN: continue
    rsc.push(trace.y[i])

    var S = 0'f32
    var n = 0'f32
    for k in max(0, i - w)..min(trace.y.high, i + w):
      if trace.y[k].classify == fcNaN: continue
      S += trace.y[k]
      n += 1
    if n > 1:
      trace.y[i] = S / n

proc smooth(by_chrom: var TableRef[string, seq[Trace]], rscs: var seq[RunningStat], w:int=9) =
  for traces in by_chrom.mvalues:
    for i, t in traces.mpairs:
      t.smooth(rscs[i], w)

proc drop(trace: var Trace, chrom: string, rsc: RunningStat, z:float64=1.3) =
  var ys = newSeq[float32]()
  var xs = newSeq[int]()

  #[
  var rsc: RunningStat

  for y in trace.y:
    if y.classify != fcNaN:
      rsc.push(y)

  var lo = rsc.mean - 1.5 * rsc.standardDeviation
  var hi = rsc.mean + 1.5 * rsc.standardDeviation
  echo lo, " ", hi
  ]#
  var lo = min(rsc.mean - z * rsc.standardDeviation, 0.8)
  var hi = max(rsc.mean + z * rsc.standardDeviation, 1.3)

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

    ys.add(min(2.2, y))
    xs.add(i)
    lastNa = false
  trace.x = xs
  trace.y = ys

proc drop(by_chrom: var TableRef[string, seq[Trace]], rscs: var seq[RunningStat], z:float64) =
  for chrom, traces in by_chrom.mpairs:
    for i, t in traces.mpairs:
      t.drop(chrom, rscs[i], z)

proc read_extracted*(path: string, cnt: var counts) =
  # read a single sample, used by versus
  var f = newFileStream(path, fmRead)
  if f == nil:
    raise newException(IOError, "couldn't open sketch file:" & path)
  var sl: uint8 = 0
  discard f.readData(sl.addr, sizeof(sl))
  doAssert sl == formatVersion, &"expected matching versions got {sl}, expected {formatVersion}"

  discard f.readData(sl.addr, sizeof(sl))
  cnt.sample_name = newString(sl)
  var n_sites: uint16
  var nx_sites: uint16
  var ny_sites: uint16


  discard f.readData(cnt.sample_name[0].addr, sl.int)
  discard f.readData(n_sites.addr, n_sites.sizeof.int)
  discard f.readData(nx_sites.addr, nx_sites.sizeof.int)
  discard f.readData(ny_sites.addr, ny_sites.sizeof.int)
  if cnt.sites.len.uint16 != n_sites:
    cnt.sites = newSeq[allele_count](n_sites)
  if cnt.x_sites.len.uint16 != nx_sites:
    cnt.x_sites = newSeq[allele_count](nx_sites)
  if cnt.y_sites.len.uint16 != ny_sites:
    cnt.y_sites = newSeq[allele_count](ny_sites)
  if nsites > 0'u16:
    doAssert n_sites.int * sizeof(cnt.sites[0]) == f.readData(cnt.sites[0].addr, nsites.int * sizeof(cnt.sites[0]))
  if nxsites > 0'u16:
    doAssert nx_sites.int * sizeof(cnt.x_sites[0]) == f.readData(cnt.x_sites[0].addr, nx_sites.int * sizeof(cnt.x_sites[0]))
  if nysites > 0'u16:
    doAssert ny_sites.int * sizeof(cnt.y_sites[0]) == f.readData(cnt.y_sites[0].addr, ny_sites.int * sizeof(cnt.y_sites[0]))
  f.close()

proc depth_main*() =

  var argv = commandLineParams()
  if argv.len == 0: argv = @["-h"]
  if argv[0] == "depthview": argv = argv[1..argv.high]

  var p = newParser("somalier depth"):
    help("depth plot on somalier-extracted data")
    option("-w", "--window", default="5", help="smoothing window")
    option("-v", "--vcf", help="sites vcf file to indicate genomic positions")
    option("-o", default="somalier-depthview.html", help="html output file")
    arg("extracted", nargs= -1, help="$sample.somalier files for each sample.")

  var opts = p.parse(argv)
  if opts.help:
    quit 0

  if opts.extracted.len == 0:
    echo p.help
    quit "send argument for extracted files"

  var
    sites = readSites(opts.vcf)
    text_by_chrom = newTable[string, seq[string]]()
    w = parseInt(opts.window)
  for s in sites:
    if s.chrom notin text_by_chrom: text_by_chrom[s.chrom] = newSeq[string]()
    text_by_chrom[s.chrom].add(&"{s.chrom}:{s.position}")

  var cnt : counts
  var by_chrom = newTable[string, seq[Trace]]()
  var samples = newSeq[string]()
  for i, f in opts.extracted:
    read_extracted(f, cnt)
    by_chrom.sample_normalize(cnt, sites)
    samples.add(cnt.sample_name)
  var rscs = newSeq[RunningStat](opts.extracted.len)

  var z = 1.3
  if samples.len > 1500:
    z = 1.5

  by_chrom.colNormalize()
  by_chrom.smooth(rscs, w)
  by_chrom.drop(rscs, z)

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





