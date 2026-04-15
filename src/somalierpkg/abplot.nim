import os
import json
import strformat
import argparse
import algorithm
import strutils
import sequtils
import ./common
import ./depthview
import ./relate

const tmpl_html = staticRead("abplot.html")

type SiteBucket = enum
  sbAutosome, sbX, sbY

type SiteLookup = object
  bucket: SiteBucket
  idx: int

proc jsonValue(x: float): JsonNode =
  if x < 0:
    result = newJNull()
  else:
    result = %x

proc toSiteLabels(sites: seq[Site]): seq[string] =
  result = newSeqOfCap[string](sites.len)
  for site in sites:
    result.add(&"{site.chrom}:{site.position + 1} {site.A_allele}>{site.B_allele}")

proc toIndices(n: int): seq[int] =
  result = newSeq[int](n)
  for i in 0..<n:
    result[i] = i

proc chromSortKey(chrom: string): tuple[group: int, num: int, label: string] =
  var c = chrom
  if c.startsWith("chr"):
    c = c[3..c.high]
  if c in ["NC_000023.10", "NC_000023.11"]:
    c = "X"
  elif c in ["NC_000024.9", "NC_000024.10"]:
    c = "Y"
  if c.len > 0 and c.allCharsInSet({'0'..'9'}):
    return (0, parseInt(c), chrom)
  if c == "X":
    return (1, 23, chrom)
  if c == "Y":
    return (1, 24, chrom)
  return (2, 0, chrom)

proc displayOrder(sites: seq[Site]): seq[int] =
  result = @[]
  for i, site in sites:
    result.add(i)
  sort(result, proc(a, b: int): int =
    let sa = sites[a]
    let sb = sites[b]
    result = cmp(chromSortKey(sa.chrom), chromSortKey(sb.chrom))
    if result == 0:
      result = cmp(sa.position, sb.position)
  )

proc chromosomeTicks(sites: seq[Site]): JsonNode =
  var tickVals = newJArray()
  var tickText = newJArray()
  var startIdx = 0
  var idx = 0
  var currentChrom = ""

  for site in sites:
    if currentChrom.len == 0:
      currentChrom = site.chrom
      startIdx = idx
    elif site.chrom != currentChrom:
      tickVals.add(%startIdx)
      tickText.add(%currentChrom)
      currentChrom = site.chrom
      startIdx = idx
    inc idx

  if currentChrom.len > 0:
    tickVals.add(%startIdx)
    tickText.add(%currentChrom)

  result = %*{
    "tickmode": "array",
    "tickvals": tickVals,
    "ticktext": tickText
  }

proc buildLookups(sites: seq[Site]): seq[SiteLookup] =
  result = newSeq[SiteLookup](sites.len)
  var ai = 0
  var xi = 0
  var yi = 0
  for i, site in sites:
    case site.chrom
    of "X", "chrX", "NC_000023.10", "NC_000023.11":
      result[i] = SiteLookup(bucket: sbX, idx: xi)
      inc xi
    of "Y", "chrY", "NC_000024.9", "NC_000024.10":
      result[i] = SiteLookup(bucket: sbY, idx: yi)
      inc yi
    else:
      result[i] = SiteLookup(bucket: sbAutosome, idx: ai)
      inc ai

proc getCount(cnt: counts, lookup: SiteLookup): allele_count =
  case lookup.bucket
  of sbAutosome:
    result = cnt.sites[lookup.idx]
  of sbX:
    result = cnt.x_sites[lookup.idx]
  of sbY:
    result = cnt.y_sites[lookup.idx]

proc addTrace(traces: var seq[JsonNode], sampleName, metricName: string,
    x: seq[int], y: JsonNode, labels: seq[string], xaxis: string, yaxis: string,
    opacity: float = 1.0) =
  traces.add(%*{
    "type": "scattergl",
    "mode": "markers",
    "name": sampleName & " " & metricName,
    "x": x,
    "y": y,
    "text": labels,
    "opacity": opacity,
    "marker": {
      "size": 3
    },
    "xaxis": xaxis,
    "yaxis": yaxis,
    "hovertemplate": "%{fullData.name}<br>index=%{x}<br>site=%{text}<br>value=%{y}<extra></extra>"
  })

proc AB_plot_main*() =
  var argv = commandLineParams()
  if argv.len == 0:
    argv = @["-h"]
  if argv[0] == "AB_plot":
    argv = argv[1..argv.high]

  var p = newParser("somalier AB_plot"):
    help("plot allele balance and depth from one or more .somalier files")
    option("-s", "--sites", help="sites vcf file used to label plotted site indexes")
    option("-o", "--output", default="somalier-AB-plot.html", help="html output file")
    option("-d", "--min-depth", default="1", help="minimum depth required to show allele balance")
    arg("extracted", nargs = -1, help="$sample.somalier files for each sample.")

  let opts = p.parse(argv)
  if opts.help:
    quit 0

  if opts.sites.len == 0:
    echo p.help
    quit "[somalier] --sites file required"
  if opts.extracted.len == 0:
    echo p.help
    quit "[somalier] send argument for extracted files"

  var sites = readSites(opts.sites)
  let lookups = buildLookups(sites)
  let order = displayOrder(sites)
  var orderedSites = newSeq[Site](order.len)
  for i, j in order:
    orderedSites[i] = sites[j]
  let siteLabels = orderedSites.toSiteLabels()
  let chromTicks = chromosomeTicks(orderedSites)
  let xs = toIndices(siteLabels.len)
  let minDepth = parseInt(opts.min_depth)
  var traces = newSeq[JsonNode]()
  var cnt: counts

  for path in opts.extracted:
    read_extracted(path, cnt)
    if cnt.sites.len != lookups.countIt(it.bucket == sbAutosome):
      quit &"[somalier] expected {lookups.countIt(it.bucket == sbAutosome)} autosomal sites from --sites but found {cnt.sites.len} in {path}"
    if cnt.x_sites.len != lookups.countIt(it.bucket == sbX):
      quit &"[somalier] expected {lookups.countIt(it.bucket == sbX)} X sites from --sites but found {cnt.x_sites.len} in {path}"
    if cnt.y_sites.len != lookups.countIt(it.bucket == sbY):
      quit &"[somalier] expected {lookups.countIt(it.bucket == sbY)} Y sites from --sites but found {cnt.y_sites.len} in {path}"

    var abValues = newJArray()
    var depthValues = newJArray()
    for j in order:
      let ac = cnt.getCount(lookups[j])
      abValues.add(jsonValue(ac.ab(minDepth)))
      depthValues.add(%(ac.depth.int))

    traces.addTrace(cnt.sample_name, "depth", xs, depthValues, siteLabels, "x", "y", 0.55)
    traces.addTrace(cnt.sample_name, "AB", xs, abValues, siteLabels, "x", "y2")

  var html = tmpl_html.replace("<INPUT_TRACES>", $(%traces))
  html = html.replace("<INPUT_XAXIS>", $chromTicks)
  var fh: File
  if not open(fh, opts.output, fmWrite):
    quit "[somalier] couldn't open output file"
  fh.write(html)
  fh.close()
