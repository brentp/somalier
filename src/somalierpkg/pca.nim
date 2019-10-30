import os
import times
import strformat
import strutils
import arraymancer
import ./relate
import ./depthview
import argparse
import sets
import arraymancer
import tables

proc read_labels(path: string): TableRef[string, string] =
  result = newTable[string, string]()
  if path == "": return
  for line in path.lines:
    var toks = line.strip().split("\t")
    result[toks[0]] = toks[1]

proc getLabelOrders(l:TableRef[string, string]): TableRef[string, int] =
  result = newTable[string, int]()
  var ancs: seq[string]
  for sample, ancestry in l:
    var idx = ancs.find(ancestry)
    if idx == -1:
      ancs.add(ancestry)
      idx = ancs.high
    result[ancestry] = idx

proc pca_main*() =

  var argv = commandLineParams()
  if argv.len == 0: argv = @["-h"]
  if argv[0] == "pca": argv = argv[1..argv.high]

  var p = newParser("somalier pca"):
    help("dimensionality reduction")
    option("--labels", help="file with ancestry labels")
    arg("extracted", nargs= -1, help="$sample.somalier files for each sample.")

  var opts = p.parse(argv)
  if opts.help:
    quit 0

  var labels = read_labels(opts.labels)

  if opts.extracted.len == 0:
    echo p.help
    quit "send argument for extracted files"

  var n_samples = opts.extracted.len
  var mat = newSeq[seq[float32]](n_samples)

  var orders = getLabelOrders(labels)
  var ys = newSeq[int]()

  var cnt : counts
  for i, f in opts.extracted:


    read_extracted(f, cnt)
    ys.add(orders[labels[cnt.sample_name]])

    if i > 0: doAssert cnt.sites.len == mat[0].len, &"bad number of sites {cnt.sites.len} in {f}"
    var vec = newSeq[float32](cnt.sites.len)
    for j, ac in cnt.sites:
      vec[j] = ac.ab(5).alts.float32
    mat[i] = vec

  var T = mat.toTensor#.transpose
  echo "y:", ys[0..20]
  let y = ys.toTensor #.astype(float32)#.unsqueeze(0).transpose
  echo T.shape, " ", y.shape



  let
    ctx = newContext Tensor[float32]
    nHidden = 60
    nOut = ys.toHashSet.len
    x = ctx.variable T

  network ctx, AncestryNet:
    layers:
      fc1: Linear(T.shape[1], nHidden)
      classifier: Linear(nHidden, nOut)

    forward x:
      x.fc1.relu.classifier

  let
    model = ctx.init(AncestryNet)
    optim = model.optimizerSGD(learning_rate = 0.01'f32)


  for epoch in 0..<100:

    let
      y_pred = model.forward(x)
      #loss = mse_loss(y_pred, y)
      loss = y_pred.sparse_softmax_cross_entropy(y)

    echo y_pred.value.data[0..20]
    echo "Epoch is: " & $epoch
    echo "Loss is:  " & $loss.value.data

    loss.backprop()
    optim.update()



  echo "T shape:", T.shape
  var t0 = cpuTime()
  var res = T.pca(5)
  echo "time for pca:", cpuTime() - t0
  echo "components shape:", res.components.shape
  echo "results shape:", res.projected.shape
  res.components.write_npy("comps.npy")

  t0 = cpuTime()
  var proj = T * res.components
  echo "time to project:", cpuTime() - t0
  proj.write_npy("proj.npy")
  echo "proj shape:", proj.shape
  #for i in 0..<proj.shape[0]:
  #  echo proj[i, _]

when isMainModule:
  pca_main()

