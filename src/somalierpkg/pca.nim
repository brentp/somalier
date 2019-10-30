import os
import times
import sequtils
import random
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
  #randomize()

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
  randomize()
  shuffle(opts.extracted)
  for i, f in opts.extracted:


    read_extracted(f, cnt)
    ys.add(orders[labels[cnt.sample_name]])

    if i > 0: doAssert cnt.sites.len == mat[0].len, &"bad number of sites {cnt.sites.len} in {f}"
    var vec = newSeq[float32](cnt.sites.len)
    for j, ac in cnt.sites:
      vec[j] = ac.ab(5).alts.float32
    mat[i] = vec


  var T = mat.toTensor()
  let Y = ys.toTensor() #.astype(float32)#.unsqueeze(0).transpose
  echo T.shape, " ", Y.shape
  var t0 = cpuTime()
  var res = T.pca(8)
  echo "time for pca:", cpuTime() - t0
  echo "components shape:", res.components.shape
  echo "results shape:", res.projected.shape
  res.components.write_npy("comps.npy")
  var R = res.projected

  let
    ctx = newContext Tensor[float32]
    nHidden = 32
    nOut = ys.toHashSet.len
    X = ctx.variable R

  network ctx, AncestryNet:
    layers:
      x: Input([1, R.shape[1]])
      fc1: Linear(R.shape[1], nHidden)
      classifier: Linear(nHidden, nOut)

    forward x:
      x.fc1.relu.classifier

  let
    model = ctx.init(AncestryNet)
    optim = model.optimizerSGD(learning_rate = 0.005'f32)
    batch_size = 32
    t00 = cpuTime()

  for epoch in 0..<10000:

    for batch_id in 0..<(X.value.shape[0] - 100) div batch_size:
      let offset = batch_id * batch_size
      let x = X[offset ..< offset + batch_size, _]
      let y = Y[offset ..< offset + batch_size]

      let
        clf = model.forward(x)
        #loss = mse_loss(y_pred, y)
        loss = clf.sparse_softmax_cross_entropy(y)


      if batch_id == 0 and epoch mod 100 == 0:
        ctx.no_grad_mode:
          let ypred = model.forward(X[2332..<2504, _]).value.softmax.argmax(axis=1).squeeze
        echo "accuracy on unseen data:", accuracy_score(Y[2332..<2504], y_pred)
        #echo "true:", ys[0..200]
        #echo "pred:", y_pred.toSeq
        echo "Epoch is: " & $epoch, " batch is:", batch_id
        echo "Loss is:  " & $loss.value.data
        echo "total time:", cpuTime() - t00

      loss.backprop()
      optim.update()


  ctx.no_grad_mode:
    let ypred = model.forward(X[2432..<2504, _]).value.softmax.argmax(axis=1).squeeze
    echo "accuracy on unseen data:", accuracy_score(Y[2432..<2504], y_pred)


  t0 = cpuTime()
  var proj = T * res.components
  echo "time to project:", cpuTime() - t0
  proj.write_npy("proj.npy")
  echo "proj shape:", proj.shape
  #for i in 0..<proj.shape[0]:
  #  echo proj[i, _]

when isMainModule:
  pca_main()

