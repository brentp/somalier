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
    if line[0] == '#': continue
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

proc split_labeled_samples(paths: seq[string]): tuple[labeled: seq[string], test:seq[string]] =
  var after = false
  for p in paths:
    if p == "++":
      after = true
      continue

    if after:
      result.test.add(p)
    else:
      result.labeled.add(p)

proc pca_main*() =

  var argv = commandLineParams()
  if argv.len == 0: argv = @["-h"]
  if argv[0] == "pca": argv = argv[1..argv.high]
  #randomize()

  var p = newParser("somalier pca"):
    help("dimensionality reduction")
    option("--labels", help="file with ancestry labels")
    option("-o", "--output-prefix", help="prefix for output files", default="somalier-ancestry")
    option("--n-pcs", help="number of principal components to use in the reduced dataset", default="12")
    option("--nn-hidden-size", help="shape of hidden layer in neural network", default="16")
    option("--nn-batch-size", help="batch size fo training neural network", default="32")
    option("--nn-test-samples", help="number of labeled samples to test for NN convergence", default="101")
    arg("extracted", nargs= -1, help="$sample.somalier files for each sample. place labelled samples first followed by '++' then *.somalier for query samples")

  var opts = p.parse(argv)
  if opts.help:
    quit 0

  var labels = read_labels(opts.labels)

  if opts.output_prefix.endswith("/"):
    opts.output_prefix &= "/somalier"
  if not opts.output_prefix.endswith("."):
    opts.output_prefix &= "."

  if opts.extracted.len == 0:
    echo p.help
    quit "send argument for extracted files"

  var (labeled_samples, query_samples) = opts.extracted.split_labeled_samples

  var train_mat = newSeq[seq[float32]](labeled_samples.len)
  var query_mat = newSeq[seq[float32]](query_samples.len)

  var orders = getLabelOrders(labels)
  var int_labels = newSeq[int]()
  var labeled_sample_names: seq[string]
  var query_sample_names: seq[string]

  var cnt : counts
  randomize()
  shuffle(labeled_samples)
  for i, f in labeled_samples:

    read_extracted(f, cnt)
    int_labels.add(orders[labels[cnt.sample_name]])
    labeled_sample_names.add(cnt.sample_name)

    if i > 0: doAssert cnt.sites.len == train_mat[0].len, &"bad number of sites {cnt.sites.len} in {f}"
    var vec = newSeq[float32](cnt.sites.len)
    for j, ac in cnt.sites:
      vec[j] = ac.ab(5).alts.float32
    train_mat[i] = vec

  for i, f in query_samples:
    read_extracted(f, cnt)
    query_sample_names.add(cnt.sample_name)
    doAssert cnt.sites.len == train_mat[0].len, &"bad number of sites {cnt.sites.len} in {f}"
    var vec = newSeq[float32](cnt.sites.len)
    for j, ac in cnt.sites:
      vec[j] = ac.ab(5).alts.float32
    query_mat[i] = vec

  var
    nPCs = parseInt(opts.n_pcs)
    T = train_mat.toTensor()
    Y = int_labels.toTensor() #.astype(float32)#.unsqueeze(0).transpose
    t0 = cpuTime()
    res = T.pca(nPCs)
    R = res.projected
  stderr.write_line &"[somalier] time for dimensionality reduction to shape {res.projected.shape}: {cpuTime() - t0:.2f} seconds"

  let
    ctx = newContext Tensor[float32]
    nHidden = parseInt(opts.nn_hidden_size)
    nOut = int_labels.toHashSet.len
    X = ctx.variable R
    nn_test_samples = parseInt(opts.nn_test_samples)

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
    batch_size = parseInt(opts.nn_batch_size)
  t0 = cpuTime()

  # train the model
  for epoch in 0..<1000:

    for batch_id in 0..<X.value.shape[0] div batch_size:

      let offset = batch_id * batch_size
      if offset > X.value.shape[0] - nn_test_samples:
        break

      let offset_stop = min(offset + batch_size,  X.value.shape[0] - nn_test_samples)
      let x = X[offset ..< offset_stop, _]
      let y = Y[offset ..< offset_stop]

      let
        clf = model.forward(x)
        loss = clf.sparse_softmax_cross_entropy(y)

      loss.backprop()
      optim.update()

    if epoch mod 100 == 0:
      ctx.no_grad_mode:
        let
          clf = model.forward(X[X.value.shape[0] - nn_test_samples..<X.value.shape[0], _])
          y_pred = clf.value.softmax.argmax(axis=1).squeeze
          y = Y[X.value.shape[0] - nn_test_samples..<X.value.shape[0]]
          loss = clf.sparse_softmax_cross_entropy(y).value.data[0]
          accuracy = accuracy_score(y, y_pred)
      stderr.write_line &"[somalier] Epoch:{epoch}. loss: {loss:.5f}. accuracy on unseen data: {accuracy:.3f}.  total-time: {cpuTime() - t0:.2f}"
      if epoch >= 200 and ((loss < 0.005 and accuracy > 0.986) or (accuracy >= 0.999 and loss < 0.025)):
        stderr.write_line &"[somalier] breaking with trained model at this accuracy and loss"
        break


  ctx.no_grad_mode:
    let t_probs = model.forward(X).value.softmax #.argmax(axis=1).squeeze
    #echo ypred
    #echo ypred.shape
    #echo ypred.argmax(axis=1).squeeze

  let
    Q = query_mat.toTensor()
    q_proj = Q * res.components
    q_var = ctx.variable q_proj
    q_probs = model.forward(q_var).value.softmax
    q_pred = q_probs.argmax(axis=1).squeeze
    t_pred = t_probs.argmax(axis=1).squeeze


  echo &"reduced query set to: {q_proj.shape}"
  echo "predictions:"
  echo q_pred
  echo "ordered:", orders

  var ancestry_fh: File
  if not open(ancestry_fh, opts.output_prefix & "tsv", fmWrite):
    quit "couldn't open output file {opts.output_prefix & \"tsv\"}"

  var header = @["#sample_id", "predicted_ancestry", "given_ancestry"]
  var inv_orders = newSeq[string](orders.len)
  # maintain order of ancestries
  header.setLen(header.len + orders.len)
  for k, v in orders:
    inv_orders[v] = k
    header[2 + v] = k & "_prob"
  for ip in 0..<nPCs:
    header.add("PC" & $(ip + 1))

  ancestry_fh.write_line(join(header, "\t"))

  for i, s in labeled_sample_names:
    var line = @[s, inv_orders[t_pred[i]], inv_orders[int_labels[i]]]
    for j in 0..<orders.len:
      line.add(formatFloat(t_probs[i, j], ffDecimal, precision=4))
    for j in 0..<nPcs:
      line.add(formatFloat(R[i, j], ffDecimal, precision=4))
    ancestry_fh.write_line(join(line, "\t"))

  for i, s in query_sample_names:
    var line = @[s, inv_orders[q_pred[i]], ""]
    for j in 0..<orders.len:
      line.add(formatFloat(q_probs[i, j], ffDecimal, precision=4))
    for j in 0..<nPcs:
      line.add(formatFloat(q_proj[i, j], ffDecimal, precision=4))
    ancestry_fh.write_line(join(line, "\t"))

  #t0 = cpuTime()
  #var proj = T * res.components
  #echo "time to project:", cpuTime() - t0
  #proj.write_npy("proj.npy")
  #echo "proj shape:", proj.shape
  #for i in 0..<proj.shape[0]:
  #  echo proj[i, _]
  ancestry_fh.close

when isMainModule:
  pca_main()

