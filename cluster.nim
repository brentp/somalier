import strutils

type centroid = ref object
  key: int
  index: int
  dataIndexes: seq[int]

proc `$`*(c:centroid): string =
  return "centroid(key:$1, index:$2, indexes: $3)" % [$c.key, $c.index, $c.dataIndexes]

type Cluster[T] = ref object
  centroids*: seq[centroid]
  index: seq[centroid]
  dists: seq[seq[T]]
  mins: seq[int]

proc mergeClosest[T](cluster: var Cluster[T], threshold: T): bool =
  var 
    minKey = 0
    minDist = T.high

  for i, c in cluster.centroids:
    var dist = cluster.dists[c.key][cluster.mins[c.key]]
    if dist < minDist:
      minKey = c.key
      minDist = dist

  if minDist > threshold:
    return false

  var
    c1 = cluster.index[minKey]
    c2 = cluster.index[cluster.mins[minKey]]

  var merged = centroid(key: c1.key, dataIndexes:c1.dataIndexes)
  merged.dataIndexes.add(c2.dataIndexes)
  cluster.centroids[c1.index] = merged
  cluster.centroids.delete(c2.index)
  cluster.index[c1.key] = merged


  # update distances
  for i, ci in cluster.centroids:
    if c1.key == ci.key:
      cluster.dists[c1.key][ci.key] = T.high
      cluster.dists[ci.key][c1.key] = T.high
      continue


    # min == single linkage, max == complete linkage
    var dist = max(cluster.dists[c1.key][ci.key], cluster.dists[c2.key][ci.key])
    cluster.dists[c1.key][ci.key] = dist
    cluster.dists[ci.key][c1.key] = dist


  # update mins
  for i, ci in cluster.centroids:
    cluster.centroids[i].index = i
    if cluster.mins[ci.key] == c1.key or cluster.mins[ci.key] == c2.key:
      var mk = ci.key
      for j, cj in cluster.centroids:
        if cluster.dists[ci.key][cj.key] < cluster.dists[ci.key][mk]:
          mk = cj.key
      cluster.mins[ci.key] = mk
  return true


proc hcluster*[T](dists: seq[seq[T]], threshold:T): seq[seq[int]] =
  var
    cluster = Cluster[T](
                      centroids : newSeq[centroid](dists.len),
                      index : newSeq[centroid](dists.len),
                      mins: newSeq[int](dists.len),
                      dists: dists)

  for i in 0..<dists.len:
    cluster.centroids[i] = centroid(index:i, key:i, dataIndexes: @[i])
    cluster.index[i] = cluster.centroids[i]
    var isset = false
    for j, d in dists[i]:
      if i == j: continue
      if d < dists[i][cluster.mins[i]] or not isset:
        cluster.mins[i] = j
        isset = true


  var merged = cluster.mergeClosest(threshold)

  while merged:
    merged = cluster.mergeClosest(threshold)

  result = newSeq[seq[int]](cluster.centroids.len)
  for i, c in cluster.centroids:
    result[i] = c.dataIndexes

