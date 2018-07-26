import strformat
import strutils
import tables

type
  Sample* = ref object
    family_id*: string
    id*: string
    paternal_id: string
    maternal_id: string
    mom*: Sample
    dad*: Sample
    sex: int
    affected: bool
    kids: seq[Sample]
    i:int

proc valueOrn9(v:string): string =
  if v == nil or v == "":
    return "-9"
  return v

proc `$`*(s:Sample): string =
  var
    family_id = s.family_id
    sample_id = s.id
    paternal_id = valueOrn9(s.paternal_id)
    maternal_id = valueOrn9(s.maternal_id)
    affected = "2"
    sex = s.sex
  if not s.affected:
    affected = "1"
  return fmt("{family_id}\t{sample_id}\t{paternal_id}\t{maternal_id}\t{affected}\t{sex}")

proc spouse(s:Sample): Sample {.inline.} =
  if s.kids == nil or s.kids.len == 0: return nil
  var k = s.kids[0]
  if k.dad == s: return k.mom
  return k.dad

proc parse_ped(path: string): seq[Sample] =
  result = new_seq_of_cap[Sample](10)

  var look = newTable[string,Sample]()

  for line in lines(path):
    if line[0] == '#': continue
    var toks = line.strip().split('\t')

    var s = Sample(family_id: toks[0], id: toks[1], kids:new_seq[Sample](), paternal_id: toks[2], maternal_id:toks[3])
    result.add(s)
    look[s.id] = s

  for s in result:
    if s.paternal_id in look:
      s.dad = look[s.paternal_id]
      s.dad.kids.add(s)
    if s.maternal_id in look:
      s.mom = look[s.maternal_id]
      s.mom.kids.add(s)

