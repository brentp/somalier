import sys
import collections

def to_exons(path):
    D = collections.defaultdict(list)
    for line in open(path):
        if line[0] == '#': continue
        toks = line.strip().split("\t")
        if toks[2] != "exon": continue

        info = dict(tuple([y.strip('"') for y in x.split(" ")]) for x in toks[8].split("; "))
        assert info["gene_id"] == info["transcript_id"], info
        D[info["gene_id"]].append((toks[0], int(toks[3]), int(toks[4])))
    E = {}

    # convert list into gene_$i
    for k, vs in D.items():
        for i, loc in enumerate(vs, start=1):
            E[k + "_" + str(i)] = loc
    return E

def to_bed(by_exon, f):
    for i, line in enumerate(open(f)):
        if i == 0:
            print("chrom\tstart\tend\t" + line)
            continue
        (n, line) = line.split("\t", 1)
        chr, start, end = by_exon[n]
        print("%s\t%d\t%d\t%s\t%s" % (chr, start, end, n, line), end="")

if __name__ == "__main__":
    by_exon = to_exons(sys.argv[1])

    to_bed(by_exon, sys.argv[2])

