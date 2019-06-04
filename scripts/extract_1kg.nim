import os
import math
import hts/vcf
import strutils
import strformat

var ivcf:VCF
var sites_vcf:VCF
var ovcf:VCF

var sites_vcf_path = paramStr(1)
var vcf_tmpl = paramStr(2)

if not open(ivcf, vcf_tmpl.replace("{}", "1")):
  quit "couldn't open file"

if not open(ovcf, "extracted-1kg.bcf", mode="w"):
  quit "couldn't open file"

ovcf.copy_header(ivcf.header)

if not open(sites_vcf, sites_vcf_path):
  quit "couldn't open sites vcf"

doAssert ovcf.write_header

type Site* = object
  A_allele*: string
  B_allele*: string
  chrom*: string
  position*: int

proc get_variant(ivcf:VCF, site:Site): Variant =
  for v in ivcf.query(&"{site.chrom}:{site.position+1}-{site.position+2}"):
    if v.start == site.position and (
      (v.REF == site.A_allele and v.ALT[0] == site.B_allele) or
      (v.REF == site.B_allele and v.ALT[0] == site.A_allele)):
      return v.copy()

var last_rid = -1
for v in sites_vcf:
  var chrom = $v.CHROM
  if chrom.startswith("chr"):
    chrom = chrom[3..chrom.high]
  if last_rid != v.rid:
    ivcf.close
    echo "opening new chrom:", $v.CHROM
    if not open(ivcf, vcf_tmpl.replace("{}", chrom), threads=2):
      quit "couldn't open input vcf file"
    last_rid = v.rid

  var site = Site(A_allele: v.REF, B_allele:v.ALT[0], chrom: $v.CHROM, position: v.start)
  var variant = ivcf.get_variant(site)
  if variant == nil:
    echo "skipping:", site
    continue
  echo "writing:", site

  doAssert ovcf.write_variant(variant)

ovcf.close
