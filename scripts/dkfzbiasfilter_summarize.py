#!/usr/bin/env python
"""Summarize filtered variants, quantifying known damage artifacts.

DKFZ bias filter looks at two types of artifacts:

- Strand bias (bSeq), imbalance between forward and reverse reads; due to sequencing errors
- PCR bias (bPCR) imbalance between read 1 and read 2, due to DNA damage:
   Oxidative (oxoG) damage -- G->T (read 1), C->A (read 2)
   FFPE deamination -- C->T, G->A


INFO fields have two values:
  - ACGTNactgnPLUS -- read pairs oriented F1 R2
  - ACGTNactgnMINUS -- read pairs oriented F2 R1
Within each there are two sets of values:
  - ACGTN -- sequencing forward strand (F1 for PLUS, F2 for MINUS)
  - acgtn -- sequencing reverse strand (R2 for PLUS, R1 for MINUS)

Usage:
   summarize_filters.py <filtered.vcf>
"""
from __future__ import print_function
import collections
import sys

import pysam

def main(filter_vcf):
    changes_pcr = collections.defaultdict(int)
    changes_seq = collections.defaultdict(int)
    depths_pcr = collections.defaultdict(int)
    depths_seq = collections.defaultdict(int)
    depths_pass = collections.defaultdict(int)
    with pysam.VariantFile(filter_vcf) as bcf_in:
        damage_filters = set([f.name for f in bcf_in.header.filters.values()
                              if f.description.startswith("Variant allele shows a bias")])
        for rec in bcf_in.fetch():
            if len(set(rec.filter) & damage_filters) == len(rec.filter):
                change = "%s->%s" % (rec.ref, rec.alts[0])
                if "bPcr" in set(rec.filter):
                    changes_pcr[change] += 1
                    depths_pcr = damage_support(rec, depths_pcr)
                if "bSeq" in set(rec.filter):
                    changes_seq[change] += 1
                    seqerror_support(rec, depths_seq)
            elif list(rec.filter) == ["PASS"]:
                depths_pass = pass_support(rec, depths_pass)
    print("* Damage")
    for change, val in changes_pcr.items():
        if val > 25:
            print(change, val)
    for depth, count in sorted(depths_pcr.items()):
        if count > 25:
            print(depth, count)
    print("* Bias")
    for change, val in changes_seq.items():
        if val > 25:
            print(change, val)
    for depth, count in sorted(depths_seq.items()):
        if count > 25:
            print(depth, count)
    print("* Pass")
    for depth, count in sorted(depths_pass.items()):
        if count > 25 and depth < 15:
            print(depth, count)

def pass_support(rec, depths):
    summarized = dict(summarize_support(rec))
    alt_depth = summarized["forward"].get(rec.alts[0], 0) + summarized["reverse"].get(rec.alts[0], 0)
    depths[alt_depth] += 1
    return depths

def seqerror_support(rec, depths):
    """Calculate support for variant based on strand bias sequencing errors.

    Summarize bias of calls between F1/F2 and R1/R2
    """
    summarized = dict(summarize_support(rec))
    alt_depth = summarized["forward"].get(rec.alts[0], 0) + summarized["reverse"].get(rec.alts[0], 0)
    depths[alt_depth] += 1
    return depths

def damage_support(rec, depths):
    """Calculate support for variant based on likely damage PCR bias.

    For damage support, summarize bias of calls between F1/R1 and F2/R2
    """
    summarized = dict(summarize_support(rec))
    alt_depth = summarized["read1"].get(rec.alts[0], 0) + summarized["read2"].get(rec.alts[0], 0)
    depths[alt_depth] += 1
    return depths

def summarize_support(rec):
    """Provide a summary of count support across first/second forward/reverse.
    """
    support = collections.defaultdict(lambda: collections.defaultdict(int))
    baseorder = "ACGTN"
    fields = list(baseorder) + [x.lower() for x in baseorder]
    for (orientation, fwd_key, rev_key) in [("PLUS", "F1", "R2"), ("MINUS", "F2", "R1")]:
        key = baseorder + baseorder.lower() + orientation
        for i, (base, val) in enumerate(zip(fields, rec.info[key])):
            read_key = fwd_key if i < len(baseorder) else rev_key
            support[read_key][base.upper()] += int(val)
    summarized = []
    for name, cur, keys in [("read1", collections.defaultdict(int), ["F1", "R1"]),
                            ("read2", collections.defaultdict(int), ["F2", "R2"]),
                            ("forward", collections.defaultdict(int), ["F1", "F2"]),
                            ("reverse", collections.defaultdict(int), ["R1", "R2"])]:
        for base in list("ACGT"):
            for k in keys:
                if support[k][base] > 0:
                    cur[base] += support[k][base]
        summarized.append((name, dict(cur)))
    return summarized

if __name__ == "__main__":
    main(*sys.argv[1:])
