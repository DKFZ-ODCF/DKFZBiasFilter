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
import optparse
import sys

import pysam
import yaml

def main(filter_vcf, opts):
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
                if len(rec.ref) == 1 and len(rec.alts[0]) == 1:
                    change = "%s->%s" % (rec.ref, rec.alts[0])
                    if "bPcr" in set(rec.filter):
                        changes_pcr[change] += 1
                        depths_pcr = damage_support(rec, depths_pcr)
                    if "bSeq" in set(rec.filter):
                        changes_seq[change] += 1
                        seqerror_support(rec, depths_seq)
            elif list(rec.filter) == ["PASS"]:
                depths_pass = pass_support(rec, depths_pass)
    changes = {}
    depths = {}
    for val, change in _organize_changes(changes_pcr, opts):
        changes["%s damage" % change] = val
    depths["damage"] = {}
    for depth, count in sorted(depths_pcr.items()):
        if depth < opts.maxdepth:
            depths["damage"][depth] = count
    for val, change in _organize_changes(changes_seq, opts):
        changes["%s bias" % change] = val
    depths["bias"] = {}
    for depth, count in sorted(depths_seq.items()):
        if depth < opts.maxdepth:
            depths["bias"][depth] = count
    depths["pass"] = {}
    for depth, count in sorted(depths_pass.items()):
        if depth < opts.maxdepth:
            depths["pass"][depth] = count
    out = {"changes": changes, "depths": depths}
    if opts.sample:
        out["sample"] = opts.sample
    out_handle = open(opts.outfile, "w") if opts.outfile else sys.stdout
    yaml.safe_dump(out, out_handle, default_flow_style=False, allow_unicode=False)

def _organize_changes(changes, opts):
    out = []
    total = 0
    for change, val in changes.items():
        total += val
        out.append((val, change))
    out.sort(reverse=True)
    return out[:opts.maxchanges] + [(total, "total")]

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
    parser = optparse.OptionParser("Summary statistics for filtered VCF output\n"
                                   "  dkfzbiasfilter_summarize.py <vcf_file>")

    # Optional Arguments
    parser.add_option('--maxchanges', type='int', default=6,
                      help='Maximum number of changes to include for each sample')
    parser.add_option('--maxdepth', type='int', default=5,
                      help='Maximum read depth for reporting changes.')
    parser.add_option('--sample', type='string',
                      help='Sample name')
    parser.add_option('--outfile', type='string',
                      help='Output file (defaults to stdout)')
    options, args = parser.parse_args()
    if len(args) != 1:
        parser.print_help()
    else:
        main(args[0], options)
