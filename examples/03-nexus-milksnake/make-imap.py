#!/usr/bin/env python3
"""Rebuild imap.txt from the taxon labels in loci/*.nex.

The published labels carry the author's own lineage assignment as their
leading token:

    getula_FTB966_NJ_Cumberland      allele 1 of voucher FTB966
    getula_a_FTB966_NJ_Cumberland    allele 2 of the same individual
    SA_FHGO3090_Ecuador_Pinchincha   a South American triangulum lineage

so the Imap is read off the data rather than invented here.  Two outgroup
genera are named with two words and are special-cased.

Usage:  python3 make-imap.py > imap.txt
"""
import glob
import os
import re
import sys

TWO_WORD = ("Arizona_elegans", "Cemophora_coccinea")


def species_of(label):
    for genus in TWO_WORD:
        if label.startswith(genus):
            return genus
    return label.split("_", 1)[0]


def taxon_labels(path):
    text = open(path, errors="replace").read()
    block = re.search(r"taxlabels(.*?);", text, re.S | re.I)
    if not block:
        sys.exit(f"{path}: no taxlabels block")
    return [t.strip() for t in block.group(1).split() if t.strip()]


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    files = sorted(glob.glob(os.path.join(here, "loci", "*.nex")))
    if not files:
        sys.exit("no loci/*.nex found")
    assignment = {}
    for path in files:
        for label in taxon_labels(path):
            assignment.setdefault(label, species_of(label))
    for label in sorted(assignment):
        print(f"{label}\t{assignment[label]}")


main()
