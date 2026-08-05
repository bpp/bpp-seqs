#!/usr/bin/env python3
"""Rebuild imap.txt from the taxon labels in loci/*.phylip.

Labels are  I<specimen>_<identification>_seq<n>, e.g.

    I25523_Botryllidae_Botrylloides_praelongus_seq1
    I25550_Botryllidae___seq1            family only, species undetermined
    I28212____seq1                       no identification recorded

The identification between the specimen number and the `_seq` suffix is the
authors' own, so it is used verbatim as the group label. Specimens the study
left undetermined are grouped as `unidentified` rather than being guessed at;
BPP will treat that as one more population, which is what an unresolved set of
specimens actually is.

Usage:  python3 make-imap.py > imap.txt
"""
import glob
import os
import re
import sys

LABEL = re.compile(r"^(I\d+)_(.*?)_seq\d+$")


def group_of(label):
    m = LABEL.match(label)
    if not m:
        return "unidentified"
    ident = re.sub(r"_+", "_", m.group(2)).strip("_")
    return ident or "unidentified"


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    files = sorted(glob.glob(os.path.join(here, "loci", "*.phylip")))
    if not files:
        sys.exit("no loci/*.phylip found")
    assignment = {}
    for path in files:
        with open(path, errors="replace") as fh:
            fh.readline()                      # "<ntax> <nchar>"
            for line in fh:
                parts = line.split()
                if parts:
                    assignment.setdefault(parts[0], group_of(parts[0]))
    for label in sorted(assignment):
        print(f"{label}\t{assignment[label]}")


main()
