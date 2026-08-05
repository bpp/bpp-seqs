#!/usr/bin/env python3
"""Build loci/*.fasta and imap.txt from the full Anastrepha Dryad deposit.

The published deposit is far too large to vendor (293 locus alignments of ~700
sequences each). This script cuts it down to a small, self-consistent example
without touching a single residue: it selects specimens belonging to the
*Anastrepha fraterculus* species group, keeps the loci that cover them best,
and copies those sequences through verbatim.

Run it only if you want to regenerate or widen the example -- `fetch-data.sh`
explains how to obtain the source archive first.

Usage:
    python3 make-subset.py DS1_aln_dir AnastrephaNames_LabCodes.csv
"""
import collections
import csv
import glob
import os
import sys

# Species from the fraterculus group with enough specimens to be interesting,
# and few enough to keep the example quick.
SPECIES = ["fraterculus", "obliqua", "suspensa", "distincta", "turpiniae"]
MAX_PER_SPECIES = 6
N_LOCI = 10


def read_labels(csv_path):
    """lab code -> species epithet, for the fraterculus group only."""
    wanted = {}
    with open(csv_path, newline="", errors="replace") as fh:
        for row in csv.reader(fh):
            if len(row) < 2:
                continue
            code, desc = row[0].strip(), row[1].strip()
            parts = desc.split("_")
            if len(parts) < 2 or parts[0] != "frat":
                continue
            if parts[1] in SPECIES:
                wanted[code] = parts[1]
    return wanted


def read_fasta(path):
    seqs, name = collections.OrderedDict(), None
    for line in open(path, errors="replace"):
        line = line.rstrip("\n")
        if line.startswith(">"):
            name = line[1:].strip()
            seqs[name] = []
        elif name:
            seqs[name].append(line.strip())
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    aln_dir, csv_path = sys.argv[1], sys.argv[2]
    here = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.join(here, "loci")
    os.makedirs(out_dir, exist_ok=True)

    code2sp = read_labels(csv_path)

    # Deterministic specimen choice: lowest lab codes first, capped per species.
    per_sp = collections.defaultdict(list)
    for code in sorted(code2sp):
        per_sp[code2sp[code]].append(code)
    chosen = {}
    for sp in SPECIES:
        for code in per_sp[sp][:MAX_PER_SPECIES]:
            chosen[code] = sp
    if not chosen:
        sys.exit("no specimens matched -- is the CSV path right?")

    # Rank loci by how many chosen specimens they contain, then by length.
    scored = []
    for path in sorted(glob.glob(os.path.join(aln_dir, "*.fasta"))):
        seqs = read_fasta(path)
        hits = [n for n in seqs if n in chosen]
        if hits:
            scored.append((len(hits), len(next(iter(seqs.values()))), path, seqs))
    scored.sort(key=lambda r: (-r[0], -r[1]))

    written = set()
    for _, _, path, seqs in scored[:N_LOCI]:
        name = os.path.basename(path)
        with open(os.path.join(out_dir, name), "w") as out:
            for label in sorted(seqs):
                if label in chosen:
                    out.write(f">{label}\n{seqs[label]}\n")
                    written.add(label)

    # Only map specimens that actually occur in the selected loci; not every
    # specimen was sequenced at every locus.
    with open(os.path.join(here, "imap.txt"), "w") as out:
        for code in sorted(written):
            out.write(f"{code}\t{chosen[code]}\n")

    print(f"wrote {min(N_LOCI, len(scored))} loci and "
          f"{len(written)} Imap entries to {out_dir}")


main()
