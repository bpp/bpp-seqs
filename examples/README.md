# bpp-seqs examples

Two small, self-contained datasets you can run end to end. Every run
**inspects** first (`--dry-run`) and then **converts** to BPP sequence format.
Add `--json` to any command for machine-readable output (the human progress
lines go to stderr, so `--json` stdout is always pure JSON).

Paths below assume you run from inside each example directory, with the
`bpp-seqs` binary two levels up (`../../bpp-seqs`). Adjust if it's on your PATH.

## 01-fasta — per-locus alignments → BPP

Two aligned loci (`locusA.fa`, `locusB.fa`) for four individuals, plus an Imap
assigning them to two populations. Already-aligned sequences need no phasing.

```
cd 01-fasta

# 1. Inspect: what does conversion still need?
../../bpp-seqs --dry-run locusA.fa locusB.fa --imap imap.txt

# 2. Convert to BPP format.
../../bpp-seqs locusA.fa locusB.fa --imap imap.txt --out demo1
```

Produces `demo1.txt` (BPP sequence file), `demo1.imap`, and per-locus stats.
Both loci pass QC (2 loci, 4 sequences). Use `phase = 0` in the control file —
the sequences are already resolved lineages.

## 02-bam — aligned reads → BPP

Four indexed BAMs (`ind1..4.bam`) aligned to `ref.fa`, a BED of five loci, and
the same Imap. bpp-seqs pileups each BED interval into a locus alignment.

```
cd 02-bam

# 1. Inspect (cross-checks BAM contigs vs. reference vs. BED vs. Imap).
../../bpp-seqs --dry-run ind1.bam ind2.bam ind3.bam ind4.bam ref.fa loci.bed --imap imap.txt

# 2. Convert.
../../bpp-seqs ind1.bam ind2.bam ind3.bam ind4.bam ref.fa loci.bed --imap imap.txt --out demo2
```

Produces `demo2.txt` + `demo2.imap` + stats. Of the five loci, four pass and one
(`locus4`) is dropped as `insufficient_snps` (no variable sites). These are
unphased diploid genotype calls, so the inspection recommends `phase = 1` in the
control file (heterozygous sites are encoded as IUPAC codes and BPP resolves the
phase analytically).

## Constructing loci when you don't have a BED

If you have a reference but no loci defined, tile it into candidate loci:

```
bpp-seqs windows ref.fa --window-size 500 --min-spacing 10000 \
                 --autosomes-only --skip-edges 100000 --out loci.bed
```

then feed the resulting BED to the BAM workflow above.
