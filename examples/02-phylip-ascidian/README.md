# 02 — botryllid ascidians: per-locus PHYLIP → BPP (`phylip2bpp`)

Real anchored-hybrid-enrichment (AHE) alignments from the botryllid ascidian
phylogenomics of Nydam et al. (2021). Each locus is a separate relaxed-PHYLIP
file, which is the layout most phylogenomic pipelines emit.

| | |
|---|---|
| **Workflow** | `phylip2bpp` |
| **Data** | 20 nuclear AHE loci, 533–5284 bp each |
| **Sequences** | 55 per locus |
| **Groups** | 17 — species where the study identified one, plus family-level and undetermined sets |
| **Source** | Dryad [doi:10.5061/dryad.3r2280gf7](https://doi.org/10.5061/dryad.3r2280gf7), file `Nydam2021_T500_FinalAlignments.zip` |
| **Paper** | Nydam M., Lemmon A., Cherry J., Kortyna M., Clancy D., Hernandez C., Cohen C.S. (2021) Phylogenomic and morphological relationships among the botryllid ascidians. *Scientific Reports* 11. [doi:10.1038/s41598-021-87255-2](https://doi.org/10.1038/s41598-021-87255-2) |
| **License** | CC0-1.0 (public domain dedication) |

The published archive holds 179 loci. `loci/` vendors the 20 with the most
complete taxon sampling (all 55 specimens present, longest alignments first) to
keep the repository small. The files themselves are unmodified — no sequence,
label, or locus has been edited. `fetch-data.sh` documents how to obtain the
full set.

## Run it

```sh
# 1. Inspect.
../../bpp-seqs --dry-run loci/*.phylip --imap imap.txt

# 2. Convert to BPP format.
../../bpp-seqs --out ascidian loci/*.phylip --imap imap.txt
```

All 20 loci pass QC, giving `ascidian.txt`, `ascidian.imap`,
`ascidian.stats.tsv` and `ascidian.loci.tsv`.

## Where the Imap comes from

`imap.txt` is generated from the taxon labels by `make-imap.py`:

```sh
python3 make-imap.py > imap.txt
```

Labels are `I<specimen>_<identification>_seq<n>`, and the identification is the
authors' own:

```
I25523_Botryllidae_Botrylloides_praelongus_seq1  → Botryllidae_Botrylloides_praelongus
I25550_Botryllidae___seq1                        → Botryllidae   (family only)
I28212____seq1                                   → unidentified
```

This dataset is a good illustration of a real-world wrinkle: a majority of the
specimens were **not** identified to species by the study, so most sequences
land in family-level or `unidentified` groups. The script keeps them rather
than guessing at assignments. For an actual BPP analysis you would replace
`imap.txt` with the delimitation you want to test — the point here is that the
grouping is a decision you make, not something the file format dictates.

## Phasing

Already-resolved sequences, so `--phasing` does not apply (it affects BAM/CRAM
and gVCF input only).
