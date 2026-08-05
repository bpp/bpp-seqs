# 01 — *Anastrepha* fruit flies: per-locus FASTA → BPP (`fasta2bpp`)

Real phylogenomic alignments from Norrbom et al. (2025), subset to the
*Anastrepha fraterculus* species group — a genuine cryptic-species complex, and
so the kind of problem BPP exists to address.

| | |
|---|---|
| **Workflow** | `fasta2bpp` |
| **Data** | 10 nuclear loci, one aligned FASTA per locus |
| **Sequences** | 24 specimens across 5 species |
| **Groups** | `fraterculus`, `obliqua`, `suspensa`, `distincta`, `turpiniae` |
| **Source** | Dryad [doi:10.5061/dryad.7d7wm3864](https://doi.org/10.5061/dryad.7d7wm3864), file `Individual_Locus_Alignments.zip` (`DS1_aln`) |
| **Paper** | Norrbom A.L., Rodriguez E.J., Steck G.J., Cassel B.K., Ruiz-Arce R., Muller A., Gangadin A., Savaris M. (2025) A new phylogeny of *Anastrepha* (Diptera: Tephritidae) based on nuclear loci obtained by phylogenomic methods. *Systematic Entomology* 51. [doi:10.1111/syen.70003](https://doi.org/10.1111/syen.70003) |
| **License** | CC0-1.0 (public domain dedication) |

## Run it

```sh
# 1. Inspect.
../../bpp-seqs --dry-run loci/*.fasta --imap imap.txt

# 2. Convert to BPP format.
../../bpp-seqs --out anastrepha loci/*.fasta --imap imap.txt
```

All 10 loci pass QC, giving `anastrepha.txt`, `anastrepha.imap`,
`anastrepha.stats.tsv` and `anastrepha.loci.tsv`.

## How this subset was made

The published deposit is 293 locus alignments of roughly 700 sequences each —
far too large to vendor. `make-subset.py` cuts it down **without altering a
single residue**: it picks specimens in the *fraterculus* group, keeps the 10
loci covering them best, and copies those sequences through verbatim. See
`fetch-data.sh` for how to get the source archive, then:

```sh
python3 make-subset.py path/to/DS1_aln path/to/AnastrephaNames_LabCodes_*.csv
```

Edit `SPECIES`, `MAX_PER_SPECIES` and `N_LOCI` at the top of the script to widen
the example.

## Where the Imap comes from

Alignment labels are bare lab codes (`a01b4523`), so the species assignment is
read from the authors' own metadata table shipped in the same deposit, whose
descriptors have the form
`speciesGroup_speciesName_localityCode_specimenID`:

```
a01b4523,pseudop_curitibana_BRA_PR_USNMENT01354523
a056630,frat_distincta_BRA_PA_USNMENT01556630
```

so `a056630` → `distincta`. Nothing here is guessed at: the mapping comes from
the published CSV.

Note that `imap.txt` lists only the specimens that actually appear in the
selected loci — not every specimen was sequenced at every locus, which is normal
for capture data and is exactly why `bpp-seqs` cross-checks the two directions.

## Phasing

These are already-resolved sequences, so `--phasing` does not apply — it affects
BAM/CRAM and gVCF input only.
