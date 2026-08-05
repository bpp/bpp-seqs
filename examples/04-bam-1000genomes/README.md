# 04 — 1000 Genomes: aligned reads → BPP (`bam2bpp`)

Real human sequencing data from the 1000 Genomes high-coverage (30×, NYGC,
GRCh38) release. Unlike examples 01–03, the data is **not committed** — it is
fetched by `setup.sh`, which needs no credentials.

| | |
|---|---|
| **Workflow** | `bam2bpp` |
| **Data** | 4 samples × `chr22:20,000,000–20,100,000` (100 kb) |
| **Groups** | `EUR` (HG00096, HG00099), `AFR` (NA19017, NA19019) |
| **Source** | [1000 Genomes / IGSR](https://www.internationalgenome.org/), NYGC 30× GRCh38 release |
| **Download** | ~50 MB reference + ~6 MB of read slices |
| **Requires** | `samtools`, `curl` |

## Run it

```sh
./setup.sh          # fetch reference + BAM slices, write imap.txt and loci.bed

# 1. Inspect: cross-checks BAM contigs vs reference vs BED vs Imap.
../../bpp-seqs --dry-run *.bam chr22.fa loci.bed --imap imap.txt

# 2. Convert.
../../bpp-seqs --out demo *.bam chr22.fa loci.bed --imap imap.txt
```

18 of the 20 loci pass QC; the other two are dropped for having no segregating
sites, which `demo.stats.tsv` records per locus.

## What this example shows that 01–03 cannot

**Loci have to be defined.** Reads are not organised into loci, so this workflow
needs a BED. `setup.sh` builds one with `bpp-seqs windows`, tiling chr22 into
1 kb windows at least 4 kb apart and masking everything outside the slice:

```sh
bpp-seqs windows chr22.fa --window-size 1000 --min-spacing 4000 \
    --include-chrom chr22 --skip-edges 0 --exclude-regions mask.bed --out loci.bed
```

**Diploid calls have to become sequences.** That is what `--phasing` controls,
and it applies here but not to the alignment examples:

```sh
../../bpp-seqs --out iupac   *.bam chr22.fa loci.bed --imap imap.txt                    # 4 seqs
../../bpp-seqs --out split   *.bam chr22.fa loci.bed --imap imap.txt --phasing split    # 8 seqs
../../bpp-seqs --out haploid *.bam chr22.fa loci.bed --imap imap.txt --phasing haploid  # 4 seqs
```

(`haploid` collapses each individual to one major-allele sequence, which removes
the within-individual heterozygosity — one further locus then falls below
`--min-snps` and 17 rather than 18 survive.)

The default `iupac` emits one sequence per individual with heterozygous sites as
ambiguity codes, and bpp-seqs then recommends `phase = 1 1` — one digit per
species, telling BPP it must resolve those heterozygotes itself. `split` instead
emits two sequences per individual, and the recommendation becomes `phase = 0 0`.

**Base-level filtering applies.** `--min-bq`, `--min-mq`, `--min-dp` and
`--het-freq` govern pileup base calling and are meaningful only here (and are
silently ignored by the alignment workflows).
