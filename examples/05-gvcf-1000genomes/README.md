# 05 — 1000 Genomes: gVCF → BPP (`gvcf2bpp`)

The same real 1000 Genomes slice as example 04, but ingested as a multi-sample
**gVCF** rather than as reads. Data is fetched by `setup.sh`, not committed.

| | |
|---|---|
| **Workflow** | `gvcf2bpp` |
| **Data** | 4 samples × `chr22:20,000,000–20,100,000`, called into one cohort gVCF |
| **Groups** | `EUR` (HG00096, HG00099), `AFR` (NA19017, NA19019) |
| **Source** | Called with `bcftools` from example 04's read slices |
| **Requires** | `bcftools`, `tabix`, `samtools`, `curl` |

`setup.sh` reuses example 04's downloads instead of fetching the 50 MB reference
a second time, and runs that example's setup first if needed.

## Run it

```sh
./setup.sh          # build per-sample gVCFs, merge into cohort.g.vcf.gz

# 1. Inspect.
../../bpp-seqs --dry-run cohort.g.vcf.gz loci.bed --imap imap.txt

# 2. Convert.
../../bpp-seqs --out demo cohort.g.vcf.gz loci.bed --imap imap.txt
```

18 of 20 loci pass QC — the same result example 04 reaches from the reads, which
is a useful check that the two ingest paths agree.

## Why a gVCF and not a plain VCF

A plain VCF lists variant sites only, so a position absent from the file is
ambiguous: it might match the reference, or it might not have been sequenced at
all. Reconstructing a locus alignment needs that distinction, which is why
bpp-seqs reports `vcf_not_recommended` for an ordinary VCF and asks for the BAMs
or a gVCF instead.

A gVCF resolves it by emitting **reference blocks** — spans of confidently
called reference bases carried as `END`/`MIN_DP` records. That is what
`bcftools call --gvcf 0` produces in `setup.sh`, and what lets bpp-seqs write a
real base rather than `N`. You can see the bands in the file itself:

```sh
bcftools view cohort.g.vcf.gz | grep -m3 END=
```

## Caveat: base-level filters do not apply here

`--min-bq`, `--min-mq`, `--min-dp` and `--het-freq` act during pileup and so
affect BAM/CRAM input only. They are accepted but ignored on this path — the
genotype calls in the gVCF have already been made by the caller, so filtering
belongs in the `bcftools` step (or in a `bcftools filter` pass) instead.
