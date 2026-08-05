# bpp-seqs examples

Five worked examples, one per conversion workflow, all built on **real
published data** with a citable source and an open licence. Every run inspects
first and then converts, so you can see what bpp-seqs checks before it commits
to anything. Add `--json` to any command for machine-readable output — the
human progress lines go to stderr, so `--json` stdout is always pure JSON.

| Example | Workflow | Data | Ships with the repo? |
|---------|----------|------|----------------------|
| [01-fasta-anastrepha](01-fasta-anastrepha/) | `fasta2bpp` | 10 nuclear loci, 24 fruit-fly specimens, 5 species | yes |
| [02-phylip-ascidian](02-phylip-ascidian/) | `phylip2bpp` | 20 AHE loci, 55 ascidian specimens | yes |
| [03-nexus-milksnake](03-nexus-milksnake/) | `nexus2bpp` | 12 loci, 248 sequences, 23 milksnake lineages | yes |
| [04-bam-1000genomes](04-bam-1000genomes/) | `bam2bpp` | 4 human samples, 100 kb of chr22 | no — `./setup.sh` |
| [05-gvcf-1000genomes](05-gvcf-1000genomes/) | `gvcf2bpp` | the same slice as a cohort gVCF | no — `./setup.sh` |

Examples 01–03 are alignment inputs: small enough to commit, so they run
straight after a clone with no network. Examples 04–05 are read-based inputs,
whose data is fetched by `setup.sh` (no credentials needed) and never committed.

Paths in each README assume you are inside that example's directory with the
binary two levels up (`../../bpp-seqs`). Adjust if it is on your `PATH`.

## The quickest look

```sh
cd 03-nexus-milksnake
../../bpp-seqs --dry-run loci/*.nex --imap imap.txt         # inspect
../../bpp-seqs --out milksnake loci/*.nex --imap imap.txt   # convert
```

## What each example is for

**01–03 — alignment inputs.** Already one alignment per locus, so they need no
BED and no `--phasing`; the interesting part is how sample names become an Imap.
Each shows a different real-world naming convention: the lineage encoded in the
label (03), an identification of varying completeness (02), or bare lab codes
that must be joined to the study's own metadata table (01).

**04–05 — read-based inputs.** Loci must be defined with a BED (example 04
builds one with `bpp-seqs windows`), and diploid genotype calls must be turned
into sequences, which is what `--phasing` controls. Example 05 covers the same
region through a gVCF instead, and explains why a plain VCF will not do.

## Provenance and licensing

The committed datasets are redistributed under **CC0-1.0** public-domain
dedications, and each example's README carries the DOI, the paper, and the exact
source file. Where a dataset was too large to vendor whole, a script in the
directory reproduces the subset from the published archive, and no sequence data
is ever modified — only selected. If you use these data in your own work, cite
the original studies, not this repository.
