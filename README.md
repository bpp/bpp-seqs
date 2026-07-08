# bpp-seqs

Inspect biological sequence-data files and convert them to the per-locus
alignment format required by **BPP** (Bayesian Phylogenetics and
Phylogeography). It detects file types from content, cross-validates the inputs,
and either converts them or tells you exactly what is still missing.

## Install

**Homebrew (macOS):**

```
brew install bpp/tap/bpp-seqs
```

**Linux (HPC):** download the static `linux-x86_64` binary from the
[latest release](https://github.com/bpp/bpp-seqs/releases) — htslib is linked
in, so there is no runtime dependency to load:

```
tar xzf bpp-seqs-*-linux-x86_64.tar.gz
./bpp-seqs-*/bpp-seqs --version
```

## Build from source

Requires [htslib](https://github.com/samtools/htslib). On macOS:
`brew install htslib`. On Linux: `conda install -c bioconda htslib` or build
from source.

```
make           # release build (dynamic htslib)
make debug     # -g + ASan/UBSan
make test      # run the integration suite
make install   # install to /usr/local/bin (override with PREFIX=...)
```

For a self-contained binary, link htslib statically (this is what the release
CI does):

```
( cd htslib && ./configure --disable-libcurl --disable-plugins && make )
make HTSLIB_A="$PWD/htslib/libhts.a" HTSLIB_INC="$PWD/htslib"
```

## How it works: inspect, then convert

Every run **inspects** first — classifies each input, samples records, and
cross-validates the set (do the BAM `@SQ` contigs match the reference and the
BED? are the Imap samples present in the data?). From that it picks a
**workflow** and decides whether it can proceed:

- If `--out PREFIX` is given **and** every required input is present, it
  **converts** and writes the BPP files.
- Otherwise (or with `--dry-run`, or without `--out`) it stops at inspection and
  reports `missing[]` — the specific pieces still needed.

```
bpp-seqs [options] file1 file2 ...
```

Files may be given in any order; type is detected from content, not extension.

## Input formats and the workflow each selects

| Inputs you provide | Workflow | Still needs |
|--------------------|----------|-------------|
| BAM/CRAM + reference FASTA + BED + Imap | `bam2bpp` | — (ready) |
| BAM/CRAM, no BED | `bam2bpp_needs_bed` | a BED (loci) |
| BAM/CRAM + BED, no reference | `bam2bpp_needs_ref` | a reference FASTA (+`.fai`) |
| gVCF + BED + Imap | `gvcf2bpp` | — (ready) |
| gVCF, no BED | `gvcf2bpp_needs_bed` | a BED (loci) and an Imap |
| Aligned FASTA (one MSA per locus) + Imap | `fasta2bpp` | — (ready) |
| PHYLIP + Imap | `phylip2bpp` | — (ready) |
| NEXUS + Imap | `nexus2bpp` | — (ready) |
| FASTQ | `needs_alignment_first` | align to a reference first (BAM+ref+BED) |
| Plain (non-gVCF) VCF | `vcf_not_recommended` | use the BAMs or a gVCF instead |
| Assembled contigs (FASTA) | `assembly_not_supported` | define loci and align first |

A convertible workflow that is missing a piece names it in `missing[]` with a
description of the exact format expected, rather than failing silently.

### Loci (BED)

Read-based inputs (BAM/CRAM, gVCF) need a **BED** defining the locus intervals:
`chrom <TAB> start <TAB> end [<TAB> name]`, 0-based start, exclusive end. Each
interval becomes one BPP locus, and its chrom names must match the data's
contigs. Alignment inputs (FASTA MSA / PHYLIP / NEXUS) are already per-locus, so
they need no BED.

### Imap

The **Imap** (sample → population mapping) is required for every conversion. It
may be passed with `--imap FILE` or positionally. `bpp-seqs` cross-checks it
against the samples actually present and flags orphans in either direction.

### Phasing (`--phasing`)

How diploid genotype calls become sequences (the `--phasing` flag — this is
*not* the BPP control-file `phase=` line, which is reported separately as
`recommended_phase`):

- `iupac` (default) — heterozygous sites as IUPAC ambiguity codes; one sequence
  per individual.
- `split` — two unphased haplotypes per diploid individual.
- `haploid` — a single major-allele sequence per individual.
- `vcf` — apply true phase from a phased VCF (`--phased-vcf FILE`).

## Options

```
General:   --out PREFIX  --dry-run  --json  --json-indent N  --quiet  --version  -h/--help
Input:     --imap FILE   --reference FILE
Phasing:   --phasing MODE   --phased-vcf FILE
Filtering: --min-bq 20  --min-mq 20  --min-dp 5  --het-freq 0.20
           --min-length 100  --max-missing 0.5  --min-snps 1  --keep-invariant
```

`--reference FILE` forces a FASTA to be treated as the reference (skipping the
content-based reference-vs-contigs classifier). Run `bpp-seqs --help` for the
authoritative list and defaults. Unknown or abbreviated long options are
rejected (exact match only), so a mistyped flag is an error, not a silent no-op.

## Output

A conversion writes four files under `PREFIX`:

| File | Contents |
|------|----------|
| `PREFIX.txt` | the BPP sequence file (per-locus alignments) |
| `PREFIX.imap` | the Imap actually used |
| `PREFIX.stats.tsv` | per-locus statistics |
| `PREFIX.loci.tsv` | the locus table |

## JSON output (`--json`)

`--json` emits a single machine-readable object (the same shape for inspection
and conversion), so `bpp-seqs` can drive a pipeline. Top-level fields:

| Field | Meaning |
|-------|---------|
| `status` | `complete` (converted / ready) or `incomplete` (something missing) |
| `workflow` | the selected workflow (table above), or `unknown` |
| `ready_to_run` | true iff a convertible workflow has no missing items |
| `files_provided[]` | per input: `type`, `n_samples`, `sample_names`, `has_phase_info`, `phased_fraction`, `has_coverage_bands`, `warnings`, … |
| `cross_validation` | contig/reference/Imap consistency checks and any `issues[]` |
| `missing[]` | each gap as `{item, description}` — e.g. `imap_file`, `bed_file`, `reference_fasta` |
| `recommended_phase` | the suggested BPP control-file `phase=` line (`"0 0 0 0"`/`"1 1 1 1"`, one per species) |
| `output_files` | `null`, or `{sequences, imap, stats, loci}` after a conversion |
| `warnings[]` | non-fatal advisories |

## Examples

```sh
# Inspect a gVCF: what is still needed before it can be converted?
bpp-seqs --json --dry-run cohort.g.vcf.gz
#   -> workflow: gvcf2bpp_needs_bed, missing: [imap_file, bed_file]

# Convert BAMs to BPP input
bpp-seqs --out run1 aln/*.bam reference.fa loci.bed samples.imap

# Convert a gVCF with explicit loci and haploid phasing
bpp-seqs --out run2 --phasing haploid cohort.g.vcf.gz loci.bed samples.imap

# Aligned FASTA loci straight to BPP
bpp-seqs --out run3 locus_*.fa samples.imap
```
