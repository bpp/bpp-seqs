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

Two subcommands sit alongside this flow: `bpp-seqs windows` builds a BED of
candidate loci when you don't have one, and `bpp-seqs extract` subsets an
already-converted BPP file. See [Subcommands](#subcommands).

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

Unknown or abbreviated long options are rejected (exact match only), so a
mistyped flag is an error, not a silent no-op. `bpp-seqs --help` prints the same
list.

### General

| Option | Meaning |
|--------|---------|
| `--out PREFIX` | Output file prefix; required for conversion |
| `--dry-run` | Inspect only; never convert |
| `--json` | Emit JSON on stdout instead of human-readable text |
| `--json-indent N` | JSON indentation width [2] |
| `--quiet` | Suppress stderr progress messages |
| `-h`, `--help` | Show help and exit |
| `--version` | Show version and exit |

### Input

| Option | Meaning |
|--------|---------|
| `--imap FILE` | Imap file (sample → population mapping) |
| `--reference FILE` | Force FILE to be treated as the reference FASTA, skipping the content-based reference-vs-contigs classifier |

### Phasing

Applies to BAM/CRAM and gVCF input; alignment inputs (FASTA MSA / PHYLIP /
NEXUS) are already resolved sequences and ignore it.

| Option | Meaning |
|--------|---------|
| `--phasing MODE` | `iupac` (default), `split`, `haploid`, or `vcf` — see [Phasing](#phasing---phasing) above |
| `--phased-vcf FILE` | Phased VCF supplying true phase for `--phasing vcf` |

### Base/read filtering

Applied while piling up reads, so these affect **BAM/CRAM input only**. They are
accepted but ignored for gVCF and alignment inputs.

| Option | Default | Meaning |
|--------|---------|---------|
| `--min-bq INT` | 20 | Minimum base quality for a read base to count |
| `--min-mq INT` | 20 | Minimum mapping quality for a read to count |
| `--min-dp INT` | 5 | Minimum depth to call a base (below this → `N`) |
| `--het-freq FLOAT` | 0.20 | Minor-allele frequency threshold for a heterozygote call |

### Locus filtering

Applied by every conversion workflow.

| Option | Default | Meaning |
|--------|---------|---------|
| `--min-length INT` | 100 | Minimum locus length in bp |
| `--max-missing FLOAT` | 0.5 | Maximum fraction of `N` per locus |
| `--min-snps INT` | 1 | Minimum segregating sites per locus |
| `--keep-invariant` | off | Keep loci with no variation |

### Exit status

`0` on success — including an inspection that reports missing items, since
`incomplete` is a usable answer, not an error. `1` on a conversion or system
error.

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

For worked examples on real published data — one per workflow, with citations
and licences — see [`examples/`](examples/). Examples 01–03 (FASTA, PHYLIP,
NEXUS) ship with the repository and run straight after a clone; 04–05 (BAM,
gVCF) fetch a 1000 Genomes slice with `./setup.sh`.

## Subcommands

Two verbs sit alongside the main inspect/convert flow. Each takes its own
options; run `bpp-seqs VERB --help` for the authoritative list.

### `bpp-seqs windows` — build a BED of candidate loci

Tiles a genome with fixed-size windows and writes them as a BED, for when you
have BAM/CRAM or gVCF input but no locus definitions yet. `INPUT` may be a
reference FASTA, a BAM/CRAM, or a VCF/gVCF — anything carrying chromosome names
and lengths.

```
bpp-seqs windows INPUT --window-size W [filters] --out OUT.bed
```

| Option | Meaning |
|--------|---------|
| `--window-size W` | Locus size in bp (required) |
| `--out FILE` | Output BED path; `-` writes to stdout (required) |
| `--step S` | Stride between window starts [window-size, i.e. non-overlapping] |
| `--include-chrom NAME[,…]` | Only consider these chromosomes |
| `--exclude-chrom NAME[,…]` | Skip these chromosomes |
| `--autosomes-only` | Skip sex chromosomes, mitochondria, and unplaced / random / alt / decoy contigs by name heuristic |
| `--skip-edges N` | Drop the first and last N bp of every chromosome |
| `--exclude-regions FILE` | Subtract the intervals in this BED before windowing |
| `--min-spacing G` | Require ≥ G bp between consecutive kept windows on a chromosome |
| `--n-loci N` | Randomly sample N windows from the survivors (omit to keep all) |
| `--seed S` | Seed for that sampling [0] |
| `--json` | Emit a JSON summary on stdout in addition to the BED |
| `-h`, `--help` | Show help and exit |

Filters apply in the order listed. Output is a 4-column BED (chrom, 0-based
start, exclusive end, `locusN`), ready to pass straight back in as loci.

```sh
# 500 windows of ~500 bp, ≥10 kb apart, autosomes only, edges skipped
bpp-seqs windows ref.fa --window-size 500 --min-spacing 10000 \
                 --n-loci 500 --autosomes-only --skip-edges 100000 \
                 --seed 42 --out loci.bed
```

### `bpp-seqs extract` — subset an existing BPP file

Reads a BPP sequence file and writes a new one holding a subset of its loci. A
sibling `<INPUT>.loci.tsv` (if present) supplies source coordinates and is
filtered to the same subset; `<INPUT>.imap` is copied through unchanged.

```
bpp-seqs extract INPUT.txt --out PREFIX [selection ...]
```

| Option | Meaning |
|--------|---------|
| `--out PREFIX` | Output file prefix (required) |
| `--loci NAME[,…]` | Keep these locus names (union with `--loci-file`) |
| `--loci-file FILE` | One locus name per line |
| `--range A-B[,C-D]` | Keep loci by 1-based index (1 = first locus) |
| `--first N` | Keep the first N loci |
| `--last N` | Keep the last N loci |
| `--chrom NAME` | Keep loci whose `source_chrom` is NAME (needs a `.loci.tsv`) |
| `--min-sites N` | Keep loci with `n_sites` ≥ N |
| `--max-sites N` | Keep loci with `n_sites` ≤ N |
| `--invert` | Keep the loci *not* matching the selection |
| `--imap FILE` | Use this Imap instead of `<INPUT>.imap` |
| `--loci-tsv FILE` | Use this provenance instead of `<INPUT>.loci.tsv` |
| `--json` | Emit a JSON summary on stdout |
| `--quiet` | Suppress stderr progress |
| `-h`, `--help` | Show help and exit |

At least one selection is required, and the different kinds compose with AND.
`--range`, `--first`, and `--last` union with each other, and index positions in
the **input** file — not positions within some other filter's result. Writes
`PREFIX.txt`, plus `PREFIX.loci.tsv` and `PREFIX.imap` when those sidecars exist.

```sh
# Every chr1 locus with at least 500 sites
bpp-seqs extract run1.txt --chrom chr1 --min-sites 500 --out subset

# A quick 50-locus test set
bpp-seqs extract run1.txt --first 50 --out smoke
```
