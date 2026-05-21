# bpp-seqs

Inspect biological sequence data files and convert them to the alignment
format required by BPP (Bayesian Phylogenetics and Phylogeography).

## Build

Requires htslib. On macOS: `brew install htslib`. On Linux:
`conda install -c bioconda htslib` or build from source.

```
make           # release build
make debug     # -g + ASan/UBSan
make test      # run integration tests
make install   # install to /usr/local/bin (override with PREFIX=...)
```

## Usage

```
bpp-seqs [options] file1 file2 ...
```

Files are provided in any order. File type is detected from content. The tool
inspects all inputs, cross-validates them, and either converts (when all
required inputs are present and `--out PREFIX` is given) or reports what is
missing.

Key options:

```
  --out PREFIX          Output file prefix (required for conversion)
  --dry-run             Inspect only, do not convert
  --json                Emit JSON instead of human-readable text
  --quiet               Suppress stderr progress messages
  --imap FILE           Sample → population map (also accepted positionally)

  --phasing MODE        iupac (default), split, haploid, vcf
  --phased-vcf FILE     Phased VCF for --phasing vcf
  --min-bq INT          Minimum base quality [20]
  --min-mq INT          Minimum mapping quality [20]
  --min-dp INT          Minimum depth to call a base [5]
  --het-freq FLOAT      Minor allele freq for het call [0.20]
  --min-length INT      Minimum locus length [100]
  --max-missing FLOAT   Max fraction of N per locus [0.5]
  --min-snps INT        Minimum segregating sites per locus [1]
  --keep-invariant      Keep loci with no variation
```

See `bpp-seqs-spec.md` for the full specification.
