# 1000 Genomes mini-test dataset

A small end-to-end test of bpp-seqs against real human sequencing data.

## Samples (4, distinct super-populations)

| Sample | Pop | Super-pop | Description |
|---|---|---|---|
| HG00096 | GBR | EUR | British in England/Scotland |
| HG01595 | KHV | EAS | Kinh in Ho Chi Minh City, Vietnam |
| NA19017 | LWK | AFR | Luhya in Webuye, Kenya |
| HG01112 | CLM | AMR | Colombian in Medellin, Colombia |

## Region

GRCh38 chr22:20000000-20100000 (100 kb on the q-arm — gene-rich, well-covered).

## Files

Produced by `setup.sh`:
- `HG00096.bam`, `HG01595.bam`, `NA19017.bam`, `HG01112.bam` (+ `.bai`)
- `chr22.fa` (+ `.fai`)
- `samples.imap`

Produced by the pipeline:
- `loci.bed` — from `bpp-seqs windows`
- `result.{txt,imap,stats.tsv,loci.tsv}` — from `bpp-seqs`

## Running

```bash
./tests/realdata/setup.sh    # downloads ~5 MB BAMs + ~50 MB chr22 ref
./tests/realdata/run.sh      # runs windows + bam2bpp pipeline
```

## Validating the output with BPP itself

If a BPP binary is available, the generated `result.txt` / `result.imap`
can be fed straight into a short MCMC run:

```bash
cd tests/realdata
~/repos/bpp/src/bpp --cfile bpp.ctl
```

`bpp.ctl` defines a 4-population, A00 (parameter estimation only) analysis
with a topology `(AFR,((EUR,AMR),EAS))` and a 200-burnin / 200-sample MCMC.
On the fixture this runs in well under one second and produces sensible
θ estimates (AFR highest, AMR lowest — matching the known out-of-Africa
diversity pattern and the founder-effect bottleneck in the Americas).
