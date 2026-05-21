#!/usr/bin/env bash
# tests/realdata/run.sh — end-to-end bpp-seqs pipeline on the 1KG mini-slice.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"
cd "$here"

BIN="$root/bpp-seqs"
[[ -x "$BIN" ]] || { (cd "$root" && make); }

# Restrict windowing to the region we actually sliced: chr22:20000000-20100000.
# windows itself tiles all chr22, so we use --exclude-regions to mask out
# everything outside the slice.
cat > mask.bed <<EOF
chr22	0	20000000
chr22	20100000	60000000
EOF

echo "=== step 1: bpp-seqs windows ==="
"$BIN" windows GRCh38_chr22.fa \
    --window-size 1000 \
    --min-spacing 5000 \
    --skip-edges 0 \
    --exclude-regions mask.bed \
    --include-chrom chr22 \
    --out loci.bed
echo
echo "  windows generated:"
wc -l < loci.bed
head -3 loci.bed
echo
echo "=== step 2: bpp-seqs convert (bam2bpp, --phasing vcf) ==="
"$BIN" --out result --keep-invariant --max-missing 0.9 \
    --phasing vcf --phased-vcf phased.vcf.gz \
    HG00096.bam HG01595.bam NA19017.bam HG01112.bam \
    GRCh38_chr22.fa loci.bed samples.imap 2>&1 | tail -25

echo
echo "=== output summary ==="
ls -la result.*
echo
echo "=== first locus in result.txt ==="
head -7 result.txt
echo
echo "=== result.loci.tsv ==="
head -4 result.loci.tsv
echo "  ..."
tail -3 result.loci.tsv
echo
echo "=== result.imap ==="
cat result.imap
