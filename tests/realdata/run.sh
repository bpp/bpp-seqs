#!/usr/bin/env bash
# tests/realdata/run.sh — exercise multiple locus-partitioning strategies
# on the 1000 Genomes mini-slice (chr22:20M-20.1M, 8 samples).
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/../.." && pwd)"
cd "$here"

BIN="$root/bpp-seqs"
[[ -x "$BIN" ]] || { (cd "$root" && make); }

# Mask out everything outside our 100 kb slice, since `windows` tiles the
# whole chromosome and we only have data in chr22:20000000-20100000.
cat > mask.bed <<EOF
chr22	0	20000000
chr22	20100000	60000000
EOF

# ─── Four partitioning strategies ─────────────────────────────────────────
declare -A STRATEGIES=(
    [adjacent]="--window-size 1000"                   # all contiguous 1 kb blocks
    [small]="--window-size 500  --min-spacing 1000"   # ~33 short loci
    [medium]="--window-size 1000 --min-spacing 5000"  # ~17 medium loci
    [large]="--window-size 2000 --min-spacing 10000"  #  ~9 long loci
)

for tag in adjacent small medium large; do
    flags=${STRATEGIES[$tag]}
    out="result_${tag}"
    bed="loci_${tag}.bed"

    echo
    echo "=========================================================="
    echo "  Strategy: ${tag}   ($flags)"
    echo "=========================================================="

    rm -f "$bed" "$out".*

    # 1) Generate the BED.
    "$BIN" windows GRCh38_chr22.fa $flags --include-chrom chr22 \
        --skip-edges 0 --exclude-regions mask.bed --out "$bed"

    echo "    BED: $(wc -l < "$bed") windows"
    head -2 "$bed" | sed 's/^/      /'

    # 2) Convert with --phasing vcf so we get 2 haplotypes per individual.
    "$BIN" --quiet --keep-invariant --max-missing 0.9 \
        --out "$out" \
        --phasing vcf --phased-vcf phased.vcf.gz \
        HG00096.bam HG00099.bam HG01595.bam HG01596.bam \
        NA19017.bam NA19019.bam HG01112.bam HG01113.bam \
        GRCh38_chr22.fa "$bed" samples.imap

    # 3) Summarise the result.
    n_loci=$(grep -cE '^[0-9]+ [0-9]+$' "$out.txt")
    n_seqs=$(awk '/^[0-9]+ [0-9]+$/ {print $1; exit}' "$out.txt")
    total_snps=$(awk -F'\t' 'NR>1 {sum+=$3} END {print sum+0}' "$out.stats.tsv")
    mean_dp=$(awk -F'\t' 'NR>1 {s+=$5; n++} END {if(n) printf "%.1f", s/n}' "$out.stats.tsv")
    rec_phase=$(awk -F'\t' 'NR>1 && /^# phase/ {print}' "$out.stats.tsv" || true)

    printf "    BPP result:  %d loci x %d sequences, total SNPs = %d, mean depth = %sx\n" \
        "$n_loci" "$n_seqs" "$total_snps" "$mean_dp"
done

# ─── gvcf2bpp on the merged real gVCF (different ingest pathway) ──────────
echo
echo "=========================================================="
echo "  gvcf2bpp on real merged gVCF (8 samples)"
echo "=========================================================="
rm -f result_gvcf.*
"$BIN" --quiet --keep-invariant --max-missing 0.5 \
    --out result_gvcf \
    merged.g.vcf.gz loci_medium.bed samples.imap

n_loci=$(grep -cE '^[0-9]+ [0-9]+$' result_gvcf.txt)
n_seqs=$(awk '/^[0-9]+ [0-9]+$/ {print $1; exit}' result_gvcf.txt)
total_snps=$(awk -F'\t' 'NR>1 {sum+=$3} END {print sum+0}' result_gvcf.stats.tsv)
printf "    BPP result:  %d loci x %d sequences, total SNPs = %d\n" "$n_loci" "$n_seqs" "$total_snps"

# ─── extract: take the first 5 loci from the medium partition ─────────────
echo
echo "=========================================================="
echo "  extract: first 5 loci of the medium partition"
echo "=========================================================="
"$BIN" extract result_medium.txt --first 5 --out result_medium_first5
n_loci_extract=$(grep -cE '^[0-9]+ [0-9]+$' result_medium_first5.txt)
echo "    extracted $n_loci_extract loci into result_medium_first5.txt"
echo
echo "All strategies done.  Compare summary files:"
ls -la result_*.txt result_*.loci.tsv | column -t
