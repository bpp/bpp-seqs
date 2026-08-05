#!/usr/bin/env bash
# examples/05-gvcf-1000genomes/setup.sh
#
# Build a real multi-sample gVCF from the same 1000 Genomes slice used by
# example 04, so the gvcf2bpp path can be exercised on genuine coverage-band
# data without needing GATK or DeepVariant.
#
# Reuses example 04's downloads rather than fetching the 50 MB reference twice,
# running that example's setup first if it has not been run already.
#
# Requires: bcftools, tabix, samtools, curl.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
bam_dir="$here/../04-bam-1000genomes"
cd "$here"

command -v bcftools >/dev/null || { echo "bcftools is required" >&2; exit 1; }
command -v tabix    >/dev/null || { echo "tabix is required"    >&2; exit 1; }

REGION="chr22:20000000-20100000"
SAMPLES=(HG00096 HG00099 NA19017 NA19019)

# 1) Make sure the reads and reference from example 04 are present.
if [[ ! -s "$bam_dir/chr22.fa" || ! -s "$bam_dir/HG00096.bam" ]]; then
    echo "  example 04 data not found; running its setup first ..."
    "$bam_dir/setup.sh"
fi

# 2) One gVCF per sample. `--gvcf 0` emits reference blocks (coverage bands)
#    rather than variant-only records, which is what gvcf2bpp needs in order to
#    tell "matches the reference" from "not sequenced".
for s in "${SAMPLES[@]}"; do
    if [[ -s "$s.g.vcf.gz" ]]; then
        echo "  skip: $s.g.vcf.gz already present"
        continue
    fi
    echo "  calling $s ..."
    bcftools mpileup -f "$bam_dir/chr22.fa" \
                     --annotate FORMAT/DP,FORMAT/AD \
                     --max-depth 200 -r "$REGION" "$bam_dir/$s.bam" 2>/dev/null \
        | bcftools call --gvcf 0 -m 2>/dev/null \
        | bgzip > "$s.g.vcf.gz"
    tabix -p vcf "$s.g.vcf.gz"
done

# 3) Merge into the single multi-sample gVCF that bpp-seqs reads.
if [[ ! -s cohort.g.vcf.gz ]]; then
    echo "  merging into cohort.g.vcf.gz ..."
    bcftools merge --gvcf "$bam_dir/chr22.fa" \
        HG00096.g.vcf.gz HG00099.g.vcf.gz NA19017.g.vcf.gz NA19019.g.vcf.gz \
        -Oz -o cohort.g.vcf.gz 2>/dev/null
    tabix -fp vcf cohort.g.vcf.gz
fi

# 4) Reuse example 04's loci and Imap.
cp "$bam_dir/loci.bed" loci.bed
cp "$bam_dir/imap.txt" imap.txt

echo
echo "  setup done:"
ls -lh cohort.g.vcf.gz loci.bed imap.txt | sed 's/^/    /'
