#!/usr/bin/env bash
# tests/realdata/setup.sh — fetch a 4-sample / 100-kb slice from the
# 1000 Genomes high-coverage (30x, NYGC, GRCh38) release.
#
# CRAMs are reference-compressed, so the chr22 reference must match the
# one used by NYGC: GRCh38 with the no-alt analysis set (GRCh38_full_analysis_set_plus_decoy_hla.fa)
# using "chr"-prefixed contig names.

set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"

REGION="chr22:20000000-20100000"

# Two individuals per super-population: gives 4 chromosomes per population
# for θ estimation when the data are phased.
# (sample, super-pop, ENA CRAM URL)
SAMPLES=(
  "HG00096:EUR:https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240114/HG00096.final.cram"
  "HG00099:EUR:https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240116/HG00099.final.cram"
  "HG01595:EAS:https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3242062/HG01595.final.cram"
  "HG01596:EAS:https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3242063/HG01596.final.cram"
  "NA19017:AFR:https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239683/NA19017.final.cram"
  "NA19019:AFR:https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239684/NA19019.final.cram"
  "HG01112:AMR:https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241828/HG01112.final.cram"
  "HG01113:AMR:https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241829/HG01113.final.cram"
)

# 1) GRCh38 chr22 reference for CRAM decoding.
#    NYGC pipeline used GRCh38_full_analysis_set_plus_decoy_hla.fa, but
#    the chr22 sequence is identical to the plain GRCh38 chr22 from NCBI.
if [[ ! -s GRCh38_chr22.fa ]]; then
    echo "  fetching GRCh38 chr22 reference (~50 MB) ..."
    curl -fsSL "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz" \
        | gunzip -c > GRCh38_chr22.fa
fi
samtools faidx GRCh38_chr22.fa

# 2) Slice each CRAM to the region. Each slice is small (~few hundred KB)
#    even though the source CRAM is 16+ GB.
for entry in "${SAMPLES[@]}"; do
    IFS=":" read -r sample pop url <<<"$entry"
    out="$sample.bam"
    if [[ -s "$out" && -s "$out.bai" ]]; then
        echo "  skip: $out already present"
        continue
    fi
    echo "  fetching $sample ($pop) slice from $REGION ..."
    samtools view -b -h -T GRCh38_chr22.fa "$url" "$REGION" > "$out.tmp"
    samtools sort -o "$out" "$out.tmp"
    samtools index "$out"
    rm -f "$out.tmp"
done

# 3) Phased VCF slice (NYGC phased panel, GRCh38, chr22).
#    Lets us exercise --phasing vcf instead of pileup-derived IUPAC.
if [[ ! -s phased.vcf.gz ]]; then
    echo "  fetching phased VCF slice ..."
    PHASED_URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    bcftools view -O z -r "$REGION" \
        -s HG00096,HG00099,HG01595,HG01596,NA19017,NA19019,HG01112,HG01113 \
        "$PHASED_URL" > phased.vcf.gz
    bcftools index -ft phased.vcf.gz
    bcftools index -fc phased.vcf.gz
fi

# 4) Per-sample bcftools-style gVCFs, merged into one multi-sample gVCF.
#    Lets us exercise the gvcf2bpp workflow on real coverage-band data
#    without needing GATK / DeepVariant.
if [[ ! -s merged.g.vcf.gz ]]; then
    echo "  generating per-sample gVCFs ..."
    for entry in "${SAMPLES[@]}"; do
        IFS=":" read -r sample _ _ <<<"$entry"
        if [[ ! -s "${sample}.g.vcf.gz" ]]; then
            bcftools mpileup -f GRCh38_chr22.fa \
                             --annotate FORMAT/DP,FORMAT/AD \
                             --max-depth 200 -r "$REGION" "${sample}.bam" 2>/dev/null \
                | bcftools call --gvcf 0 -m 2>/dev/null \
                | bgzip > "${sample}.g.vcf.gz"
            tabix -p vcf "${sample}.g.vcf.gz"
        fi
    done
    echo "  merging into one multi-sample gVCF ..."
    bcftools merge --gvcf GRCh38_chr22.fa \
        HG00096.g.vcf.gz HG00099.g.vcf.gz HG01595.g.vcf.gz HG01596.g.vcf.gz \
        NA19017.g.vcf.gz NA19019.g.vcf.gz HG01112.g.vcf.gz HG01113.g.vcf.gz \
        -Oz -o merged.g.vcf.gz 2>/dev/null
    tabix -fp vcf merged.g.vcf.gz
fi

# 5) Imap (sample -> super-population for population structure).
cat > samples.imap <<'EOF'
HG00096	EUR
HG00099	EUR
HG01595	EAS
HG01596	EAS
NA19017	AFR
NA19019	AFR
HG01112	AMR
HG01113	AMR
EOF

echo
echo "  setup done — files in $here:"
ls -lh
