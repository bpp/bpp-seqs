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

# (sample, super-pop, ENA CRAM URL)
SAMPLES=(
  "HG00096:EUR:https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240114/HG00096.final.cram"
  "HG01595:EAS:https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3242062/HG01595.final.cram"
  "NA19017:AFR:https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239683/NA19017.final.cram"
  "HG01112:AMR:https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3241828/HG01112.final.cram"
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

# 3) Imap (sample -> super-population for population structure).
cat > samples.imap <<'EOF'
HG00096	EUR
HG01595	EAS
NA19017	AFR
HG01112	AMR
EOF

echo
echo "  setup done — files in $here:"
ls -lh
