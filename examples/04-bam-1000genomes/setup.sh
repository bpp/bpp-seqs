#!/usr/bin/env bash
# examples/04-bam-1000genomes/setup.sh
#
# Fetch a small slice of real human sequencing data from the 1000 Genomes
# high-coverage (30x, NYGC, GRCh38) release. Nothing here needs credentials.
#
# Downloads ~50 MB (chr22 reference) plus ~1 MB of read slices; the source
# CRAMs are 16+ GB each but only the requested region is transferred.
#
# Requires: samtools, curl.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
cd "$here"

command -v samtools >/dev/null || { echo "samtools is required" >&2; exit 1; }

REGION="chr22:20000000-20100000"     # 100 kb on the gene-rich q-arm

# Two individuals from each of two super-populations: enough for BPP to have
# four chromosomes per population once the diploid calls are split.
SAMPLES=(
  "HG00096:EUR:https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240114/HG00096.final.cram"
  "HG00099:EUR:https://ftp.sra.ebi.ac.uk/vol1/run/ERR324/ERR3240116/HG00099.final.cram"
  "NA19017:AFR:https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239683/NA19017.final.cram"
  "NA19019:AFR:https://ftp.sra.ebi.ac.uk/vol1/run/ERR323/ERR3239684/NA19019.final.cram"
)

# 1) chr22 reference. The CRAMs are reference-compressed, so this is needed to
#    decode them as well as to call bases later.
if [[ ! -s chr22.fa ]]; then
    echo "  fetching GRCh38 chr22 reference (~50 MB) ..."
    curl -fsSL "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz" \
        | gunzip -c > chr22.fa
fi
[[ -s chr22.fa.fai ]] || samtools faidx chr22.fa

# 2) Slice each CRAM to the region and store it as an indexed BAM.
for entry in "${SAMPLES[@]}"; do
    IFS=":" read -r sample pop url <<<"$entry"
    if [[ -s "$sample.bam" && -s "$sample.bam.bai" ]]; then
        echo "  skip: $sample.bam already present"
        continue
    fi
    echo "  fetching $sample ($pop) over $REGION ..."
    samtools view -b -h -T chr22.fa "$url" "$REGION" > "$sample.tmp.bam"
    samtools sort -o "$sample.bam" "$sample.tmp.bam"
    samtools index "$sample.bam"
    rm -f "$sample.tmp.bam"
done

# 3) Imap: sample -> super-population.
cat > imap.txt <<'EOF'
HG00096	EUR
HG00099	EUR
NA19017	AFR
NA19019	AFR
EOF

# 4) Loci. Read-based input needs a BED, and `bpp-seqs windows` builds one.
#    The slice only covers 100 kb, so mask out the rest of the chromosome.
cat > mask.bed <<'EOF'
chr22	0	20000000
chr22	20100000	60000000
EOF

BIN="$here/../../bpp-seqs"
[[ -x "$BIN" ]] || { echo "  building bpp-seqs ..."; (cd "$here/../.." && make >/dev/null); }

"$BIN" windows chr22.fa \
    --window-size 1000 --min-spacing 4000 \
    --include-chrom chr22 --skip-edges 0 \
    --exclude-regions mask.bed --out loci.bed

echo
echo "  setup done: $(wc -l < loci.bed) loci in loci.bed, 4 samples"
ls -lh chr22.fa *.bam imap.txt loci.bed 2>/dev/null | sed 's/^/    /'
