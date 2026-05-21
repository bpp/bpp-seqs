#!/usr/bin/env bash
# tests/run_tests.sh — integration tests for bpp-seqs
#
# Each test is named, runs the binary, and asserts something about the
# output (exit code, file contents, or JSON shape).  Failures print a
# diagnostic and the script exits non-zero.

set -uo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
root="$(cd "$here/.." && pwd)"
bin="$root/bpp-seqs"
data="$here/data"
tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT

pass=0
fail=0

run() {
    local name="$1"; shift
    if "$@"; then
        echo "  PASS  $name"
        pass=$((pass + 1))
    else
        echo "  FAIL  $name"
        fail=$((fail + 1))
    fi
}

# ── Setup: synthesize test reference if absent ─────────────────────────────
if [[ ! -f "$data/test_ref.fa" ]]; then
    python3 -c "
seq = 'ACGT' * 625
with open('$data/test_ref.fa','w') as f:
    f.write('>chr1\n')
    for i in range(0, len(seq), 60): f.write(seq[i:i+60] + '\n')
"
    samtools faidx "$data/test_ref.fa"
fi

# ── Scenario 1: dry-run, missing imap, expect incomplete (exit 0) ──────────
t1() {
    local out
    out=$("$bin" --dry-run \
              "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
              "$data/test_ref.fa" "$data/loci.bed" 2>/dev/null) || return 1
    echo "$out" | grep -q "Workflow: bam2bpp" && \
    echo "$out" | grep -q "imap_file"
}
run "1. dry-run without imap reports missing imap" t1

# ── Scenario 2: with imap, status ready, --dry-run still no output files ──
t2() {
    local out
    out=$("$bin" --dry-run \
              "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
              "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" 2>/dev/null) || return 1
    echo "$out" | grep -q "Cross-validation: OK" && \
    ! echo "$out" | grep -q "Missing:"
}
run "2. dry-run with imap reports ready, no missing" t2

# ── Scenario 3: full conversion, output identical to expected ──────────────
t3() {
    "$bin" --quiet --out "$tmp/m" \
        "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
        "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    diff -q "$tmp/m.txt"  "$data/test_out.txt"  >/dev/null || return 1
    diff -q "$tmp/m.imap" "$data/test_out.imap" >/dev/null
}
run "3. bam2bpp conversion matches test_out.txt and test_out.imap" t3

# ── Scenario 4: --json --dry-run without imap → status incomplete ──────────
t4() {
    local out
    out=$("$bin" --json --dry-run \
              "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
              "$data/test_ref.fa" "$data/loci.bed" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
d = json.loads(sys.argv[1])
assert d["status"] == "incomplete", d["status"]
assert d["workflow"] == "bam2bpp"
assert d["ready_to_run"] is False
assert len(d["missing"]) == 1 and d["missing"][0]["item"] == "imap_file"
PY
}
run "4. --json --dry-run without imap returns status=incomplete" t4

# ── Scenario 5: --json full conversion → status complete + files exist ────
t5() {
    local out
    out=$("$bin" --json --quiet --out "$tmp/jm" \
              "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
              "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" 2>/dev/null) || return 1
    python3 - "$out" <<PY
import json, os, sys
d = json.loads(sys.argv[1])
assert d["status"] == "complete", d["status"]
of = d["output_files"]
for key in ("sequences","imap","stats"):
    assert os.path.isfile(of[key]), of[key]
PY
}
run "5. --json full conversion returns complete and writes files" t5

# ── Scenario 6: fasta2bpp converter ────────────────────────────────────────
t6() {
    "$bin" --quiet --out "$tmp/msa" \
        "$data/msa1.fa" "$data/msa2.fa" "$data/imap.txt" >/dev/null 2>&1 || return 1
    [[ -f "$tmp/msa.txt" && -f "$tmp/msa.imap" ]] && \
        grep -qE '^4 [0-9]+$' "$tmp/msa.txt" && \
        grep -qE '^\^ind1' "$tmp/msa.txt"
}
run "6. fasta2bpp converts MSA fixtures" t6

# ── Scenario 7: phylip2bpp converter ───────────────────────────────────────
t7() {
    "$bin" --quiet --out "$tmp/phy" \
        "$data/aln.phy" "$data/imap.txt" >/dev/null 2>&1 || return 1
    [[ -f "$tmp/phy.txt" && -f "$tmp/phy.imap" ]] && \
        grep -qE '^4 [0-9]+$' "$tmp/phy.txt" && \
        grep -qE '^\^ind1' "$tmp/phy.txt"
}
run "7. phylip2bpp converts PHYLIP fixture" t7

# ── Scenario 8: nexus2bpp converter ────────────────────────────────────────
t8() {
    "$bin" --quiet --min-length 50 --keep-invariant --out "$tmp/nex" \
        "$data/aln.nex" "$data/imap.txt" >/dev/null 2>&1 || return 1
    [[ -f "$tmp/nex.txt" && -f "$tmp/nex.imap" ]] && \
        # two loci → two "n_seqs n_sites" header lines
        [[ $(grep -cE '^4 [0-9]+$' "$tmp/nex.txt") -eq 2 ]] && \
        grep -qE '^\^ind1' "$tmp/nex.txt"
}
run "8. nexus2bpp converts NEXUS with charsets" t8

# ── Scenario 9: multilocus PHYLIP (multiple n_seqs n_sites blocks per file) ──
t_multi_phy() {
    [[ -f "$data/multi.phy" ]] || return 0  # skip if absent
    "$bin" --quiet --out "$tmp/mp" \
        "$data/multi.phy" "$data/imap.txt" >/dev/null 2>&1 || return 1
    # File should contain three "4 <length>" headers, one per locus
    [[ $(grep -cE '^4 [0-9]+$' "$tmp/mp.txt") -eq 3 ]]
}
run "9. phylip2bpp parses multilocus PHYLIP (3 loci, varying lengths)" t_multi_phy

# ── Scenarios 10–13: sanity checks ─────────────────────────────────────────

# Duplicate individual id within a locus
t_sanity_dup() {
    python3 -c "
print('>ind1'); print('ACGT'*30 + 'A')
print('>ind1'); print('ACGT'*30 + 'C')
print('>ind2'); print('ACGT'*30 + 'G')
print('>ind3'); print('ACGT'*30 + 'T')
" > "$tmp/dup.fa"
    "$bin" --quiet --keep-invariant --out "$tmp/od" "$tmp/dup.fa" "$data/imap.txt" 2>"$tmp/err" >/dev/null
    grep -q "duplicate individual id" "$tmp/err"
}
run "10. sanity: duplicate individual id is detected and rejected" t_sanity_dup

# Invalid characters sanitized to N
t_sanity_chars() {
    python3 -c "
print('>ind1'); print('ACGTACGT' * 14 + 'ACGZZZ?Q')
print('>ind2'); print('ACGTACGT' * 14 + 'ACGTACGT')
print('>ind3'); print('ACGTACGT' * 14 + 'ACGTACGT')
print('>ind4'); print('ACGTACGT' * 14 + 'ACGTACGT')
" > "$tmp/chars.fa"
    "$bin" --quiet --keep-invariant --out "$tmp/oc" "$tmp/chars.fa" "$data/imap.txt" 2>"$tmp/err" >/dev/null
    grep -q "sanitized to 'N'" "$tmp/err"
}
run "11. sanity: invalid characters sanitized with warning" t_sanity_chars

# Imap orphan: sample missing from Imap
t_sanity_orphan_seq() {
    python3 -c "
print('>ind1');     print('ACGTACGT' * 15)
print('>ind2');     print('ACGTACGT' * 15)
print('>ind3');     print('ACGTACGT' * 15)
print('>missing1'); print('ACGTACGT' * 15)
" > "$tmp/o1.fa"
    "$bin" --quiet --keep-invariant --out "$tmp/oo1" "$tmp/o1.fa" "$data/imap.txt" 2>"$tmp/err" >/dev/null
    grep -q "'missing1' present in sequence file but not in Imap" "$tmp/err"
}
run "12. sanity: sequence id missing from Imap is reported" t_sanity_orphan_seq

# Imap orphan: Imap row not used by any locus
t_sanity_orphan_imap() {
    cat > "$tmp/imap_orph.txt" <<EOF
ind1	popA
ind2	popA
ind3	popB
ind4	popB
extra	popC
EOF
    "$bin" --quiet --keep-invariant --out "$tmp/oo2" "$data/msa1.fa" "$tmp/imap_orph.txt" 2>"$tmp/err" >/dev/null
    grep -q "'extra' has no matching sequence" "$tmp/err"
}
run "13. sanity: unused Imap entry is reported" t_sanity_orphan_imap

# ── Scenarios 14–15: declared-dimension mismatch warnings ────────────────

# PHYLIP declares fewer sites than the data actually provides
t_phy_excess_sites() {
    python3 -c "
print('   4   100')
for s in ['ind1','ind2','ind3','ind4']:
    print(f'{s:<10}' + 'ACGT'*30)
print()
" > "$tmp/over.phy"
    "$bin" --quiet --keep-invariant --min-length 50 --out "$tmp/op" \
        "$tmp/over.phy" "$data/imap.txt" 2>"$tmp/err" >/dev/null
    grep -q "characters beyond declared n_sites" "$tmp/err"
}
run "14. PHYLIP excess characters beyond declared n_sites trigger warning" t_phy_excess_sites

# NEXUS declares fewer nchar than the data provides
t_nex_excess_chars() {
    cat > "$tmp/over.nex" <<EOF
#NEXUS
BEGIN DATA;
    DIMENSIONS NTAX=4 NCHAR=60;
    FORMAT DATATYPE=DNA MISSING=N GAP=-;
    MATRIX
        ind1   $(python3 -c "print('ACGT'*30)")
        ind2   $(python3 -c "print('ACGT'*30)")
        ind3   $(python3 -c "print('ACGT'*30)")
        ind4   $(python3 -c "print('ACGT'*30)")
    ;
END;
EOF
    "$bin" --quiet --keep-invariant --min-length 30 --out "$tmp/on" \
        "$tmp/over.nex" "$data/imap.txt" 2>"$tmp/err" >/dev/null
    grep -q "characters beyond declared nchar" "$tmp/err"
}
run "15. NEXUS excess characters beyond declared nchar trigger warning" t_nex_excess_chars

# ── Scenario 16: per-locus provenance file (.loci.tsv) ─────────────────────
t_loci_tsv() {
    "$bin" --quiet --out "$tmp/lt" \
        "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
        "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    [[ -f "$tmp/lt.loci.tsv" ]] || return 1
    # Header line
    head -1 "$tmp/lt.loci.tsv" | grep -q '^locus_name	source_kind	source_file	source_chrom	source_start	source_end	source_stride	length	n_seqs$' || return 1
    # One data row per passing locus, BED provenance
    [[ $(awk 'NR>1 && $2=="BED" && $4=="chr1"' "$tmp/lt.loci.tsv" | wc -l) -eq 4 ]]
}
run "16. bam2bpp emits .loci.tsv with BED provenance for each passing locus" t_loci_tsv

# ── Spec-compliance and feature scenarios ──────────────────────────────────

# Partial-BAM workflow: BAM + ref + Imap, no BED → missing bed_file entry
t_needs_bed() {
    local out
    out=$("$bin" --json --dry-run \
              "$data/ind1.bam" "$data/test_ref.fa" "$data/imap.txt" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
d = json.loads(sys.argv[1])
assert d["workflow"] == "bam2bpp_needs_bed", d["workflow"]
items = [m["item"] for m in d["missing"]]
assert "bed_file" in items, items
PY
}
run "17. bam2bpp_needs_bed adds bed_file to missing[]" t_needs_bed

# Partial-BAM workflow: BAM + BED + Imap, no ref → missing reference_fasta
t_needs_ref() {
    local out
    out=$("$bin" --json --dry-run \
              "$data/ind1.bam" "$data/loci.bed" "$data/imap.txt" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
d = json.loads(sys.argv[1])
assert d["workflow"] == "bam2bpp_needs_ref", d["workflow"]
items = [m["item"] for m in d["missing"]]
assert "reference_fasta" in items, items
PY
}
run "18. bam2bpp_needs_ref adds reference_fasta to missing[]" t_needs_ref

# --quiet really silences bam2bpp writer
t_quiet() {
    local err
    err=$("$bin" --quiet --out "$tmp/q" \
            "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
            "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" 2>&1 >/dev/null)
    ! echo "$err" | grep -q "Wrote BPP"
}
run "19. --quiet suppresses bam2bpp Wrote-* stderr lines" t_quiet

# gVCF↔BED chrom mismatch is detected
t_gvcf_chrom_mismatch() {
    echo -e "chrZZZ\t0\t400\tlocus1" > "$tmp/wrong.bed"
    "$bin" --dry-run "$data/tiny.g.vcf.gz" "$tmp/wrong.bed" "$data/imap.txt" 2>&1 | \
        grep -qE 'CHROMOSOME_MISMATCH.*##contig'
}
run "20. gVCF↔BED chrom mismatch surfaces CHROMOSOME_MISMATCH" t_gvcf_chrom_mismatch

# gvcf2bpp coverage-band fill: positions inside a NON_REF block are not all N
t_gvcf_block_fill() {
    "$bin" --quiet --keep-invariant --max-missing 1.0 --min-length 50 \
        --out "$tmp/gv" \
        "$data/tiny.g.vcf.gz" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    # In locus1, the first NON_REF block is chr1:100-199 (REF=A). After fill,
    # ind1 should have a run of A's, not all N's.
    awk '/^4 / {lc++} lc==1 && /^\^ind1/ {print}' "$tmp/gv.txt" | grep -q 'AAAAAAAAAA'
}
run "21. gvcf2bpp fills non-variant coverage blocks across full span" t_gvcf_block_fill

# gvcf2bpp --phasing split produces 2*n samples with _1/_2 ids
t_gvcf_split() {
    "$bin" --quiet --phasing split --keep-invariant --max-missing 1.0 \
        --min-length 50 --out "$tmp/sp" \
        "$data/tiny.g.vcf.gz" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    # First locus header should report 8 sequences (4 samples * 2 haplotypes)
    awk '/^[0-9]+ [0-9]+$/{print $1; exit}' "$tmp/sp.txt" | grep -q '^8$' && \
    grep -q '^\^ind1_1' "$tmp/sp.txt" && grep -q '^\^ind1_2' "$tmp/sp.txt"
}
run "22. gvcf2bpp --phasing split emits 2 sequences per sample" t_gvcf_split

# gvcf2bpp --phasing haploid keeps single sequence per sample
t_gvcf_haploid() {
    "$bin" --quiet --phasing haploid --keep-invariant --max-missing 1.0 \
        --min-length 50 --out "$tmp/hp" \
        "$data/tiny.g.vcf.gz" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    awk '/^[0-9]+ [0-9]+$/{print $1; exit}' "$tmp/hp.txt" | grep -q '^4$'
}
run "23. gvcf2bpp --phasing haploid keeps one sequence per sample" t_gvcf_haploid

# CRAM round-trips identically to BAM
t_cram() {
    [[ -f "$data/ind1.cram" ]] || return 0
    "$bin" --quiet --out "$tmp/cr" \
        "$data/ind1.cram" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
        "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    diff -q "$tmp/cr.txt" "$data/test_out.txt" >/dev/null
}
run "24. CRAM input round-trips byte-identical to BAM" t_cram

# Workflow advisory for VCF-not-recommended
t_vcf_advisory() {
    # Re-use the gVCF fixture but trick the workflow by giving VCF (not gVCF).
    # Easiest path: hand a VCF without coverage bands.
    python3 -c "
import gzip
content = '''##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##contig=<ID=chr1,length=2500>
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tind1\tind2\tind3\tind4
chr1\t200\t.\tA\tT\t100\t.\t.\tGT\t0/0\t0/1\t1/1\t0/0
'''
open('$tmp/plain.vcf','w').write(content)
"
    bgzip -f "$tmp/plain.vcf" && tabix -p vcf "$tmp/plain.vcf.gz" >/dev/null
    local out
    out=$("$bin" --json --dry-run "$tmp/plain.vcf.gz" "$data/loci.bed" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
d = json.loads(sys.argv[1])
assert d["workflow"] == "vcf_not_recommended", d["workflow"]
assert "advisory" in d and "gVCF" in d["advisory"], d.get("advisory")
PY
}
run "25. VCF (non-gVCF) workflow reports vcf_not_recommended + advisory" t_vcf_advisory

# ── Scenario 26: gvcf2bpp converter ────────────────────────────────────────
t9() {
    [[ -f "$data/tiny.g.vcf.gz" && -f "$data/tiny.g.vcf.gz.tbi" ]] || return 0  # skip if absent
    "$bin" --quiet --min-length 50 --keep-invariant --max-missing 1.0 \
        --out "$tmp/gv" \
        "$data/tiny.g.vcf.gz" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    [[ -f "$tmp/gv.txt" ]]
}
run "26. gvcf2bpp converts tiny gVCF fixture" t9

# ── Summary ───────────────────────────────────────────────────────────────
echo
echo "Tests: $pass passed, $fail failed"
[[ $fail -eq 0 ]]
