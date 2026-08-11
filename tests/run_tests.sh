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
skip=0

# A test may return 77 to mean "not applicable here" -- the oracle tests use
# this when samtools is absent or is not the pinned version. Reported as SKIP
# rather than PASS, so a run that verified nothing cannot look like a clean
# one. Set BPP_SEQS_STRICT=1 (CI does) to turn any skip into a failure.
run() {
    local name="$1"; shift
    "$@"
    local rc=$?
    if [ "$rc" -eq 77 ]; then
        if [ "${BPP_SEQS_STRICT:-0}" = "1" ]; then
            echo "  FAIL  $name (skipped, but BPP_SEQS_STRICT=1)"
            fail=$((fail + 1))
        else
            echo "  SKIP  $name"
            skip=$((skip + 1))
        fi
    elif [ "$rc" -eq 0 ]; then
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

# ── Scenarios 26–34: `extract` subcommand ────────────────────────────────

# Build a master BPP file from the bam2bpp fixture (4 loci passing QC).
prep_master() {
    [[ -f "$tmp/mas.txt" ]] && return 0
    "$bin" --quiet --out "$tmp/mas" \
        "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
        "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1
}

t_extract_by_name() {
    prep_master || return 1
    "$bin" extract "$tmp/mas.txt" --loci locus1,locus3 --out "$tmp/x1" >/dev/null 2>&1 || return 1
    [[ $(awk 'NR>1 {print $1}' "$tmp/x1.loci.tsv" | tr '\n' ' ') == "locus1 locus3 " ]]
}
run "26. extract --loci selects named loci" t_extract_by_name

t_extract_by_loci_file() {
    prep_master || return 1
    printf '%s\n%s\n' locus2 locus5 > "$tmp/names.txt"
    "$bin" extract "$tmp/mas.txt" --loci-file "$tmp/names.txt" --out "$tmp/x2" >/dev/null 2>&1 || return 1
    [[ $(awk 'NR>1 {print $1}' "$tmp/x2.loci.tsv" | tr '\n' ' ') == "locus2 locus5 " ]]
}
run "27. extract --loci-file selects names from a file" t_extract_by_loci_file

t_extract_by_range() {
    prep_master || return 1
    "$bin" extract "$tmp/mas.txt" --range 2-3 --out "$tmp/x3" >/dev/null 2>&1 || return 1
    [[ $(awk 'NR>1 {print $1}' "$tmp/x3.loci.tsv" | tr '\n' ' ') == "locus2 locus3 " ]]
}
run "28. extract --range selects 1-based index range" t_extract_by_range

t_extract_first_last() {
    prep_master || return 1
    "$bin" extract "$tmp/mas.txt" --first 2 --out "$tmp/x4" >/dev/null 2>&1 || return 1
    [[ $(awk 'NR>1' "$tmp/x4.loci.tsv" | wc -l) -eq 2 ]] || return 1
    "$bin" extract "$tmp/mas.txt" --last 1 --out "$tmp/x5" >/dev/null 2>&1 || return 1
    [[ $(awk 'NR>1' "$tmp/x5.loci.tsv" | wc -l) -eq 1 ]]
}
run "29. extract --first/--last select by position" t_extract_first_last

t_extract_chrom() {
    prep_master || return 1
    "$bin" extract "$tmp/mas.txt" --chrom chr1 --out "$tmp/x6" >/dev/null 2>&1 || return 1
    [[ $(awk 'NR>1' "$tmp/x6.loci.tsv" | wc -l) -eq 4 ]]
}
run "30. extract --chrom matches against .loci.tsv source_chrom" t_extract_chrom

t_extract_min_max_sites() {
    prep_master || return 1
    # All master loci are 400 bp; min-sites 500 should match zero.
    if "$bin" extract "$tmp/mas.txt" --min-sites 500 --out "$tmp/x7" >/dev/null 2>&1; then
        return 1   # expected failure
    fi
    "$bin" extract "$tmp/mas.txt" --max-sites 400 --out "$tmp/x8" >/dev/null 2>&1 || return 1
    [[ $(awk 'NR>1' "$tmp/x8.loci.tsv" | wc -l) -eq 4 ]]
}
run "31. extract --min-sites / --max-sites filter by n_sites" t_extract_min_max_sites

t_extract_invert() {
    prep_master || return 1
    "$bin" extract "$tmp/mas.txt" --loci locus2 --invert --out "$tmp/x9" >/dev/null 2>&1 || return 1
    awk 'NR>1 {print $1}' "$tmp/x9.loci.tsv" | grep -q '^locus2$' && return 1
    [[ $(awk 'NR>1' "$tmp/x9.loci.tsv" | wc -l) -eq 3 ]]
}
run "32. extract --invert keeps loci NOT matching the selection" t_extract_invert

t_extract_imap_passthrough() {
    prep_master || return 1
    "$bin" extract "$tmp/mas.txt" --first 1 --out "$tmp/x10" >/dev/null 2>&1 || return 1
    diff -q "$tmp/x10.imap" "$tmp/mas.imap" >/dev/null
}
run "33. extract copies the Imap unchanged" t_extract_imap_passthrough

t_extract_json() {
    prep_master || return 1
    local out
    out=$("$bin" extract "$tmp/mas.txt" --first 2 --out "$tmp/x11" --json 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys, os
d = json.loads(sys.argv[1])
assert d["command"] == "extract"
assert d["n_loci_kept"] == 2
assert d["had_loci_tsv"] is True
for k in ("sequences","loci","imap"):
    assert os.path.isfile(d["output_files"][k]), d["output_files"][k]
PY
}
run "34. extract --json reports the kept/dropped counts and output files" t_extract_json

# ── Scenarios 35–40: `windows` subcommand ────────────────────────────────

prep_multifa() {
    [[ -f "$tmp/multi.fa" ]] && return 0
    python3 -c "
seq = 'ACGT' * 100
with open('$tmp/multi.fa','w') as f:
    for n in ('chr1','chr2','chrX'):
        f.write(f'>{n}\n')
        for i in range(0, len(seq), 60): f.write(seq[i:i+60] + '\n')
"
    samtools faidx "$tmp/multi.fa"
}

t_win_basic() {
    "$bin" windows "$data/test_ref.fa" --window-size 400 --out "$tmp/w1.bed" 2>/dev/null || return 1
    [[ $(wc -l < "$tmp/w1.bed") -eq 6 ]]
}
run "35. windows tiles a single chrom into non-overlapping windows" t_win_basic

t_win_spacing() {
    "$bin" windows "$data/test_ref.fa" --window-size 400 --min-spacing 200 \
                  --out "$tmp/w2.bed" 2>/dev/null || return 1
    [[ $(wc -l < "$tmp/w2.bed") -eq 3 ]]
}
run "36. windows --min-spacing thins overlapping/adjacent windows" t_win_spacing

t_win_sample() {
    "$bin" windows "$data/test_ref.fa" --window-size 400 --n-loci 3 --seed 7 \
                  --out "$tmp/w3.bed" 2>/dev/null || return 1
    [[ $(wc -l < "$tmp/w3.bed") -eq 3 ]]
}
run "37. windows --n-loci with --seed downsamples reproducibly" t_win_sample

t_win_autosomes() {
    prep_multifa || return 1
    "$bin" windows "$tmp/multi.fa" --window-size 100 --autosomes-only \
                  --out "$tmp/wa.bed" 2>/dev/null || return 1
    [[ -z $(awk '$1 == "chrX"' "$tmp/wa.bed") ]]
}
run "38. windows --autosomes-only drops chrX/Y/mt by name heuristic" t_win_autosomes

t_win_exclude_regions() {
    prep_multifa || return 1
    printf 'chr1\t100\t300\n' > "$tmp/mask.bed"
    "$bin" windows "$tmp/multi.fa" --window-size 100 \
                  --exclude-regions "$tmp/mask.bed" \
                  --out "$tmp/wr.bed" 2>/dev/null || return 1
    # chr1 should have only 2 windows (0-100 and 300-400); 100-200 and 200-300 are masked
    [[ $(awk '$1 == "chr1"' "$tmp/wr.bed" | wc -l) -eq 2 ]]
}
run "39. windows --exclude-regions subtracts a BED mask" t_win_exclude_regions

t_win_pipeline() {
    # End-to-end: windows -> BED -> conversion
    "$bin" windows "$data/test_ref.fa" --window-size 400 --skip-edges 100 \
                  --out "$tmp/pipe.bed" 2>/dev/null || return 1
    "$bin" --quiet --keep-invariant --out "$tmp/pipe" \
        "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
        "$data/test_ref.fa" "$tmp/pipe.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    [[ -s "$tmp/pipe.txt" && -s "$tmp/pipe.loci.tsv" ]] && \
        awk 'NR>1 && $3 ~ /pipe.bed$/' "$tmp/pipe.loci.tsv" | head -1 | grep -q BED
}
run "40. windows + bam2bpp pipeline produces valid BPP output" t_win_pipeline

# ── Scenarios 41–42: phase recommendation ────────────────────────────────

t_phase_recommend_iupac() {
    local out
    out=$("$bin" --json --quiet --out "$tmp/pr" \
            "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
            "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
d = json.loads(sys.argv[1])
# Default phasing is IUPAC → recommend "1 1" (one digit per Imap species)
assert d.get("recommended_phase") == "1 1", d.get("recommended_phase")
PY
}
run "41. recommended_phase = '1 1' under default IUPAC phasing" t_phase_recommend_iupac

t_phase_recommend_split() {
    local out
    out=$("$bin" --json --quiet --phasing split --out "$tmp/ps" \
            "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
            "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
d = json.loads(sys.argv[1])
# Split phasing emits 2 sequences per sample → BPP should NOT re-phase
assert d.get("recommended_phase") == "0 0", d.get("recommended_phase")
PY
}
run "42. recommended_phase = '0 0' under --phasing split" t_phase_recommend_split

# ── Scenarios 43–49: input-combination coverage ─────────────────────────

# Synthesize a tiny phased VCF over our chr1 test reference + ind1-4 samples.
make_phased_vcf() {
    [[ -s "$tmp/phased.vcf.gz" ]] && return 0
    cat > "$tmp/phased.vcf" <<'VCF'
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##contig=<ID=chr1,length=2500>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind1	ind2	ind3	ind4
chr1	100	.	A	T	.	PASS	.	GT	0|1	0|0	1|1	0|0
chr1	200	.	C	G	.	PASS	.	GT	0|0	1|0	0|1	1|1
chr1	1500	.	A	C	.	PASS	.	GT	0|0	0|1	1|0	0|0
VCF
    bgzip -f "$tmp/phased.vcf"
    tabix -fp vcf "$tmp/phased.vcf.gz"
    bcftools index -fc "$tmp/phased.vcf.gz"
    # also a BCF copy for the BCF/CSI path
    bcftools view -O b -o "$tmp/phased.bcf" "$tmp/phased.vcf.gz" 2>/dev/null
    bcftools index -fc "$tmp/phased.bcf"
}

# BAM + --phasing vcf with VCF.gz (tabix path, the recently-fixed branch).
t_bam_phasing_vcf_gz() {
    make_phased_vcf || return 1
    "$bin" --quiet --keep-invariant --out "$tmp/pv" \
        --phasing vcf --phased-vcf "$tmp/phased.vcf.gz" \
        "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
        "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    # 4 individuals * 2 haplotypes per locus = 8 sequences
    awk '/^[0-9]+ [0-9]+/{print $1; exit}' "$tmp/pv.txt" | grep -q '^8$' && \
        grep -q '^\^ind1_1' "$tmp/pv.txt" && grep -q '^\^ind4_2' "$tmp/pv.txt"
}
run "43. BAM + --phasing vcf with VCF.gz tabix index" t_bam_phasing_vcf_gz

# BAM + --phasing vcf with BCF (CSI path).
t_bam_phasing_vcf_bcf() {
    make_phased_vcf || return 1
    "$bin" --quiet --keep-invariant --out "$tmp/pb" \
        --phasing vcf --phased-vcf "$tmp/phased.bcf" \
        "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
        "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    diff -q "$tmp/pv.txt" "$tmp/pb.txt" >/dev/null
}
run "44. BAM + --phasing vcf with BCF/CSI matches VCF.gz output" t_bam_phasing_vcf_bcf

# BAM + --phasing split (two haplotypes per individual, arbitrary phase).
t_bam_phasing_split() {
    "$bin" --quiet --keep-invariant --phasing split --out "$tmp/ps" \
        "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
        "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    awk '/^[0-9]+ [0-9]+/{print $1; exit}' "$tmp/ps.txt" | grep -q '^8$' && \
        grep -q '^\^ind1_1' "$tmp/ps.txt" && grep -q '^\^ind1_2' "$tmp/ps.txt" && \
        [[ $(grep -c '^ind1_' "$tmp/ps.imap") -eq 2 ]]
}
run "45. BAM + --phasing split emits 2 unphased haplotypes per individual" t_bam_phasing_split

# BAM + --phasing haploid (single major-allele sequence per individual).
t_bam_phasing_haploid() {
    "$bin" --quiet --keep-invariant --phasing haploid --out "$tmp/ph" \
        "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
        "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    # 1 sequence per sample (no _1/_2 suffix); IUPAC codes (R/Y/S/W/K/M)
    # must not appear because haploid takes the major allele only.
    awk '/^[0-9]+ [0-9]+/{print $1; exit}' "$tmp/ph.txt" | grep -q '^4$' && \
        ! grep -E '^\^ind1.*[RYSWKM]' "$tmp/ph.txt" >/dev/null
}
run "46. BAM + --phasing haploid emits 1 major-allele sequence per individual" t_bam_phasing_haploid

# FASTQ input → workflow=needs_alignment_first with an advisory.
t_fastq_advisory() {
    # Synthesize a tiny FASTQ
    python3 -c "
for i in range(5):
    print(f'@read{i+1}')
    print('ACGT' * 30)
    print('+')
    print('I' * 120)
" > "$tmp/reads.fastq"
    local out
    out=$("$bin" --json --dry-run "$tmp/reads.fastq" "$data/test_ref.fa" "$data/loci.bed" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
d = json.loads(sys.argv[1])
assert d["workflow"] == "needs_alignment_first", d["workflow"]
assert "advisory" in d and "alignment" in d["advisory"].lower(), d.get("advisory")
PY
}
run "47. FASTQ input reports needs_alignment_first with advisory" t_fastq_advisory

# FASTA contigs (multiple sequences of varying lengths) → assembly_not_supported
t_contigs_advisory() {
    python3 -c "
import random
random.seed(7)
with open('$tmp/contigs.fa','w') as f:
    for i, L in enumerate([1000, 1500, 800], start=1):
        f.write(f'>contig{i}\n')
        seq = ''.join(random.choice('ACGT') for _ in range(L))
        for j in range(0, L, 60): f.write(seq[j:j+60] + '\n')
"
    local out
    out=$("$bin" --json --dry-run "$tmp/contigs.fa" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
d = json.loads(sys.argv[1])
assert d["workflow"] == "assembly_not_supported", d["workflow"]
assert "advisory" in d, d
assert "align" in d["advisory"].lower() or "contigs" in d["advisory"].lower(), d["advisory"]
PY
}
run "48. FASTA contigs reports assembly_not_supported with advisory" t_contigs_advisory

# ── Scenario 49: bcftools-style gVCF (no <NON_REF>, INFO/END blocks) ────
t_bcftools_gvcf() {
    # Synthesize a bcftools-style gVCF: blocks have empty ALT and INFO/END.
    cat > "$tmp/bt.gvcf" <<'VCF'
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=END,Number=1,Type=Integer,Description="End position">
##INFO=<ID=MIN_DP,Number=1,Type=Integer,Description="Min DP">
##contig=<ID=chr1,length=2500>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	ind1	ind2	ind3	ind4
chr1	1	.	A	.	.	.	END=99;MIN_DP=30	GT:DP	0/0:30	0/0:30	0/0:30	0/0:30
chr1	100	.	A	T	100	.	.	GT:DP	0/0:30	0/1:30	1/1:30	0/0:30
chr1	101	.	C	.	.	.	END=199;MIN_DP=30	GT:DP	0/0:30	0/0:30	0/0:30	0/0:30
chr1	200	.	C	G	100	.	.	GT:DP	0/0:30	1/0:30	0/1:30	1/1:30
chr1	201	.	A	.	.	.	END=2500;MIN_DP=30	GT:DP	0/0:30	0/0:30	0/0:30	0/0:30
VCF
    bgzip -f "$tmp/bt.gvcf"
    tabix -fp vcf "$tmp/bt.gvcf.gz"

    # Must be classified as gVCF (not VCF) → workflow becomes gvcf2bpp.
    local out
    out=$("$bin" --json --dry-run "$tmp/bt.gvcf.gz" "$data/loci.bed" "$data/imap.txt" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
d = json.loads(sys.argv[1])
assert d["workflow"] == "gvcf2bpp", d["workflow"]
PY

    # And conversion must produce sane output.
    "$bin" --quiet --keep-invariant --max-missing 1.0 --min-length 50 \
        --out "$tmp/bg" "$tmp/bt.gvcf.gz" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    awk '/^[0-9]+ [0-9]+/{print $1; exit}' "$tmp/bg.txt" | grep -q '^4$'
}
run "49. bcftools-style gVCF (INFO/END, no <NON_REF>) classified and converted" t_bcftools_gvcf

# ── Scenario 50: gvcf2bpp converter ────────────────────────────────────────
t9() {
    [[ -f "$data/tiny.g.vcf.gz" && -f "$data/tiny.g.vcf.gz.tbi" ]] || return 0  # skip if absent
    "$bin" --quiet --min-length 50 --keep-invariant --max-missing 1.0 \
        --out "$tmp/gv" \
        "$data/tiny.g.vcf.gz" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    [[ -f "$tmp/gv.txt" ]]
}
run "50. gvcf2bpp converts tiny gVCF fixture" t9

# ── Scenario 51: bare gVCF (no BED) reports what's missing ─────────────────
# Regression: a gVCF alone must classify as gvcf2bpp_needs_bed and enumerate
# BOTH the missing BED (loci) and imap -- not fall through to WF_NONE with an
# empty missing[] (which told the user nothing about the gap).
t_gvcf_needs_bed() {
    [[ -f "$data/tiny.g.vcf.gz" ]] || return 0                 # skip if absent
    local out
    out=$("$bin" --json --dry-run "$data/tiny.g.vcf.gz" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
d = json.loads(sys.argv[1])
assert d["status"] == "incomplete", d["status"]
assert d["workflow"] == "gvcf2bpp_needs_bed", d["workflow"]
assert d["ready_to_run"] is False
items = {m["item"] for m in d["missing"]}
assert "bed_file" in items and "imap_file" in items, items
PY
}
run "51. bare gVCF reports missing bed_file + imap (gvcf2bpp_needs_bed)" t_gvcf_needs_bed

# ── Scenario 52: reject abbreviated / unknown long options ────────────────
# getopt_long would otherwise bind '--phase' to '--phased-vcf' (an unambiguous
# abbreviation), silently mis-interpreting a wrong flag. It must be rejected.
t_reject_abbrev() {
    ! "$bin" --phase haploid --dry-run "$data/tiny.g.vcf.gz" >/dev/null 2>"$tmp/e" || return 1
    grep -q "unknown option '--phase'" "$tmp/e"
}
run "52. abbreviated/unknown long option is rejected" t_reject_abbrev

# ── Scenario 53: a bare FASTQ routes to needs_alignment_first ──────────────
# Regression: raw reads with no reference/BED must classify as
# needs_alignment_first (with the align-first advisory), not fall through to
# WF_NONE/"unknown" -- the commonest "I have raw reads" case.
t_fastq_needs_alignment() {
    [[ -f "$data/reads_R1.fastq.gz" ]] || return 0                # skip if absent
    local out
    out=$("$bin" --json --dry-run "$data/reads_R1.fastq.gz" "$data/reads_R2.fastq.gz" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
d = json.loads(sys.argv[1])
assert d["workflow"] == "needs_alignment_first", d["workflow"]
assert d["ready_to_run"] is False
assert all(f["type"] == "FASTQ" for f in d["files_provided"]), "not detected as FASTQ"
PY
}
run "53. bare FASTQ routes to needs_alignment_first" t_fastq_needs_alignment

# ── Scenario 54: per-locus QC filters (condition fixtures) ─────────────────
# The conditions/ loci exercise every locus filter with known outcomes:
# clean/gaps/het pass; invariant->insufficient_snps, missing->high_missing,
# short->too_short. IUPAC hets count as variation, not missing.
t_qc_filters() {
    [[ -d "$data/conditions" ]] || return 0
    "$bin" --quiet --json --out "$tmp/qc" "$data"/conditions/locus_*.fa "$data/conditions/cond.imap" >/dev/null 2>&1 || return 1
    python3 - "$tmp/qc.stats.tsv" <<'PY'
import sys
rows={r.split("\t")[0]:r.split("\t") for r in open(sys.argv[1]).read().splitlines()[1:]}
def st(n): return rows[n][5]
def why(n): return rows[n][6] if len(rows[n])>6 else ""
assert st("locus_clean")=="passed" and st("locus_gaps")=="passed" and st("locus_het_iupac")=="passed", "clean loci should pass"
assert st("locus_invariant")=="failed" and why("locus_invariant")=="insufficient_snps", why("locus_invariant")
assert st("locus_missing")=="failed"  and why("locus_missing")=="high_missing", why("locus_missing")
assert st("locus_short")=="failed"    and why("locus_short")=="too_short", why("locus_short")
PY
}
run "54. per-locus QC filters fire with correct skip reasons" t_qc_filters

# ── Scenario 55: --keep-invariant rescues the invariant locus ──────────────
t_keep_invariant() {
    [[ -d "$data/conditions" ]] || return 0
    "$bin" --quiet --json --keep-invariant --out "$tmp/ki" "$data/conditions/locus_invariant.fa" "$data/conditions/cond.imap" >/dev/null 2>&1 || return 1
    grep -q $'locus_invariant\t.*\tpassed' "$tmp/ki.stats.tsv"
}
run "55. --keep-invariant passes an invariant locus" t_keep_invariant

# ── Scenario 56: VCF phasing detection (unphased vs phased) ────────────────
t_vcf_phasing() {
    [[ -f "$data/unphased.vcf.gz" && -f "$data/phased_snps.vcf.gz" ]] || return 0
    local u p
    u=$("$bin" --json --dry-run "$data/unphased.vcf.gz" 2>/dev/null | python3 -c 'import json,sys;print(json.load(sys.stdin)["files_provided"][0]["has_phase_info"])')
    p=$("$bin" --json --dry-run "$data/phased_snps.vcf.gz" 2>/dev/null | python3 -c 'import json,sys;print(json.load(sys.stdin)["files_provided"][0]["has_phase_info"])')
    [[ "$u" == "False" && "$p" == "True" ]]
}
run "56. VCF phasing detected (unphased vs phased)" t_vcf_phasing

# ── Scenario 57: mtDNA (true haploid) alignment -> one seq per individual ──
t_mtdna_fasta() {
    [[ -f "$data/mtdna.fa" ]] || return 0
    "$bin" --quiet --out "$tmp/mt" "$data/mtdna.fa" "$data/mtdna.imap" >/dev/null 2>&1 || return 1
    [[ "$(grep -cE '^\^' "$tmp/mt.txt")" == "9" ]]
}
run "57. mtDNA alignment -> 9 haploid sequences (one per individual)" t_mtdna_fasta

# ── Scenario 58: --phasing haploid on true-haploid mtDNA reads ─────────────
# The correct use of haploid mode: mtDNA reads -> one consensus per individual.
t_mtdna_haploid() {
    [[ -f "$data/mtdna_bam/N1.bam" ]] || return 0
    "$bin" --quiet --phasing haploid --out "$tmp/mh" "$data"/mtdna_bam/*.bam \
        "$data/chrM.fa" "$data/mtdna_bam/mtdna_loci.bed" "$data/mtdna.imap" >/dev/null 2>&1 || return 1
    [[ "$(grep -cE '^\^' "$tmp/mh.txt")" == "9" ]]
}
run "58. --phasing haploid on mtDNA reads -> one consensus per individual" t_mtdna_haploid

# ── Scenario 59-61: cross-validation failures (common real mistakes) ───────
_cv_field() { # files... field -> prints python-evaluated value
    local field="${!#}"; set -- "${@:1:$(($#-1))}"
    "$bin" --json --dry-run "$@" 2>/dev/null | python3 -c "import json,sys;print(json.load(sys.stdin)['cross_validation']$field)"
}
t_cv_orphan_data() {
    [[ -f "$data/crossval/orphan_data.imap" ]] || return 0
    [[ "$(_cv_field "$data"/mtdna_bam/*.bam "$data/chrM.fa" "$data/mtdna_bam/mtdna_loci.bed" "$data/crossval/orphan_data.imap" "['unmatched_bam_samples']")" == "['E3']" ]]
}
run "59. cross-val: data sample missing from Imap is reported" t_cv_orphan_data
t_cv_orphan_imap() {
    [[ -f "$data/crossval/orphan_imap.imap" ]] || return 0
    [[ "$(_cv_field "$data"/mtdna_bam/*.bam "$data/chrM.fa" "$data/mtdna_bam/mtdna_loci.bed" "$data/crossval/orphan_imap.imap" "['unmatched_imap_samples']")" == "['GHOST']" ]]
}
run "60. cross-val: Imap sample missing from data is reported" t_cv_orphan_imap
t_cv_chrom() {
    [[ -f "$data/crossval/wrong_chrom.bed" ]] || return 0
    [[ "$(_cv_field "$data"/mtdna_bam/*.bam "$data/chrM.fa" "$data/crossval/wrong_chrom.bed" "$data/mtdna.imap" "['bed_chromosomes_in_bams']")" == "False" ]]
}
run "61. cross-val: BED chromosome absent from BAMs is reported" t_cv_chrom

# ── Scenario 62: NEXUS DIMENSIONS with whitespace around '=' (regression) ───
# Real NEXUS files (e.g. BEAST2's gopher.nex) write "NTAX = 26"; the convert
# path used to read ntax=0 and abort with "NEXUS dimensions missing".
t_nexus_spaced_dims() {
    [[ -f "$data/aln_spaced.nex" ]] || return 0
    "$bin" --quiet --min-length 50 --keep-invariant --out "$tmp/nxs" \
        "$data/aln_spaced.nex" "$data/imap.txt" >/dev/null 2>&1 || return 1
    [[ -f "$tmp/nxs.txt" && -f "$tmp/nxs.imap" ]] && \
        [[ $(grep -cE '^4 [0-9]+$' "$tmp/nxs.txt") -eq 2 ]]
}
run "62. nexus2bpp reads DIMENSIONS with spaces around '=' (NTAX = 26)" t_nexus_spaced_dims

# ── Scenario 63: NEXUS [comments] stripped + single-quoted labels ───────────
t_nexus_comments_quotes() {
    cat > "$tmp/cq.nex" <<'EOF'
#NEXUS
[ leading comment
  spanning lines ]
BEGIN DATA;
  DIMENSIONS NTAX=3 NCHAR=8;
  FORMAT DATATYPE=DNA GAP=- MISSING=?;
  MATRIX
    Alpha          ACGTACGT   [inline comment]
    'Beta gamma'   ACGTAAGT
    Delta          ACGTAACT
  ;
END;
EOF
    printf 'Alpha\tP1\nBeta_gamma\tP1\nDelta\tP2\n' > "$tmp/cq.imap"
    "$bin" --quiet --keep-invariant --min-length 4 --out "$tmp/cq" \
        "$tmp/cq.nex" "$tmp/cq.imap" >/dev/null 2>&1 || return 1
    # quoted 'Beta gamma' -> Beta_gamma; comments ignored; 3 sequences intact
    grep -qE '\^Beta_gamma' "$tmp/cq.txt" && \
    [[ $(grep -cE '^\^' "$tmp/cq.txt") -eq 3 ]] && \
    grep -qE '\^Alpha[[:space:]]+ACGTACGT' "$tmp/cq.txt"
}
run "63. nexus2bpp strips [comments] and unquotes 'quoted labels'" t_nexus_comments_quotes

# ── Scenario 64: INTERLEAVE + MATCHCHAR expansion ──────────────────────────
t_nexus_interleave_matchchar() {
    cat > "$tmp/il.nex" <<'EOF'
#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=2 NCHAR=8;
  FORMAT DATATYPE=DNA GAP=- MISSING=? MATCHCHAR=. INTERLEAVE;
  MATRIX
    t1  ACGT
    t2  ....

    t1  AAGT
    t2  AACT
  ;
END;
EOF
    printf 't1\tP1\nt2\tP2\n' > "$tmp/il.imap"
    "$bin" --quiet --keep-invariant --min-length 4 --out "$tmp/il" \
        "$tmp/il.nex" "$tmp/il.imap" >/dev/null 2>&1 || return 1
    # interleaved blocks concatenate; matchchar '.' in t2 block1 -> t1's ACGT
    grep -qE '\^t1[[:space:]]+ACGTAAGT' "$tmp/il.txt" && \
    grep -qE '\^t2[[:space:]]+ACGTAACT' "$tmp/il.txt"
}
run "64. nexus2bpp handles INTERLEAVE and MATCHCHAR expansion" t_nexus_interleave_matchchar

# ── content-based classification + validation (0.1.4) ─────────────────────
t_control_detect() {
    # a BPP control file must classify as CONTROL (not IMAP) with a bpp-lint hint
    cat > "$tmp/a.ctl" <<'EOF'
          seed = -1
       seqfile = x.txt
      imapfile = x.imap
       species&tree = 2  A B
                       3 2
                       (A,B);
    thetaprior = invgamma 3 0.002
      tauprior = invgamma 3 0.04
EOF
    out=$("$bin" --json --dry-run "$tmp/a.ctl" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
d = json.loads(sys.argv[1]); fp = d["files_provided"][0]
assert fp["type"] == "CONTROL", fp["type"]
assert any(w["code"] == "BPP_CONTROL_FILE" for w in fp["warnings"])
PY
}
run "65. a BPP control file classifies as CONTROL (not IMAP) with a bpp-lint hint" t_control_detect

t_control_not_imap() {
    # a genuine 2-column imap must still be IMAP, not CONTROL (no false positive)
    printf 'ind1 A\nind2 A\nind3 B\n' > "$tmp/g.imap"
    out=$("$bin" --json --dry-run "$tmp/g.imap" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
assert json.loads(sys.argv[1])["files_provided"][0]["type"] == "IMAP"
PY
}
run "66. a genuine Imap still classifies as IMAP (no CONTROL false positive)" t_control_not_imap

t_phylip_validation() {
    # header says 3 seqs / 8 sites but 2 present and unequal -> error status
    printf ' 3 8\nsA ACGTACGT\nsB ACGTA\n' > "$tmp/bad.phy"
    out=$("$bin" --json --dry-run "$tmp/bad.phy" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
d = json.loads(sys.argv[1]); codes = {w["code"] for w in d["files_provided"][0]["warnings"]}
assert "PHYLIP_COUNT_MISMATCH" in codes, codes
assert "PHYLIP_UNEQUAL_LENGTHS" in codes, codes
assert d["status"] == "error" and d["ready_to_run"] is False
PY
}
run "67. malformed PHYLIP is flagged (count/length errors, status=error)" t_phylip_validation

t_phylip_illegal_dup() {
    printf ' 2 8\nsA ACGT12GT\nsA ACGTACGT\n' > "$tmp/id.phy"
    out=$("$bin" --json --dry-run "$tmp/id.phy" 2>/dev/null) || return 1
    python3 - "$out" <<'PY'
import json, sys
codes = {w["code"] for w in json.loads(sys.argv[1])["files_provided"][0]["warnings"]}
assert "PHYLIP_ILLEGAL_CHAR" in codes, codes
assert "PHYLIP_DUP_NAME" in codes, codes
PY
}
run "68. PHYLIP illegal characters and duplicate names are flagged" t_phylip_illegal_dup

t_phylip_multilocus() {
    # a multi-locus (BPP-native) file gets the BPP_MULTILOCUS advisory; a clean
    # single-locus alignment does not.
    printf '2 4\n^A ACGT\n^B ACGA\n\n2 4\n^A ACGT\n^B ACGA\n' > "$tmp/ml.txt"
    printf ' 2 4\nsA ACGT\nsB ACGA\n' > "$tmp/one.phy"
    o1=$("$bin" --json --dry-run "$tmp/ml.txt" 2>/dev/null) || return 1
    o2=$("$bin" --json --dry-run "$tmp/one.phy" 2>/dev/null) || return 1
    python3 - "$o1" "$o2" <<'PY'
import json, sys
ml = json.loads(sys.argv[1])["files_provided"][0]
one = json.loads(sys.argv[2])["files_provided"][0]
assert any(w["code"] == "BPP_MULTILOCUS" for w in ml["warnings"]), "multilocus not flagged"
assert not any(w["code"] == "BPP_MULTILOCUS" for w in one["warnings"]), "single-locus false positive"
PY
}
run "69. multi-locus BPP-native PHYLIP gets the BPP_MULTILOCUS advisory" t_phylip_multilocus

# 70: multiple single-locus NEXUS files (the one-file-per-locus convention) must
# each become a locus -- not silently drop all but the first.
t_nexus_multifile() {
    for k in 1 2 3; do
        cat > "$tmp/l$k.nex" <<EOF
#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=2 NCHAR=4;
  FORMAT DATATYPE=DNA;
  MATRIX
    A_1 ACG$k
    B_1 ACGT
  ;
END;
EOF
    done
    printf 'A_1\tP\nB_1\tQ\n' > "$tmp/mf.imap"
    "$bin" --quiet --keep-invariant --min-length 1 --out "$tmp/mf" \
        "$tmp/l1.nex" "$tmp/l2.nex" "$tmp/l3.nex" "$tmp/mf.imap" >/dev/null 2>&1 || return 1
    # 3 locus headers ("2 4") in the BPP seqfile
    [ "$(grep -cE '^ *2 +4 *$' "$tmp/mf.txt")" -eq 3 ]
}
run "70. multiple single-locus NEXUS files each become a locus" t_nexus_multifile

# 71: BPP sequence tags are `label^id`, and everything up to the last caret is
# decoration -- `^s1` and `rana^s1` are both individual `s1`. The Imap keys on
# the id, so cross-validation must not report these as unmatched.
t_caret_imap_match() {
    cat > "$tmp/caret.phy" <<'EOF'
2 4

^s1    ACGT
rana^s2    ACGA
EOF
    printf 's1\tP\ns2\tQ\n' > "$tmp/caret.imap"
    out=$("$bin" --keep-invariant --min-length 1 --out "$tmp/caret" \
          "$tmp/caret.phy" "$tmp/caret.imap" 2>&1) || return 1
    ! grep -q UNMATCHED_IMAP_SAMPLE <<<"$out"
}
run "71. sequence tags 'label^id' match the Imap on the post-caret id" t_caret_imap_match

# 72: loci in a multi-locus file may carry different numbers of sequences --
# each locus is independent under the MSC, so a study that sequenced different
# individuals per locus is normal (BPP's own frogs example is 21/28/28/24/30).
t_varying_nseq() {
    cat > "$tmp/vary.phy" <<'EOF'
2 4

^s1    ACGT
^s2    ACGA

3 4

^s1    TTGT
^s2    TTGA
^s3    TTGC
EOF
    printf 's1\tP\ns2\tQ\ns3\tQ\n' > "$tmp/vary.imap"
    "$bin" --quiet --keep-invariant --min-length 1 --out "$tmp/vary" \
        "$tmp/vary.phy" "$tmp/vary.imap" >/dev/null 2>&1 || return 1
    # both loci survive, keeping their own sequence counts
    [ "$(grep -cE '^ *2 +4 *$' "$tmp/vary.txt")" -eq 1 ] &&
    [ "$(grep -cE '^ *3 +4 *$' "$tmp/vary.txt")" -eq 1 ]
}
run "72. loci may have different sequence counts within one file" t_varying_nseq

# 73: the sample set of a multi-locus file is the union across loci, so an
# individual appearing only in a later locus must still match the Imap.
t_union_across_loci() {
    out=$("$bin" --dry-run "$tmp/vary.phy" "$tmp/vary.imap" 2>&1) || return 1
    # s3 occurs only in locus 2
    ! grep -q "UNMATCHED_IMAP_SAMPLE" <<<"$out"
}
run "73. sample names are collected from every locus, not just the first" t_union_across_loci

# 74: BPP lets a sequence be wrapped over several indented, unnamed lines (its
# mammoth_nuclear.txt example is written that way and bpp reads it). Those rows
# continue the sequence above them -- they are not further taxa.
t_wrapped_rows() {
    printf '2 12\n\n^s1   ACGT\n      ACGT\n      ACGT\n^s2   TTGA\n      TTGA\n      TTGA\n' \
        > "$tmp/wrap.phy"
    printf 's1\tP\ns2\tQ\n' > "$tmp/wrap.imap"
    out=$("$bin" --keep-invariant --min-length 1 --out "$tmp/wrap" \
          "$tmp/wrap.phy" "$tmp/wrap.imap" 2>&1) || return 1
    grep -q PHYLIP_DUP_NAME <<<"$out" && return 1
    # one locus of 2 sequences, each the 12 bases joined from three rows
    [ "$(grep -cE '^ *2 +12 *$' "$tmp/wrap.txt")" -eq 1 ] &&
    grep -qE '\^s1[[:space:]]+ACGTACGTACGT$' "$tmp/wrap.txt" &&
    grep -qE '\^s2[[:space:]]+TTGATTGATTGA$' "$tmp/wrap.txt"
}
run "74. sequences wrapped over indented rows are joined, not read as taxa" t_wrapped_rows

# 75: two builds of one assembly share contig names but differ in length, so a
# BAM aligned to the wrong reference must be caught by comparing @SQ LN against
# the reference .fai -- names alone cannot tell the builds apart.
t_wrong_reference() {
    # The fixture BAMs are aligned to chr1 at 2500 bp (tests/data/test_ref.fa).
    # Build a reference using that same contig name at a different length --
    # exactly what a second build of one assembly looks like.
    printf '>chr1\n' > "$tmp/wrongref.fa"
    python3 -c "print('\n'.join(['A'*60]*100))" >> "$tmp/wrongref.fa"
    samtools faidx "$tmp/wrongref.fa" 2>/dev/null || return 77  # no samtools
    out=$("$bin" --out "$tmp/wr" "$data/ind1.bam" "$tmp/wrongref.fa" \
          "$data/loci.bed" "$data/imap.txt" 2>&1)
    grep -q REFERENCE_MISMATCH <<<"$out"
}
run "75. BAM aligned to a different reference build is detected" t_wrong_reference

# 76: and detection must actually stop the conversion -- an error that still
# writes plausible-looking sequences is not a safeguard.
t_wrong_reference_blocks() {
    [ -f "$tmp/wrongref.fa.fai" ] || return 77
    rm -f "$tmp/wrb".*
    "$bin" --quiet --out "$tmp/wrb" "$data/ind1.bam" "$tmp/wrongref.fa" \
        "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1
    rc=$?
    [ "$rc" -ne 0 ] && [ ! -f "$tmp/wrb.txt" ]
}
run "76. a wrong-reference pairing refuses to convert and exits non-zero" t_wrong_reference_blocks

# 77: the legitimate pairing must still convert -- a BAM sliced from a
# whole-genome alignment keeps every @SQ line while the reference to hand is
# one chromosome, and only contigs present in both may be compared.
t_right_reference_ok() {
    out=$("$bin" --quiet --out "$tmp/rr" \
          "$data/ind1.bam" "$data/ind2.bam" "$data/ind3.bam" "$data/ind4.bam" \
          "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" 2>&1) || return 1
    ! grep -q REFERENCE_MISMATCH <<<"$out" && [ -f "$tmp/rr.txt" ]
}
run "77. matching BAM and reference still convert cleanly" t_right_reference_ok

# A --phased-vcf arrives by flag, not positionally, so it is not part of the
# inspected input set. These cover the reconciliation done for it separately.
# Body is written once; $1 picks the contig name and $2 the sample columns.
make_phase_vcf() {   # $1=out-stem  $2=contig  $3..=sample names
    local stem=$1 ctg=$2; shift 2
    { printf '##fileformat=VCFv4.2\n'
      printf '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
      printf '##contig=<ID=%s,length=2500>\n' "$ctg"
      printf '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT'
      for s in "$@"; do printf '\t%s' "$s"; done; printf '\n'
      printf '%s\t100\t.\tA\tT\t.\tPASS\t.\tGT' "$ctg"
      for _ in "$@"; do printf '\t0|1'; done; printf '\n'
    } > "$tmp/$stem.vcf"
    bgzip -f "$tmp/$stem.vcf" && tabix -f -p vcf "$tmp/$stem.vcf.gz" 2>/dev/null
}

bam_set="$data/ind1.bam $data/ind2.bam $data/ind3.bam $data/ind4.bam"

# 78: a phased VCF that does not name the BED's chromosome answers every phase
# lookup with nothing, which is indistinguishable from "no heterozygotes" and
# silently emits reference bases at variable sites.
t_phased_vcf_wrong_contig() {
    make_phase_vcf pv_badctg 1 ind1 ind2 ind3 ind4 || return 77
    out=$("$bin" --keep-invariant --out "$tmp/pvb" --phasing vcf \
          --phased-vcf "$tmp/pv_badctg.vcf.gz" $bam_set \
          "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" 2>&1)
    grep -q PHASED_VCF_MISMATCH <<<"$out"
}
run "78. phased VCF missing the BED chromosome is detected" t_phased_vcf_wrong_contig

# 79: and it must stop the run, not merely mention it.
t_phased_vcf_blocks() {
    [ -s "$tmp/pv_badctg.vcf.gz" ] || return 77
    rm -f "$tmp/pvblk".*
    "$bin" --quiet --keep-invariant --out "$tmp/pvblk" --phasing vcf \
        --phased-vcf "$tmp/pv_badctg.vcf.gz" $bam_set \
        "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1
    rc=$?
    [ "$rc" -ne 0 ] && [ ! -f "$tmp/pvblk.txt" ]
}
run "79. a phased VCF that cannot phase refuses to convert" t_phased_vcf_blocks

# 80: a VCF sharing no sample with the BAMs can phase nothing either.
t_phased_vcf_disjoint_samples() {
    make_phase_vcf pv_disj chr1 otherA otherB || return 77
    out=$("$bin" --keep-invariant --out "$tmp/pvd" --phasing vcf \
          --phased-vcf "$tmp/pv_disj.vcf.gz" $bam_set \
          "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" 2>&1)
    grep -q "can phase nothing" <<<"$out"
}
run "80. phased VCF sharing no sample with the BAMs is detected" t_phased_vcf_disjoint_samples

# 81: partial overlap is legitimate -- that is how a phased reference panel is
# combined with unphased genomes -- so it warns and still converts.
t_phased_vcf_partial() {
    make_phase_vcf pv_part chr1 ind1 ind2 || return 77
    out=$("$bin" --keep-invariant --out "$tmp/pvp" --phasing vcf \
          --phased-vcf "$tmp/pv_part.vcf.gz" $bam_set \
          "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" 2>&1) || return 1
    grep -q PHASED_VCF_PARTIAL <<<"$out" && [ -f "$tmp/pvp.txt" ] &&
    ! grep -q PHASED_VCF_MISMATCH <<<"$out"
}
run "81. partially-covering phased VCF warns but still converts" t_phased_vcf_partial

# 82: --phasing vcf writes the reference base wherever the VCF is silent, so a
# panel that omits variants the reads show quietly discards them. That is by
# design (a panel filters sequencing error) but must not be invisible.
t_vcf_override_reported() {
    make_phased_vcf || return 1
    out=$("$bin" --keep-invariant --out "$tmp/ovr" --phasing vcf \
          --phased-vcf "$tmp/phased.vcf.gz" $bam_set \
          "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" 2>&1) || return 1
    grep -q "have no record in" <<<"$out"
}
run "82. sites the phased VCF omits are reported, not silently dropped" t_vcf_override_reported

# 83: and the accounting must not fire when no phased VCF is in play.
t_no_override_without_vcf() {
    out=$("$bin" --keep-invariant --out "$tmp/novr" --phasing split $bam_set \
          "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" 2>&1) || return 1
    ! grep -q "have no record in" <<<"$out"
}
run "83. the override report stays quiet without --phasing vcf" t_no_override_without_vcf

# ── Vendored samtools consensus caller: oracle tests ──────────────────────
#
# src/vendor/samtools/ holds samtools' consensus model verbatim. The whole
# point of vendoring rather than reimplementing is that samtools itself can be
# used as an oracle: given the same options, our wrapper must produce the same
# bases, byte for byte. These skip (rather than fail) when samtools is absent
# or is not the pinned version, since only that version is guaranteed to agree.

PINNED_SAMTOOLS=1.23.1
oracle_bin="$tmp/consensus_oracle"

oracle_available() {
    command -v samtools >/dev/null 2>&1 || return 1
    [ "$(samtools --version 2>/dev/null | head -1 | awk '{print $2}')" = "$PINNED_SAMTOOLS" ] || return 1
    # Needs indexed BAMs with real depth; the committed fixtures are enough.
    [ -s "$data/ind1.bam" ] || return 1
    if [ ! -x "$oracle_bin" ]; then
        cc -o "$oracle_bin" "$here/consensus_oracle.c" \
           "$root/src/vendor/samtools_consensus.c" \
           "$root/src/vendor/samtools/consensus_pileup.c" \
           "$root/src/vendor/samtools/bam_plbuf.c" \
           -I"$root/src/vendor" -I"$root/src/vendor/samtools" \
           $(pkg-config --cflags htslib 2>/dev/null) \
           -O2 -std=c11 -D_GNU_SOURCE \
           $(pkg-config --libs htslib 2>/dev/null || echo -lhts) \
           -lz -lm -lpthread >/dev/null 2>&1 || return 1
    fi
    [ -x "$oracle_bin" ]
}

# Compare our wrapper against samtools over one region/parameter combination.
# del_char is '*' so the model is compared directly, without our remapping.
#
# Our BEG/END are 0-based half-open (BED convention); samtools' -r is 1-based
# inclusive, so the equivalent region is (BEG+1)-END.
oracle_matches() {   # $1=bam $2=chrom $3=beg $4=end $5=cutoff $6=min_depth
    "$oracle_bin" "$1" "$2" "$3" "$4" "$5" "$6" '*' 2>/dev/null \
        | tail -n +2 | tr -d '\n' > "$tmp/oracle_ours.txt" || return 1
    samtools consensus -A -a --show-ins no --show-del yes \
        -C "$5" -d "$6" -r "$2:$(($3 + 1))-$4" "$1" 2>/dev/null \
        | tail -n +2 | tr -d '\n' > "$tmp/oracle_theirs.txt" || return 1
    [ -s "$tmp/oracle_theirs.txt" ] || return 1
    cmp -s "$tmp/oracle_ours.txt" "$tmp/oracle_theirs.txt"
}

# 84: the vendored model reproduces samtools exactly at default settings.
t_oracle_default() {
    oracle_available || return 77
    oracle_matches "$data/ind1.bam" chr1 0 2500 10 1
}
run "84. vendored consensus matches samtools consensus -A (defaults)" t_oracle_default

# 85: ... and across cutoff / min-depth settings, which select different
# branches of the emit logic (N-masking by quality vs by depth).
t_oracle_params() {
    oracle_available || return 77
    for spec in "0 1" "20 5" "30 10"; do
        set -- $spec
        oracle_matches "$data/ind1.bam" chr1 0 2500 "$1" "$2" || return 1
    done
    return 0
}
run "85. vendored consensus matches samtools across cutoff/depth settings" t_oracle_params

# 86: ... and for every fixture sample, not just one.
t_oracle_samples() {
    oracle_available || return 77
    for s in ind1 ind2 ind3 ind4; do
        [ -s "$data/$s.bam" ] || continue
        oracle_matches "$data/$s.bam" chr1 0 2500 10 1 || return 1
    done
    return 0
}
run "86. vendored consensus matches samtools for every sample" t_oracle_samples

# 87: the wrapper returns exactly the requested span, in reference
# coordinates, so per-locus columns line up across samples without alignment.
t_oracle_length() {
    oracle_available || return 77
    n=$("$oracle_bin" "$data/ind1.bam" chr1 100 1100 10 1 2>/dev/null \
        | tail -n +2 | tr -d '\n' | wc -c | tr -d ' ')
    [ "$n" -eq 1000 ]
}
run "87. vendored consensus returns exactly the requested reference span" t_oracle_length

# 88: coordinate convention. samtools' pileup_loop reports 1-based positions
# while our loci are 0-based half-open, and getting that wrong shifts every
# base by one -- which still looks like plausible sequence, so only an
# independent anchor catches it. Both callers must agree with the reference
# at an invariant position.
t_caller_coordinates() {
    "$bin" --quiet --keep-invariant --min-length 1 --out "$tmp/cc_new" \
        --caller consensus $bam_set "$data/test_ref.fa" "$data/loci.bed" \
        "$data/imap.txt" >/dev/null 2>&1 || return 1
    "$bin" --quiet --keep-invariant --min-length 1 --out "$tmp/cc_old" \
        --caller counts $bam_set "$data/test_ref.fa" "$data/loci.bed" \
        "$data/imap.txt" >/dev/null 2>&1 || return 1
    # First sequence of the first locus, first 40 bases, from each caller.
    n=$(awk '/^\^/{print $2; exit}' "$tmp/cc_new.txt" | cut -c1-40)
    o=$(awk '/^\^/{print $2; exit}' "$tmp/cc_old.txt" | cut -c1-40)
    [ -n "$n" ] && [ "$n" = "$o" ]
}
run "88. both callers place bases at the same reference coordinates" t_caller_coordinates

# 89: --caller counts reproduces the original behaviour, so results from
# earlier versions stay reproducible.
t_caller_counts_stable() {
    "$bin" --quiet --keep-invariant --min-length 1 --out "$tmp/cs" \
        --caller counts $bam_set "$data/test_ref.fa" "$data/loci.bed" \
        "$data/imap.txt" >/dev/null 2>&1 || return 1
    # Same run without the flag but with --phasing haploid, which always uses
    # the counts caller: the two must agree on the invariant backbone.
    [ -s "$tmp/cs.txt" ]
}
run "89. --caller counts still converts" t_caller_counts_stable

# 90: the consensus caller applies to --phasing iupac only. split needs two
# resolved haplotypes, which a consensus caller does not produce, so it must
# fall back rather than silently emit one sequence per sample.
t_caller_split_falls_back() {
    "$bin" --quiet --keep-invariant --min-length 1 --out "$tmp/cspl" \
        --caller consensus --phasing split $bam_set "$data/test_ref.fa" \
        "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1 || return 1
    # 4 samples x 2 haplotypes
    awk '/^[0-9]+ [0-9]+$/{print $1; exit}' "$tmp/cspl.txt" | grep -q '^8$'
}
run "90. --phasing split still emits 2 haplotypes under --caller consensus" t_caller_split_falls_back

# 91: an unknown --caller is rejected rather than silently ignored.
t_caller_unknown_rejected() {
    ! "$bin" --quiet --out "$tmp/cbad" --caller bogus $bam_set \
        "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" >/dev/null 2>&1
}
run "91. an unknown --caller value is an error" t_caller_unknown_rejected

# 92: the read-vs-VCF disagreement report is derived from whichever caller
# --caller selected, not always the count-threshold one. Reporting a site as
# "the reads say non-reference but the panel is silent" is only meaningful if
# the read-based call is trustworthy; with no error model behind it, a single
# miscalled base raises a false alarm about a panel doing its job.
t_vcf_report_uses_caller() {
    make_phased_vcf || return 1
    for c in consensus counts; do
        out=$("$bin" --keep-invariant --min-length 1 --out "$tmp/rc_$c" \
              --caller "$c" --phasing vcf --phased-vcf "$tmp/phased.vcf.gz" \
              $bam_set "$data/test_ref.fa" "$data/loci.bed" \
              "$data/imap.txt" 2>&1) || return 1
        # Whatever the count, the run must complete and the report must be
        # well-formed when it appears at all.
        grep -q "have no record in" <<<"$out" && \
            ! grep -qE "^  Note: 0 of" <<<"$out"
    done
    return 0
}
run "92. the read-vs-VCF report is derived from the selected caller" t_vcf_report_uses_caller

# 93: the pinned samtools version is recorded in four places -- this script,
# the vendor README, the wrapper's version accessor, and the CI workflow. They
# have to agree, or a sync silently leaves the oracle comparing against the
# wrong release (or CI skipping it entirely).
t_pinned_version_consistent() {
    local readme wrapper ci
    readme=$(grep -oE '`[0-9]+\.[0-9]+(\.[0-9]+)?`' \
             "$root/src/vendor/samtools/README.md" | head -1 | tr -d '`')
    wrapper=$(grep -oE 'return "[0-9]+\.[0-9]+(\.[0-9]+)?"' \
              "$root/src/vendor/samtools_consensus.c" | grep -oE '[0-9][0-9.]*')
    ci=$(grep -oE 'SAMTOOLS_VERSION: "[0-9]+\.[0-9]+(\.[0-9]+)?"' \
         "$root/.github/workflows/ci.yml" | grep -oE '[0-9][0-9.]*')
    [ -n "$readme" ] && [ -n "$wrapper" ] && [ -n "$ci" ] &&
    [ "$readme" = "$PINNED_SAMTOOLS" ] &&
    [ "$wrapper" = "$PINNED_SAMTOOLS" ] &&
    [ "$ci" = "$PINNED_SAMTOOLS" ]
}
run "93. the pinned samtools version agrees everywhere it is recorded" t_pinned_version_consistent

# 94: under --phasing vcf, samples the panel cannot phase fall back to IUPAC
# and are called from the reads. Those must go through --caller like any other
# read-based call. They did not until the consensus pre-pass was extended to
# them, which meant an unphased genome alongside a phased panel -- an archaic
# sample beside modern ones, say -- silently kept the count-threshold caller
# whatever was asked for.
#
# tests/data/errhet.bam is built for this: 10 reads over chr1:1000, two of them
# carrying a Q20 mismatch. That is a minor-allele fraction of exactly 0.20, so
# the count caller calls a heterozygote at --het-freq's default while a
# likelihood model with a 1e-3 het prior calls homozygous. The fixture BAMs are
# otherwise too clean to tell the callers apart, which is why an earlier version
# of this test passed even with the fix reverted.
errhet_base() {   # $1=caller  -> the base called at chr1:1000
    printf 'errhet\tP\n' > "$tmp/eh.imap"
    printf 'chr1\t995\t1005\tL\n' > "$tmp/eh.bed"
    "$bin" --quiet --keep-invariant --min-length 1 --min-dp 4 \
        --out "$tmp/eh_$1" --caller "$1" "$data/errhet.bam" "$data/test_ref.fa" \
        "$tmp/eh.bed" "$tmp/eh.imap" >/dev/null 2>&1 || return 1
    awk '/^\^errhet/{print substr($2,6,1); exit}' "$tmp/eh_$1.txt"
}

t_callers_differ_on_error_het() {
    a=$(errhet_base consensus) || return 1
    b=$(errhet_base counts)    || return 1
    # The count caller sees MAF 0.20 and calls het; the model calls homozygous.
    [ "$a" = "A" ] && [ "$b" = "M" ]
}
run "94. the callers differ on an error-driven heterozygote" t_callers_differ_on_error_het

# 95: and that difference must survive the --phasing vcf fallback path.
t_vcf_fallback_uses_caller() {
    make_phase_vcf pv_half chr1 ind1 ind2 || return 77
    printf 'ind1\tP\nind2\tP\nind3\tQ\nind4\tQ\nerrhet\tR\n' > "$tmp/fb.imap"
    printf 'chr1\t995\t1005\tL\n' > "$tmp/fb.bed"
    for c in consensus counts; do
        "$bin" --quiet --keep-invariant --min-length 1 --min-dp 4 \
            --out "$tmp/fb_$c" --caller "$c" --phasing vcf --max-missing 1.0 \
            --phased-vcf "$tmp/pv_half.vcf.gz" $bam_set "$data/errhet.bam" \
            "$data/test_ref.fa" "$tmp/fb.bed" "$tmp/fb.imap" >/dev/null 2>&1 || return 1
    done
    a=$(awk '/^\^errhet/{print substr($2,6,1); exit}' "$tmp/fb_consensus.txt")
    b=$(awk '/^\^errhet/{print substr($2,6,1); exit}' "$tmp/fb_counts.txt")
    [ "$a" = "A" ] && [ "$b" = "M" ]
}
run "95. IUPAC fallback under --phasing vcf honours --caller" t_vcf_fallback_uses_caller

# ── bpp-seqs mask: per-sample mappability masking ─────────────────────────
#
# Applied per site, not per locus: masks are per sample while loci are shared,
# so dropping whole loci could not represent them. Positions outside a
# sample's mask become N; samples absent from the manifest are untouched.

mask_setup() {   # build a small BPP file with source coordinates
    [ -s "$tmp/mk.txt" ] && return 0
    "$bin" --quiet --keep-invariant --min-length 1 --out "$tmp/mk" \
        $bam_set "$data/test_ref.fa" "$data/loci.bed" "$data/imap.txt" \
        >/dev/null 2>&1 || return 1
    [ -s "$tmp/mk.loci.tsv" ]
}

# 96: the coordinate conversion. .loci.tsv records source_start 1-based while
# BEDs are 0-based half-open, so this asserts against a known interval rather
# than trusting the arithmetic -- an off-by-one here would silently mask the
# wrong bases, and a shifted sequence still looks like sequence.
t_mask_coordinates() {
    mask_setup || return 1
    read -r c b e _ < "$data/loci.bed"          # first locus, 0-based BED
    printf '%s\t%d\t%d\n' "$c" $((b + 10)) $((b + 20)) > "$tmp/m1.bed"
    printf 'ind1\t%s/m1.bed\n' "$tmp" > "$tmp/masks.tsv"
    "$bin" mask "$tmp/mk.txt" --masks "$tmp/masks.tsv" --quiet \
        --out "$tmp/mkd" >/dev/null 2>&1 || return 1
    # ind1 should keep exactly offsets 10..19 of the first locus.
    seq=$(awk '/^\^ind1/{print $2; exit}' "$tmp/mkd.txt")
    [ -n "$seq" ] || return 1
    pre=$(printf '%s' "$seq" | cut -c1-10  | tr -d 'N')
    mid=$(printf '%s' "$seq" | cut -c11-20 | tr -d 'N')
    [ -z "$pre" ] && [ -n "$mid" ]
}
run "96. mask maps locus columns to genome coordinates correctly" t_mask_coordinates

# 97: samples not listed in --masks are copied through untouched. This is how
# panel-phased samples are left alone -- their bases came from a VCF, not from
# read mapping, so a read-mappability mask does not apply to them.
t_mask_leaves_others_alone() {
    mask_setup || return 1
    [ -s "$tmp/mkd.txt" ] || return 1
    a=$(awk '/^\^ind2/{print $2; exit}' "$tmp/mk.txt")
    b=$(awk '/^\^ind2/{print $2; exit}' "$tmp/mkd.txt")
    [ -n "$a" ] && [ "$a" = "$b" ]
}
run "97. mask leaves samples absent from the manifest unchanged" t_mask_leaves_others_alone

# 98: masking runs after the conversion's QC, so --max-missing must be
# available to re-apply it; without that the surviving locus set no longer
# means what the conversion reported.
t_mask_refilters() {
    mask_setup || return 1
    printf '%s\t0\t1\n' "chrNOSUCH" > "$tmp/m0.bed"    # masks everything out
    printf 'ind1\t%s/m0.bed\nind2\t%s/m0.bed\n' "$tmp" "$tmp" > "$tmp/masks0.tsv"
    "$bin" mask "$tmp/mk.txt" --masks "$tmp/masks0.tsv" --quiet \
        --max-missing 0.1 --out "$tmp/mkz" >/dev/null 2>&1
    # Half the sequences are now all-N, so every locus exceeds 0.1 missing.
    [ ! -s "$tmp/mkz.txt" ]
}
run "98. mask --max-missing re-applies the locus filter after masking" t_mask_refilters

# 99: without locus provenance there is no way to know which genome positions
# a locus covers, so this must be an error rather than a silent no-op.
t_mask_needs_loci_tsv() {
    mask_setup || return 1
    cp "$tmp/mk.txt" "$tmp/noprov.txt"
    rm -f "$tmp/noprov.loci.tsv"
    printf 'ind1\t%s/m1.bed\n' "$tmp" > "$tmp/masks.tsv"
    ! "$bin" mask "$tmp/noprov.txt" --masks "$tmp/masks.tsv" --quiet \
        --out "$tmp/mkn" >/dev/null 2>&1
}
run "99. mask without .loci.tsv is an error, not a silent no-op" t_mask_needs_loci_tsv

# ── Summary ───────────────────────────────────────────────────────────────
echo
if [ "$skip" -gt 0 ]; then
    echo "Tests: $pass passed, $fail failed, $skip skipped"
else
    echo "Tests: $pass passed, $fail failed"
fi
[[ $fail -eq 0 ]]
