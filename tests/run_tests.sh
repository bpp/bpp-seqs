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
    samtools faidx "$tmp/wrongref.fa" 2>/dev/null || return 0   # skip if no samtools
    out=$("$bin" --out "$tmp/wr" "$data/ind1.bam" "$tmp/wrongref.fa" \
          "$data/loci.bed" "$data/imap.txt" 2>&1)
    grep -q REFERENCE_MISMATCH <<<"$out"
}
run "75. BAM aligned to a different reference build is detected" t_wrong_reference

# 76: and detection must actually stop the conversion -- an error that still
# writes plausible-looking sequences is not a safeguard.
t_wrong_reference_blocks() {
    [ -f "$tmp/wrongref.fa.fai" ] || return 0
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

# ── Summary ───────────────────────────────────────────────────────────────
echo
echo "Tests: $pass passed, $fail failed"
[[ $fail -eq 0 ]]
