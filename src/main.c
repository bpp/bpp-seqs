/* main.c — bpp-seqs entry point. */

#define _POSIX_C_SOURCE 200809L

#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <htslib/faidx.h>

#include "inspect.h"
#include "cross_validate.h"
#include "workflow.h"
#include "output.h"
#include "sanity.h"
#include "cmd_extract.h"
#include "cmd_windows.h"
#include "converters/converters.h"
#include "converters/aln_writer.h"
#include "bam2bpp/bam2bpp.h"
#include "bam2bpp/vcf_phase.h"

#define BPP_SEQS_VERSION "0.1.6"

/* ────────────────────────────────────────────────────────────────────────
 * CLI options
 * ───────────────────────────────────────────────────────────────────── */

typedef struct {
    char **inputs;          /* positional files + values from --imap */
    int    n_inputs;

    char  *out_prefix;
    int    dry_run;
    int    json;
    int    json_indent;
    int    quiet;

    char  *imap_path_arg;   /* from --imap (may overlap with positional) */
    char  *reference_path;  /* from --reference: force this file to BS_FASTA_REFERENCE */

    /* Filtering / conversion options (passed to bam2bpp / converters) */
    int    min_bq;
    int    min_mq;
    int    min_dp;
    double het_freq;
    int    min_length;
    double max_missing;
    int    min_snps;
    int    keep_invariant;
    Phasing phasing;
    char  *phased_vcf;
} CLI;

static void print_usage(FILE *fp, const char *prog)
{
    fprintf(fp,
"bpp-seqs %s — convert sequence data files to BPP input format.\n"
"\n"
"Usage: %s [options] file1 file2 ...\n"
"\n"
"Files are provided in any order; type is detected from content.\n"
"With --dry-run or when an Imap is missing, only inspection is performed.\n"
"\n"
"Subcommands (run `%s VERB --help` for their options):\n"
"  extract INPUT.txt     Write a new BPP sequence file holding a subset of\n"
"                        an existing one's loci (by name, index, chrom, size).\n"
"  windows INPUT         Tile a genome into fixed-size candidate loci and\n"
"                        write them as a BED for the conversion flow below.\n"
"\n"
"General:\n"
"  --out PREFIX          Output file prefix for BPP files (required for conversion)\n"
"  --dry-run             Inspect files only; do not convert\n"
"  --json                Output JSON instead of human-readable text\n"
"  --json-indent N       JSON indentation width [2]\n"
"  --quiet               Suppress stderr progress messages\n"
"  -h, --help            Show help and exit\n"
"  --version             Show version and exit\n"
"\n"
"Input:\n"
"  --imap FILE           Imap file (sample → population mapping)\n"
"  --reference FILE      Designate FILE as the reference FASTA (skips the\n"
"                        content-based REFERENCE vs CONTIGS classifier).\n"
"\n"
"Phasing (how diploid calls become sequences; BAM/CRAM and gVCF only):\n"
"  --phasing MODE        iupac (default), split, haploid, vcf\n"
"                          iupac   one sequence per individual, hets as\n"
"                                  IUPAC ambiguity codes\n"
"                          split   two unphased haplotypes per individual\n"
"                          haploid one major-allele sequence per individual\n"
"                          vcf     apply true phase from --phased-vcf\n"
"  --phased-vcf FILE     Phased VCF for --phasing vcf mode\n"
"\n"
"Base/read filtering (BAM/CRAM pileups only; ignored for other inputs):\n"
"  --min-bq INT          Minimum base quality [20]\n"
"  --min-mq INT          Minimum mapping quality [20]\n"
"  --min-dp INT          Minimum depth to call a base [5]\n"
"  --het-freq FLOAT      Minor-allele frequency threshold for het call [0.20]\n"
"\n"
"Locus filtering (all workflows):\n"
"  --min-length INT      Minimum locus length in bp [100]\n"
"  --max-missing FLOAT   Maximum fraction of N per locus [0.5]\n"
"  --min-snps INT        Minimum segregating sites per locus [1]\n"
"  --keep-invariant      Keep loci with no variation\n"
"\n"
"Outputs (written when --out PREFIX is given and nothing is missing):\n"
"  PREFIX.txt            BPP sequence file (per-locus alignments)\n"
"  PREFIX.imap           the Imap actually used\n"
"  PREFIX.stats.tsv      per-locus statistics\n"
"  PREFIX.loci.tsv       the locus table\n"
"\n"
"Exit status: 0 on success (including an inspection reporting missing items),\n"
"1 on a conversion or system error. Long options must match exactly;\n"
"abbreviations are rejected.\n",
        BPP_SEQS_VERSION, prog, prog);
}

enum {
    OPT_OUT = 256,
    OPT_DRY_RUN,
    OPT_JSON,
    OPT_JSON_INDENT,
    OPT_QUIET,
    OPT_VERSION,
    OPT_IMAP,
    OPT_REFERENCE,
    OPT_PHASING,
    OPT_PHASED_VCF,
    OPT_MIN_BQ,
    OPT_MIN_MQ,
    OPT_MIN_DP,
    OPT_HET_FREQ,
    OPT_MIN_LENGTH,
    OPT_MAX_MISSING,
    OPT_MIN_SNPS,
    OPT_KEEP_INV
};

static const struct option LONG_OPTS[] = {
    {"out",            required_argument, NULL, OPT_OUT},
    {"dry-run",        no_argument,       NULL, OPT_DRY_RUN},
    {"json",           no_argument,       NULL, OPT_JSON},
    {"json-indent",    required_argument, NULL, OPT_JSON_INDENT},
    {"quiet",          no_argument,       NULL, OPT_QUIET},
    {"version",        no_argument,       NULL, OPT_VERSION},
    {"imap",           required_argument, NULL, OPT_IMAP},
    {"reference",      required_argument, NULL, OPT_REFERENCE},
    {"phasing",        required_argument, NULL, OPT_PHASING},
    {"phased-vcf",     required_argument, NULL, OPT_PHASED_VCF},
    {"min-bq",         required_argument, NULL, OPT_MIN_BQ},
    {"min-mq",         required_argument, NULL, OPT_MIN_MQ},
    {"min-dp",         required_argument, NULL, OPT_MIN_DP},
    {"het-freq",       required_argument, NULL, OPT_HET_FREQ},
    {"min-length",     required_argument, NULL, OPT_MIN_LENGTH},
    {"max-missing",    required_argument, NULL, OPT_MAX_MISSING},
    {"min-snps",       required_argument, NULL, OPT_MIN_SNPS},
    {"keep-invariant", no_argument,       NULL, OPT_KEEP_INV},
    {"help",           no_argument,       NULL, 'h'},
    {NULL, 0, NULL, 0}
};

static void cli_init(CLI *c)
{
    memset(c, 0, sizeof(*c));
    c->json_indent  = 2;
    c->min_bq       = DEFAULT_MIN_BQ;
    c->min_mq       = DEFAULT_MIN_MQ;
    c->min_dp       = DEFAULT_MIN_DP;
    c->het_freq     = DEFAULT_HET_FREQ;
    c->min_length   = DEFAULT_MIN_LEN;
    c->max_missing  = DEFAULT_MAX_MISS;
    c->min_snps     = DEFAULT_MIN_SNPS;
    c->phasing      = PHASE_IUPAC;
}

static void cli_free(CLI *c)
{
    for (int i = 0; i < c->n_inputs; i++) free(c->inputs[i]);
    free(c->inputs);
    free(c->out_prefix);
    free(c->imap_path_arg);
    free(c->reference_path);
    free(c->phased_vcf);
}

static void cli_push_input(CLI *c, const char *path)
{
    /* dedupe by path */
    for (int i = 0; i < c->n_inputs; i++) {
        if (strcmp(c->inputs[i], path) == 0) return;
    }
    c->inputs = (char **)realloc(c->inputs, sizeof(char *) * (size_t)(c->n_inputs + 1));
    c->inputs[c->n_inputs++] = strdup(path);
}

/* getopt_long silently accepts unambiguous ABBREVIATIONS of long options
 * (e.g. '--phase' binds to '--phased-vcf'), so a mistyped or wrong flag can be
 * misinterpreted instead of rejected. Require every '--name[=value]' token to
 * match a known option EXACTLY; a bare '--' (end-of-options) is left alone. */
static int reject_option_abbreviations(int argc, char **argv)
{
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if (a[0] != '-' || a[1] != '-' || a[2] == '\0') continue;
        size_t len = strcspn(a + 2, "=");
        int exact = 0;
        for (const struct option *o = LONG_OPTS; o->name; o++) {
            if (strlen(o->name) == len && strncmp(a + 2, o->name, len) == 0) { exact = 1; break; }
        }
        if (!exact) {
            fprintf(stderr, "Error: unknown option '--%.*s'. Run with --help for usage.\n",
                    (int)len, a + 2);
            return -1;
        }
    }
    return 0;
}

static int parse_cli(int argc, char **argv, CLI *c)
{
    int opt;
    if (reject_option_abbreviations(argc, argv) != 0) return -1;
    while ((opt = getopt_long(argc, argv, "h", LONG_OPTS, NULL)) != -1) {
        switch (opt) {
            case OPT_OUT:        c->out_prefix = strdup(optarg); break;
            case OPT_DRY_RUN:    c->dry_run = 1; break;
            case OPT_JSON:       c->json = 1; break;
            case OPT_JSON_INDENT: c->json_indent = atoi(optarg); break;
            case OPT_QUIET:      c->quiet = 1; break;
            case OPT_VERSION:    printf("bpp-seqs %s\n", BPP_SEQS_VERSION); exit(0);
            case OPT_IMAP:
                free(c->imap_path_arg);
                c->imap_path_arg = strdup(optarg);
                cli_push_input(c, optarg);
                break;
            case OPT_REFERENCE:
                free(c->reference_path);
                c->reference_path = strdup(optarg);
                cli_push_input(c, optarg);
                break;
            case OPT_PHASING:
                if      (strcmp(optarg, "iupac")   == 0) c->phasing = PHASE_IUPAC;
                else if (strcmp(optarg, "split")   == 0) c->phasing = PHASE_SPLIT;
                else if (strcmp(optarg, "haploid") == 0) c->phasing = PHASE_HAPLOID;
                else if (strcmp(optarg, "vcf")     == 0) c->phasing = PHASE_VCF;
                else { fprintf(stderr, "Error: unknown --phasing '%s'\n", optarg); return -1; }
                break;
            case OPT_PHASED_VCF: c->phased_vcf  = strdup(optarg); break;
            case OPT_MIN_BQ:     c->min_bq      = atoi(optarg); break;
            case OPT_MIN_MQ:     c->min_mq      = atoi(optarg); break;
            case OPT_MIN_DP:     c->min_dp      = atoi(optarg); break;
            case OPT_HET_FREQ:   c->het_freq    = atof(optarg); break;
            case OPT_MIN_LENGTH: c->min_length  = atoi(optarg); break;
            case OPT_MAX_MISSING: c->max_missing = atof(optarg); break;
            case OPT_MIN_SNPS:   c->min_snps    = atoi(optarg); break;
            case OPT_KEEP_INV:   c->keep_invariant = 1; break;
            case 'h':            print_usage(stdout, argv[0]); exit(0);
            default:             return -1;
        }
    }
    for (int i = optind; i < argc; i++) cli_push_input(c, argv[i]);
    if (c->n_inputs == 0) {
        fprintf(stderr, "Error: no input files. Run with --help for usage.\n");
        return -1;
    }
    return 0;
}

/* Build a "recommended command" string for status=ready. */
static char *build_recommended_command(const CLI *c)
{
    /* "bpp-seqs --out mydata <files>" */
    const char *prefix = c->out_prefix ? c->out_prefix : "mydata";
    size_t cap = 256 + strlen(prefix);
    for (int i = 0; i < c->n_inputs; i++) cap += strlen(c->inputs[i]) + 1;
    char *s = (char *)malloc(cap);
    if (!s) return NULL;
    int n = snprintf(s, cap, "bpp-seqs --out %s", prefix);
    for (int i = 0; i < c->n_inputs; i++) {
        n += snprintf(s + n, cap - n, " %s", c->inputs[i]);
    }
    return s;
}

static FileInfo *find_by_type(FileInfo **files, int n, FileType t)
{
    for (int i = 0; i < n; i++) if (files[i]->ft == t) return files[i];
    return NULL;
}

static const char *human_skip_reason(const char *r)
{
    if (!r) return NULL;
    return r;
}

/* ────────────────────────────────────────────────────────────────────────
 * bam2bpp conversion (mirrors logic from the bam2bpp standalone main.c)
 * ───────────────────────────────────────────────────────────────────── */

static int run_bam2bpp(const CLI *c, FileInfo **files, int n_files,
                       ConversionResult *cr)
{
    FileInfo *ref_fi  = find_by_type(files, n_files, BS_FASTA_REFERENCE);
    FileInfo *bed_fi  = find_by_type(files, n_files, BS_BED);
    FileInfo *imap_fi = find_by_type(files, n_files, BS_IMAP);
    if (!ref_fi || !bed_fi || !imap_fi) {
        fprintf(stderr, "Error: bam2bpp requires a reference FASTA, BED, and Imap.\n");
        return 1;
    }

    /* Collect BAM paths */
    char **bam_paths = NULL;
    int n_bams = 0;
    for (int i = 0; i < n_files; i++) {
        if (files[i]->ft == BS_BAM || files[i]->ft == BS_CRAM) {
            bam_paths = realloc(bam_paths, sizeof(char *) * (size_t)(n_bams + 1));
            bam_paths[n_bams++] = strdup(files[i]->path);
        }
    }
    if (n_bams == 0) { free(bam_paths); return 1; }

    /* Build a bam2bpp Args */
    Args a;
    memset(&a, 0, sizeof(a));
    a.bam_paths     = bam_paths;
    a.n_bams        = n_bams;
    a.ref_path      = strdup(ref_fi->path);
    a.bed_path      = strdup(bed_fi->path);
    a.imap_path     = strdup(imap_fi->path);
    a.out_prefix    = strdup(c->out_prefix);
    a.min_bq        = c->min_bq;
    a.min_mq        = c->min_mq;
    a.min_dp        = c->min_dp;
    a.het_freq      = c->het_freq;
    a.phasing       = c->phasing;
    a.phased_vcf    = c->phased_vcf ? strdup(c->phased_vcf) : NULL;
    a.min_length    = c->min_length;
    a.max_missing   = c->max_missing;
    a.min_snps      = c->min_snps;
    a.keep_invariant = c->keep_invariant;
    a.max_depth     = DEFAULT_MAX_DEPTH;
    a.unphased_policy = UNPHASED_MISSING;

    /* --- Open inputs (mirrors bam2bpp/main.c) --- */
    faidx_t *fai = fai_load(a.ref_path);
    if (!fai) {
        fprintf(stderr, "Error: cannot load FASTA index for '%s'\n", a.ref_path);
        goto fail;
    }
    int n_loci;
    Locus *loci = load_bed(a.bed_path, &n_loci);
    if (!loci || n_loci == 0) {
        fprintf(stderr, "Error: no usable loci in '%s'\n", a.bed_path);
        fai_destroy(fai); goto fail;
    }
    int n_imap;
    ImapEntry *imap = load_imap(a.imap_path, &n_imap);
    if (!imap) { fai_destroy(fai); goto fail; }

    BamFile **bams = open_bams(a.bam_paths, a.n_bams, a.max_depth);
    if (!bams) { fai_destroy(fai); goto fail; }

    /* Initial per-sample phasing mirrors the global mode; the VCF
     * classification below may override per sample so phased samples emit
     * 2 haplotypes while unphased samples emit 1 IUPAC sequence. */
    for (int i = 0; i < a.n_bams; i++) bams[i]->sample_phasing = a.phasing;

    VcfPhase *vcf = NULL;
    if (a.phasing == PHASE_VCF) {
        vcf = open_vcf_phase(a.phased_vcf);
        if (!vcf) { close_bams(bams, a.n_bams); fai_destroy(fai); goto fail; }

        /* Per-sample classification: PHASE_VCF for fully-phased samples
         * present in the VCF, PHASE_IUPAC for everything else (samples not
         * in the VCF, or with any unphased heterozygote in the locus set).
         * This lets a single run mix phased reference panels with unphased
         * single-individual samples (e.g. ancient genomes). */
        VcfPhaseSampleStat *ps = vcf_phase_classify(vcf, bams, a.n_bams,
                                                    loci, n_loci);
        int n_phased = 0, n_iupac = 0;
        for (int i = 0; i < a.n_bams; i++) {
            if (ps[i].not_in_vcf || ps[i].n_unphased > 0) {
                bams[i]->sample_phasing = PHASE_IUPAC;
                n_iupac++;
            } else {
                bams[i]->sample_phasing = PHASE_VCF;
                n_phased++;
            }
        }
        if (!c->quiet) {
            fprintf(stderr,
                "  --phasing vcf: %d sample%s fully phased (split into "
                "haplotypes), %d sample%s unphased (IUPAC).\n",
                n_phased, n_phased == 1 ? "" : "s",
                n_iupac,  n_iupac  == 1 ? "" : "s");
            for (int i = 0; i < a.n_bams; i++) {
                if (ps[i].not_in_vcf) {
                    fprintf(stderr, "         %-20s NOT in VCF → IUPAC\n",
                            ps[i].sample);
                } else if (ps[i].n_unphased > 0) {
                    double frac = ps[i].n_het ?
                        100.0 * ps[i].n_unphased / ps[i].n_het : 0.0;
                    fprintf(stderr,
                        "         %-20s %d/%d unphased (%.1f%%) → IUPAC\n",
                        ps[i].sample, ps[i].n_unphased, ps[i].n_het, frac);
                }
            }
            if (n_iupac > 0)
                fprintf(stderr,
                    "         BPP control file: phase = 0 for phased "
                    "species, 1 for unphased.\n");
        }
        vcf_phase_stats_free(ps, a.n_bams);

        /* If nothing got classified as PHASE_VCF, close the VCF — pass 2
         * is not needed; pass 1 still does the right thing per sample. */
        if (n_phased == 0) {
            close_vcf_phase(vcf);
            vcf = NULL;
            a.phasing = PHASE_IUPAC;
        }
    }

    bam2bpp_writer_set_quiet(c->quiet);

    if (!c->quiet) {
        fprintf(stderr, "bpp-seqs: %d loci, %d samples\n", n_loci, a.n_bams);
    }

    LocusResult *results = calloc((size_t)n_loci, sizeof(LocusResult));
    LocusStat   *stats   = calloc((size_t)n_loci, sizeof(LocusStat));
    int n_results = 0;
    int n_too_short = 0, n_high_missing = 0, n_insufficient_snps = 0;

    for (int li = 0; li < n_loci; li++) {
        Locus *loc = &loci[li];
        stats[li].locus_name = strdup(loc->name);
        stats[li].locus_len  = (int)(loc->end - loc->start);

        if ((int)(loc->end - loc->start) < a.min_length) {
            stats[li].skip_reason = "too_short";
            n_too_short++;
            continue;
        }

        char region[1024];
        snprintf(region, sizeof(region),
                 "%s:%" PRId32 "-%" PRId32,
                 loc->chrom, loc->start + 1, loc->end);
        int ref_len = 0;
        char *ref_seq = fai_fetch(fai, region, &ref_len);
        if (!ref_seq || ref_len != (int)(loc->end - loc->start)) {
            fprintf(stderr, "Warning: cannot fetch reference for '%s'; skipping\n", loc->name);
            free(ref_seq);
            stats[li].skip_reason = "ref_error";
            continue;
        }
        for (int j = 0; j < ref_len; j++) ref_seq[j] = (char)toupper((unsigned char)ref_seq[j]);

        LocusResult result; memset(&result, 0, sizeof(result));
        double mean_dp = 0.0;
        int rc;
        if (a.phasing == PHASE_VCF) {
            rc = process_locus_vcf(bams, a.n_bams, vcf, loc, ref_seq, &a, &result, &mean_dp);
        } else {
            rc = process_locus(bams, a.n_bams, loc, &a, &result, &mean_dp);
        }
        free(ref_seq);
        if (rc != 0) { stats[li].skip_reason = "pileup_error"; continue; }

        double miss_frac = compute_missing_frac(result.seqs, result.n_seqs, result.locus_len);
        int    n_snps    = count_snps          (result.seqs, result.n_seqs, result.locus_len);
        stats[li].missing_frac = miss_frac;
        stats[li].mean_depth   = mean_dp;
        stats[li].n_snps       = n_snps;

        if (miss_frac > a.max_missing) {
            stats[li].skip_reason = "high_missing"; n_high_missing++;
            free_locus_result(&result); continue;
        }
        if (!a.keep_invariant && n_snps < a.min_snps) {
            stats[li].skip_reason = "insufficient_snps"; n_insufficient_snps++;
            free_locus_result(&result); continue;
        }

        /* For SPLIT/VCF phasing, internal seq_names use '<sample>^<n>' to
         * mark haplotypes.  The bam2bpp writer later transforms '^' to '_'
         * so the per-locus BPP ids are unique across samples.  Apply the
         * same transform here so the sanity check sees the same ids the
         * writer will emit (otherwise it would flag every locus as having
         * duplicate id '1' / '2'). */
        if (a.phasing == PHASE_SPLIT || a.phasing == PHASE_VCF) {
            for (int s = 0; s < result.n_seqs; s++) {
                for (char *p = result.seq_names[s]; *p; p++)
                    if (*p == '^') *p = '_';
            }
        }

        /* Structural sanity check (lengths, ids, char set) just before
         * accepting this locus for output. */
        {
            SanityReport sr; sanity_report_init(&sr);
            int ok = sanity_check_locus(result.name,
                                        result.seq_names, result.seqs,
                                        result.n_seqs, result.locus_len,
                                        &sr, stderr);
            if (ok != 0) {
                stats[li].skip_reason = "sanity_failed";
                free_locus_result(&result);
                continue;
            }
        }

        results[n_results++] = result;
        if (!c->quiet && ((li + 1) % 500 == 0 || li == n_loci - 1)) {
            fprintf(stderr, "  %d / %d loci processed (%d passed)\n", li + 1, n_loci, n_results);
        }
    }

    if (n_results == 0) {
        fprintf(stderr, "Error: no loci passed QC.\n");
        goto fail2;
    }

    /* Per-locus sequence count comes from results[].n_seqs (locus-level),
     * but for the summary we report the maximum across loci — which equals
     * the sum of per-sample emission counts (1 or 2 each). */
    int n_seqs = 0;
    for (int i = 0; i < a.n_bams; i++) {
        n_seqs += (bams[i]->sample_phasing == PHASE_SPLIT ||
                   bams[i]->sample_phasing == PHASE_VCF) ? 2 : 1;
    }

    write_bpp_file  (a.out_prefix, results, n_results, n_seqs);
    write_imap_file (a.out_prefix, bams, a.n_bams, imap, n_imap, a.phasing);
    write_stats_file(a.out_prefix, stats, n_loci);

    /* Imap coverage check: every individual id present in the BAM-derived
     * sequence output should have an Imap row, and every Imap row should
     * be used by some sequence. */
    {
        char **ids = NULL; int n_ids = 0;
        for (int i = 0; i < n_results; i++) {
            for (int j = 0; j < results[i].n_seqs; j++) {
                char *id = sanity_individual_id(results[i].seq_names[j]);
                if (!id) continue;
                /* For SPLIT and VCF phasing, seq_names are like "ind1^1";
                 * the post-caret id "1" won't match the Imap. Use the full
                 * sample name in that case, transformed to "ind1_1" to
                 * mirror the writer. */
                if (a.phasing == PHASE_SPLIT || a.phasing == PHASE_VCF) {
                    free(id);
                    const char *raw = results[i].seq_names[j];
                    id = strdup(raw);
                    for (char *p = id; *p; p++) if (*p == '^') *p = '_';
                }
                int dup = 0;
                for (int k = 0; k < n_ids; k++) if (strcmp(ids[k], id) == 0) { dup = 1; break; }
                if (dup) { free(id); continue; }
                ids = (char **)realloc(ids, sizeof(char *) * (size_t)(n_ids + 1));
                ids[n_ids++] = id;
            }
        }
        /* Build the Imap sample list as it will appear in the output file.
         * Per-sample: PHASE_SPLIT and PHASE_VCF samples contribute
         * <sample>_1 and <sample>_2; everything else contributes <sample>. */
        int n_imap_eff = 0;
        for (int i = 0; i < a.n_bams; i++) {
            n_imap_eff += (bams[i]->sample_phasing == PHASE_SPLIT ||
                           bams[i]->sample_phasing == PHASE_VCF) ? 2 : 1;
        }
        char **imap_eff = (char **)calloc((size_t)n_imap_eff, sizeof(char *));
        int k = 0;
        for (int i = 0; i < a.n_bams; i++) {
            if (bams[i]->sample_phasing == PHASE_SPLIT ||
                bams[i]->sample_phasing == PHASE_VCF) {
                char buf[256];
                snprintf(buf, sizeof(buf), "%s_1", bams[i]->sample);
                imap_eff[k++] = strdup(buf);
                snprintf(buf, sizeof(buf), "%s_2", bams[i]->sample);
                imap_eff[k++] = strdup(buf);
            } else {
                imap_eff[k++] = strdup(bams[i]->sample);
            }
        }
        sanity_check_imap_coverage(ids, n_ids, imap_eff, n_imap_eff, stderr);
        for (int i = 0; i < n_ids; i++) free(ids[i]);
        free(ids);
        for (int i = 0; i < n_imap_eff; i++) free(imap_eff[i]);
        free(imap_eff);
    }

    /* Populate ConversionResult */
    {
        char p[1024];
        snprintf(p, sizeof(p), "%s.txt",        a.out_prefix); cr->out_sequences = strdup(p);
        snprintf(p, sizeof(p), "%s.imap",       a.out_prefix); cr->out_imap      = strdup(p);
        snprintf(p, sizeof(p), "%s.stats.tsv",  a.out_prefix); cr->out_stats     = strdup(p);
    }
    cr->n_loci_input        = n_loci;
    cr->n_loci_passed       = n_results;
    cr->n_loci_failed       = n_loci - n_results;
    cr->n_sequences         = n_seqs;
    cr->n_too_short         = n_too_short;
    cr->n_high_missing      = n_high_missing;
    cr->n_insufficient_snps = n_insufficient_snps;
    cr->has_results         = 1;

    /* Recommended `phase = ...` line — one digit per Imap species.
     *   PHASE_VCF   → 0 (data already phased; BPP must not re-phase)
     *   PHASE_SPLIT → 0 (two arbitrary-phase sequences per sample)
     *   PHASE_HAPLOID → 0 (one allele per sample, no phasing needed)
     *   PHASE_IUPAC → 1 (each IUPAC het = unphased het, BPP must phase)
     */
    {
        /* count distinct populations in the Imap */
        int n_species = 0;
        char **seen = NULL;
        for (int i = 0; i < n_imap; i++) {
            int dup = 0;
            for (int j = 0; j < n_species; j++)
                if (strcmp(seen[j], imap[i].population) == 0) { dup = 1; break; }
            if (!dup) {
                seen = realloc(seen, sizeof(char *) * (size_t)(n_species + 1));
                seen[n_species++] = imap[i].population;
            }
        }
        if (n_species > 0) {
            /* Per-species flag.  A species is "phased" (flag=0) only if every
             * BAM sample mapping to that population is PHASE_VCF / SPLIT /
             * HAPLOID.  If any sample is PHASE_IUPAC, the species needs BPP
             * to do phase sampling → flag=1. */
            size_t cap = (size_t)n_species * 2 + 1;
            char *buf = (char *)malloc(cap);
            int off = 0;
            for (int i = 0; i < n_species; i++) {
                int flag = 0;   /* default: phased */
                int saw_sample = 0;
                for (int b = 0; b < a.n_bams; b++) {
                    const char *bpop = imap_lookup(imap, n_imap, bams[b]->sample);
                    if (!bpop || strcmp(bpop, seen[i]) != 0) continue;
                    saw_sample = 1;
                    if (bams[b]->sample_phasing == PHASE_IUPAC) {
                        flag = 1;
                        break;
                    }
                }
                if (!saw_sample) flag = 1;  /* unreached population — assume unphased */
                off += snprintf(buf + off, cap - (size_t)off,
                                "%s%d", i ? " " : "", flag);
            }
            cr->recommended_phase = buf;
        }
        free(seen);
    }
    for (int li = 0; li < n_loci; li++) {
        const char *status = stats[li].skip_reason ? "failed" : "passed";
        conversion_result_add_locus(cr,
            stats[li].locus_name, stats[li].locus_len, stats[li].n_snps,
            stats[li].missing_frac, stats[li].mean_depth,
            status, human_skip_reason(stats[li].skip_reason));
        conversion_result_set_locus_source(cr,
            "BED", a.bed_path, loci[li].chrom,
            loci[li].start + 1, loci[li].end, 1);
    }

    /* Write <prefix>.loci.tsv for the bam2bpp path (passing loci only). */
    {
        LocusProv *items = (LocusProv *)calloc((size_t)n_results, sizeof(LocusProv));
        int k = 0;
        for (int li = 0; li < n_loci && k < n_results; li++) {
            if (stats[li].skip_reason) continue;
            items[k].name   = loci[li].name;
            items[k].kind   = "BED";
            items[k].file   = a.bed_path;
            items[k].chrom  = loci[li].chrom;
            items[k].start  = loci[li].start + 1;
            items[k].end    = loci[li].end;
            items[k].stride = 1;
            items[k].length = (int)(loci[li].end - loci[li].start);
            items[k].n_seqs = n_seqs;
            k++;
        }
        if (write_loci_tsv(a.out_prefix, items, n_results) == 0) {
            char p[1024];
            snprintf(p, sizeof(p), "%s.loci.tsv", a.out_prefix);
            cr->out_loci = strdup(p);
        }
        free(items);
    }

    /* Cleanup */
    for (int i = 0; i < n_results; i++) free_locus_result(&results[i]);
    free(results);
    for (int i = 0; i < n_loci; i++) free(stats[i].locus_name);
    free(stats);
    close_bams(bams, a.n_bams);
    if (vcf) close_vcf_phase(vcf);
    free_loci(loci, n_loci);
    free_imap(imap, n_imap);
    fai_destroy(fai);

    for (int i = 0; i < a.n_bams; i++) free(a.bam_paths[i]);
    free(a.bam_paths);
    free(a.ref_path); free(a.bed_path); free(a.imap_path);
    free(a.out_prefix); free(a.phased_vcf);
    return 0;

fail2:
    for (int i = 0; i < n_results; i++) free_locus_result(&results[i]);
    free(results);
    for (int i = 0; i < n_loci; i++) free(stats[i].locus_name);
    free(stats);
    close_bams(bams, a.n_bams);
    if (vcf) close_vcf_phase(vcf);
    free_loci(loci, n_loci);
    free_imap(imap, n_imap);
    fai_destroy(fai);
fail:
    for (int i = 0; i < a.n_bams; i++) free(a.bam_paths[i]);
    free(a.bam_paths);
    free(a.ref_path); free(a.bed_path); free(a.imap_path);
    free(a.out_prefix); free(a.phased_vcf);
    return 1;
}

/* ────────────────────────────────────────────────────────────────────────
 * main
 * ───────────────────────────────────────────────────────────────────── */

int main(int argc, char **argv)
{
    if (argc < 2) { print_usage(stderr, argv[0]); return 1; }

    /* Subcommand dispatch: if argv[1] is a recognised verb, hand off
     * with argv[1..] (the verb becomes argv[0] for the subcommand). */
    if (strcmp(argv[1], "extract") == 0) {
        return cmd_extract(argc - 1, argv + 1);
    }
    if (strcmp(argv[1], "windows") == 0) {
        return cmd_windows(argc - 1, argv + 1);
    }
    /* No verb matched → fall through to the existing inspect/convert flow. */

    CLI cli;
    cli_init(&cli);
    if (parse_cli(argc, argv, &cli) != 0) { cli_free(&cli); return 1; }

    /* --- Inspect every input --- */
    FileInfo **files = (FileInfo **)calloc((size_t)cli.n_inputs, sizeof(FileInfo *));
    for (int i = 0; i < cli.n_inputs; i++) {
        files[i] = inspect_file(cli.inputs[i]);
        if (!files[i]) {
            fprintf(stderr, "Error: cannot open or read '%s'\n", cli.inputs[i]);
            for (int j = 0; j < i; j++) file_info_free(files[j]);
            free(files); cli_free(&cli);
            return 1;
        }
    }

    /* --- Apply --reference override (force a multi-seq FASTA to REFERENCE) --- */
    if (cli.reference_path) {
        FileInfo *ref_fi = NULL;
        for (int i = 0; i < cli.n_inputs; i++) {
            if (strcmp(cli.inputs[i], cli.reference_path) == 0) { ref_fi = files[i]; break; }
        }
        if (!ref_fi) {
            fprintf(stderr, "Error: --reference '%s' was not inspected (not in input list?)\n",
                    cli.reference_path);
        } else if (file_info_force_reference(ref_fi) != 0) {
            fprintf(stderr, "Error: --reference '%s' is not a FASTA (detected: %s)\n",
                    cli.reference_path, file_type_name(ref_fi->ft));
        }
    }

    /* --- Cross-validate --- */
    CrossValidation *cv = cross_validate(files, cli.n_inputs);

    /* --- Decide workflow --- */
    WorkflowDecision *d = workflow_decide(files, cli.n_inputs);

    /* --- Conversion (if not dry-run and ready) --- */
    ConversionResult cr; conversion_result_init(&cr);
    int conversion_attempted = 0;
    int conversion_rc = 0;
    if (!cli.dry_run && d->ready_to_run) {
        if (!cli.out_prefix) {
            fprintf(stderr, "Error: --out PREFIX is required for conversion.\n");
            conversion_rc = 1;
        } else {
            FileInfo *imap_fi = find_by_type(files, cli.n_inputs, BS_IMAP);
            ConvertOpts copts = {
                .min_length    = cli.min_length,
                .max_missing   = cli.max_missing,
                .min_snps      = cli.min_snps,
                .keep_invariant = cli.keep_invariant,
                .quiet         = cli.quiet,
                .phasing       = cli.phasing
            };
            switch (d->workflow) {
                case WF_BAM2BPP:
                    conversion_attempted = 1;
                    conversion_rc = run_bam2bpp(&cli, files, cli.n_inputs, &cr);
                    break;
                case WF_FASTA2BPP:
                    conversion_attempted = 1;
                    conversion_rc = convert_fasta_msa(files, cli.n_inputs, imap_fi,
                                                      cli.out_prefix, &copts, &cr);
                    break;
                case WF_PHYLIP2BPP:
                    conversion_attempted = 1;
                    conversion_rc = convert_phylip(files, cli.n_inputs, imap_fi,
                                                   cli.out_prefix, &copts, &cr);
                    break;
                case WF_NEXUS2BPP:
                    conversion_attempted = 1;
                    conversion_rc = convert_nexus(files, cli.n_inputs, imap_fi,
                                                  cli.out_prefix, &copts, &cr);
                    break;
                case WF_GVCF2BPP:
                    conversion_attempted = 1;
                    conversion_rc = convert_gvcf(files, cli.n_inputs, imap_fi,
                                                 cli.out_prefix, &copts, &cr);
                    break;
                default:
                    break;
            }
        }
    }

    /* --- Emit --- */
    char *rec = NULL;
    if (d->ready_to_run && !cr.has_results) {
        rec = build_recommended_command(&cli);
    }

    if (cli.json) {
        print_json(stdout, cli.json_indent,
                   files, cli.n_inputs, cv, d, &cr,
                   rec, BPP_SEQS_VERSION);
    } else {
        print_human(stdout, files, cli.n_inputs, cv, d, &cr, BPP_SEQS_VERSION);
    }

    /* --- Cleanup --- */
    free(rec);
    conversion_result_free(&cr);
    workflow_decision_free(d);
    cross_validation_free(cv);
    for (int i = 0; i < cli.n_inputs; i++) file_info_free(files[i]);
    free(files);
    cli_free(&cli);

    /* Exit code policy:
     *   0 — success (inspect-only, ready, or completed conversion)
     *   1 — conversion or system error
     * "incomplete" (missing Imap) is not an error: the agent uses it. */
    if (conversion_attempted && conversion_rc != 0) return 1;
    return 0;
}
