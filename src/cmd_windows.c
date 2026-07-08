#define _POSIX_C_SOURCE 200809L
#define _DARWIN_C_SOURCE 1   /* for strcasestr on macOS */
#define _GNU_SOURCE      1   /* for strcasestr on glibc/Linux */

#include "cmd_windows.h"
#include "inspect.h"
#include "json_writer.h"

#include <ctype.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sys/stat.h>

#include <htslib/faidx.h>

static void print_help(FILE *fp)
{
    fprintf(fp,
"Usage: bpp-seqs windows INPUT --window-size W [filters] --out OUT.bed\n"
"\n"
"Generate a BED file of candidate BPP loci by tiling a genome with\n"
"fixed-size windows.  INPUT may be a FASTA reference, a BAM/CRAM, or a\n"
"VCF/gVCF — anything that contains chromosome names and lengths.\n"
"\n"
"Required:\n"
"  INPUT                    reference FASTA (or BAM/CRAM or VCF/gVCF)\n"
"  --window-size W          locus size in bp\n"
"  --out FILE               output BED path (use \"-\" for stdout)\n"
"\n"
"Filters (applied in this order):\n"
"  --step S                 stride between window starts (default: window-size,\n"
"                           giving non-overlapping windows)\n"
"  --include-chrom NAME[,...]\n"
"                           only consider these chromosomes\n"
"  --exclude-chrom NAME[,...]\n"
"                           skip these chromosomes\n"
"  --autosomes-only         skip sex chromosomes, mitochondria, unplaced /\n"
"                           random / alt / decoy contigs by name heuristic\n"
"  --skip-edges N           drop the first and last N bp of every chromosome\n"
"  --exclude-regions FILE   subtract intervals in FILE.bed before windowing\n"
"  --min-spacing G          require >= G bp between consecutive kept windows\n"
"                           on the same chromosome\n"
"  --n-loci N --seed S      randomly sample N windows from the survivors\n"
"\n"
"Other:\n"
"  --json                   emit a JSON summary on stdout in addition to BED\n"
"  -h, --help               this message\n"
"\n"
"Example: 500 ~500-bp loci, ≥10 kb apart, on autosomes, edges skipped\n"
"  bpp-seqs windows ref.fa --window-size 500 --min-spacing 10000 \\\n"
"                  --n-loci 500 --autosomes-only --skip-edges 100000 \\\n"
"                  --seed 42 --out loci.bed\n");
}

/* ────────────────────────────────────────────────────────────────────────
 * Per-chromosome length list
 * ───────────────────────────────────────────────────────────────────── */

typedef struct {
    char    *name;
    int64_t  length;
} Chrom;

static void free_chroms(Chrom *c, int n)
{
    if (!c) return;
    for (int i = 0; i < n; i++) free(c[i].name);
    free(c);
}

/* Get chrom names + lengths from a FASTA (.fai) using htslib. */
static Chrom *chroms_from_fasta(const char *path, int *n_out)
{
    faidx_t *fai = fai_load(path);
    if (!fai) return NULL;
    int n = faidx_nseq(fai);
    Chrom *c = (Chrom *)calloc((size_t)n, sizeof(Chrom));
    for (int i = 0; i < n; i++) {
        const char *nm = faidx_iseq(fai, i);
        c[i].name   = strdup(nm ? nm : "");
        c[i].length = faidx_seq_len(fai, nm);
    }
    fai_destroy(fai);
    *n_out = n;
    return c;
}

/* From a FileInfo populated by inspect_file.  For any FASTA variant we
 * use the .fai (creating it if missing) so the chrom list is exact and
 * fast.  For BAM/CRAM we use @SQ; for VCF/gVCF we use ##contig. */
static Chrom *chroms_from_file_info(const FileInfo *fi, int *n_out)
{
    if (!fi) return NULL;
    if (fi->ft == BS_FASTA_REFERENCE || fi->ft == BS_FASTA_MSA ||
        fi->ft == BS_FASTA_CONTIGS) {
        return chroms_from_fasta(fi->path, n_out);
    }
    if (fi->ft == BS_BAM || fi->ft == BS_CRAM) {
        if (fi->n_seq_refs <= 0) return NULL;
        Chrom *c = (Chrom *)calloc((size_t)fi->n_seq_refs, sizeof(Chrom));
        for (int i = 0; i < fi->n_seq_refs; i++) {
            c[i].name   = strdup(fi->seq_refs[i].name);
            c[i].length = fi->seq_refs[i].length;
        }
        *n_out = fi->n_seq_refs;
        return c;
    }
    if (fi->ft == BS_VCF || fi->ft == BS_GVCF) {
        if (fi->n_vcf_contigs <= 0) return NULL;
        Chrom *c = (Chrom *)calloc((size_t)fi->n_vcf_contigs, sizeof(Chrom));
        for (int i = 0; i < fi->n_vcf_contigs; i++) {
            c[i].name   = strdup(fi->vcf_contigs[i].name);
            c[i].length = fi->vcf_contigs[i].length;
        }
        *n_out = fi->n_vcf_contigs;
        return c;
    }
    return NULL;
}

/* ────────────────────────────────────────────────────────────────────────
 * Filters
 * ───────────────────────────────────────────────────────────────────── */

static int str_in_csv(const char *s, const char *csv)
{
    if (!csv) return 0;
    size_t n = strlen(s);
    const char *p = csv;
    while (*p) {
        const char *e = p;
        while (*e && *e != ',') e++;
        size_t m = (size_t)(e - p);
        if (m == n && strncmp(p, s, m) == 0) return 1;
        p = (*e == ',') ? e + 1 : e;
    }
    return 0;
}

/* Heuristic for "non-autosomal" contigs: sex/mito/unplaced/alt/decoy. */
static int is_non_autosomal(const char *n)
{
    if (!n || !*n) return 1;
    /* Sex / mito */
    if (strcasecmp(n, "x") == 0 || strcasecmp(n, "y") == 0 ||
        strcasecmp(n, "w") == 0 || strcasecmp(n, "z") == 0 ||
        strcasecmp(n, "mt") == 0 || strcasecmp(n, "m") == 0) return 1;
    if (strcasecmp(n, "chrx") == 0 || strcasecmp(n, "chry") == 0 ||
        strcasecmp(n, "chrm") == 0 || strcasecmp(n, "chrmt") == 0) return 1;
    /* Substring markers in unplaced / alt / decoy / random / patches */
    const char *needles[] = {
        "_random", "Un_", "_alt", "_decoy", "_hap", "_fix",
        "EBV", "decoy", "hap", "patch", NULL
    };
    for (int i = 0; needles[i]; i++) {
        if (strcasestr(n, needles[i]) != NULL) return 1;
    }
    return 0;
}

/* Load a BED file of intervals to mask out. Returns array of {chrom, beg, end}. */
typedef struct { char *chrom; int64_t beg, end; } Region;

static Region *load_mask_bed(const char *path, int *n_out)
{
    FILE *f = fopen(path, "r");
    if (!f) return NULL;
    Region *r = NULL; int n = 0, cap = 0;
    char line[8192];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        char *p = line; char *fields[4] = {0}; int nf = 0;
        fields[nf++] = p;
        while (*p && nf < 4) {
            if (*p == '\t') { *p = '\0'; fields[nf++] = p + 1; p++; continue; }
            if (*p == '\n' || *p == '\r') { *p = '\0'; break; }
            p++;
        }
        if (nf < 3) continue;
        if (n >= cap) { cap = cap ? cap * 2 : 64;
            r = (Region *)realloc(r, sizeof(Region) * (size_t)cap); }
        r[n].chrom = strdup(fields[0]);
        r[n].beg   = strtoll(fields[1], NULL, 10);
        r[n].end   = strtoll(fields[2], NULL, 10);
        n++;
    }
    fclose(f);
    *n_out = n;
    return r;
}

static void free_mask(Region *r, int n)
{
    if (!r) return;
    for (int i = 0; i < n; i++) free(r[i].chrom);
    free(r);
}

/* Return 1 if [w_beg, w_end) overlaps any mask interval on this chrom. */
static int window_masked(const char *chrom, int64_t w_beg, int64_t w_end,
                         const Region *mask, int n_mask)
{
    for (int i = 0; i < n_mask; i++) {
        if (strcmp(chrom, mask[i].chrom) != 0) continue;
        if (w_beg < mask[i].end && w_end > mask[i].beg) return 1;
    }
    return 0;
}

/* ────────────────────────────────────────────────────────────────────────
 * Output
 * ───────────────────────────────────────────────────────────────────── */

typedef struct {
    char    *chrom;
    int64_t  beg;       /* 0-based inclusive */
    int64_t  end;       /* 0-based exclusive */
    char    *name;      /* e.g. "locus42" */
} Window;

/* ────────────────────────────────────────────────────────────────────────
 * Entry point
 * ───────────────────────────────────────────────────────────────────── */

enum {
    OPT_WINDOW = 256, OPT_STEP, OPT_INCLUDE, OPT_EXCLUDE, OPT_AUTOSOMES,
    OPT_SKIP_EDGES, OPT_EXCLUDE_REGIONS, OPT_MIN_SPACING,
    OPT_N_LOCI, OPT_SEED, OPT_OUT, OPT_JSON
};

static const struct option LONG_OPTS[] = {
    {"window-size",     required_argument, NULL, OPT_WINDOW},
    {"step",            required_argument, NULL, OPT_STEP},
    {"include-chrom",   required_argument, NULL, OPT_INCLUDE},
    {"exclude-chrom",   required_argument, NULL, OPT_EXCLUDE},
    {"autosomes-only",  no_argument,       NULL, OPT_AUTOSOMES},
    {"skip-edges",      required_argument, NULL, OPT_SKIP_EDGES},
    {"exclude-regions", required_argument, NULL, OPT_EXCLUDE_REGIONS},
    {"min-spacing",     required_argument, NULL, OPT_MIN_SPACING},
    {"n-loci",          required_argument, NULL, OPT_N_LOCI},
    {"seed",            required_argument, NULL, OPT_SEED},
    {"out",             required_argument, NULL, OPT_OUT},
    {"json",            no_argument,       NULL, OPT_JSON},
    {"help",            no_argument,       NULL, 'h'},
    {NULL, 0, NULL, 0}
};

int cmd_windows(int argc, char **argv)
{
    int64_t window_size = 0;
    int64_t step        = 0;          /* 0 → defaults to window_size */
    int64_t min_spacing = 0;
    int64_t skip_edges  = 0;
    int     autosomes_only = 0;
    char   *include_chrom = NULL;
    char   *exclude_chrom = NULL;
    char   *exclude_regions = NULL;
    int64_t n_loci = 0;               /* 0 = no downsampling */
    unsigned int seed = 0;
    int     seed_set = 0;
    char   *out_path = NULL;
    int     json_mode = 0;

    optind = 1;
    int opt;
    while ((opt = getopt_long(argc, argv, "h", LONG_OPTS, NULL)) != -1) {
        switch (opt) {
            case OPT_WINDOW:           window_size = strtoll(optarg, NULL, 10); break;
            case OPT_STEP:             step        = strtoll(optarg, NULL, 10); break;
            case OPT_INCLUDE:          include_chrom = strdup(optarg); break;
            case OPT_EXCLUDE:          exclude_chrom = strdup(optarg); break;
            case OPT_AUTOSOMES:        autosomes_only = 1; break;
            case OPT_SKIP_EDGES:       skip_edges = strtoll(optarg, NULL, 10); break;
            case OPT_EXCLUDE_REGIONS:  exclude_regions = strdup(optarg); break;
            case OPT_MIN_SPACING:      min_spacing = strtoll(optarg, NULL, 10); break;
            case OPT_N_LOCI:           n_loci = strtoll(optarg, NULL, 10); break;
            case OPT_SEED:             seed = (unsigned)strtoul(optarg, NULL, 10); seed_set = 1; break;
            case OPT_OUT:              out_path = strdup(optarg); break;
            case OPT_JSON:             json_mode = 1; break;
            case 'h':                  print_help(stdout); return 0;
            default:                   print_help(stderr); return 1;
        }
    }
    if (optind >= argc) {
        fprintf(stderr, "Error: missing INPUT (reference FASTA / BAM / VCF).\n");
        print_help(stderr); return 1;
    }
    const char *input = argv[optind];
    if (window_size <= 0) { fprintf(stderr, "Error: --window-size W (>0) is required.\n"); return 1; }
    if (!out_path)        { fprintf(stderr, "Error: --out FILE is required.\n");           return 1; }
    if (step <= 0) step = window_size;
    if (n_loci > 0 && !seed_set) seed = 0;

    /* Inspect input to discover chromosomes */
    FileInfo *fi = inspect_file(input);
    if (!fi || fi->ft == BS_UNKNOWN) {
        fprintf(stderr, "Error: could not determine type of '%s'\n", input);
        file_info_free(fi); return 1;
    }
    int n_chroms = 0;
    Chrom *chroms = chroms_from_file_info(fi, &n_chroms);
    file_info_free(fi);
    if (!chroms || n_chroms == 0) {
        fprintf(stderr, "Error: no chromosomes found in '%s'\n", input);
        return 1;
    }

    /* Load mask, if any */
    int n_mask = 0;
    Region *mask = NULL;
    if (exclude_regions) {
        mask = load_mask_bed(exclude_regions, &n_mask);
        if (!mask) {
            fprintf(stderr, "Error: cannot read --exclude-regions '%s'\n", exclude_regions);
            free_chroms(chroms, n_chroms);
            return 1;
        }
    }

    /* Generate candidate windows */
    Window *wins = NULL;
    int n_wins = 0, cap = 0;
    int64_t total_genome_bp_considered = 0;
    int64_t total_chroms_used = 0;

    for (int ci = 0; ci < n_chroms; ci++) {
        const char *nm = chroms[ci].name;
        int64_t len = chroms[ci].length;

        /* Chromosome filters */
        if (include_chrom && !str_in_csv(nm, include_chrom)) continue;
        if (exclude_chrom &&  str_in_csv(nm, exclude_chrom)) continue;
        if (autosomes_only && is_non_autosomal(nm)) continue;

        int64_t lo = skip_edges;
        int64_t hi = len - skip_edges;
        if (hi - lo < window_size) continue;

        total_chroms_used++;
        total_genome_bp_considered += (hi - lo);

        for (int64_t s = lo; s + window_size <= hi; s += step) {
            int64_t e = s + window_size;
            if (mask && window_masked(nm, s, e, mask, n_mask)) continue;
            if (n_wins >= cap) {
                cap = cap ? cap * 2 : 1024;
                wins = (Window *)realloc(wins, sizeof(Window) * (size_t)cap);
            }
            wins[n_wins].chrom = strdup(nm);
            wins[n_wins].beg   = s;
            wins[n_wins].end   = e;
            wins[n_wins].name  = NULL;
            n_wins++;
        }
    }

    /* min-spacing: greedy walk within each chromosome */
    int n_after_spacing = n_wins;
    if (min_spacing > 0 && n_wins > 0) {
        int *keep = (int *)calloc((size_t)n_wins, sizeof(int));
        int64_t last_end = INT64_MIN;
        const char *last_chrom = NULL;
        for (int i = 0; i < n_wins; i++) {
            if (!last_chrom || strcmp(wins[i].chrom, last_chrom) != 0) {
                keep[i] = 1; last_chrom = wins[i].chrom; last_end = wins[i].end;
                continue;
            }
            if (wins[i].beg - last_end >= min_spacing) {
                keep[i] = 1; last_end = wins[i].end;
            }
        }
        int j = 0;
        for (int i = 0; i < n_wins; i++) {
            if (keep[i]) wins[j++] = wins[i];
            else { free(wins[i].chrom); free(wins[i].name); }
        }
        n_wins = j;
        free(keep);
        n_after_spacing = n_wins;
    }

    /* random downsample to N */
    int n_after_sample = n_wins;
    if (n_loci > 0 && n_wins > (int)n_loci) {
        srand(seed);
        /* Fisher–Yates shuffle then truncate */
        for (int i = n_wins - 1; i > 0; i--) {
            int j = rand() % (i + 1);
            Window t = wins[i]; wins[i] = wins[j]; wins[j] = t;
        }
        for (int i = (int)n_loci; i < n_wins; i++) {
            free(wins[i].chrom); free(wins[i].name);
        }
        n_wins = (int)n_loci;
        /* Re-sort by chrom, beg for tidy BED output */
        /* simple insertion sort suffices for the modest counts BPP uses */
        for (int i = 1; i < n_wins; i++) {
            Window key = wins[i]; int j = i - 1;
            while (j >= 0 && (strcmp(wins[j].chrom, key.chrom) > 0 ||
                              (strcmp(wins[j].chrom, key.chrom) == 0 && wins[j].beg > key.beg))) {
                wins[j+1] = wins[j]; j--;
            }
            wins[j+1] = key;
        }
        n_after_sample = n_wins;
    }

    /* Assign locus names */
    for (int i = 0; i < n_wins; i++) {
        char buf[32]; snprintf(buf, sizeof(buf), "locus%d", i + 1);
        wins[i].name = strdup(buf);
    }

    /* Write BED */
    FILE *out = (strcmp(out_path, "-") == 0) ? stdout : fopen(out_path, "w");
    if (!out) { perror(out_path); free_chroms(chroms, n_chroms); free_mask(mask, n_mask);
        for (int i = 0; i < n_wins; i++) { free(wins[i].chrom); free(wins[i].name); }
        free(wins); free(out_path); free(include_chrom); free(exclude_chrom); free(exclude_regions);
        return 1;
    }
    for (int i = 0; i < n_wins; i++) {
        fprintf(out, "%s\t%lld\t%lld\t%s\n",
                wins[i].chrom, (long long)wins[i].beg, (long long)wins[i].end, wins[i].name);
    }
    if (out != stdout) fclose(out);

    if (json_mode) {
        JsonWriter w; jw_init(&w, stdout, 2);
        jw_obj_open(&w);
        jw_kv_str (&w, "command", "windows");
        jw_kv_str (&w, "input", input);
        jw_kv_str (&w, "out", out_path);
        jw_kv_int (&w, "window_size", window_size);
        jw_kv_int (&w, "step", step);
        jw_kv_int (&w, "min_spacing", min_spacing);
        jw_kv_int (&w, "skip_edges", skip_edges);
        jw_kv_bool(&w, "autosomes_only", autosomes_only);
        jw_kv_int (&w, "n_chroms_considered", total_chroms_used);
        jw_kv_int (&w, "bp_considered", total_genome_bp_considered);
        jw_kv_int (&w, "n_windows_initial", n_after_spacing == n_wins && n_after_sample == n_wins
                                              ? n_wins : (n_after_spacing - 0));
        jw_kv_int (&w, "n_windows_after_spacing", n_after_spacing);
        jw_kv_int (&w, "n_windows_after_sample",  n_after_sample);
        jw_kv_int (&w, "n_windows_emitted", n_wins);
        if (n_loci > 0) jw_kv_int(&w, "seed", seed);
        jw_obj_close(&w);
        jw_finish(&w);
    } else {
        fprintf(stderr,
            "bpp-seqs windows: %d windows on %lld chrom%s "
            "(%lld bp considered, window=%lld, step=%lld",
            n_wins, (long long)total_chroms_used,
            total_chroms_used == 1 ? "" : "s",
            (long long)total_genome_bp_considered,
            (long long)window_size, (long long)step);
        if (min_spacing > 0) fprintf(stderr, ", spacing>=%lld", (long long)min_spacing);
        if (n_loci > 0)      fprintf(stderr, ", sampled to %lld with seed=%u",
                                     (long long)n_loci, seed);
        fprintf(stderr, ")\n");
    }

    /* Cleanup */
    for (int i = 0; i < n_wins; i++) { free(wins[i].chrom); free(wins[i].name); }
    free(wins);
    free_chroms(chroms, n_chroms);
    free_mask(mask, n_mask);
    free(out_path); free(include_chrom); free(exclude_chrom); free(exclude_regions);
    return 0;
}
