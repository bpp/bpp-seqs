#define _POSIX_C_SOURCE 200809L

#include "cmd_extract.h"
#include "bpp_parser.h"
#include "converters/aln_writer.h"
#include "json_writer.h"

#include <ctype.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static void print_help(FILE *fp)
{
    fprintf(fp,
"Usage: bpp-seqs extract INPUT.txt --out PREFIX [selection ...]\n"
"\n"
"Read a BPP-format sequence file and write a new one containing a\n"
"subset of its loci.  A sibling <INPUT>.loci.tsv (if present) is used\n"
"to match by source coordinates and is filtered to the same subset.\n"
"<INPUT>.imap is copied through unchanged.\n"
"\n"
"Required:\n"
"  INPUT.txt                BPP sequence file produced by `bpp-seqs`\n"
"  --out PREFIX             output file prefix\n"
"\n"
"Selection (at least one; filters compose with AND):\n"
"  --loci NAME[,NAME...]    keep these locus names (union with --loci-file)\n"
"  --loci-file FILE         one locus name per line\n"
"  --range A-B[,C-D]        keep loci by 1-based index (1 = first locus)\n"
"  --first N                keep the first N loci\n"
"  --last N                 keep the last N loci\n"
"  --chrom NAME             keep loci whose source_chrom == NAME\n"
"                           (requires <INPUT>.loci.tsv or --loci-tsv)\n"
"  --min-sites N            keep loci with n_sites >= N\n"
"  --max-sites N            keep loci with n_sites <= N\n"
"  --invert                 keep loci NOT matching the selection\n"
"\n"
"Sidecar overrides (auto-detected by default):\n"
"  --imap FILE              use this Imap instead of <INPUT>.imap\n"
"  --loci-tsv FILE          use this provenance instead of <INPUT>.loci.tsv\n"
"\n"
"Other:\n"
"  --json                   emit JSON summary on stdout\n"
"  --quiet                  suppress stderr progress\n"
"  -h, --help               this message\n");
}

/* ────────────────────────────────────────────────────────────────────────
 * Selection sets
 * ───────────────────────────────────────────────────────────────────── */

typedef struct {
    /* Name-based filter: union of --loci and --loci-file */
    char **names;
    int    n_names;
    int    have_name_filter;

    /* Index-based filter: union of --range, --first, --last */
    int   *index_keep;       /* 1-based index → 1 if kept */
    int    index_cap;
    int    have_index_filter;
    int    first_n;          /* 0 = no --first */
    int    last_n;           /* 0 = no --last */

    /* Single-value filters */
    char  *chrom;            /* NULL = none */
    int    min_sites;        /* -1 = none */
    int    max_sites;        /* -1 = none */

    int    invert;
} Selection;

static void sel_init(Selection *s) { memset(s, 0, sizeof(*s)); s->min_sites = s->max_sites = -1; }

static void sel_free(Selection *s)
{
    for (int i = 0; i < s->n_names; i++) free(s->names[i]);
    free(s->names);
    free(s->index_keep);
    free(s->chrom);
}

static void sel_add_name(Selection *s, const char *name)
{
    s->names = (char **)realloc(s->names, sizeof(char *) * (size_t)(s->n_names + 1));
    s->names[s->n_names++] = strdup(name);
    s->have_name_filter = 1;
}

static void sel_load_names_file(Selection *s, const char *path)
{
    FILE *fp = fopen(path, "r");
    if (!fp) { fprintf(stderr, "Error: cannot open --loci-file '%s'\n", path); return; }
    char line[4096];
    while (fgets(line, sizeof(line), fp)) {
        char *p = line;
        while (*p && isspace((unsigned char)*p)) p++;
        char *e = p;
        while (*e && !isspace((unsigned char)*e)) e++;
        size_t n = (size_t)(e - p);
        if (n > 0) {
            char tmp[1024];
            if (n >= sizeof(tmp)) n = sizeof(tmp) - 1;
            memcpy(tmp, p, n); tmp[n] = '\0';
            sel_add_name(s, tmp);
        }
    }
    fclose(fp);
}

static void sel_add_range(Selection *s, const char *spec, int max_loci)
{
    /* spec is "A-B" or "A-B,C-D,..." */
    if (!s->index_keep) {
        s->index_cap = max_loci + 1;
        s->index_keep = (int *)calloc((size_t)s->index_cap, sizeof(int));
    }
    s->have_index_filter = 1;
    char *buf = strdup(spec);
    for (char *tok = strtok(buf, ","); tok; tok = strtok(NULL, ",")) {
        int a = 0, b = 0;
        if (sscanf(tok, "%d-%d", &a, &b) == 2) {
            if (a < 1) a = 1; if (b > max_loci) b = max_loci;
            for (int i = a; i <= b; i++) if (i < s->index_cap) s->index_keep[i] = 1;
        } else if (sscanf(tok, "%d", &a) == 1) {
            if (a >= 1 && a < s->index_cap) s->index_keep[a] = 1;
        }
    }
    free(buf);
}

static void sel_apply_first_last(Selection *s, int n_total)
{
    if (s->first_n <= 0 && s->last_n <= 0) return;
    if (!s->index_keep) {
        s->index_cap = n_total + 1;
        s->index_keep = (int *)calloc((size_t)s->index_cap, sizeof(int));
    }
    s->have_index_filter = 1;
    if (s->first_n > 0) {
        int hi = s->first_n > n_total ? n_total : s->first_n;
        for (int i = 1; i <= hi; i++) s->index_keep[i] = 1;
    }
    if (s->last_n > 0) {
        int lo = n_total - s->last_n + 1;
        if (lo < 1) lo = 1;
        for (int i = lo; i <= n_total; i++) s->index_keep[i] = 1;
    }
}

/* Returns 1 if this locus should be kept (before --invert). */
static int sel_matches(const Selection *s, const BppLocus *l, int idx_1based)
{
    if (s->have_name_filter) {
        int hit = 0;
        for (int i = 0; i < s->n_names; i++) {
            if (l->name && strcmp(l->name, s->names[i]) == 0) { hit = 1; break; }
        }
        if (!hit) return 0;
    }
    if (s->have_index_filter) {
        if (idx_1based >= s->index_cap || !s->index_keep[idx_1based]) return 0;
    }
    if (s->chrom) {
        if (!l->source_chrom || strcmp(l->source_chrom, s->chrom) != 0) return 0;
    }
    if (s->min_sites > 0 && l->n_sites < s->min_sites) return 0;
    if (s->max_sites > 0 && l->n_sites > s->max_sites) return 0;
    return 1;
}

/* ────────────────────────────────────────────────────────────────────────
 * I/O helpers
 * ───────────────────────────────────────────────────────────────────── */

static int file_exists(const char *p) { struct stat st; return p && stat(p, &st) == 0; }

static char *replace_ext(const char *path, const char *new_suffix)
{
    /* Drop a single ".txt" suffix if present, then append new_suffix. */
    size_t n = strlen(path);
    const char *dot = strrchr(path, '.');
    size_t base = (dot && strcmp(dot, ".txt") == 0) ? (size_t)(dot - path) : n;
    char *r = (char *)malloc(base + strlen(new_suffix) + 1);
    memcpy(r, path, base);
    strcpy(r + base, new_suffix);
    return r;
}

/* Copy src → dst as bytes.  Returns 0 on success. */
static int copy_file(const char *src, const char *dst)
{
    FILE *in = fopen(src, "rb");
    if (!in) return -1;
    FILE *out = fopen(dst, "wb");
    if (!out) { fclose(in); return -1; }
    char buf[8192];
    size_t got;
    while ((got = fread(buf, 1, sizeof(buf), in)) > 0) fwrite(buf, 1, got, out);
    fclose(in); fclose(out);
    return 0;
}

/* Write the selected loci to <prefix>.txt in BPP format. */
static int write_bpp_subset(const char *prefix, BppLocus *loci, int n_loci, int *keep)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s.txt", prefix);
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); return -1; }
    int first = 1;
    for (int i = 0; i < n_loci; i++) {
        if (!keep[i]) continue;
        if (!first) fputc('\n', f);
        first = 0;
        fprintf(f, "%d %d\n\n", loci[i].n_seqs, loci[i].n_sites);
        int max_name = 0;
        for (int j = 0; j < loci[i].n_seqs; j++) {
            int l = (int)strlen(loci[i].seq_names[j]);
            if (l > max_name) max_name = l;
        }
        int col = max_name + 4;
        for (int j = 0; j < loci[i].n_seqs; j++) {
            int wrote = fprintf(f, "%s", loci[i].seq_names[j]);
            int pad = col - wrote; if (pad < 1) pad = 1;
            for (int p = 0; p < pad; p++) fputc(' ', f);
            fprintf(f, "%s\n", loci[i].seqs[j]);
        }
    }
    fclose(f);
    return 0;
}

/* ────────────────────────────────────────────────────────────────────────
 * Entry point
 * ───────────────────────────────────────────────────────────────────── */

enum {
    OPT_OUT = 256, OPT_LOCI, OPT_LOCI_FILE, OPT_RANGE,
    OPT_FIRST, OPT_LAST, OPT_CHROM, OPT_MIN_SITES, OPT_MAX_SITES,
    OPT_INVERT, OPT_IMAP, OPT_LOCI_TSV, OPT_JSON, OPT_QUIET
};

static const struct option LONG_OPTS[] = {
    {"out",        required_argument, NULL, OPT_OUT},
    {"loci",       required_argument, NULL, OPT_LOCI},
    {"loci-file",  required_argument, NULL, OPT_LOCI_FILE},
    {"range",      required_argument, NULL, OPT_RANGE},
    {"first",      required_argument, NULL, OPT_FIRST},
    {"last",       required_argument, NULL, OPT_LAST},
    {"chrom",      required_argument, NULL, OPT_CHROM},
    {"min-sites",  required_argument, NULL, OPT_MIN_SITES},
    {"max-sites",  required_argument, NULL, OPT_MAX_SITES},
    {"invert",     no_argument,       NULL, OPT_INVERT},
    {"imap",       required_argument, NULL, OPT_IMAP},
    {"loci-tsv",   required_argument, NULL, OPT_LOCI_TSV},
    {"json",       no_argument,       NULL, OPT_JSON},
    {"quiet",      no_argument,       NULL, OPT_QUIET},
    {"help",       no_argument,       NULL, 'h'},
    {NULL, 0, NULL, 0}
};

int cmd_extract(int argc, char **argv)
{
    Selection sel; sel_init(&sel);
    char *out_prefix = NULL;
    char *range_spec = NULL;
    char *imap_override = NULL;
    char *loci_tsv_override = NULL;
    int json_mode = 0;
    int quiet = 0;

    optind = 1;  /* called as a subcommand; argv[1] is the input filename */

    /* We need to scan positional args before/around getopt. The expected
     * shape is: "bpp-seqs extract INPUT.txt --out ... [--flag ...]".
     * Caller's argv already has "extract" stripped; argv[0] = "extract",
     * argv[1..] = the rest. */

    int opt;
    while ((opt = getopt_long(argc, argv, "h", LONG_OPTS, NULL)) != -1) {
        switch (opt) {
            case OPT_OUT:        out_prefix = strdup(optarg); break;
            case OPT_LOCI: {
                char *buf = strdup(optarg);
                for (char *t = strtok(buf, ","); t; t = strtok(NULL, ",")) sel_add_name(&sel, t);
                free(buf);
                break;
            }
            case OPT_LOCI_FILE:  sel_load_names_file(&sel, optarg); break;
            case OPT_RANGE:      range_spec = strdup(optarg); break;
            case OPT_FIRST:      sel.first_n = atoi(optarg); break;
            case OPT_LAST:       sel.last_n  = atoi(optarg); break;
            case OPT_CHROM:      free(sel.chrom); sel.chrom = strdup(optarg); break;
            case OPT_MIN_SITES:  sel.min_sites = atoi(optarg); break;
            case OPT_MAX_SITES:  sel.max_sites = atoi(optarg); break;
            case OPT_INVERT:     sel.invert = 1; break;
            case OPT_IMAP:       imap_override = strdup(optarg); break;
            case OPT_LOCI_TSV:   loci_tsv_override = strdup(optarg); break;
            case OPT_JSON:       json_mode = 1; break;
            case OPT_QUIET:      quiet = 1; break;
            case 'h':            print_help(stdout); sel_free(&sel); free(out_prefix); return 0;
            default:             print_help(stderr); sel_free(&sel); free(out_prefix); return 1;
        }
    }

    /* Remaining positional argument is the input file. */
    if (optind >= argc) {
        fprintf(stderr, "Error: missing INPUT.txt argument.\n");
        print_help(stderr);
        sel_free(&sel); free(out_prefix);
        return 1;
    }
    const char *input_path = argv[optind];

    if (!out_prefix) {
        fprintf(stderr, "Error: --out PREFIX is required.\n");
        sel_free(&sel);
        return 1;
    }

    int have_any_selector =
        sel.have_name_filter || sel.have_index_filter || range_spec ||
        sel.first_n > 0 || sel.last_n > 0 || sel.chrom ||
        sel.min_sites > 0 || sel.max_sites > 0;
    if (!have_any_selector) {
        fprintf(stderr, "Error: at least one selection flag is required.\n"
                        "       (--loci, --loci-file, --range, --first, --last,\n"
                        "        --chrom, --min-sites, or --max-sites)\n");
        sel_free(&sel); free(out_prefix);
        return 1;
    }

    /* Parse input */
    int n_loci = 0;
    BppLocus *loci = bpp_parse_file(input_path, &n_loci);
    if (!loci || n_loci == 0) {
        fprintf(stderr, "Error: no loci parsed from '%s'\n", input_path);
        sel_free(&sel); free(out_prefix); free(range_spec);
        free(imap_override); free(loci_tsv_override);
        return 1;
    }

    /* Locate sidecars */
    char *tsv_path = loci_tsv_override
        ? strdup(loci_tsv_override)
        : replace_ext(input_path, ".loci.tsv");
    char *imap_path = imap_override
        ? strdup(imap_override)
        : replace_ext(input_path, ".imap");
    int have_tsv  = file_exists(tsv_path);
    int have_imap = file_exists(imap_path);
    if (have_tsv) bpp_attach_loci_tsv(loci, n_loci, tsv_path);

    /* Now that we know n_loci, materialize index-based selectors. */
    if (range_spec) sel_add_range(&sel, range_spec, n_loci);
    sel_apply_first_last(&sel, n_loci);

    /* If --chrom was requested but no provenance is loaded, that's an error. */
    if (sel.chrom && !have_tsv) {
        fprintf(stderr, "Error: --chrom requires <INPUT>.loci.tsv (or --loci-tsv).\n");
        bpp_loci_free(loci, n_loci);
        sel_free(&sel); free(out_prefix); free(range_spec);
        free(tsv_path); free(imap_path);
        free(imap_override); free(loci_tsv_override);
        return 1;
    }

    /* Decide which loci to keep. */
    int *keep = (int *)calloc((size_t)n_loci, sizeof(int));
    int n_kept = 0;
    for (int i = 0; i < n_loci; i++) {
        int m = sel_matches(&sel, &loci[i], i + 1);
        keep[i] = sel.invert ? !m : m;
        if (keep[i]) n_kept++;
    }

    if (n_kept == 0) {
        fprintf(stderr, "Error: selection matched no loci.\n");
        free(keep);
        bpp_loci_free(loci, n_loci);
        sel_free(&sel); free(out_prefix); free(range_spec);
        free(tsv_path); free(imap_path);
        free(imap_override); free(loci_tsv_override);
        return 1;
    }

    /* Write outputs. */
    if (write_bpp_subset(out_prefix, loci, n_loci, keep) != 0) {
        free(keep); bpp_loci_free(loci, n_loci);
        sel_free(&sel); free(out_prefix); free(range_spec);
        free(tsv_path); free(imap_path);
        free(imap_override); free(loci_tsv_override);
        return 1;
    }

    char *out_txt   = replace_ext(out_prefix, ".txt");
    char *out_imap  = NULL;
    char *out_tsv   = NULL;
    (void)out_txt; /* path already known to be <prefix>.txt */

    /* Filtered loci.tsv */
    if (have_tsv) {
        LocusProv *items = (LocusProv *)calloc((size_t)n_kept, sizeof(LocusProv));
        int k = 0;
        for (int i = 0; i < n_loci; i++) {
            if (!keep[i]) continue;
            items[k].name   = loci[i].name;
            items[k].kind   = loci[i].source_kind;
            items[k].file   = loci[i].source_file;
            items[k].chrom  = loci[i].source_chrom;
            items[k].start  = loci[i].source_start;
            items[k].end    = loci[i].source_end;
            items[k].stride = loci[i].source_stride > 0 ? loci[i].source_stride : 1;
            items[k].length = loci[i].n_sites;
            items[k].n_seqs = loci[i].n_seqs;
            k++;
        }
        write_loci_tsv(out_prefix, items, n_kept);
        char p[1024]; snprintf(p, sizeof(p), "%s.loci.tsv", out_prefix);
        out_tsv = strdup(p);
        free(items);
    }

    /* Pass-through imap */
    if (have_imap) {
        char p[1024]; snprintf(p, sizeof(p), "%s.imap", out_prefix);
        copy_file(imap_path, p);
        out_imap = strdup(p);
    }

    if (json_mode) {
        JsonWriter w; jw_init(&w, stdout, 2);
        jw_obj_open(&w);
        jw_kv_str (&w, "command", "extract");
        jw_kv_str (&w, "input",   input_path);
        jw_kv_int (&w, "n_loci_input",  n_loci);
        jw_kv_int (&w, "n_loci_kept",   n_kept);
        jw_kv_int (&w, "n_loci_dropped", n_loci - n_kept);
        jw_kv_bool(&w, "had_loci_tsv",  have_tsv);
        jw_kv_bool(&w, "had_imap",      have_imap);
        jw_key(&w, "output_files"); jw_obj_open(&w);
        char p[1024]; snprintf(p, sizeof(p), "%s.txt", out_prefix);
        jw_kv_str(&w, "sequences", p);
        jw_kv_str(&w, "loci",      out_tsv);
        jw_kv_str(&w, "imap",      out_imap);
        jw_obj_close(&w);
        jw_obj_close(&w);
        jw_finish(&w);
    } else if (!quiet) {
        char p[1024]; snprintf(p, sizeof(p), "%s.txt", out_prefix);
        fprintf(stdout,
            "bpp-seqs extract: kept %d / %d loci from %s\n"
            "  %s  BPP sequence file\n", n_kept, n_loci, input_path, p);
        if (out_tsv)  fprintf(stdout, "  %s  per-locus provenance\n", out_tsv);
        if (out_imap) fprintf(stdout, "  %s  Imap (copied unchanged)\n", out_imap);
    }

    free(out_txt); free(out_imap); free(out_tsv);
    free(keep);
    bpp_loci_free(loci, n_loci);
    sel_free(&sel);
    free(out_prefix); free(range_spec);
    free(tsv_path); free(imap_path);
    free(imap_override); free(loci_tsv_override);
    return 0;
}
