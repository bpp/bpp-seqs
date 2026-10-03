#define _POSIX_C_SOURCE 200809L

#include "cmd_check.h"
#include "inspect.h"
#include "lineio.h"
#include "json_writer.h"
#include "nexus.h"

#include <ctype.h>
#include <getopt.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static void print_help(FILE *fp)
{
    fprintf(fp,
"Usage: bpp-seqs check FILE [--nloci N] [--json]\n"
"\n"
"Identify the format of a sequence file and check that its contents agree\n"
"with what it declares: the number of loci, and, for every locus, the number\n"
"of sequences and the number of sites in each sequence. Nothing is written.\n"
"\n"
"Formats (detected from content, gzip allowed):\n"
"  BPP / PHYLIP   each locus starts with an `<n_seqs> <n_sites>` header; every\n"
"                 locus is checked, in sequential, wrapped (indented\n"
"                 continuation rows) or interleaved layout\n"
"  FASTA          one locus; all records must have the same length\n"
"  NEXUS          one locus, or one per CHARSET; the MATRIX must match\n"
"                 DIMENSIONS ntax/nchar and each CHARSET must lie within nchar\n"
"\n"
"Options:\n"
"  --nloci N                also require exactly N loci (e.g. the control\n"
"                           file's nloci)\n"
"  --json                   emit the report as JSON on stdout\n"
"  -h, --help               this message\n"
"\n"
"Lines of any length are read whole, so loci longer than 64 kb on one line\n"
"are counted correctly.\n"
"\n"
"Exit status: 0 if the file passes, 2 if it has errors (listed in the\n"
"report), 1 on a usage error or a file that cannot be read or is not a\n"
"sequence alignment.\n"
"\n"
"Example\n"
"  bpp-seqs check run1.txt --nloci 349\n");
}

/* ────────────────────────────────────────────────────────────────────────
 * Report
 * ───────────────────────────────────────────────────────────────────── */

typedef struct {
    int   locus;          /* 1-based; 0 = file-level */
    long  line;           /* 1-based; 0 = not tied to a line */
    char  code[24];
    char *msg;
} CheckErr;

typedef struct {
    long        line;            /* line the locus starts on */
    int         declared_seqs;   /* -1 when the format declares none */
    int         declared_sites;  /* -1 when the format declares none */
    int         seqs;            /* sequences found */
    long        min_sites, max_sites;
    const char *layout;          /* NULL when not applicable */
    int         ok;
} LocusRec;

typedef struct {
    const char *format;
    CheckErr   *errs;
    int         n_errs, cap_errs;
    LocusRec   *loci;
    int         n_loci, cap_loci;
} Report;

static void add_err(Report *r, int locus, long line, const char *code,
                    const char *fmt, ...)
{
    if (r->n_errs >= r->cap_errs) {
        r->cap_errs = r->cap_errs ? r->cap_errs * 2 : 16;
        r->errs = (CheckErr *)realloc(r->errs, sizeof(CheckErr) * (size_t)r->cap_errs);
    }
    CheckErr *e = &r->errs[r->n_errs++];
    e->locus = locus;
    e->line  = line;
    snprintf(e->code, sizeof e->code, "%s", code);
    char buf[512];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    e->msg = strdup(buf);
}

static LocusRec *add_locus(Report *r)
{
    if (r->n_loci >= r->cap_loci) {
        r->cap_loci = r->cap_loci ? r->cap_loci * 2 : 64;
        r->loci = (LocusRec *)realloc(r->loci, sizeof(LocusRec) * (size_t)r->cap_loci);
    }
    LocusRec *l = &r->loci[r->n_loci++];
    memset(l, 0, sizeof(*l));
    l->declared_seqs = l->declared_sites = -1;
    l->ok = 1;
    return l;
}

static int is_blank(const char *s)
{
    for (; *s; s++) if (!isspace((unsigned char)*s)) return 0;
    return 1;
}

/* A locus header is exactly two positive integers on a line. */
static int parse_header(const char *line, int *ns, int *nc)
{
    int a = 0, b = 0, used = 0;
    if (sscanf(line, " %d %d%n", &a, &b, &used) != 2 || a <= 0 || b <= 0) return 0;
    if (!is_blank(line + used)) return 0;
    *ns = a; *nc = b;
    return 1;
}

static long count_seq_chars(const char *s)
{
    long n = 0;
    for (; *s; s++) if (!isspace((unsigned char)*s)) n++;
    return n;
}

/* ────────────────────────────────────────────────────────────────────────
 * BPP / PHYLIP
 *
 * Row layouts follow the rest of bpp-seqs (inspect.c, bpp_parser.c):
 *   sequential   one named row per sequence
 *   wrapped      a named row followed by indented, nameless continuation rows
 *   interleaved  named rows for the first block, then nameless blocks that
 *                cycle through the sequences in order; detected, as in
 *                bpp_parser.c, by the first sequence still being short once
 *                every sequence has been named
 * ───────────────────────────────────────────────────────────────────── */

typedef struct {
    int    idx;
    long   line;
    int    nseq, nsites;
    int    named;
    long   extra;         /* named rows beyond the declared count */
    long   extra_line;
    char **names;
    long  *name_lines;
    long  *lens;
    int    interleaved, wrapped;
    int    ilv;
} PhyLocus;

static void phy_begin(PhyLocus *L, int idx, long line, int nseq, int nsites)
{
    memset(L, 0, sizeof(*L));
    L->idx = idx; L->line = line; L->nseq = nseq; L->nsites = nsites;
    L->names      = (char **)calloc((size_t)nseq, sizeof(char *));
    L->name_lines = (long *) calloc((size_t)nseq, sizeof(long));
    L->lens       = (long *) calloc((size_t)nseq, sizeof(long));
}

static void phy_finish(PhyLocus *L, Report *r)
{
    LocusRec *rec = add_locus(r);
    rec->line           = L->line;
    rec->declared_seqs  = L->nseq;
    rec->declared_sites = L->nsites;
    rec->seqs           = L->named + (int)L->extra;
    rec->layout = L->interleaved ? "interleaved" : L->wrapped ? "wrapped" : "sequential";

    int errs_before = r->n_errs;
    if (L->named < L->nseq)
        add_err(r, L->idx, L->line, "SEQ_COUNT",
                "locus %d declares %d sequences but has %d",
                L->idx, L->nseq, L->named);
    else if (L->extra > 0)
        add_err(r, L->idx, L->extra_line, "SEQ_COUNT",
                "locus %d declares %d sequences but has %ld",
                L->idx, L->nseq, (long)L->nseq + L->extra);

    long mn = -1, mx = -1;
    int n_bad = 0;
    for (int i = 0; i < L->named; i++) {
        if (mn < 0 || L->lens[i] < mn) mn = L->lens[i];
        if (L->lens[i] > mx) mx = L->lens[i];
        if (L->lens[i] != L->nsites) n_bad++;
    }
    rec->min_sites = mn < 0 ? 0 : mn;
    rec->max_sites = mx < 0 ? 0 : mx;

    if (n_bad > 0 && n_bad == L->named && mn == mx) {
        /* Every sequence agrees with every other: the header is the odd one out. */
        add_err(r, L->idx, L->line, "SITE_COUNT",
                "locus %d declares %d sites but all %d sequences have %ld",
                L->idx, L->nsites, L->named, mn);
    } else if (n_bad > 0) {
        for (int i = 0; i < L->named; i++)
            if (L->lens[i] != L->nsites)
                add_err(r, L->idx, L->name_lines[i], "SITE_COUNT",
                        "locus %d: sequence '%s' has %ld sites; header declares %d",
                        L->idx, L->names[i], L->lens[i], L->nsites);
    }
    rec->ok = (r->n_errs == errs_before);

    for (int i = 0; i < L->named; i++) free(L->names[i]);
    free(L->names); free(L->name_lines); free(L->lens);
    memset(L, 0, sizeof(*L));
}

static void phy_row(PhyLocus *L, const char *line, long lineno)
{
    int indented = (line[0] == ' ' || line[0] == '\t');
    const char *p = line;
    int target;

    if (!L->interleaved && L->named == L->nseq && L->nseq > 0 &&
        L->lens[0] < L->nsites) {
        L->interleaved = 1;
        L->ilv = 0;
    }

    if (L->interleaved) {
        target = L->ilv % L->nseq;
        L->ilv++;
    } else if (indented && L->named > 0) {
        target = L->named - 1;
        L->wrapped = 1;
    } else if (L->named < L->nseq) {
        while (*p && isspace((unsigned char)*p)) p++;
        const char *q = p;
        while (*q && !isspace((unsigned char)*q)) q++;
        target = L->named++;
        L->names[target] = strndup(p, (size_t)(q - p));
        L->name_lines[target] = lineno;
        p = q;
    } else {
        if (L->extra++ == 0) L->extra_line = lineno;
        return;
    }
    L->lens[target] += count_seq_chars(p);
}

static int check_phylip(const char *path, Report *r)
{
    gzFile gz = gzopen(path, "rb");
    if (!gz) { fprintf(stderr, "Error: cannot open '%s'\n", path); return -1; }

    char  *buf = NULL;
    size_t cap = 0;
    long   lineno = 0, len;
    int    in_locus = 0, n = 0;
    PhyLocus L;
    memset(&L, 0, sizeof(L));

    while ((len = gz_getline(gz, &buf, &cap)) >= 0) {
        lineno++;
        if (is_blank(buf)) {
            if (in_locus && L.interleaved) L.ilv = 0;
            continue;
        }
        int ns, nc;
        if (parse_header(buf, &ns, &nc)) {
            if (in_locus) phy_finish(&L, r);
            phy_begin(&L, ++n, lineno, ns, nc);
            in_locus = 1;
            continue;
        }
        /* Detection guarantees the first line is a header, so in_locus holds. */
        if (in_locus) phy_row(&L, buf, lineno);
    }
    if (in_locus) phy_finish(&L, r);
    gzclose(gz);
    free(buf);

    r->format = n > 1 ? "BPP" : "PHYLIP";
    return 0;
}

/* ────────────────────────────────────────────────────────────────────────
 * FASTA: one locus, no declared counts; records must share a length
 * ───────────────────────────────────────────────────────────────────── */

static int check_fasta(const char *path, Report *r)
{
    gzFile gz = gzopen(path, "rb");
    if (!gz) { fprintf(stderr, "Error: cannot open '%s'\n", path); return -1; }

    char  *buf = NULL;
    size_t cap = 0;
    long   lineno = 0, len;
    int    n = 0, ncap = 0;
    char **names = NULL;
    long  *lens = NULL, *lines = NULL;

    while ((len = gz_getline(gz, &buf, &cap)) >= 0) {
        lineno++;
        if (buf[0] == '>') {
            if (n >= ncap) {
                ncap = ncap ? ncap * 2 : 64;
                names = (char **)realloc(names, sizeof(char *) * (size_t)ncap);
                lens  = (long *) realloc(lens,  sizeof(long)   * (size_t)ncap);
                lines = (long *) realloc(lines, sizeof(long)   * (size_t)ncap);
            }
            const char *q = buf + 1;
            while (*q && !isspace((unsigned char)*q)) q++;
            names[n] = strndup(buf + 1, (size_t)(q - buf - 1));
            lens[n]  = 0;
            lines[n] = lineno;
            n++;
        } else if (n > 0) {
            lens[n - 1] += count_seq_chars(buf);
        }
    }
    gzclose(gz);
    free(buf);

    r->format = "FASTA";
    if (n == 0) { free(names); free(lens); free(lines); return 0; }

    LocusRec *rec = add_locus(r);
    rec->line = lines[0];
    rec->seqs = n;
    rec->min_sites = rec->max_sites = lens[0];
    for (int i = 1; i < n; i++) {
        if (lens[i] < rec->min_sites) rec->min_sites = lens[i];
        if (lens[i] > rec->max_sites) rec->max_sites = lens[i];
    }
    for (int i = 0; i < n; i++) {
        if (lens[i] == 0)
            add_err(r, 1, lines[i], "EMPTY_SEQ", "sequence '%s' is empty", names[i]);
        else if (lens[i] != lens[0])
            add_err(r, 1, lines[i], "SITE_COUNT",
                    "sequence '%s' has %ld sites; the first ('%s') has %ld",
                    names[i], lens[i], names[0], lens[0]);
    }
    rec->ok = (r->n_errs == 0);

    for (int i = 0; i < n; i++) free(names[i]);
    free(names); free(lens); free(lines);
    return 0;
}

/* ────────────────────────────────────────────────────────────────────────
 * NEXUS: the shared reader enforces ntax/nchar; loci are its CHARSETs
 * ───────────────────────────────────────────────────────────────────── */

static int check_nexus(const char *path, Report *r)
{
    NexusDoc d;
    memset(&d, 0, sizeof(d));
    r->format = "NEXUS";
    if (nexus_parse(path, &d) != 0) {
        add_err(r, 0, 0, "NEXUS_PARSE", "%s", d.err[0] ? d.err : "could not parse NEXUS file");
        nexus_free(&d);
        return 0;
    }
    if (d.n_excess > 0)
        add_err(r, 0, 0, "SITE_COUNT",
                "MATRIX has %ld character%s beyond the declared nchar=%d",
                d.n_excess, d.n_excess == 1 ? "" : "s", d.nchar);

    int nloci = d.n_charsets ? d.n_charsets : 1;
    for (int i = 0; i < nloci; i++) {
        int errs_before = r->n_errs;
        LocusRec *rec = add_locus(r);
        rec->declared_seqs  = d.ntax;
        rec->seqs           = d.ntax;
        if (!d.n_charsets) {
            rec->declared_sites = d.nchar;
            rec->min_sites = rec->max_sites = d.nchar;
        } else {
            int s = d.cs_start[i], e = d.cs_end[i] < 0 ? d.nchar : d.cs_end[i];
            int st = d.cs_stride[i] > 0 ? d.cs_stride[i] : 1;
            if (s < 1 || e < s || e > d.nchar)
                add_err(r, i + 1, 0, "CHARSET_RANGE",
                        "CHARSET '%s' (%d-%d) lies outside the %d declared sites",
                        d.cs_name[i], s, e, d.nchar);
            long nsites = (s >= 1 && e >= s) ? (long)(e - s) / st + 1 : 0;
            rec->min_sites = rec->max_sites = nsites;
        }
        rec->ok = (r->n_errs == errs_before);
    }
    nexus_free(&d);
    return 0;
}

/* ────────────────────────────────────────────────────────────────────────
 * Output
 * ───────────────────────────────────────────────────────────────────── */

#define TEXT_MAX_ERRORS 20

static void print_text(const char *path, const Report *r, int expected)
{
    printf("File:    %s\n", path);
    printf("Format:  %s%s\n", r->format,
           strcmp(r->format, "BPP") == 0 ? " (multi-locus sequence file)" : "");
    if (expected >= 0) printf("Loci:    %d (expected %d)\n", r->n_loci, expected);
    else               printf("Loci:    %d\n", r->n_loci);

    if (r->n_loci > 0) {
        int  smin = r->loci[0].seqs, smax = smin;
        long lmin = r->loci[0].min_sites, lmax = r->loci[0].max_sites, total = 0;
        int  n_seq = 0, n_wrap = 0, n_ilv = 0;
        for (int i = 0; i < r->n_loci; i++) {
            const LocusRec *l = &r->loci[i];
            if (l->seqs < smin) smin = l->seqs;
            if (l->seqs > smax) smax = l->seqs;
            if (l->min_sites < lmin) lmin = l->min_sites;
            if (l->max_sites > lmax) lmax = l->max_sites;
            total += l->declared_sites >= 0 ? l->declared_sites : l->max_sites;
            if (l->layout) {
                if (strcmp(l->layout, "interleaved") == 0) n_ilv++;
                else if (strcmp(l->layout, "wrapped") == 0) n_wrap++;
                else n_seq++;
            }
        }
        if (smin == smax) printf("Sequences per locus: %d\n", smin);
        else              printf("Sequences per locus: %d-%d\n", smin, smax);
        if (lmin == lmax) printf("Sites per locus:     %ld (total %ld)\n", lmin, total);
        else              printf("Sites per locus:     %ld-%ld (total %ld)\n", lmin, lmax, total);
        if (n_seq + n_wrap + n_ilv > 0) {
            if (n_wrap == 0 && n_ilv == 0)      printf("Layout:  sequential\n");
            else if (n_seq == 0 && n_ilv == 0)  printf("Layout:  wrapped\n");
            else if (n_seq == 0 && n_wrap == 0) printf("Layout:  interleaved\n");
            else printf("Layout:  mixed (%d sequential, %d wrapped, %d interleaved)\n",
                        n_seq, n_wrap, n_ilv);
        }
    }

    if (r->n_errs == 0) { printf("\nResult:  OK\n"); return; }
    printf("\nResult:  %d error%s\n", r->n_errs, r->n_errs == 1 ? "" : "s");
    for (int i = 0; i < r->n_errs && i < TEXT_MAX_ERRORS; i++) {
        const CheckErr *e = &r->errs[i];
        if (e->line > 0) printf("  line %ld: [%s] %s\n", e->line, e->code, e->msg);
        else             printf("  [%s] %s\n", e->code, e->msg);
    }
    if (r->n_errs > TEXT_MAX_ERRORS)
        printf("  ... and %d more (use --json for the full list)\n",
               r->n_errs - TEXT_MAX_ERRORS);
}

static void print_json(const char *path, const char *detected, const Report *r,
                       int expected)
{
    JsonWriter w; jw_init(&w, stdout, 2);
    jw_obj_open(&w);
    jw_kv_str(&w, "file", path);
    jw_kv_str(&w, "detected_type", detected);
    jw_kv_str(&w, "format", r->format);
    jw_kv_bool(&w, "ok", r->n_errs == 0);
    jw_kv_int(&w, "n_loci", r->n_loci);
    if (expected >= 0) jw_kv_int(&w, "expected_nloci", expected);
    else               jw_kv_null(&w, "expected_nloci");
    jw_kv_int(&w, "n_errors", r->n_errs);

    jw_key(&w, "errors");
    jw_arr_open(&w);
    for (int i = 0; i < r->n_errs; i++) {
        const CheckErr *e = &r->errs[i];
        jw_obj_open(&w);
        jw_kv_str(&w, "code", e->code);
        if (e->locus > 0) jw_kv_int(&w, "locus", e->locus); else jw_kv_null(&w, "locus");
        if (e->line > 0)  jw_kv_int(&w, "line", e->line);   else jw_kv_null(&w, "line");
        jw_kv_str(&w, "message", e->msg);
        jw_obj_close(&w);
    }
    jw_arr_close(&w);

    jw_key(&w, "loci");
    jw_arr_open(&w);
    for (int i = 0; i < r->n_loci; i++) {
        const LocusRec *l = &r->loci[i];
        jw_obj_open(&w);
        jw_kv_int(&w, "index", i + 1);
        if (l->line > 0) jw_kv_int(&w, "line", l->line); else jw_kv_null(&w, "line");
        if (l->declared_seqs >= 0)  jw_kv_int(&w, "declared_seqs", l->declared_seqs);
        else                        jw_kv_null(&w, "declared_seqs");
        if (l->declared_sites >= 0) jw_kv_int(&w, "declared_sites", l->declared_sites);
        else                        jw_kv_null(&w, "declared_sites");
        jw_kv_int(&w, "seqs", l->seqs);
        jw_kv_int(&w, "min_sites", l->min_sites);
        jw_kv_int(&w, "max_sites", l->max_sites);
        if (l->layout) jw_kv_str(&w, "layout", l->layout); else jw_kv_null(&w, "layout");
        jw_kv_bool(&w, "ok", l->ok);
        jw_obj_close(&w);
    }
    jw_arr_close(&w);
    jw_obj_close(&w);
    jw_finish(&w);
}

/* ────────────────────────────────────────────────────────────────────────
 * Entry point
 * ───────────────────────────────────────────────────────────────────── */

enum { OPT_NLOCI = 256, OPT_JSON };

static const struct option LONG_OPTS[] = {
    {"nloci", required_argument, NULL, OPT_NLOCI},
    {"json",  no_argument,       NULL, OPT_JSON},
    {"help",  no_argument,       NULL, 'h'},
    {NULL, 0, NULL, 0}
};

int cmd_check(int argc, char **argv)
{
    int json_mode = 0, expected = -1, opt;
    while ((opt = getopt_long(argc, argv, "h", LONG_OPTS, NULL)) != -1) {
        switch (opt) {
            case OPT_NLOCI: {
                char *end = NULL;
                long v = strtol(optarg, &end, 10);
                if (!end || *end || v < 0) {
                    fprintf(stderr, "Error: --nloci needs a non-negative integer, got '%s'\n", optarg);
                    return 1;
                }
                expected = (int)v;
                break;
            }
            case OPT_JSON: json_mode = 1; break;
            case 'h':      print_help(stdout); return 0;
            default:       print_help(stderr); return 1;
        }
    }
    if (optind != argc - 1) {
        fprintf(stderr, "Error: check takes exactly one FILE\n\n");
        print_help(stderr);
        return 1;
    }
    const char *path = argv[optind];

    FileType ft = detect_file_type(path);
    Report r;
    memset(&r, 0, sizeof(r));
    int rc;
    switch (ft) {
        case BS_PHYLIP:          rc = check_phylip(path, &r); break;
        case BS_FASTA_MSA:
        case BS_FASTA_CONTIGS:
        case BS_FASTA_REFERENCE: rc = check_fasta(path, &r);  break;
        case BS_NEXUS:           rc = check_nexus(path, &r);  break;
        default:
            fprintf(stderr, "Error: '%s' is %s, not a sequence alignment "
                    "(BPP/PHYLIP, FASTA or NEXUS)\n", path,
                    ft == BS_UNKNOWN ? "of an unrecognised or unreadable type"
                                     : file_type_name(ft));
            return 1;
    }
    if (rc != 0) return 1;

    if (r.n_loci == 0 && r.n_errs == 0)
        add_err(&r, 0, 0, "NO_LOCI", "no loci found");
    if (expected >= 0 && r.n_loci != expected)
        add_err(&r, 0, 0, "NLOCI", "file has %d loc%s; expected %d",
                r.n_loci, r.n_loci == 1 ? "us" : "i", expected);

    if (json_mode) print_json(path, file_type_name(ft), &r, expected);
    else           print_text(path, &r, expected);

    int failed = r.n_errs > 0;
    for (int i = 0; i < r.n_errs; i++) free(r.errs[i].msg);
    free(r.errs);
    free(r.loci);
    return failed ? 2 : 0;
}
