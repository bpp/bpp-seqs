#define _POSIX_C_SOURCE 200809L

#include "cmd_mask.h"
#include "bpp_parser.h"
#include "json_writer.h"
#include "sanity.h"

#include <ctype.h>
#include <getopt.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

static void print_help(FILE *fp)
{
    fprintf(fp,
"Usage: bpp-seqs mask INPUT.txt --masks FILE --out PREFIX [options]\n"
"\n"
"Replace positions outside a sample's mappability mask with 'N' in an\n"
"existing BPP sequence file. Short reads cannot be placed uniquely across\n"
"repeats and segmental duplications; reads from a paralogous locus map to\n"
"the wrong place and carry the differences between the paralogs with them,\n"
"which reads as a heterozygote no base caller can reject -- the mapper made\n"
"the error, not the sequencer. Masking removes those positions.\n"
"\n"
"Masked per site, not per locus. Masks are per sample (each genome has its\n"
"own coverage and so its own callable set), and loci are shared across\n"
"samples, so dropping whole loci could not represent them; writing 'N' lets\n"
"one sample be missing where another has data. BPP reads 'N' as missing for\n"
"that sequence at that site provided the control file sets cleandata = 0.\n"
"\n"
"Required:\n"
"  INPUT.txt                BPP sequence file produced by `bpp-seqs`\n"
"  --masks FILE             two columns, <sample> <TAB> <path to BED>\n"
"  --out PREFIX             output file prefix\n"
"\n"
"The BED files list positions that PASS (the convention the 1000 Genomes\n"
"accessibility masks and the archaic genome filter beds are distributed in).\n"
"Anything outside them becomes 'N'.\n"
"\n"
"Only samples named in --masks are touched; every other sequence is copied\n"
"through unchanged. That is how to leave panel-phased samples alone, whose\n"
"bases came from a VCF rather than from read mapping.\n"
"\n"
"Options:\n"
"  --max-missing FLOAT      after masking, drop loci whose missing fraction\n"
"                           exceeds this. Off by default: masking is applied\n"
"                           after the conversion already ran its QC, so\n"
"                           without this the surviving locus set no longer\n"
"                           means what the conversion reported.\n"
"  --loci-tsv FILE          use this provenance instead of <INPUT>.loci.tsv\n"
"  --imap FILE              use this Imap instead of <INPUT>.imap\n"
"  --json                   emit a JSON summary on stdout\n"
"  --quiet                  suppress stderr progress\n"
"  -h, --help               this message\n"
"\n"
"Outputs:\n"
"  PREFIX.txt               the masked BPP sequence file\n"
"  PREFIX.loci.tsv          provenance for the kept loci (if input had one)\n"
"  PREFIX.imap              copied from the input's Imap (if one was found)\n"
"\n"
"A <INPUT>.loci.tsv is required: without the source coordinates of each\n"
"locus there is no way to know which genome positions its columns describe.\n"
"\n"
"Example\n"
"  bpp-seqs mask run1.txt --masks masks.tsv --max-missing 0.5 --out run1.masked\n");
}

/* ────────────────────────────────────────────────────────────────────────
 * Mask manifest and interval sets
 * ───────────────────────────────────────────────────────────────────── */

typedef struct { int64_t beg, end; } Iv;      /* 0-based half-open */

typedef struct {
    char   *chrom;
    Iv     *iv;
    int     n, cap;
} ChromIv;

typedef struct {
    char    *sample;
    char    *path;
    ChromIv *chroms;
    int      n_chroms;
    int64_t  n_pass;      /* total bp marked passing, for reporting */
} SampleMask;

static int iv_cmp(const void *a, const void *b)
{
    const Iv *x = (const Iv *)a, *y = (const Iv *)b;
    if (x->beg < y->beg) return -1;
    if (x->beg > y->beg) return 1;
    return 0;
}

static ChromIv *mask_chrom(SampleMask *m, const char *chrom, int create)
{
    for (int i = 0; i < m->n_chroms; i++)
        if (strcmp(m->chroms[i].chrom, chrom) == 0) return &m->chroms[i];
    if (!create) return NULL;
    m->chroms = (ChromIv *)realloc(m->chroms,
                                   sizeof(ChromIv) * (size_t)(m->n_chroms + 1));
    ChromIv *c = &m->chroms[m->n_chroms++];
    memset(c, 0, sizeof(*c));
    c->chrom = strdup(chrom);
    return c;
}

/* Contig naming differs between sources ("chr22" vs "22"); compare on the
 * part after any "chr" prefix so a mask and a reference need not agree. */
static const char *bare_chrom(const char *c)
{
    if (!c) return c;
    if (strncmp(c, "chr", 3) == 0 || strncmp(c, "CHR", 3) == 0) return c + 3;
    return c;
}

static int load_mask_bed(SampleMask *m)
{
    FILE *f = fopen(m->path, "r");
    if (!f) {
        fprintf(stderr, "Error: cannot open mask BED '%s' for sample '%s'\n",
                m->path, m->sample);
        return -1;
    }
    char line[8192];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char c[512];
        long long b, e;
        if (sscanf(line, "%511s %lld %lld", c, &b, &e) != 3) continue;
        if (e <= b) continue;
        ChromIv *ci = mask_chrom(m, bare_chrom(c), 1);
        if (ci->n >= ci->cap) {
            ci->cap = ci->cap ? ci->cap * 2 : 1024;
            ci->iv  = (Iv *)realloc(ci->iv, sizeof(Iv) * (size_t)ci->cap);
        }
        ci->iv[ci->n].beg = b;
        ci->iv[ci->n].end = e;
        ci->n++;
        m->n_pass += e - b;
    }
    fclose(f);
    for (int i = 0; i < m->n_chroms; i++)
        qsort(m->chroms[i].iv, (size_t)m->chroms[i].n, sizeof(Iv), iv_cmp);
    return 0;
}

/* Mark keep[0..len) for genome span [beg, beg+len) on `chrom`. Intervals are
 * sorted, so walk from the first that could overlap. */
static void mark_pass(const SampleMask *m, const char *chrom,
                      int64_t beg, int len, unsigned char *keep)
{
    ChromIv *ci = mask_chrom((SampleMask *)m, bare_chrom(chrom), 0);
    if (!ci) return;                       /* no mask for this contig: all N */
    int64_t end = beg + len;
    int lo = 0, hi = ci->n;
    while (lo < hi) {                      /* first interval with end > beg */
        int mid = (lo + hi) / 2;
        if (ci->iv[mid].end <= beg) lo = mid + 1; else hi = mid;
    }
    for (int i = lo; i < ci->n && ci->iv[i].beg < end; i++) {
        int64_t s = ci->iv[i].beg > beg ? ci->iv[i].beg : beg;
        int64_t e = ci->iv[i].end < end ? ci->iv[i].end : end;
        for (int64_t p = s; p < e; p++) keep[p - beg] = 1;
    }
}

static SampleMask *load_manifest(const char *path, int *n_out, int quiet)
{
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "Error: cannot open --masks '%s'\n", path); return NULL; }
    SampleMask *ms = NULL;
    int n = 0;
    char line[8192];
    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n') continue;
        char s[1024], p[4096];
        if (sscanf(line, "%1023s %4095s", s, p) != 2) continue;
        ms = (SampleMask *)realloc(ms, sizeof(SampleMask) * (size_t)(n + 1));
        memset(&ms[n], 0, sizeof(SampleMask));
        ms[n].sample = strdup(s);
        ms[n].path   = strdup(p);
        if (load_mask_bed(&ms[n]) != 0) { fclose(f); free(ms); return NULL; }
        if (!quiet)
            fprintf(stderr, "  mask %-24s %.2f Mb passing\n",
                    ms[n].sample, (double)ms[n].n_pass / 1e6);
        n++;
    }
    fclose(f);
    if (n == 0) fprintf(stderr, "Error: --masks '%s' lists no samples\n", path);
    *n_out = n;
    return ms;
}

static void masks_free(SampleMask *ms, int n)
{
    for (int i = 0; i < n; i++) {
        for (int c = 0; c < ms[i].n_chroms; c++) {
            free(ms[i].chroms[c].chrom);
            free(ms[i].chroms[c].iv);
        }
        free(ms[i].chroms);
        free(ms[i].sample);
        free(ms[i].path);
    }
    free(ms);
}

/* Which mask applies to a sequence, or NULL.
 *
 * Sequence ids are the part after the last '^'. Under --phasing split/vcf a
 * diploid sample becomes "<sample>_1" and "<sample>_2", so an exact miss is
 * retried against the id with that suffix removed -- one mask covers both
 * haplotypes of an individual. */
static SampleMask *mask_for(SampleMask *ms, int n, const char *seq_name)
{
    const char *id = sanity_seq_id(seq_name);
    for (int i = 0; i < n; i++)
        if (strcmp(ms[i].sample, id) == 0) return &ms[i];

    size_t L = strlen(id);
    if (L > 2 && id[L-2] == '_' && (id[L-1] == '1' || id[L-1] == '2')) {
        for (int i = 0; i < n; i++)
            if (strncmp(ms[i].sample, id, L - 2) == 0 &&
                ms[i].sample[L-2] == '\0') return &ms[i];
    }
    return NULL;
}

/* ────────────────────────────────────────────────────────────────────────
 * I/O helpers (mirroring cmd_extract.c)
 * ───────────────────────────────────────────────────────────────────── */

static int file_exists(const char *p) { struct stat st; return p && stat(p, &st) == 0; }

static char *replace_ext(const char *path, const char *suffix)
{
    size_t n = strlen(path);
    const char *dot = strrchr(path, '.');
    size_t base = (dot && strcmp(dot, ".txt") == 0) ? (size_t)(dot - path) : n;
    char *r = (char *)malloc(base + strlen(suffix) + 1);
    memcpy(r, path, base);
    strcpy(r + base, suffix);
    return r;
}

static int copy_file(const char *src, const char *dst)
{
    FILE *in = fopen(src, "rb");
    if (!in) return -1;
    FILE *out = fopen(dst, "wb");
    if (!out) { fclose(in); return -1; }
    char buf[8192]; size_t got;
    while ((got = fread(buf, 1, sizeof(buf), in)) > 0) fwrite(buf, 1, got, out);
    fclose(in); fclose(out);
    return 0;
}

static double missing_frac(const BppLocus *l)
{
    int64_t tot = 0, miss = 0;
    for (int i = 0; i < l->n_seqs; i++)
        for (int j = 0; j < l->n_sites; j++) {
            tot++;
            char c = l->seqs[i][j];
            if (c == 'N' || c == 'n' || c == '?' || c == '-') miss++;
        }
    return tot ? (double)miss / (double)tot : 0.0;
}

/* ────────────────────────────────────────────────────────────────────────
 * main
 * ───────────────────────────────────────────────────────────────────── */

enum { OPT_MASKS = 256, OPT_OUT, OPT_MAXMISS, OPT_LOCI_TSV, OPT_IMAP,
       OPT_JSON, OPT_QUIET };

static const struct option LONG_OPTS[] = {
    {"masks",       required_argument, NULL, OPT_MASKS},
    {"out",         required_argument, NULL, OPT_OUT},
    {"max-missing", required_argument, NULL, OPT_MAXMISS},
    {"loci-tsv",    required_argument, NULL, OPT_LOCI_TSV},
    {"imap",        required_argument, NULL, OPT_IMAP},
    {"json",        no_argument,       NULL, OPT_JSON},
    {"quiet",       no_argument,       NULL, OPT_QUIET},
    {"help",        no_argument,       NULL, 'h'},
    {NULL, 0, NULL, 0}
};

int cmd_mask(int argc, char **argv)
{
    char  *masks_path = NULL, *out_prefix = NULL;
    char  *loci_tsv = NULL, *imap_override = NULL;
    double max_missing = -1.0;          /* <0 = do not re-filter */
    int    json_mode = 0, quiet = 0, opt;

    while ((opt = getopt_long(argc, argv, "h", LONG_OPTS, NULL)) != -1) {
        switch (opt) {
            case OPT_MASKS:    masks_path    = strdup(optarg); break;
            case OPT_OUT:      out_prefix    = strdup(optarg); break;
            case OPT_MAXMISS:  max_missing   = atof(optarg);   break;
            case OPT_LOCI_TSV: loci_tsv      = strdup(optarg); break;
            case OPT_IMAP:     imap_override = strdup(optarg); break;
            case OPT_JSON:     json_mode = 1; break;
            case OPT_QUIET:    quiet = 1; break;
            case 'h':          print_help(stdout); return 0;
            default:           print_help(stderr); return 1;
        }
    }
    if (optind >= argc) {
        fprintf(stderr, "Error: no input BPP sequence file.\n");
        print_help(stderr); return 1;
    }
    const char *input = argv[optind];
    if (!masks_path || !out_prefix) {
        fprintf(stderr, "Error: --masks and --out are both required.\n");
        print_help(stderr); return 1;
    }

    int n_loci = 0;
    BppLocus *loci = bpp_parse_file(input, &n_loci);
    if (!loci || n_loci == 0) {
        fprintf(stderr, "Error: no loci parsed from '%s'\n", input);
        return 1;
    }

    char *tsv = loci_tsv ? strdup(loci_tsv) : replace_ext(input, ".loci.tsv");
    if (!file_exists(tsv)) {
        fprintf(stderr,
            "Error: no locus provenance found ('%s').\n"
            "       Masking needs each locus's source coordinates to know which\n"
            "       genome positions its columns describe. Supply --loci-tsv.\n",
            tsv);
        free(tsv); bpp_loci_free(loci, n_loci); return 1;
    }
    bpp_attach_loci_tsv(loci, n_loci, tsv);

    int n_masks = 0;
    SampleMask *ms = load_manifest(masks_path, &n_masks, quiet);
    if (!ms) { free(tsv); bpp_loci_free(loci, n_loci); return 1; }

    /* ── apply ─────────────────────────────────────────────────────────── */

    int64_t n_masked = 0, n_considered = 0;
    int n_seq_masked = 0, n_seq_untouched = 0, n_no_coords = 0;
    unsigned char *keep = NULL;
    int keep_cap = 0;

    for (int i = 0; i < n_loci; i++) {
        BppLocus *l = &loci[i];
        if (!l->source_chrom || l->source_start < 0) { n_no_coords++; continue; }
        /* .loci.tsv records source_start 1-based; masks are 0-based half-open. */
        int64_t beg = (int64_t)l->source_start - 1;

        if (l->n_sites > keep_cap) {
            keep_cap = l->n_sites;
            keep = (unsigned char *)realloc(keep, (size_t)keep_cap);
        }

        for (int s = 0; s < l->n_seqs; s++) {
            SampleMask *m = mask_for(ms, n_masks, l->seq_names[s]);
            if (!m) { if (i == 0) n_seq_untouched++; continue; }
            if (i == 0) n_seq_masked++;

            memset(keep, 0, (size_t)l->n_sites);
            mark_pass(m, l->source_chrom, beg, l->n_sites, keep);
            for (int j = 0; j < l->n_sites; j++) {
                n_considered++;
                if (keep[j]) continue;
                if (l->seqs[s][j] != 'N') n_masked++;
                l->seqs[s][j] = 'N';
            }
        }
    }
    free(keep);

    /* ── optional re-filter ────────────────────────────────────────────── */

    int *pass = (int *)malloc(sizeof(int) * (size_t)n_loci);
    int n_kept = 0;
    for (int i = 0; i < n_loci; i++) {
        pass[i] = (max_missing < 0.0) || (missing_frac(&loci[i]) <= max_missing);
        if (pass[i]) n_kept++;
    }

    if (n_kept == 0) {
        fprintf(stderr, "Error: no loci survive --max-missing %.3f after masking.\n",
                max_missing);
        free(pass); free(tsv); masks_free(ms, n_masks);
        bpp_loci_free(loci, n_loci);
        return 1;
    }

    /* ── write ─────────────────────────────────────────────────────────── */

    char path[4096];
    snprintf(path, sizeof(path), "%s.txt", out_prefix);
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); free(pass); return 1; }
    int first = 1;
    for (int i = 0; i < n_loci; i++) {
        if (!pass[i]) continue;
        BppLocus *l = &loci[i];
        if (!first) fputc('\n', f);
        first = 0;
        fprintf(f, "%d %d\n\n", l->n_seqs, l->n_sites);
        int w = 0;
        for (int s = 0; s < l->n_seqs; s++) {
            int len = (int)strlen(l->seq_names[s]);
            if (len > w) w = len;
        }
        for (int s = 0; s < l->n_seqs; s++) {
            int wrote = fprintf(f, "%s", l->seq_names[s]);
            int pad = w + 4 - wrote; if (pad < 1) pad = 1;
            for (int p = 0; p < pad; p++) fputc(' ', f);
            fprintf(f, "%s\n", l->seqs[s]);
        }
    }
    fclose(f);

    /* Sidecars: Imap through unchanged, provenance filtered to kept loci. */
    char *in_imap = imap_override ? strdup(imap_override) : replace_ext(input, ".imap");
    char out_imap[4096];
    snprintf(out_imap, sizeof(out_imap), "%s.imap", out_prefix);
    int have_imap = file_exists(in_imap) && copy_file(in_imap, out_imap) == 0;

    char out_tsv[4096];
    snprintf(out_tsv, sizeof(out_tsv), "%s.loci.tsv", out_prefix);
    FILE *ti = fopen(tsv, "r"), *to = fopen(out_tsv, "w");
    int have_tsv = 0;
    if (ti && to) {
        char line[8192];
        if (fgets(line, sizeof(line), ti)) fputs(line, to);   /* header */
        for (int i = 0; i < n_loci && fgets(line, sizeof(line), ti); i++)
            if (pass[i]) fputs(line, to);
        have_tsv = 1;
    }
    if (ti) fclose(ti);
    if (to) fclose(to);

    if (!quiet) {
        fprintf(stderr,
            "\n  %d sequence%s masked, %d left unchanged (not in --masks)\n",
            n_seq_masked, n_seq_masked == 1 ? "" : "s", n_seq_untouched);
        fprintf(stderr,
            "  %" PRId64 " of %" PRId64 " masked-sample positions set to N (%.1f%%)\n",
            n_masked, n_considered,
            n_considered ? 100.0 * (double)n_masked / (double)n_considered : 0.0);
        if (n_no_coords)
            fprintf(stderr,
                "  Warning: %d loci had no source coordinates and were left unmasked.\n",
                n_no_coords);
        if (max_missing >= 0.0)
            fprintf(stderr, "  %d / %d loci kept at --max-missing %.3f\n",
                    n_kept, n_loci, max_missing);
        else
            fprintf(stderr,
                "  %d loci written, none re-filtered. Masking ran after the\n"
                "  conversion's own QC, so pass --max-missing to re-apply it.\n",
                n_kept);
        fprintf(stderr, "\n  %s.txt\n", out_prefix);
        if (have_imap) fprintf(stderr, "  %s.imap\n", out_prefix);
        if (have_tsv)  fprintf(stderr, "  %s.loci.tsv\n", out_prefix);
    }

    if (json_mode) {
        JsonWriter w; jw_init(&w, stdout, 2);
        jw_obj_open(&w);
        jw_kv_str (&w, "command", "mask");
        jw_kv_str (&w, "input", input);
        jw_kv_int (&w, "n_loci_in", n_loci);
        jw_kv_int (&w, "n_loci_out", n_kept);
        jw_kv_int (&w, "n_sequences_masked", n_seq_masked);
        jw_kv_int (&w, "n_sequences_untouched", n_seq_untouched);
        jw_kv_int (&w, "n_positions_masked", (int)n_masked);
        jw_kv_int (&w, "n_positions_considered", (int)n_considered);
        jw_kv_str (&w, "sequences", path);
        jw_obj_close(&w);
        fputc('\n', stdout);
    }

    free(pass); free(tsv); free(in_imap);
    free(masks_path); free(out_prefix); free(loci_tsv); free(imap_override);
    masks_free(ms, n_masks);
    bpp_loci_free(loci, n_loci);
    return 0;
}
