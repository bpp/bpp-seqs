/* inspect.c — content-based file type detection and per-type inspection. */

#define _POSIX_C_SOURCE 200809L

#include "inspect.h"

#include <ctype.h>
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <zlib.h>

#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include <htslib/kstring.h>

/* ────────────────────────────────────────────────────────────────────────
 * small utilities
 * ───────────────────────────────────────────────────────────────────── */

static char *xstrdup(const char *s)
{
    if (!s) return NULL;
    size_t n = strlen(s);
    char *r = (char *)malloc(n + 1);
    if (!r) return NULL;
    memcpy(r, s, n + 1);
    return r;
}

static int starts_with(const char *s, const char *prefix)
{
    return strncmp(s, prefix, strlen(prefix)) == 0;
}

static int istarts_with(const char *s, const char *prefix)
{
    size_t n = strlen(prefix);
    for (size_t i = 0; i < n; i++) {
        if (!s[i]) return 0;
        if (tolower((unsigned char)s[i]) != tolower((unsigned char)prefix[i])) return 0;
    }
    return 1;
}

static char *rtrim(char *s)
{
    if (!s) return s;
    size_t n = strlen(s);
    while (n > 0 && (s[n-1] == '\n' || s[n-1] == '\r' || s[n-1] == ' ' || s[n-1] == '\t')) {
        s[--n] = '\0';
    }
    return s;
}

static int file_exists(const char *path)
{
    struct stat st;
    return path && stat(path, &st) == 0;
}

static char *path_with_suffix(const char *path, const char *suffix)
{
    size_t a = strlen(path), b = strlen(suffix);
    char *r = (char *)malloc(a + b + 1);
    if (!r) return NULL;
    memcpy(r, path, a); memcpy(r + a, suffix, b); r[a + b] = '\0';
    return r;
}

static int has_companion_index(const char *path, const char *suf1, const char *suf2)
{
    int ok = 0;
    char *p = path_with_suffix(path, suf1);
    if (p && file_exists(p)) ok = 1;
    free(p);
    if (ok) return 1;
    if (suf2) {
        p = path_with_suffix(path, suf2);
        if (p && file_exists(p)) ok = 1;
        free(p);
    }
    return ok;
}

/* ────────────────────────────────────────────────────────────────────────
 * FileInfo lifecycle
 * ───────────────────────────────────────────────────────────────────── */

const char *file_type_name(FileType t)
{
    switch (t) {
        case BS_BAM:             return "BAM";
        case BS_CRAM:            return "CRAM";
        case BS_FASTQ:           return "FASTQ";
        case BS_FASTA_REFERENCE: return "FASTA_REFERENCE";
        case BS_FASTA_MSA:       return "FASTA_MSA";
        case BS_FASTA_CONTIGS:   return "FASTA_CONTIGS";
        case BS_VCF:             return "VCF";
        case BS_GVCF:            return "GVCF";
        case BS_BED:             return "BED";
        case BS_IMAP:            return "IMAP";
        case BS_PHYLIP:          return "PHYLIP";
        case BS_NEXUS:           return "NEXUS";
        case BS_UNKNOWN:
        default:                 return "UNKNOWN";
    }
}

static FileInfo *file_info_new(const char *path)
{
    FileInfo *fi = (FileInfo *)calloc(1, sizeof(FileInfo));
    if (!fi) return NULL;
    fi->path = xstrdup(path);
    fi->ft = BS_UNKNOWN;
    return fi;
}

void file_info_add_warning(FileInfo *fi, const char *code, const char *severity, const char *message)
{
    Warning *w = (Warning *)realloc(fi->warnings, sizeof(Warning) * (fi->n_warnings + 1));
    if (!w) return;
    fi->warnings = w;
    fi->warnings[fi->n_warnings].code     = xstrdup(code);
    fi->warnings[fi->n_warnings].severity = xstrdup(severity);
    fi->warnings[fi->n_warnings].message  = xstrdup(message);
    fi->n_warnings++;
}

static void free_str_array(char **a, int n)
{
    if (!a) return;
    for (int i = 0; i < n; i++) free(a[i]);
    free(a);
}

void file_info_free(FileInfo *fi)
{
    if (!fi) return;
    free(fi->path);
    for (int i = 0; i < fi->n_warnings; i++) {
        free(fi->warnings[i].code);
        free(fi->warnings[i].severity);
        free(fi->warnings[i].message);
    }
    free(fi->warnings);

    free(fi->aligner);
    free(fi->aligner_version);
    free(fi->aligner_command);
    free(fi->sample_name);
    free(fi->platform);
    free(fi->read_type);
    if (fi->seq_refs) {
        for (int i = 0; i < fi->n_seq_refs; i++) free(fi->seq_refs[i].name);
        free(fi->seq_refs);
    }
    free(fi->reference_path_in_header);
    free(fi->platform_guess);
    free(fi->recommended_aligner);
    free_str_array(fi->sequence_names, fi->n_sequence_names);
    free_str_array(fi->sample_names, fi->n_sample_names);
    free(fi->vcf_reference_in_header);
    if (fi->vcf_contigs) {
        for (int i = 0; i < fi->n_vcf_contigs; i++) free(fi->vcf_contigs[i].name);
        free(fi->vcf_contigs);
    }
    free_str_array(fi->chromosomes, fi->n_chromosomes);
    free_str_array(fi->imap_sample_names, fi->imap_n_sample_names);
    free_str_array(fi->imap_population_names, fi->imap_n_population_names);
    free_str_array(fi->imap_samples, fi->imap_n_entries);
    free_str_array(fi->imap_pops,    fi->imap_n_entries);
    free(fi->phylip_format);
    free_str_array(fi->nexus_charset_names, fi->nexus_n_charset_names);
    free(fi->nexus_charset_starts);
    free(fi->nexus_charset_ends);
    free(fi->nexus_charset_strides);
    free(fi);
}

/* ────────────────────────────────────────────────────────────────────────
 * Content detection
 *
 * Strategy: read first ~64 bytes raw. If gzip magic, read the first ~4 KB
 * decompressed via zlib and detect on that. Otherwise detect on the raw
 * head buffer.
 * ───────────────────────────────────────────────────────────────────── */

static int read_head_raw(const char *path, unsigned char *buf, size_t n)
{
    FILE *fp = fopen(path, "rb");
    if (!fp) return -1;
    size_t got = fread(buf, 1, n, fp);
    fclose(fp);
    return (int)got;
}

static int read_head_text(const char *path, char *buf, size_t cap)
{
    /* gzopen handles both plain and gzipped transparently. */
    gzFile gz = gzopen(path, "rb");
    if (!gz) return -1;
    int got = gzread(gz, buf, (unsigned)(cap - 1));
    gzclose(gz);
    if (got < 0) return -1;
    buf[got] = '\0';
    return got;
}

static int is_gzip(const unsigned char *buf, int n)
{
    return n >= 2 && buf[0] == 0x1f && buf[1] == 0x8b;
}

/* CRAM magic = "CRAM" 4 bytes */
static int is_cram(const unsigned char *buf, int n)
{
    return n >= 4 && buf[0] == 'C' && buf[1] == 'R' && buf[2] == 'A' && buf[3] == 'M';
}

/* Check that line begins with two whitespace-separated non-negative ints. */
static int line_is_two_ints(const char *line)
{
    const char *p = line;
    while (*p == ' ' || *p == '\t') p++;
    if (!isdigit((unsigned char)*p)) return 0;
    while (isdigit((unsigned char)*p)) p++;
    if (*p != ' ' && *p != '\t') return 0;
    while (*p == ' ' || *p == '\t') p++;
    if (!isdigit((unsigned char)*p)) return 0;
    while (isdigit((unsigned char)*p)) p++;
    while (*p == ' ' || *p == '\t') p++;
    /* Trailing junk → still accept (PHYLIP variants sometimes have extra) */
    return 1;
}

/* Three+ tab-separated fields, fields 2 and 3 are integers, field 1 is a
 * non-empty token. */
static int line_is_bed(const char *line)
{
    const char *p = line;
    if (!*p || *p == '#' || *p == '\n' || *p == '\r') return 0;
    /* field 1 */
    int seen_tab = 0;
    while (*p && *p != '\t' && *p != '\n' && *p != '\r') p++;
    if (*p != '\t') return 0;
    seen_tab++;
    p++;
    /* field 2 must be int */
    if (!isdigit((unsigned char)*p)) return 0;
    while (isdigit((unsigned char)*p)) p++;
    if (*p != '\t') return 0;
    seen_tab++;
    p++;
    /* field 3 must be int */
    if (!isdigit((unsigned char)*p)) return 0;
    while (isdigit((unsigned char)*p)) p++;
    if (*p != '\t' && *p != '\n' && *p != '\r' && *p != '\0') return 0;
    (void)seen_tab;
    return 1;
}

/* Two non-empty fields per line separated by whitespace (tab or spaces). */
static int line_is_imap(const char *line)
{
    const char *p = line;
    if (!*p || *p == '#' || *p == '\n' || *p == '\r') return 0;
    while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') p++;
    if (*p != ' ' && *p != '\t') return 0;
    while (*p == ' ' || *p == '\t') p++;
    if (!*p || *p == '\n' || *p == '\r') return 0;
    /* must have something here */
    while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') p++;
    /* trailing whitespace OK; another field would still be ambiguous, but
     * BED is checked first so a 3-tab-field line goes to BED. */
    return 1;
}

/* Heuristic: looks like a nucleotide sequence (ACGTUN-) */
static int line_looks_like_sequence(const char *s)
{
    int n = 0, m = 0;
    while (*s && *s != '\n' && *s != '\r' && n < 200) {
        char c = (char)toupper((unsigned char)*s);
        if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'U' || c == 'N' || c == '-') m++;
        n++;
        s++;
    }
    if (n == 0) return 0;
    return m * 4 >= n * 3;  /* ≥75% nucleotide characters */
}

static FileType detect_type_from_text(const char *head, int n)
{
    /* Look for VCF header */
    if (starts_with(head, "##fileformat=VCFv")) {
        return BS_VCF;  /* may be promoted to GVCF in inspector */
    }
    /* NEXUS */
    if (istarts_with(head, "#NEXUS") || istarts_with(head, "#nexus")) {
        return BS_NEXUS;
    }
    /* SAM (text BAM equivalent) — route to BAM inspector */
    if (n >= 3 && head[0] == '@' &&
        (starts_with(head, "@HD") || starts_with(head, "@SQ") ||
         starts_with(head, "@RG") || starts_with(head, "@PG") ||
         starts_with(head, "@CO"))) {
        return BS_BAM;
    }
    /* FASTQ: first line starts with '@' (any token), second line is sequence */
    if (n >= 2 && head[0] == '@') {
        const char *nl = strchr(head, '\n');
        if (nl && line_looks_like_sequence(nl + 1)) return BS_FASTQ;
    }
    /* FASTA */
    if (n >= 1 && head[0] == '>') return BS_FASTA_MSA;  /* refined later */
    /* PHYLIP: first line is two ints */
    {
        char line1[512];
        int i = 0;
        while (i < n && i < (int)sizeof(line1) - 1 && head[i] != '\n' && head[i] != '\r') {
            line1[i] = head[i]; i++;
        }
        line1[i] = '\0';
        if (line_is_two_ints(line1) && i < n) {
            /* must have a follow-up line with non-empty content */
            return BS_PHYLIP;
        }
    }
    /* BED before IMAP (BED is more specific) */
    {
        char line1[512];
        int i = 0;
        while (i < n && i < (int)sizeof(line1) - 1 && head[i] != '\n' && head[i] != '\r') {
            line1[i] = head[i]; i++;
        }
        line1[i] = '\0';
        if (line_is_bed(line1))  return BS_BED;
        if (line_is_imap(line1)) return BS_IMAP;
    }
    return BS_UNKNOWN;
}

static FileType detect_type(const char *path)
{
    unsigned char raw[64];
    int nr = read_head_raw(path, raw, sizeof(raw));
    if (nr <= 0) return BS_UNKNOWN;

    if (is_cram(raw, nr)) return BS_CRAM;

    if (is_gzip(raw, nr)) {
        /* Read first 4 KB decompressed. Inside, BAM has magic "BAM\1". */
        char buf[4096];
        int n = read_head_text(path, buf, sizeof(buf));
        if (n >= 4 && buf[0] == 'B' && buf[1] == 'A' && buf[2] == 'M' && buf[3] == 1) {
            return BS_BAM;
        }
        if (n > 0) return detect_type_from_text(buf, n);
        return BS_UNKNOWN;
    }

    /* Plain text: read more as text */
    char buf[4096];
    int n = read_head_text(path, buf, sizeof(buf));
    if (n <= 0) return BS_UNKNOWN;
    return detect_type_from_text(buf, n);
}

/* ────────────────────────────────────────────────────────────────────────
 * BAM / CRAM inspector
 * ───────────────────────────────────────────────────────────────────── */

/* Find and copy the value of "TAG:" up to next tab/newline within a header
 * line text. Returns a freshly allocated string or NULL. */
static char *header_tag_value(const char *line, const char *tag)
{
    /* tag is e.g. "SM:". Locate at start-of-token (preceded by tab). */
    size_t tlen = strlen(tag);
    const char *p = line;
    while ((p = strstr(p, tag)) != NULL) {
        /* Must be preceded by '\t' or be the first non-prefix token */
        if (p == line || *(p - 1) == '\t') {
            p += tlen;
            const char *end = p;
            while (*end && *end != '\t' && *end != '\n' && *end != '\r') end++;
            size_t n = (size_t)(end - p);
            char *r = (char *)malloc(n + 1);
            if (!r) return NULL;
            memcpy(r, p, n); r[n] = '\0';
            return r;
        }
        p += tlen;
    }
    return NULL;
}

static void inspect_bam_header(sam_hdr_t *hdr, FileInfo *fi, int is_cram)
{
    (void)is_cram;
    const char *htext = sam_hdr_str(hdr);
    if (!htext) return;

    /* iterate lines */
    const char *p = htext;
    int first_rg_sample_seen = 0;
    char *first_rg_sample = NULL;
    int multi_sm = 0;
    int rg_seen = 0;
    int first_pg_seen = 0;

    while (*p) {
        const char *nl = strchr(p, '\n');
        size_t llen = nl ? (size_t)(nl - p) : strlen(p);
        char *line = (char *)malloc(llen + 1);
        if (!line) break;
        memcpy(line, p, llen); line[llen] = '\0';

        if (starts_with(line, "@RG")) {
            rg_seen++;
            char *sm = header_tag_value(line, "SM:");
            char *pl = header_tag_value(line, "PL:");
            if (sm) {
                if (!first_rg_sample_seen) {
                    first_rg_sample = xstrdup(sm);
                    first_rg_sample_seen = 1;
                } else if (first_rg_sample && strcmp(first_rg_sample, sm) != 0) {
                    multi_sm = 1;
                }
                free(sm);
            }
            if (pl && !fi->platform) { fi->platform = pl; } else { free(pl); }
        } else if (starts_with(line, "@PG") && !first_pg_seen) {
            char *pn = header_tag_value(line, "PN:");
            char *vn = header_tag_value(line, "VN:");
            char *cl = header_tag_value(line, "CL:");
            if (pn) { free(fi->aligner);         fi->aligner = pn; first_pg_seen = 1; }
            if (vn) { free(fi->aligner_version); fi->aligner_version = vn; }
            if (cl) {
                free(fi->aligner_command); fi->aligner_command = cl;
                /* try to extract a reference path: heuristic — look for a
                 * token ending in .fa or .fasta or .fna or .fa.gz */
                const char *q = cl;
                while ((q = strchr(q, '.'))) {
                    if (strncmp(q, ".fa", 3) == 0 &&
                        (q[3] == '\0' || q[3] == ' ' || q[3] == '\t' ||
                         strncmp(q, ".fasta", 6) == 0 || strncmp(q, ".fna", 4) == 0 ||
                         strncmp(q, ".fa.gz", 6) == 0)) {
                        /* walk back to whitespace */
                        const char *start = q;
                        while (start > cl && *(start - 1) != ' ' && *(start - 1) != '\t') start--;
                        const char *end = q;
                        while (*end && *end != ' ' && *end != '\t') end++;
                        size_t n = (size_t)(end - start);
                        free(fi->reference_path_in_header);
                        fi->reference_path_in_header = (char *)malloc(n + 1);
                        if (fi->reference_path_in_header) {
                            memcpy(fi->reference_path_in_header, start, n);
                            fi->reference_path_in_header[n] = '\0';
                        }
                        break;
                    }
                    q++;
                }
            }
        }
        free(line);
        if (!nl) break;
        p = nl + 1;
    }

    if (!rg_seen) {
        file_info_add_warning(fi, "MISSING_RG", "warning",
            "No @RG lines in BAM header; sample name unavailable.");
    }
    if (multi_sm) {
        file_info_add_warning(fi, "MULTIPLE_RG_SAMPLES", "warning",
            "Multiple @RG SM values present; using the first.");
    }
    fi->sample_name = first_rg_sample;

    /* @SQ list via htslib targets */
    int nt = sam_hdr_nref(hdr);
    if (nt > 0) {
        fi->seq_refs = (SeqRef *)calloc((size_t)nt, sizeof(SeqRef));
        fi->n_seq_refs = nt;
        for (int i = 0; i < nt; i++) {
            const char *name = sam_hdr_tid2name(hdr, i);
            int64_t   len  = sam_hdr_tid2len (hdr, i);
            fi->seq_refs[i].name   = xstrdup(name ? name : "");
            fi->seq_refs[i].length = len;
        }
    }
}

static void inspect_bam_reads(samFile *sf, sam_hdr_t *hdr, FileInfo *fi)
{
    bam1_t *b = bam_init1();
    int n_sampled = 0;
    int n_paired = 0;
    int n_mapped = 0;
    int64_t total_len = 0;
    const int LIMIT = 10000;

    int rc;
    while (n_sampled < LIMIT && (rc = sam_read1(sf, hdr, b)) >= 0) {
        n_sampled++;
        int flag = b->core.flag;
        if (flag & BAM_FPAIRED) n_paired++;
        if (!(flag & BAM_FUNMAP)) {
            n_mapped++;
        }
        total_len += b->core.l_qseq;
    }
    bam_destroy1(b);

    fi->n_reads_sampled = n_sampled;
    fi->read_length_mean = n_sampled ? (int)(total_len / n_sampled) : 0;
    fi->mapping_rate = n_sampled ? (double)n_mapped / (double)n_sampled : 0.0;
    fi->read_type = xstrdup((n_paired * 2 > n_sampled) ? "paired" : "single");
    fi->unaligned = (n_sampled > 0 && n_mapped * 2 < n_sampled);

    /* depth estimate */
    int64_t total_ref = 0;
    for (int i = 0; i < fi->n_seq_refs; i++) total_ref += fi->seq_refs[i].length;
    if (total_ref > 0 && n_sampled > 0 && fi->read_length_mean > 0) {
        fi->mean_depth_estimate =
            (double)n_sampled * (double)fi->read_length_mean / (double)total_ref;
    }

    if (fi->unaligned) {
        file_info_add_warning(fi, "UNALIGNED_BAM", "warning",
            "More than 50% of sampled reads are unmapped. Treat this BAM as FASTQ-equivalent.");
    } else if (fi->mapping_rate > 0 && fi->mapping_rate < 0.90) {
        char msg[256];
        snprintf(msg, sizeof(msg),
            "Mapping rate %.1f%% is below 90%%. Check for species mismatch or contamination.",
            fi->mapping_rate * 100.0);
        file_info_add_warning(fi, "LOW_MAPPING_RATE", "warning", msg);
    }
}

static void inspect_bam(const char *path, FileInfo *fi, int is_cram)
{
    samFile *sf = sam_open(path, "r");
    if (!sf) {
        file_info_add_warning(fi, "OPEN_FAILED", "error", "Could not open BAM/CRAM file with htslib.");
        return;
    }
    sam_hdr_t *hdr = sam_hdr_read(sf);
    if (!hdr) {
        file_info_add_warning(fi, "HEADER_READ_FAILED", "error", "Could not read BAM/CRAM header.");
        sam_close(sf);
        return;
    }
    inspect_bam_header(hdr, fi, is_cram);
    inspect_bam_reads(sf, hdr, fi);

    if (is_cram) {
        fi->indexed = has_companion_index(path, ".crai", ".csi");
        if (!fi->indexed) {
            file_info_add_warning(fi, "MISSING_INDEX", "warning",
                "No .crai/.csi index file found alongside this CRAM.");
        }
    } else {
        fi->indexed = has_companion_index(path, ".bai", ".csi");
        if (!fi->indexed) {
            file_info_add_warning(fi, "MISSING_INDEX", "warning",
                "No .bai/.csi index file found alongside this BAM.");
        }
    }

    sam_hdr_destroy(hdr);
    sam_close(sf);
}

/* ────────────────────────────────────────────────────────────────────────
 * VCF / GVCF inspector
 * ───────────────────────────────────────────────────────────────────── */

static void inspect_vcf(const char *path, FileInfo *fi)
{
    htsFile *fp = bcf_open(path, "r");
    if (!fp) {
        file_info_add_warning(fi, "OPEN_FAILED", "error", "Could not open VCF file with htslib.");
        return;
    }
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) {
        file_info_add_warning(fi, "HEADER_READ_FAILED", "error", "Could not read VCF header.");
        hts_close(fp);
        return;
    }

    int ns = bcf_hdr_nsamples(hdr);
    fi->n_samples = ns;
    if (ns > 0) {
        fi->sample_names = (char **)calloc((size_t)ns, sizeof(char *));
        for (int i = 0; i < ns; i++) {
            fi->sample_names[i] = xstrdup(hdr->samples[i]);
        }
        fi->n_sample_names = ns;
    }

    /* extract reference from header lines */
    int nrec = hdr->nhrec;
    for (int i = 0; i < nrec; i++) {
        bcf_hrec_t *h = hdr->hrec[i];
        if (!h || !h->key) continue;
        if (strcmp(h->key, "reference") == 0 && h->value &&
            !fi->vcf_reference_in_header) {
            fi->vcf_reference_in_header = xstrdup(h->value);
        }
    }

    /* Use the parsed-contig list htslib maintains for us. */
    int n_ctg = 0;
    const char **ctg_names = bcf_hdr_seqnames(hdr, &n_ctg);
    if (ctg_names && n_ctg > 0) {
        fi->vcf_contigs = (SeqRef *)calloc((size_t)n_ctg, sizeof(SeqRef));
        fi->n_vcf_contigs = n_ctg;
        for (int i = 0; i < n_ctg; i++) {
            fi->vcf_contigs[i].name = xstrdup(ctg_names[i]);
            /* length lives under bcf_hdr_id2int/contig records; look it up */
            bcf_idpair_t *ip = hdr->id[BCF_DT_CTG];
            int64_t len = 0;
            if (ip) {
                int idx = bcf_hdr_name2id(hdr, ctg_names[i]);
                if (idx >= 0 && idx < hdr->n[BCF_DT_CTG]) {
                    len = (int64_t)hdr->id[BCF_DT_CTG][idx].val->info[0];
                }
            }
            fi->vcf_contigs[i].length = len;
        }
    }
    free(ctg_names);

    /* scan records */
    bcf1_t *rec = bcf_init();
    int n_phased = 0, n_gts = 0, n_sampled = 0;
    int is_gvcf = 0;
    const int LIMIT = 10000;
    int32_t *gt_arr = NULL;
    int gt_arr_n = 0;

    while (n_sampled < LIMIT && bcf_read(fp, hdr, rec) == 0) {
        n_sampled++;
        bcf_unpack(rec, BCF_UN_ALL);

        /* check ALT for <NON_REF> */
        for (int a = 1; a < rec->n_allele; a++) {
            if (rec->d.allele[a] && strcmp(rec->d.allele[a], "<NON_REF>") == 0) {
                is_gvcf = 1;
            }
        }
        /* check INFO for END (used in gVCF coverage bands too) */
        if (bcf_get_info_int32(hdr, rec, "END", &gt_arr, &gt_arr_n) > 0) {
            /* presence of END info doesn't itself mean gVCF, but the
             * <NON_REF> check above is the canonical marker. */
        }

        /* phase fraction */
        int ngt = bcf_get_genotypes(hdr, rec, &gt_arr, &gt_arr_n);
        if (ngt > 0 && ns > 0) {
            int ploidy = ngt / ns;
            for (int s = 0; s < ns; s++) {
                int32_t *p = gt_arr + s * ploidy;
                if (ploidy >= 2 && p[0] != bcf_int32_vector_end && p[1] != bcf_int32_vector_end) {
                    n_gts++;
                    if (bcf_gt_is_phased(p[1])) n_phased++;
                }
            }
        }
    }
    free(gt_arr);
    bcf_destroy(rec);

    /* Also scan unstructured header for GVCFBlock INFO lines, since they
     * mark gVCF even when no <NON_REF> ALTs were seen in the first 10k. */
    if (!is_gvcf) {
        for (int i = 0; i < nrec; i++) {
            bcf_hrec_t *h = hdr->hrec[i];
            if (!h) continue;
            for (int j = 0; j < h->nkeys; j++) {
                if (h->keys[j] && h->vals[j] &&
                    strstr(h->vals[j], "GVCFBlock") != NULL) {
                    is_gvcf = 1; break;
                }
            }
            if (is_gvcf) break;
        }
    }

    fi->ft = is_gvcf ? BS_GVCF : BS_VCF;
    fi->n_records_sampled = n_sampled;
    fi->phased_fraction   = n_gts ? (double)n_phased / (double)n_gts : 0.0;
    fi->has_phase_info    = (n_phased > 0);
    fi->has_coverage_bands = is_gvcf;

    if (!is_gvcf) {
        file_info_add_warning(fi, "VCF_CONSTANT_SITE_PROBLEM", "critical",
            "Standard VCF does not record coverage at non-variant sites. "
            "It is impossible to distinguish 'invariant and covered' from "
            "'not covered' at these positions. BPP requires this distinction "
            "for accurate inference. Consider using the BAM files this VCF "
            "was called from instead.");
    }

    bcf_hdr_destroy(hdr);
    hts_close(fp);
}

/* ────────────────────────────────────────────────────────────────────────
 * FASTA inspector  (REFERENCE / MSA / CONTIGS)
 * ───────────────────────────────────────────────────────────────────── */

static int int64_cmp(const void *a, const void *b)
{
    int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
    return (x > y) - (x < y);
}

static void inspect_fasta(const char *path, FileInfo *fi)
{
    gzFile gz = gzopen(path, "rb");
    if (!gz) {
        file_info_add_warning(fi, "OPEN_FAILED", "error", "Could not open FASTA file.");
        return;
    }

    char line[8192];
    int n_seqs = 0, name_cap = 16;
    char **names = (char **)malloc((size_t)name_cap * sizeof(char *));
    int64_t *lens = (int64_t *)malloc((size_t)name_cap * sizeof(int64_t));
    int len_cap = name_cap;
    int64_t total_len = 0;
    int has_gap = 0;
    int64_t missing_count = 0;
    int has_n_runs = 0;
    int cur_n_run = 0;

    int64_t cur_len = 0;

    while (gzgets(gz, line, sizeof(line)) != NULL) {
        if (line[0] == '>') {
            if (n_seqs > 0) {
                lens[n_seqs - 1] = cur_len;
                total_len += cur_len;
            }
            cur_len = 0;
            cur_n_run = 0;
            /* capture name */
            if (n_seqs >= name_cap) {
                name_cap *= 2;
                names = (char **)realloc(names, (size_t)name_cap * sizeof(char *));
                lens  = (int64_t *)realloc(lens, (size_t)name_cap * sizeof(int64_t));
                len_cap = name_cap;
            }
            char *p = line + 1;
            char *end = p;
            while (*end && *end != ' ' && *end != '\t' && *end != '\n' && *end != '\r') end++;
            size_t nlen = (size_t)(end - p);
            names[n_seqs] = (char *)malloc(nlen + 1);
            if (names[n_seqs]) { memcpy(names[n_seqs], p, nlen); names[n_seqs][nlen] = '\0'; }
            n_seqs++;
        } else {
            for (char *q = line; *q && *q != '\n' && *q != '\r'; q++) {
                char c = (char)toupper((unsigned char)*q);
                cur_len++;
                if (c == '-' || c == '.') { has_gap = 1; missing_count++; }
                else if (c == 'N' || c == '?') {
                    missing_count++;
                    cur_n_run++;
                    if (cur_n_run >= 10) has_n_runs = 1;
                } else {
                    cur_n_run = 0;
                }
            }
        }
    }
    if (n_seqs > 0) {
        lens[n_seqs - 1] = cur_len;
        total_len += cur_len;
    }
    gzclose(gz);

    fi->n_sequences  = n_seqs;
    fi->total_length = total_len;
    if (n_seqs > 0 && missing_count > 0 && total_len > 0) {
        fi->fasta_missing_fraction = (double)missing_count / (double)total_len;
    }
    fi->fasta_has_gaps = has_gap;

    /* Truncate sequence_names array to first 10 retained, but keep all
     * names in a single array for downstream consumers when needed. We
     * keep all names — the formatter truncates the display. */
    fi->sequence_names   = names;
    fi->n_sequence_names = n_seqs;
    (void)len_cap;

    /* Classification */
    if (n_seqs == 0) {
        fi->ft = BS_UNKNOWN;
        file_info_add_warning(fi, "EMPTY_FASTA", "error", "FASTA file contains no sequences.");
        free(lens);
        return;
    }
    int all_equal = 1;
    int64_t first_len = lens[0];
    int64_t min_len = lens[0], max_len = lens[0];
    for (int i = 1; i < n_seqs; i++) {
        if (lens[i] != first_len) all_equal = 0;
        if (lens[i] < min_len) min_len = lens[i];
        if (lens[i] > max_len) max_len = lens[i];
    }

    if (n_seqs == 1 && total_len > 1000000LL) {
        fi->ft = BS_FASTA_REFERENCE;
    } else if (n_seqs > 1 && all_equal) {
        fi->ft = BS_FASTA_MSA;
        fi->alignment_length = (int)first_len;
    } else if (n_seqs == 1) {
        /* Single short sequence — treat as REFERENCE if name implies a chr,
         * else as CONTIGS-of-one. We choose REFERENCE since BAMs may
         * reference small refs (e.g. tests). */
        fi->ft = BS_FASTA_REFERENCE;
    } else {
        fi->ft = BS_FASTA_CONTIGS;
        fi->length_min = min_len;
        fi->length_max = max_len;
        fi->length_mean = total_len / n_seqs;
        /* N50 */
        int64_t *sorted = (int64_t *)malloc(sizeof(int64_t) * (size_t)n_seqs);
        memcpy(sorted, lens, sizeof(int64_t) * (size_t)n_seqs);
        qsort(sorted, (size_t)n_seqs, sizeof(int64_t), int64_cmp);
        int64_t half = total_len / 2;
        int64_t acc = 0;
        int64_t n50 = 0;
        for (int i = n_seqs - 1; i >= 0; i--) {
            acc += sorted[i];
            if (acc >= half) { n50 = sorted[i]; break; }
        }
        free(sorted);
        fi->length_n50 = n50;
        fi->has_n_runs = has_n_runs;
        file_info_add_warning(fi, "ASSEMBLY_NOT_ALIGNABLE", "warning",
            "Assembled contigs cannot be used directly as BPP input. "
            "Loci must be defined and sequences aligned across individuals first.");
    }

    if (fi->ft == BS_FASTA_REFERENCE) {
        fi->indexed = has_companion_index(path, ".fai", NULL);
        if (!fi->indexed) {
            file_info_add_warning(fi, "MISSING_FAI", "warning",
                "No .fai index found. Run `samtools faidx` on this FASTA.");
        }
    }

    free(lens);
}

/* ────────────────────────────────────────────────────────────────────────
 * FASTQ inspector
 * ───────────────────────────────────────────────────────────────────── */

static void inspect_fastq(const char *path, FileInfo *fi)
{
    gzFile gz = gzopen(path, "rb");
    if (!gz) {
        file_info_add_warning(fi, "OPEN_FAILED", "error", "Could not open FASTQ file.");
        return;
    }
    char line[1 << 16];
    int n = 0;
    int64_t total = 0;
    int lmin = INT32_MAX, lmax = 0;
    const int LIMIT = 10000;
    while (n < LIMIT && gzgets(gz, line, sizeof(line)) != NULL) {
        if (line[0] != '@') break;  /* malformed */
        /* sequence on next line */
        if (gzgets(gz, line, sizeof(line)) == NULL) break;
        int slen = 0;
        for (char *p = line; *p && *p != '\n' && *p != '\r'; p++) slen++;
        total += slen;
        if (slen < lmin) lmin = slen;
        if (slen > lmax) lmax = slen;
        /* + line */
        if (gzgets(gz, line, sizeof(line)) == NULL) break;
        /* qual line */
        if (gzgets(gz, line, sizeof(line)) == NULL) break;
        n++;
    }
    gzclose(gz);
    fi->n_reads_sampled  = n;
    fi->read_length_mean = n ? (int)(total / n) : 0;
    fi->read_length_min  = (lmin == INT32_MAX) ? 0 : lmin;
    fi->read_length_max  = lmax;
    if (fi->read_length_mean > 0 && fi->read_length_mean <= 300) {
        fi->platform_guess     = xstrdup("ILLUMINA");
        fi->recommended_aligner = xstrdup("bwa mem");
    } else {
        fi->platform_guess     = xstrdup("PACBIO_HIFI or NANOPORE");
        fi->recommended_aligner = xstrdup("minimap2");
    }
    fi->read_type = xstrdup("unknown");
    file_info_add_warning(fi, "NEEDS_ALIGNMENT", "warning",
        "FASTQ input must be aligned to a reference before bpp-seqs can produce BPP output.");
}

/* ────────────────────────────────────────────────────────────────────────
 * BED inspector
 * ───────────────────────────────────────────────────────────────────── */

static int int64_cmp_asc(const void *a, const void *b)
{
    int64_t x = *(const int64_t *)a, y = *(const int64_t *)b;
    return (x > y) - (x < y);
}

static int strcmp_ptr(const void *a, const void *b)
{
    return strcmp(*(const char *const *)a, *(const char *const *)b);
}

static int find_or_add(char ***arr, int *n, const char *s)
{
    for (int i = 0; i < *n; i++) {
        if (strcmp((*arr)[i], s) == 0) return i;
    }
    *arr = (char **)realloc(*arr, (size_t)(*n + 1) * sizeof(char *));
    (*arr)[*n] = xstrdup(s);
    (*n)++;
    return *n - 1;
}

static void inspect_bed(const char *path, FileInfo *fi)
{
    gzFile gz = gzopen(path, "rb");
    if (!gz) {
        file_info_add_warning(fi, "OPEN_FAILED", "error", "Could not open BED file.");
        return;
    }
    char line[8192];
    int n_loci = 0, has_names = 0;
    int64_t *lens = NULL;
    int lens_cap = 0;
    char **chroms = NULL;
    int n_chroms = 0;

    while (gzgets(gz, line, sizeof(line)) != NULL) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        rtrim(line);
        if (!*line) continue;
        char *fields[8] = {0};
        int nf = 0;
        char *p = line;
        fields[nf++] = p;
        while (*p && nf < 8) {
            if (*p == '\t') { *p = '\0'; fields[nf++] = p + 1; }
            p++;
        }
        if (nf < 3) continue;
        int64_t start = strtoll(fields[1], NULL, 10);
        int64_t end   = strtoll(fields[2], NULL, 10);
        if (end <= start) continue;
        if (nf >= 4 && fields[3] && *fields[3]) has_names = 1;
        if (n_loci >= lens_cap) {
            lens_cap = lens_cap ? lens_cap * 2 : 64;
            lens = (int64_t *)realloc(lens, sizeof(int64_t) * (size_t)lens_cap);
        }
        lens[n_loci++] = end - start;
        find_or_add(&chroms, &n_chroms, fields[0]);
    }
    gzclose(gz);

    fi->n_loci = n_loci;
    fi->bed_has_names = has_names;
    fi->chromosomes = chroms;
    fi->n_chromosomes = n_chroms;
    if (n_loci > 0) {
        int64_t mn = lens[0], mx = lens[0], sum = 0;
        for (int i = 0; i < n_loci; i++) {
            if (lens[i] < mn) mn = lens[i];
            if (lens[i] > mx) mx = lens[i];
            sum += lens[i];
        }
        fi->bed_length_min  = mn;
        fi->bed_length_max  = mx;
        fi->bed_length_mean = sum / n_loci;
        qsort(lens, (size_t)n_loci, sizeof(int64_t), int64_cmp_asc);
        fi->bed_length_median = lens[n_loci / 2];
    }
    free(lens);
}

/* ────────────────────────────────────────────────────────────────────────
 * IMAP inspector
 * ───────────────────────────────────────────────────────────────────── */

static void inspect_imap(const char *path, FileInfo *fi)
{
    gzFile gz = gzopen(path, "rb");
    if (!gz) {
        file_info_add_warning(fi, "OPEN_FAILED", "error", "Could not open Imap file.");
        return;
    }
    char line[4096];
    int n = 0;
    char **samps = NULL;
    char **pops  = NULL;
    char **uniq_pops    = NULL;
    int    n_uniq_pops  = 0;
    while (gzgets(gz, line, sizeof(line)) != NULL) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;
        rtrim(line);
        if (!*line) continue;
        char *p = line;
        char *sample = p;
        while (*p && *p != ' ' && *p != '\t') p++;
        if (!*p) continue;
        *p = '\0'; p++;
        while (*p == ' ' || *p == '\t') p++;
        if (!*p) continue;
        char *pop = p;
        while (*p && *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') p++;
        *p = '\0';
        samps = (char **)realloc(samps, sizeof(char *) * (size_t)(n + 1));
        pops  = (char **)realloc(pops,  sizeof(char *) * (size_t)(n + 1));
        samps[n] = xstrdup(sample);
        pops[n]  = xstrdup(pop);
        n++;
        find_or_add(&uniq_pops, &n_uniq_pops, pop);
    }
    gzclose(gz);

    fi->imap_n_entries        = n;
    fi->imap_samples          = samps;
    fi->imap_pops             = pops;
    fi->imap_n_samples        = n;
    fi->imap_population_names = uniq_pops;
    fi->imap_n_population_names = n_uniq_pops;

    /* duplicate sample names into imap_sample_names for the JSON list */
    if (n > 0) {
        fi->imap_sample_names = (char **)malloc(sizeof(char *) * (size_t)n);
        for (int i = 0; i < n; i++) fi->imap_sample_names[i] = xstrdup(samps[i]);
        fi->imap_n_sample_names = n;
    }
}

/* ────────────────────────────────────────────────────────────────────────
 * PHYLIP inspector
 * ───────────────────────────────────────────────────────────────────── */

static void inspect_phylip(const char *path, FileInfo *fi)
{
    gzFile gz = gzopen(path, "rb");
    if (!gz) {
        file_info_add_warning(fi, "OPEN_FAILED", "error", "Could not open PHYLIP file.");
        return;
    }
    char line[1 << 16];
    if (gzgets(gz, line, sizeof(line)) == NULL) {
        gzclose(gz);
        return;
    }
    int nseq = 0, nsites = 0;
    if (sscanf(line, "%d %d", &nseq, &nsites) != 2 || nseq <= 0 || nsites <= 0) {
        gzclose(gz);
        return;
    }
    fi->phylip_n_sequences = nseq;
    fi->phylip_n_sites     = nsites;

    /* Read the next non-blank line to peek at format + capture first name */
    long mark_after_header = gztell(gz);
    if (gzgets(gz, line, sizeof(line)) == NULL) {
        gzclose(gz);
        fi->phylip_format = xstrdup("sequential");
        return;
    }
    rtrim(line);
    /* count sequence-like characters after first whitespace-delimited name */
    char *p = line;
    while (*p && !isspace((unsigned char)*p)) p++;
    /* count nucleotide-like chars in remainder of this line + continuation */
    int seq_chars = 0;
    int gap_chars = 0, n_chars = 0, total = 0;
    for (char *q = p; *q; q++) {
        if (isspace((unsigned char)*q)) continue;
        seq_chars++;
        total++;
        if (*q == '-' || *q == '?') gap_chars++;
        if (*q == 'N' || *q == 'n') n_chars++;
    }
    fi->phylip_format = xstrdup(seq_chars >= nsites ? "sequential" : "interleaved");

    /* approximate missing fraction over a small scan of the file */
    int64_t miss = 0, tot = 0;
    miss += gap_chars + n_chars;
    tot  += total;
    int lines_scanned = 1;
    while (lines_scanned < 200 && gzgets(gz, line, sizeof(line)) != NULL) {
        rtrim(line);
        if (!*line) continue;
        p = line;
        while (*p && !isspace((unsigned char)*p)) p++;
        for (char *q = p; *q; q++) {
            if (isspace((unsigned char)*q)) continue;
            tot++;
            if (*q == '-' || *q == '?' || *q == 'N' || *q == 'n') miss++;
        }
        lines_scanned++;
    }
    if (tot > 0) fi->phylip_missing_fraction = (double)miss / (double)tot;
    gzclose(gz);

    /* Second pass: capture the first nseq sample names so cross-validate
     * can match against the Imap. */
    gz = gzopen(path, "rb");
    if (gz) {
        gzgets(gz, line, sizeof(line));  /* header */
        char **names = (char **)calloc((size_t)nseq, sizeof(char *));
        int got = 0;
        while (got < nseq && gzgets(gz, line, sizeof(line)) != NULL) {
            int empty = 1;
            for (char *c = line; *c; c++) if (!isspace((unsigned char)*c)) { empty = 0; break; }
            if (empty) continue;
            char *q = line;
            while (*q && !isspace((unsigned char)*q)) q++;
            size_t nlen = (size_t)(q - line);
            if (nlen == 0) continue;
            names[got] = (char *)malloc(nlen + 1);
            memcpy(names[got], line, nlen);
            names[got][nlen] = '\0';
            got++;
        }
        fi->sequence_names = names;
        fi->n_sequence_names = got;
        gzclose(gz);
    }
    (void)mark_after_header;
}

/* ────────────────────────────────────────────────────────────────────────
 * NEXUS inspector
 * ───────────────────────────────────────────────────────────────────── */

static void parse_charset_range(const char *expr, int *start, int *end, int *stride)
{
    *start = -1; *end = -1; *stride = 1;
    /* Skip leading whitespace */
    while (*expr == ' ' || *expr == '\t') expr++;
    int a = 0, b = 0, s = 1;
    if (sscanf(expr, "%d-%d\\%d", &a, &b, &s) == 3) {
        *start = a; *end = b; *stride = s; return;
    }
    if (sscanf(expr, "%d-%d", &a, &b) == 2) {
        *start = a; *end = b; *stride = 1; return;
    }
    if (sscanf(expr, "%d", &a) == 1) {
        *start = a; *end = a; *stride = 1; return;
    }
}

static void inspect_nexus(const char *path, FileInfo *fi)
{
    gzFile gz = gzopen(path, "rb");
    if (!gz) {
        file_info_add_warning(fi, "OPEN_FAILED", "error", "Could not open NEXUS file.");
        return;
    }
    char line[1 << 16];
    int ntax = 0, nchar = 0;
    int n_charsets = 0;
    int charset_cap = 0;
    char **names = NULL;
    int *starts = NULL, *ends = NULL, *strides = NULL;
    int64_t miss = 0, tot = 0;
    int in_matrix = 0;

    while (gzgets(gz, line, sizeof(line)) != NULL) {
        char low[1024];
        int i; for (i = 0; i < (int)sizeof(low) - 1 && line[i]; i++) low[i] = (char)tolower((unsigned char)line[i]);
        low[i] = '\0';

        /* dimensions */
        char *p;
        if ((p = strstr(low, "ntax="))) sscanf(p + 5, "%d", &ntax);
        if ((p = strstr(low, "nchar="))) sscanf(p + 6, "%d", &nchar);

        /* matrix data scan for missing fraction */
        if (strstr(low, "matrix") && !in_matrix) {
            in_matrix = 1;
            continue;
        }
        if (in_matrix) {
            if (strchr(line, ';')) in_matrix = 0;
            /* capture taxon name (first non-blank token) */
            char *t = line;
            while (*t && isspace((unsigned char)*t)) t++;
            if (*t && *t != ';' && fi->n_sequence_names < (ntax ? ntax : 1024)) {
                char *te = t;
                while (*te && !isspace((unsigned char)*te) && *te != ';') te++;
                size_t nlen = (size_t)(te - t);
                if (nlen > 0) {
                    /* dedupe: append only if not already present */
                    int dup = 0;
                    for (int j = 0; j < fi->n_sequence_names; j++) {
                        if (strlen(fi->sequence_names[j]) == nlen &&
                            strncmp(fi->sequence_names[j], t, nlen) == 0) { dup = 1; break; }
                    }
                    if (!dup) {
                        fi->sequence_names = (char **)realloc(
                            fi->sequence_names,
                            sizeof(char *) * (size_t)(fi->n_sequence_names + 1));
                        char *cp = (char *)malloc(nlen + 1);
                        memcpy(cp, t, nlen); cp[nlen] = '\0';
                        fi->sequence_names[fi->n_sequence_names++] = cp;
                    }
                }
            }
            /* scan sequence-like chars */
            char *q = strchr(line, '\t');
            if (!q) q = strchr(line, ' ');
            if (q) {
                for (; *q; q++) {
                    if (isspace((unsigned char)*q) || *q == ';') continue;
                    tot++;
                    if (*q == '-' || *q == '?' || *q == 'N' || *q == 'n') miss++;
                }
            }
        }

        /* charset NAME = expr ; */
        if ((p = strstr(low, "charset"))) {
            const char *src = line + (p - low);
            src += 7;
            while (*src == ' ' || *src == '\t') src++;
            /* capture name up to '=' or whitespace */
            char name[128]; int ni = 0;
            while (*src && *src != '=' && *src != ' ' && *src != '\t' && ni < (int)sizeof(name) - 1) {
                name[ni++] = *src++;
            }
            name[ni] = '\0';
            const char *eq = strchr(src, '=');
            if (!eq) continue;
            eq++;
            /* range may span to ';' */
            char expr[256]; int ei = 0;
            while (*eq && *eq != ';' && ei < (int)sizeof(expr) - 1) expr[ei++] = *eq++;
            expr[ei] = '\0';

            int a, b, s; parse_charset_range(expr, &a, &b, &s);
            if (n_charsets >= charset_cap) {
                charset_cap = charset_cap ? charset_cap * 2 : 8;
                names   = (char **)realloc(names,   sizeof(char *) * (size_t)charset_cap);
                starts  = (int *)  realloc(starts,  sizeof(int)    * (size_t)charset_cap);
                ends    = (int *)  realloc(ends,    sizeof(int)    * (size_t)charset_cap);
                strides = (int *)  realloc(strides, sizeof(int)    * (size_t)charset_cap);
            }
            names[n_charsets]   = xstrdup(name);
            starts[n_charsets]  = a;
            ends[n_charsets]    = b;
            strides[n_charsets] = s;
            n_charsets++;
        }
    }
    gzclose(gz);

    fi->nexus_n_sequences      = ntax;
    fi->nexus_n_sites          = nchar;
    fi->nexus_n_loci           = n_charsets ? n_charsets : 1;
    fi->nexus_has_charsets     = (n_charsets > 0);
    fi->nexus_charset_names    = names;
    fi->nexus_n_charset_names  = n_charsets;
    fi->nexus_charset_starts   = starts;
    fi->nexus_charset_ends     = ends;
    fi->nexus_charset_strides  = strides;
    if (tot > 0) fi->nexus_missing_fraction = (double)miss / (double)tot;
}

/* ────────────────────────────────────────────────────────────────────────
 * Dispatcher
 * ───────────────────────────────────────────────────────────────────── */

FileInfo *inspect_file(const char *path)
{
    if (!path) return NULL;
    FileInfo *fi = file_info_new(path);
    if (!fi) return NULL;

    FileType t = detect_type(path);
    fi->ft = t;

    switch (t) {
        case BS_BAM:   inspect_bam(path, fi, 0); break;
        case BS_CRAM:  inspect_bam(path, fi, 1); break;
        case BS_VCF:   inspect_vcf(path, fi); break;
        /* BS_GVCF is never set by detector; inspect_vcf upgrades */
        case BS_FASTA_MSA: /* refined by inspect_fasta */
            inspect_fasta(path, fi);
            break;
        case BS_FASTQ:   inspect_fastq(path, fi); break;
        case BS_BED:     inspect_bed(path, fi); break;
        case BS_IMAP:    inspect_imap(path, fi); break;
        case BS_PHYLIP:  inspect_phylip(path, fi); break;
        case BS_NEXUS:   inspect_nexus(path, fi); break;
        case BS_GVCF:
        case BS_FASTA_REFERENCE:
        case BS_FASTA_CONTIGS:
        case BS_UNKNOWN:
        default:
            file_info_add_warning(fi, "UNKNOWN_TYPE", "warning",
                "Could not determine file type from content.");
            break;
    }
    return fi;
}

/* Unused helpers kept for future use */
static __attribute__((unused)) int _unused_int64_cmp(const void *a, const void *b) {
    return int64_cmp(a, b);
}
static __attribute__((unused)) int _unused_strcmp_ptr(const void *a, const void *b) {
    return strcmp_ptr(a, b);
}
