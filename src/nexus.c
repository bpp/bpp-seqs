#define _POSIX_C_SOURCE 200809L

#include "nexus.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>   /* strncasecmp */
#include <zlib.h>

static int is_word(int c) { return isalnum(c) || c == '_'; }

/* Read a whole (possibly gzipped) file into a NUL-terminated buffer. */
static char *slurp(const char *path)
{
    gzFile gz = gzopen(path, "rb");
    if (!gz) return NULL;
    size_t cap = 65536, len = 0;
    char *buf = (char *)malloc(cap);
    if (!buf) { gzclose(gz); return NULL; }
    int got;
    while ((got = gzread(gz, buf + len, (unsigned)(cap - len - 1))) > 0) {
        len += (size_t)got;
        if (len + 1 >= cap) {
            cap *= 2;
            char *nb = (char *)realloc(buf, cap);
            if (!nb) { free(buf); gzclose(gz); return NULL; }
            buf = nb;
        }
    }
    buf[len] = '\0';
    gzclose(gz);
    return buf;
}

/* Remove NEXUS `[ ... ]` comments (nestable), preserving single-quoted tokens
 * (inside which '[' and ']' are literal). Newlines outside comments survive,
 * so INTERLEAVE row structure is retained. Returns a new buffer; caller frees. */
static char *strip_comments(const char *s)
{
    size_t n = strlen(s);
    char *out = (char *)malloc(n + 1);
    if (!out) return NULL;
    size_t o = 0;
    int in_quote = 0;
    for (size_t i = 0; i < n; ) {
        char c = s[i];
        if (in_quote) {
            if (c == '\'') {
                if (s[i + 1] == '\'') { out[o++] = '\''; out[o++] = '\''; i += 2; continue; }
                out[o++] = c; in_quote = 0; i++; continue;
            }
            out[o++] = c; i++; continue;
        }
        if (c == '\'') { out[o++] = c; in_quote = 1; i++; continue; }
        if (c == '[') {
            int depth = 1; i++;
            while (i < n && depth > 0) { if (s[i] == '[') depth++; else if (s[i] == ']') depth--; i++; }
            continue;
        }
        out[o++] = c; i++;
    }
    out[o] = '\0';
    return out;
}

/* Find keyword `kw` (case-insensitive) as a whole word within [from, to). */
static const char *find_kw(const char *from, const char *to, const char *kw)
{
    size_t kl = strlen(kw);
    for (const char *p = from; p + kl <= to; p++) {
        if (strncasecmp(p, kw, kl) != 0) continue;
        char before = (p == from) ? ' ' : p[-1];
        char after  = (p + kl < to) ? p[kl] : ' ';
        if (!is_word((unsigned char)before) && !is_word((unsigned char)after)) return p;
    }
    return NULL;
}

/* Within [start, end), find `key`, then (whitespace) '=' (whitespace), then a
 * value token. Copies the value to out; returns 1 if found. Tolerant of any
 * whitespace around '=' (so "NTAX=26", "NTAX = 26", "NTAX =26" all work). */
static int kv_value(const char *start, const char *end, const char *key,
                    char *out, size_t cap)
{
    const char *k = find_kw(start, end, key);
    if (!k) return 0;
    const char *p = k + strlen(key);
    while (p < end && isspace((unsigned char)*p)) p++;
    if (p >= end || *p != '=') return 0;
    p++;
    while (p < end && isspace((unsigned char)*p)) p++;
    size_t o = 0;
    if (p < end && *p == '"') {              /* quoted value, e.g. SYMBOLS="01" */
        p++;
        while (p < end && *p != '"' && o < cap - 1) out[o++] = *p++;
    } else {
        while (p < end && !isspace((unsigned char)*p) &&
               *p != ';' && *p != '=' && o < cap - 1) out[o++] = *p++;
    }
    out[o] = '\0';
    return o > 0;
}

static int kv_int(const char *s, const char *e, const char *key)
{
    char b[32];
    if (!kv_value(s, e, key, b, sizeof b)) return -1;
    if (!isdigit((unsigned char)b[0])) return -1;
    return atoi(b);
}

/* FORMAT ... INTERLEAVE / INTERLEAVE=YES|NO */
static int format_interleaved(const char *s, const char *e)
{
    const char *k = find_kw(s, e, "interleave");
    if (!k) return 0;
    const char *p = k + strlen("interleave");
    while (p < e && isspace((unsigned char)*p)) p++;
    if (p < e && *p == '=') {
        p++;
        while (p < e && isspace((unsigned char)*p)) p++;
        if (p + 1 < e && strncasecmp(p, "no", 2) == 0) return 0;
    }
    return 1;
}

/* Parse a single charset range expression: "a", "a-b", "a-b\s", "a-." */
static void parse_range(const char *expr, const char *end, int *a, int *b, int *s)
{
    *a = -1; *b = -1; *s = 1;
    char buf[64]; int bi;
    while (expr < end && isspace((unsigned char)*expr)) expr++;
    bi = 0;
    while (expr < end && isdigit((unsigned char)*expr) && bi < 63) buf[bi++] = *expr++;
    buf[bi] = '\0';
    if (bi == 0) return;
    *a = atoi(buf); *b = *a;
    while (expr < end && isspace((unsigned char)*expr)) expr++;
    if (expr < end && *expr == '-') {
        expr++;
        while (expr < end && isspace((unsigned char)*expr)) expr++;
        if (expr < end && *expr == '.') { *b = -1; expr++; }   /* '.' = last site */
        else {
            bi = 0;
            while (expr < end && isdigit((unsigned char)*expr) && bi < 63) buf[bi++] = *expr++;
            buf[bi] = '\0';
            if (bi) *b = atoi(buf);
        }
    }
    while (expr < end && isspace((unsigned char)*expr)) expr++;
    if (expr < end && *expr == '\\') {
        expr++;
        while (expr < end && isspace((unsigned char)*expr)) expr++;
        bi = 0;
        while (expr < end && isdigit((unsigned char)*expr) && bi < 63) buf[bi++] = *expr++;
        buf[bi] = '\0';
        if (bi) *s = atoi(buf);
    }
}

/* Read a taxon label at *pp (single-quote aware; interior space -> '_'). */
static void read_label(const char **pp, const char *end, char *out, size_t cap)
{
    const char *p = *pp;
    size_t o = 0;
    if (p < end && *p == '\'') {
        p++;
        while (p < end) {
            if (*p == '\'') {
                if (p + 1 < end && p[1] == '\'') { if (o < cap - 1) out[o++] = '\''; p += 2; continue; }
                p++; break;
            }
            char c = *p++;
            if (c == ' ' || c == '\t') c = '_';
            if (o < cap - 1) out[o++] = c;
        }
    } else {
        while (p < end && !isspace((unsigned char)*p) && o < cap - 1) out[o++] = *p++;
    }
    out[o] = '\0';
    *pp = p;
}

/* Read the MATRIX body [start, end) into doc->taxa / doc->seq.
 *
 * Line-based: each physical line is "label seq-fragment"; the fragment is keyed
 * to the label's taxon and appended. This handles BOTH the sequential layout
 * (each taxon once, on its own line) and the interleaved layout (labels repeat
 * across blocks), so the FORMAT INTERLEAVE flag isn't needed here. Characters
 * arriving after a taxon already has nchar sites are counted and reported (a
 * header/data mismatch), matching the declared-vs-actual contract. */
static int read_matrix(const char *start, const char *end, int ntax, int nchar,
                       int interleaved, char matchchar, NexusDoc *doc)
{
    (void)interleaved;
    doc->taxa = (char **)calloc((size_t)ntax, sizeof(char *));
    doc->seq  = (char **)calloc((size_t)ntax, sizeof(char *));
    int *lens = (int *)calloc((size_t)ntax, sizeof(int));
    if (!doc->taxa || !doc->seq || !lens) { free(lens); return -1; }
    for (int i = 0; i < ntax; i++) { doc->seq[i] = (char *)malloc((size_t)nchar + 1); doc->seq[i][0] = '\0'; }

    const char *p = start;
    int nseen = 0;
    long excess = 0;
    while (p < end) {
        while (p < end && isspace((unsigned char)*p)) p++;   /* to next row */
        if (p >= end || *p == ';') break;

        char label[256];
        read_label(&p, end, label, sizeof label);
        if (label[0] == '\0') break;

        int idx = -1;
        for (int i = 0; i < nseen; i++) if (strcmp(doc->taxa[i], label) == 0) { idx = i; break; }
        if (idx < 0) {
            if (nseen < ntax) { idx = nseen; doc->taxa[nseen++] = strdup(label); }
            else idx = -1;   /* more distinct labels than ntax -> ignore data */
        }

        /* sequence fragment: to end of this physical line */
        while (p < end && *p != '\n' && *p != ';') {
            char c = *p;
            if (isspace((unsigned char)c)) { p++; continue; }
            if (c == '{' || c == '(') {                 /* ambiguity group = one site */
                char close = (c == '{') ? '}' : ')';
                p++;
                while (p < end && *p != close) p++;
                if (p < end) p++;
                if (idx >= 0) { if (lens[idx] < nchar) doc->seq[idx][lens[idx]++] = 'N'; else excess++; }
                continue;
            }
            p++;
            if (idx >= 0) {
                if (lens[idx] < nchar) doc->seq[idx][lens[idx]++] = (char)toupper((unsigned char)c);
                else excess++;
            }
        }
    }
    if (excess > 0)
        fprintf(stderr, "WARNING: NEXUS MATRIX: %ld character%s beyond declared "
                "nchar=%d were dropped (header/data mismatch).\n",
                excess, excess == 1 ? "" : "s", nchar);

    for (int i = 0; i < ntax; i++) {
        int L = lens[i]; if (L > nchar) L = nchar;
        if (doc->seq[i]) doc->seq[i][L] = '\0';
    }
    if (matchchar) {
        char mc = (char)toupper((unsigned char)matchchar);
        for (int i = 1; i < ntax; i++)
            for (int j = 0; j < nchar && doc->seq[i][j]; j++)
                if (doc->seq[i][j] == mc && doc->seq[0][j]) doc->seq[i][j] = doc->seq[0][j];
    }

    int ok = (nseen == ntax);
    for (int i = 0; i < ntax && ok; i++)
        if (!doc->taxa[i] || (int)strlen(doc->seq[i]) != nchar) ok = 0;
    free(lens);
    return ok ? 0 : -1;
}

static void parse_charsets(const char *clean, const char *end,
                           const char *mstart, const char *mend, NexusDoc *doc)
{
    int cap = 0;
    const char *p = clean;
    for (;;) {
        const char *c = find_kw(p, end, "charset");
        if (!c) break;
        if (mstart && c >= mstart && c < mend) { p = c + 7; continue; }  /* inside MATRIX */
        const char *q = c + strlen("charset");
        while (q < end && isspace((unsigned char)*q)) q++;
        char name[128]; int ni = 0;
        while (q < end && *q != '=' && !isspace((unsigned char)*q) && ni < 127) name[ni++] = *q++;
        name[ni] = '\0';
        const char *eq = (const char *)memchr(q, '=', (size_t)(end - q));
        if (!eq) { p = c + 7; continue; }
        const char *semi = (const char *)memchr(eq, ';', (size_t)(end - eq));
        const char *rend = semi ? semi : end;
        int a, b, s;
        parse_range(eq + 1, rend, &a, &b, &s);
        if (a >= 1) {
            if (doc->n_charsets >= cap) {
                cap = cap ? cap * 2 : 8;
                doc->cs_name   = (char **)realloc(doc->cs_name,   sizeof(char *) * (size_t)cap);
                doc->cs_start  = (int *)  realloc(doc->cs_start,  sizeof(int)    * (size_t)cap);
                doc->cs_end    = (int *)  realloc(doc->cs_end,    sizeof(int)    * (size_t)cap);
                doc->cs_stride = (int *)  realloc(doc->cs_stride, sizeof(int)    * (size_t)cap);
            }
            doc->cs_name[doc->n_charsets]   = strdup(name);
            doc->cs_start[doc->n_charsets]  = a;
            doc->cs_end[doc->n_charsets]    = b;
            doc->cs_stride[doc->n_charsets] = s;
            doc->n_charsets++;
        }
        p = semi ? semi + 1 : end;
    }
}

int nexus_parse(const char *path, NexusDoc *doc)
{
    memset(doc, 0, sizeof *doc);

    char *raw = slurp(path);
    if (!raw) { snprintf(doc->err, sizeof doc->err, "cannot read '%s'", path); return -1; }
    char *clean = strip_comments(raw);
    free(raw);
    if (!clean) { snprintf(doc->err, sizeof doc->err, "out of memory"); return -1; }

    const char *end = clean + strlen(clean);
    const char *mstart = find_kw(clean, end, "matrix");
    const char *head_end = mstart ? mstart : end;   /* DIMENSIONS/FORMAT precede MATRIX */

    int ntax  = kv_int(clean, head_end, "ntax");
    int nchar = kv_int(clean, head_end, "nchar");
    if (ntax <= 0 || nchar <= 0) {
        snprintf(doc->err, sizeof doc->err,
                 "missing or invalid DIMENSIONS (ntax=%d nchar=%d)", ntax, nchar);
        free(clean); return -1;
    }
    doc->ntax = ntax; doc->nchar = nchar;

    char matchchar = 0;
    const char *fmt = find_kw(clean, head_end, "format");
    if (fmt) {
        const char *fe = (const char *)memchr(fmt, ';', (size_t)(end - fmt));
        if (!fe) fe = head_end;
        char v[16];
        if (kv_value(fmt, fe, "datatype", v, sizeof v))
            for (int i = 0; v[i] && i < 15; i++) doc->datatype[i] = (char)tolower((unsigned char)v[i]);
        doc->interleaved = format_interleaved(fmt, fe);
        char mc[8];
        if (kv_value(fmt, fe, "matchchar", mc, sizeof mc)) matchchar = mc[0];
    }

    if (!mstart) { snprintf(doc->err, sizeof doc->err, "no MATRIX command"); free(clean); return -1; }
    const char *mbody = mstart + strlen("matrix");
    const char *mend  = (const char *)memchr(mbody, ';', (size_t)(end - mbody));
    if (!mend) mend = end;

    if (read_matrix(mbody, mend, ntax, nchar, doc->interleaved, matchchar, doc) != 0) {
        snprintf(doc->err, sizeof doc->err,
                 "MATRIX does not match declared ntax=%d nchar=%d", ntax, nchar);
        free(clean); return -1;
    }

    long miss = 0, tot = 0;
    for (int i = 0; i < ntax; i++)
        for (int j = 0; j < nchar; j++) {
            char c = doc->seq[i][j];
            tot++;
            if (c == '-' || c == '?' || c == 'N') miss++;
        }
    doc->missing_fraction = tot ? (double)miss / (double)tot : 0.0;

    parse_charsets(clean, end, mstart, mend, doc);

    free(clean);
    return 0;
}

void nexus_free(NexusDoc *doc)
{
    if (!doc) return;
    if (doc->taxa) { for (int i = 0; i < doc->ntax; i++) free(doc->taxa[i]); free(doc->taxa); }
    if (doc->seq)  { for (int i = 0; i < doc->ntax; i++) free(doc->seq[i]);  free(doc->seq);  }
    if (doc->cs_name) { for (int i = 0; i < doc->n_charsets; i++) free(doc->cs_name[i]); free(doc->cs_name); }
    free(doc->cs_start); free(doc->cs_end); free(doc->cs_stride);
    memset(doc, 0, sizeof *doc);
}
