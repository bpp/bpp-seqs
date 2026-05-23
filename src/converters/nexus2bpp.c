#define _POSIX_C_SOURCE 200809L

#include "aln_writer.h"
#include "converters.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>    /* strncasecmp (POSIX, declared here not in <string.h>) */
#include <zlib.h>

/* Read whole file (possibly gzipped) into memory.  Returns NUL-terminated
 * string; caller frees. */
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

/* Parse the NEXUS DATA matrix from `text`. Fills names and seqs (each of
 * length nchar). Returns 0 on success. */
static int parse_nexus_matrix(const char *text, int ntax, int nchar,
                              char ***names_out, char ***seqs_out)
{
    /* Find "matrix" (case-insensitive) after "data" / "characters" block. */
    const char *p = text;
    const char *matrix_kw = NULL;
    while (*p) {
        if ((p[0] == 'm' || p[0] == 'M') &&
            strncasecmp(p, "matrix", 6) == 0) { matrix_kw = p + 6; break; }
        p++;
    }
    if (!matrix_kw) return -1;

    char **names = (char **)calloc((size_t)ntax, sizeof(char *));
    char **seqs  = (char **)calloc((size_t)ntax, sizeof(char *));
    int   *lens  = (int *)  calloc((size_t)ntax, sizeof(int));
    for (int i = 0; i < ntax; i++) {
        seqs[i] = (char *)malloc((size_t)nchar + 1);
        seqs[i][0] = '\0';
    }

    /* Iterate "name <seq fragment>" tokens until ';'. NEXUS allows
     * interleaving and blank lines. */
    int row = 0;
    int done_first_block = 0;
    int excess_chars = 0;
    int excess_rows  = 0;
    const char *q = matrix_kw;
    while (*q && *q != ';') {
        /* skip whitespace */
        while (*q && (isspace((unsigned char)*q)) && *q != '\n') q++;
        if (*q == '\n') {
            /* end of row → next row */
            q++;
            continue;
        }
        if (!*q || *q == ';') break;

        /* Capture the row's "name + seq fragment" up to newline or ';' */
        const char *line_end = q;
        while (*line_end && *line_end != '\n' && *line_end != ';') line_end++;
        /* Parse name */
        const char *t = q;
        while (t < line_end && !isspace((unsigned char)*t)) t++;
        size_t nlen = (size_t)(t - q);

        int target;
        if (!done_first_block && row < ntax && nlen > 0) {
            /* Capture name on first occurrence */
            if (!names[row]) {
                names[row] = (char *)malloc(nlen + 1);
                memcpy(names[row], q, nlen);
                names[row][nlen] = '\0';
            }
            target = row;
        } else {
            target = row;
        }

        /* Detect a row arriving after all samples are full and the locus
         * has been completely populated — i.e., an extra taxon beyond ntax. */
        int all_full = 1;
        for (int i = 0; i < ntax; i++) if (lens[i] < nchar) { all_full = 0; break; }
        if (all_full && nlen > 0) excess_rows++;

        /* Append sequence characters from t..line_end */
        for (const char *c = t; c < line_end; c++) {
            if (isspace((unsigned char)*c)) continue;
            if (target < ntax && lens[target] < nchar) {
                seqs[target][lens[target]++] = (char)toupper((unsigned char)*c);
            } else if (target < ntax) {
                excess_chars++;
            }
        }
        row++;
        if (row >= ntax) { row = 0; done_first_block = 1; }
        q = (*line_end == '\n') ? line_end + 1 : line_end;
    }

    if (excess_chars > 0) {
        fprintf(stderr,
            "WARNING: NEXUS MATRIX: %d character%s beyond declared "
            "nchar=%d were dropped (header/data mismatch).\n",
            excess_chars, excess_chars == 1 ? "" : "s", nchar);
    }
    if (excess_rows > 0) {
        fprintf(stderr,
            "WARNING: NEXUS MATRIX: %d unexpected data row%s after "
            "ntax=%d samples were filled (extra taxa or stray content).\n",
            excess_rows, excess_rows == 1 ? "" : "s", ntax);
    }

    int ok = 1;
    for (int i = 0; i < ntax; i++) {
        seqs[i][lens[i]] = '\0';
        if (lens[i] != nchar || !names[i]) ok = 0;
    }
    free(lens);
    if (!ok) {
        for (int i = 0; i < ntax; i++) { free(names[i]); free(seqs[i]); }
        free(names); free(seqs);
        return -2;
    }
    *names_out = names;
    *seqs_out  = seqs;
    return 0;
}

int convert_nexus(FileInfo **files, int n_files,
                  FileInfo *imap_fi,
                  const char *out_prefix,
                  const ConvertOpts *opts,
                  ConversionResult *cr)
{
    /* Single NEXUS file expected (the spec doesn't require multi-file). */
    FileInfo *nx = NULL;
    for (int i = 0; i < n_files; i++) if (files[i]->ft == BS_NEXUS) { nx = files[i]; break; }
    if (!nx || !imap_fi) {
        fprintf(stderr, "Error: nexus2bpp requires a NEXUS file and an Imap.\n");
        return 1;
    }

    char *text = slurp(nx->path);
    if (!text) { fprintf(stderr, "Error: cannot read '%s'\n", nx->path); return 1; }

    int ntax = nx->nexus_n_sequences, nchar = nx->nexus_n_sites;
    if (ntax <= 0 || nchar <= 0) {
        fprintf(stderr, "Error: NEXUS dimensions missing (ntax=%d nchar=%d).\n", ntax, nchar);
        free(text); return 1;
    }

    char **names = NULL, **seqs = NULL;
    if (parse_nexus_matrix(text, ntax, nchar, &names, &seqs) != 0) {
        fprintf(stderr, "Error: failed to parse NEXUS data matrix in '%s'\n", nx->path);
        free(text); return 1;
    }
    free(text);

    /* Decide loci: one per charset, or one over all sites if none. */
    int n_loci_out = nx->nexus_has_charsets ? nx->nexus_n_charset_names : 1;
    LocusAln *loci = (LocusAln *)calloc((size_t)n_loci_out, sizeof(LocusAln));

    if (!nx->nexus_has_charsets) {
        loci[0].name = strdup("locus1");
        loci[0].length = nchar;
        loci[0].n_seqs = ntax;
        loci[0].seq_names = names;
        loci[0].seqs      = seqs;
        loci[0].source_kind   = strdup("CHARSET");
        loci[0].source_file   = strdup(nx->path);
        loci[0].source_chrom  = NULL;
        loci[0].source_start  = 1;
        loci[0].source_end    = nchar;
        loci[0].source_stride = 1;
    } else {
        for (int c = 0; c < nx->nexus_n_charset_names; c++) {
            int s = nx->nexus_charset_starts[c];   /* 1-based inclusive */
            int e = nx->nexus_charset_ends[c];
            int st = nx->nexus_charset_strides[c] ? nx->nexus_charset_strides[c] : 1;
            if (s < 1) s = 1;
            if (e < 1 || e > nchar) e = nchar;
            int len = 0;
            for (int x = s; x <= e; x += st) len++;
            loci[c].name = strdup(nx->nexus_charset_names[c]);
            loci[c].length = len;
            loci[c].n_seqs = ntax;
            loci[c].seq_names = (char **)calloc((size_t)ntax, sizeof(char *));
            loci[c].seqs      = (char **)calloc((size_t)ntax, sizeof(char *));
            loci[c].source_kind   = strdup("CHARSET");
            loci[c].source_file   = strdup(nx->path);
            loci[c].source_chrom  = NULL;
            loci[c].source_start  = s;
            loci[c].source_end    = e;
            loci[c].source_stride = st;
            for (int t = 0; t < ntax; t++) {
                loci[c].seq_names[t] = strdup(names[t]);
                loci[c].seqs[t]      = (char *)malloc((size_t)len + 1);
                int k = 0;
                for (int x = s; x <= e; x += st) {
                    loci[c].seqs[t][k++] = seqs[t][x - 1];
                }
                loci[c].seqs[t][len] = '\0';
            }
        }
        /* Free the master matrix once we've sliced it. */
        for (int i = 0; i < ntax; i++) { free(names[i]); free(seqs[i]); }
        free(names); free(seqs);
    }

    return write_alignment_outputs(out_prefix, loci, n_loci_out, imap_fi,
                                   opts->min_length, opts->max_missing,
                                   opts->min_snps, opts->keep_invariant, cr);
}
