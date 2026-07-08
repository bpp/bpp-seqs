#define _POSIX_C_SOURCE 200809L

#include "bpp_parser.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *xdup(const char *s)
{
    if (!s) return NULL;
    size_t n = strlen(s);
    char *r = (char *)malloc(n + 1);
    if (!r) return NULL;
    memcpy(r, s, n + 1);
    return r;
}

/* True if the line begins with two whitespace-separated positive ints. */
static int looks_like_header(const char *line, int *ns, int *nc)
{
    int a = 0, b = 0;
    char trailer = 0;
    int got = sscanf(line, " %d %d %c", &a, &b, &trailer);
    if (got >= 2 && a > 0 && b > 0) {
        if (got == 2 || trailer == '\0' || trailer == '\n' ||
            trailer == '\r' || isspace((unsigned char)trailer)) {
            if (ns) *ns = a;
            if (nc) *nc = b;
            return 1;
        }
    }
    return 0;
}

static int is_blank(const char *line)
{
    for (const char *p = line; *p; p++) if (!isspace((unsigned char)*p)) return 0;
    return 1;
}

/* Extract the post-caret id (or the whole name if no caret). */
static char *post_caret_id(const char *name)
{
    const char *c = strrchr(name, '^');
    return xdup(c ? c + 1 : name);
}

/* Read one locus starting from `header_line` (already buffered) using
 * `fp` for continuation lines. Returns 0 on success, -1 on EOF, -2 on
 * malformed.  On success, populates *out and consumes lines up to and
 * including the locus's data (but does not consume the next locus's
 * header — that is left in *carry_line if present). */
static int read_one_locus(FILE *fp, char *line, size_t cap,
                          BppLocus *out)
{
    int nseq = 0, nsites = 0;
    if (!looks_like_header(line, &nseq, &nsites)) return -2;

    out->n_seqs    = nseq;
    out->n_sites   = nsites;
    out->seq_names = (char **)calloc((size_t)nseq, sizeof(char *));
    out->ids       = (char **)calloc((size_t)nseq, sizeof(char *));
    out->seqs      = (char **)calloc((size_t)nseq, sizeof(char *));
    int *lens      = (int *)  calloc((size_t)nseq, sizeof(int));
    for (int i = 0; i < nseq; i++) {
        out->seqs[i] = (char *)malloc((size_t)nsites + 1);
        if (out->seqs[i]) out->seqs[i][0] = '\0';
    }

    int row = 0;
    int finished_first_block = 0;

    while (fgets(line, (int)cap, fp) != NULL) {
        if (is_blank(line)) {
            if (row > 0) finished_first_block = 1;
            row = 0;
            continue;
        }

        /* If we're between rows and the line looks like a new header AND
         * every sample is already full, stash it and stop. */
        if (row == 0 && finished_first_block) {
            int peek_a, peek_b;
            if (looks_like_header(line, &peek_a, &peek_b)) {
                int all_full = 1;
                for (int i = 0; i < nseq; i++) if (lens[i] < nsites) { all_full = 0; break; }
                if (all_full) {
                    /* Rewind so caller's next fgets() sees this header. */
                    fseek(fp, -(long)strlen(line), SEEK_CUR);
                    break;
                }
            }
        }

        char *p = line;
        const char *seq_start = p;
        int target = row;
        int got_name = 0;
        if (!finished_first_block) {
            char *name_end = p;
            while (*name_end && !isspace((unsigned char)*name_end)) name_end++;
            size_t nlen = (size_t)(name_end - p);
            if (nlen > 0 && row < nseq) {
                out->seq_names[row] = (char *)malloc(nlen + 1);
                memcpy(out->seq_names[row], p, nlen);
                out->seq_names[row][nlen] = '\0';
                out->ids[row] = post_caret_id(out->seq_names[row]);
                seq_start = name_end;
                target = row;
                got_name = 1;
            } else if (row > 0) {
                /* No name token in the first block → wrap continuation of
                 * the previously named sample (per-sample multi-line layout,
                 * e.g. the mammoth_nuclear.txt fixture BPP itself accepts).
                 * The row cursor stays put until the next named line. */
                seq_start = p;
                target = (row > nseq ? nseq : row) - 1;
            }
            /* else: row==0 with no name in the first block — malformed
             *       locus body; fall through with seq_names[0]==NULL so
             *       the validation at the end of read_one_locus returns -2. */
        }

        for (const char *q = seq_start; *q && *q != '\n' && *q != '\r'; q++) {
            if (isspace((unsigned char)*q)) continue;
            if (target < nseq && lens[target] < nsites) {
                out->seqs[target][lens[target]++] = *q;
            }
        }
        if (got_name || finished_first_block) {
            row++;
            if (row >= nseq && !finished_first_block) {
                /* End of the first batch of named lines. Disambiguate two
                 * layouts by checking whether sample 0 has been wrapped to
                 * completion:
                 *   sample 0 full → per-sample wrap; keep row at nseq so
                 *                   further nameless lines continue the
                 *                   most-recently-named sample (nseq-1).
                 *   sample 0 not full → PHYLIP-interleaved; flip into
                 *                       subsequent-block mode so the next
                 *                       nameless line targets sample 0.
                 * A blank-line separator (the canonical interleaved form)
                 * independently triggers the flip via the is_blank handler
                 * above. */
                if (lens[0] < nsites) {
                    row = 0;
                    finished_first_block = 1;
                }
            }
        }

        int all_full = 1;
        for (int i = 0; i < nseq; i++) if (lens[i] < nsites) { all_full = 0; break; }
        if (all_full) break;
    }

    int ok = 1;
    for (int i = 0; i < nseq; i++) {
        if (out->seqs[i]) out->seqs[i][lens[i]] = '\0';
        if (lens[i] != nsites) ok = 0;
        if (!out->seq_names[i]) ok = 0;
    }
    free(lens);
    return ok ? 0 : -2;
}

BppLocus *bpp_parse_file(const char *path, int *n_out)
{
    FILE *fp = fopen(path, "r");
    if (!fp) return NULL;

    BppLocus *arr = NULL;
    int n = 0, cap = 0;
    char line[1 << 16];

    /* Find the first header line */
    while (fgets(line, sizeof(line), fp) != NULL) {
        if (is_blank(line)) continue;
        int a, b;
        if (looks_like_header(line, &a, &b)) {
            for (;;) {
                if (n >= cap) {
                    cap = cap ? cap * 2 : 16;
                    arr = (BppLocus *)realloc(arr, sizeof(BppLocus) * (size_t)cap);
                }
                BppLocus *cur = &arr[n];
                memset(cur, 0, sizeof(*cur));
                cur->source_start = cur->source_end = -1;
                cur->source_stride = 1;

                int rc = read_one_locus(fp, line, sizeof(line), cur);
                if (rc != 0) {
                    /* malformed locus or EOF before finishing — best-effort drop */
                    if (cur->seqs) {
                        for (int i = 0; i < cur->n_seqs; i++) {
                            free(cur->seq_names ? cur->seq_names[i] : NULL);
                            free(cur->ids       ? cur->ids[i]       : NULL);
                            free(cur->seqs      ? cur->seqs[i]      : NULL);
                        }
                        free(cur->seq_names); free(cur->ids); free(cur->seqs);
                    }
                    break;
                }

                /* Default name if no metadata attached */
                char nm[32]; snprintf(nm, sizeof(nm), "locus%d", n + 1);
                cur->name = xdup(nm);

                n++;

                /* Look for the next header */
                int found = 0;
                while (fgets(line, sizeof(line), fp) != NULL) {
                    if (is_blank(line)) continue;
                    if (looks_like_header(line, &a, &b)) { found = 1; break; }
                }
                if (!found) break;
            }
            break;
        }
    }
    fclose(fp);

    *n_out = n;
    if (n == 0) {
        free(arr);
        return NULL;
    }
    return arr;
}

/* Parse a single TSV line into ≤ max_fields fields. Mutates the buffer. */
static int split_tsv(char *line, char **fields, int max_fields)
{
    int n = 0;
    char *p = line;
    if (max_fields <= 0) return 0;
    fields[n++] = p;
    while (*p && n < max_fields) {
        if (*p == '\t') {
            *p = '\0';
            p++;
            fields[n++] = p;
            continue;
        }
        if (*p == '\n' || *p == '\r') { *p = '\0'; break; }
        p++;
    }
    /* Trim trailing \r/\n on the last field */
    if (n > 0) {
        char *t = fields[n - 1];
        while (*t) t++;
        while (t > fields[n - 1] && (*(t - 1) == '\n' || *(t - 1) == '\r')) *--t = '\0';
    }
    return n;
}

int bpp_attach_loci_tsv(BppLocus *loci, int n_loci, const char *tsv_path)
{
    FILE *fp = fopen(tsv_path, "r");
    if (!fp) return 0;
    char line[1 << 14];
    int header_seen = 0;
    int attached = 0;
    int row_index = 0;

    while (fgets(line, sizeof(line), fp) != NULL) {
        if (is_blank(line)) continue;
        if (!header_seen) {
            header_seen = 1;
            /* Expect the standard header.  Don't enforce — just skip. */
            if (strncmp(line, "locus_name", 10) == 0) continue;
            /* If first row isn't a header, fall through to treat it as data. */
        }

        char *fields[16] = {0};
        int nf = split_tsv(line, fields, 16);
        if (nf < 1) continue;

        /* Match by name first; fall back to positional. */
        const char *name = fields[0];
        BppLocus *target = NULL;
        for (int i = 0; i < n_loci; i++) {
            if (loci[i].name && strcmp(loci[i].name, name) == 0) { target = &loci[i]; break; }
        }
        if (!target && row_index < n_loci) target = &loci[row_index];
        row_index++;
        if (!target) continue;

        /* Override the synthesized name with the metadata one. */
        if (name && *name && strcmp(name, ".") != 0) {
            free(target->name);
            target->name = xdup(name);
        }
        if (nf > 1 && strcmp(fields[1], ".") != 0) {
            free(target->source_kind);   target->source_kind   = xdup(fields[1]);
        }
        if (nf > 2 && strcmp(fields[2], ".") != 0) {
            free(target->source_file);   target->source_file   = xdup(fields[2]);
        }
        if (nf > 3 && strcmp(fields[3], ".") != 0) {
            free(target->source_chrom);  target->source_chrom  = xdup(fields[3]);
        }
        if (nf > 4 && strcmp(fields[4], ".") != 0) target->source_start  = atoi(fields[4]);
        if (nf > 5 && strcmp(fields[5], ".") != 0) target->source_end    = atoi(fields[5]);
        if (nf > 6 && strcmp(fields[6], ".") != 0) target->source_stride = atoi(fields[6]);
        attached++;
    }
    fclose(fp);
    return attached;
}

void bpp_loci_free(BppLocus *arr, int n)
{
    if (!arr) return;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < arr[i].n_seqs; j++) {
            if (arr[i].seq_names) free(arr[i].seq_names[j]);
            if (arr[i].ids)       free(arr[i].ids[j]);
            if (arr[i].seqs)      free(arr[i].seqs[j]);
        }
        free(arr[i].seq_names);
        free(arr[i].ids);
        free(arr[i].seqs);
        free(arr[i].name);
        free(arr[i].source_kind);
        free(arr[i].source_file);
        free(arr[i].source_chrom);
    }
    free(arr);
}
