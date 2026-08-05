#include "aln_writer.h"
#include "converters.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

/* Parse one PHYLIP locus block from the open gz stream.
 *
 *   Input:  *header_line is either the already-read header line
 *           ("nseq nsites"), or NULL to read it from the stream.
 *
 *   Output: out is filled.  Returns:
 *       0   on success, with the trailing block separator consumed
 *      -1   end of file (no more loci)
 *      -2   on malformed content
 *
 * Supports both sequential and interleaved layouts: data rows for the
 * current locus are accumulated round-robin until each sample has
 * collected exactly `nsites` characters.
 */
static int parse_one_locus(gzFile gz, char *carry_line, size_t carry_cap,
                           LocusAln *out, int locus_index)
{
    char line[1 << 16];

    /* Obtain the header line: either passed in via carry, or read fresh. */
    if (carry_line[0] == '\0') {
        do {
            if (gzgets(gz, line, sizeof(line)) == NULL) return -1;
            /* skip blank lines */
            int empty = 1;
            for (char *p = line; *p; p++) if (!isspace((unsigned char)*p)) { empty = 0; break; }
            if (!empty) break;
        } while (1);
    } else {
        strncpy(line, carry_line, sizeof(line) - 1);
        line[sizeof(line) - 1] = '\0';
        carry_line[0] = '\0';
    }

    int nseq = 0, nsites = 0;
    if (sscanf(line, "%d %d", &nseq, &nsites) != 2 || nseq <= 0 || nsites <= 0) {
        return -2;
    }

    char **names = (char **)calloc((size_t)nseq, sizeof(char *));
    char **seqs  = (char **)calloc((size_t)nseq, sizeof(char *));
    int   *lens  = (int *)  calloc((size_t)nseq, sizeof(int));
    if (!names || !seqs || !lens) { free(names); free(seqs); free(lens); return -2; }
    for (int i = 0; i < nseq; i++) {
        seqs[i] = (char *)malloc((size_t)nsites + 1);
        if (!seqs[i]) return -2;
        seqs[i][0] = '\0';
    }

    int row = 0;
    int cur_named = -1;          /* sample last introduced by an unindented row */
    int blank_since_name = 0;    /* a blank line has closed that row's block */
    int finished_first_block = 0;
    int all_full = 0;
    int excess_chars = 0;
    int excess_lines = 0;  /* non-blank lines arriving after all samples are full */

    while (!all_full) {
        if (gzgets(gz, line, sizeof(line)) == NULL) break;

        int empty = 1;
        for (char *p = line; *p; p++) if (!isspace((unsigned char)*p)) { empty = 0; break; }
        if (empty) {
            if (row > 0) finished_first_block = 1;
            row = 0;
            /* A blank line closes an interleaved block: indented rows after it
             * start a new block and cycle through the samples again, rather
             * than continuing the sample named above. */
            blank_since_name = 1;
            continue;
        }

        /* If the line "looks like" a new locus header (two ints, and no
         * tokens that look like sequence), AND we are between loci (row==0
         * and all samples filled), then this is a new locus header. Stash
         * it in carry_line and stop. */
        if (row == 0 && finished_first_block) {
            int peek_nseq, peek_nsites; char rest = '\0';
            int parsed = sscanf(line, " %d %d %c", &peek_nseq, &peek_nsites, &rest);
            int is_header = (parsed >= 2) && peek_nseq > 0 && peek_nsites > 0 &&
                            (parsed == 2 || rest == '\n' || rest == '\r' ||
                             isspace((unsigned char)rest) || rest == '\0');
            if (is_header) {
                int filled = 1;
                for (int i = 0; i < nseq; i++) if (lens[i] < nsites) { filled = 0; break; }
                if (filled) {
                    snprintf(carry_line, carry_cap, "%s", line);
                    break;
                }
            }
        }

        const char *seq_start;
        int target;
        int advance = 1;
        if (phylip_row_is_continuation(line) && cur_named >= 0 && !blank_since_name) {
            /* Wrapped sequential layout: an indented row that follows its
             * named row directly keeps filling that same sample, and does not
             * consume a taxon slot. (Indented rows *after* a blank line are an
             * interleaved block instead, and fall through below.) */
            seq_start = line;
            target    = cur_named;
            advance   = 0;
        } else if (!finished_first_block && row < nseq) {
            char *p = line;
            char *name_end = p;
            while (*name_end && !isspace((unsigned char)*name_end)) name_end++;
            size_t nlen = (size_t)(name_end - p);
            if (nlen > 0) {
                names[row] = (char *)malloc(nlen + 1);
                memcpy(names[row], p, nlen); names[row][nlen] = '\0';
            } else {
                names[row] = strdup("");
            }
            seq_start = name_end;
            target = row;
            cur_named = row;
            blank_since_name = 0;
        } else {
            /* Interleaved layout: unindented rows after the first block cycle
             * through the samples in order. */
            seq_start = line;
            target = row;
        }

        for (const char *q = seq_start; *q && *q != '\n' && *q != '\r'; q++) {
            if (isspace((unsigned char)*q)) continue;
            if (target < nseq && lens[target] < nsites) {
                seqs[target][lens[target]++] = (char)toupper((unsigned char)*q);
            } else if (target < nseq) {
                /* Sample already at declared n_sites — silently dropping
                 * these characters would mask a header/data mismatch. */
                excess_chars++;
            }
        }
        if (advance) {
            row++;
            if (row >= nseq) { row = 0; finished_first_block = 1; }
        }

        /* Check if every sample is full */
        all_full = 1;
        for (int i = 0; i < nseq; i++) if (lens[i] < nsites) { all_full = 0; break; }
    }

    /* After all samples are full, peek for non-blank, non-header lines —
     * those would indicate extra taxa rows beyond declared n_seqs. */
    while (gzgets(gz, line, sizeof(line)) != NULL) {
        int empty = 1;
        for (char *p = line; *p; p++) if (!isspace((unsigned char)*p)) { empty = 0; break; }
        if (empty) continue;
        int peek_nseq, peek_nsites; char rest = '\0';
        int parsed = sscanf(line, " %d %d %c", &peek_nseq, &peek_nsites, &rest);
        int is_header = (parsed >= 2) && peek_nseq > 0 && peek_nsites > 0 &&
                        (parsed == 2 || rest == '\n' || rest == '\r' ||
                         isspace((unsigned char)rest) || rest == '\0');
        if (is_header) {
            snprintf(carry_line, carry_cap, "%s", line);
            break;
        }
        excess_lines++;
    }

    if (excess_chars > 0) {
        fprintf(stderr,
            "WARNING: PHYLIP locus #%d: %d character%s beyond declared "
            "n_sites=%d were dropped (header/data mismatch).\n",
            locus_index, excess_chars, excess_chars == 1 ? "" : "s", nsites);
    }
    if (excess_lines > 0) {
        fprintf(stderr,
            "WARNING: PHYLIP locus #%d: %d unexpected data line%s after "
            "n_seqs=%d samples were filled (extra taxa or stray content).\n",
            locus_index, excess_lines, excess_lines == 1 ? "" : "s", nseq);
    }

    /* Validate */
    int ok = 1;
    for (int i = 0; i < nseq; i++) {
        if (!names[i]) names[i] = strdup("");
        seqs[i][lens[i]] = '\0';
        if (lens[i] != nsites) ok = 0;
    }
    free(lens);
    if (!ok) {
        for (int i = 0; i < nseq; i++) { free(names[i]); free(seqs[i]); }
        free(names); free(seqs);
        return -2;
    }

    char nm[64];
    snprintf(nm, sizeof(nm), "locus%d", locus_index);
    out->name      = strdup(nm);
    out->length    = nsites;
    out->n_seqs    = nseq;
    out->seq_names = names;
    out->seqs      = seqs;
    out->source_kind   = strdup("PHYLIP");
    out->source_file   = NULL;   /* set by caller (knows the path) */
    out->source_chrom  = NULL;
    out->source_start  = -1;
    out->source_end    = -1;
    out->source_stride = 1;
    return 0;
}

/* Load all loci from a single PHYLIP file. Caller frees returned array
 * by freeing each LocusAln's contents and then the array itself. */
static int load_phylip_multi(const char *path, LocusAln **out, int *n_out,
                             int starting_index)
{
    gzFile gz = gzopen(path, "rb");
    if (!gz) return -1;

    LocusAln *arr = NULL;
    int n = 0, cap = 0;
    char carry[1 << 16];
    carry[0] = '\0';
    int idx = starting_index;

    for (;;) {
        LocusAln loc; memset(&loc, 0, sizeof(loc));
        int rc = parse_one_locus(gz, carry, sizeof(carry), &loc, idx);
        if (rc == -1) break;       /* EOF */
        if (rc == -2) { gzclose(gz); return -2; }
        if (n >= cap) { cap = cap ? cap * 2 : 8;
            arr = (LocusAln *)realloc(arr, sizeof(LocusAln) * (size_t)cap); }
        arr[n++] = loc;
        idx++;
    }
    gzclose(gz);
    *out = arr;
    *n_out = n;
    return 0;
}

int convert_phylip(FileInfo **files, int n_files,
                   FileInfo *imap_fi,
                   const char *out_prefix,
                   const ConvertOpts *opts,
                   ConversionResult *cr)
{
    int n_phy = 0;
    for (int i = 0; i < n_files; i++) if (files[i]->ft == BS_PHYLIP) n_phy++;
    if (n_phy == 0 || !imap_fi) {
        fprintf(stderr, "Error: phylip2bpp requires at least one PHYLIP file and an Imap.\n");
        return 1;
    }

    LocusAln *all = NULL;
    int n_all = 0;
    int idx = 1;

    for (int i = 0; i < n_files; i++) {
        if (files[i]->ft != BS_PHYLIP) continue;
        LocusAln *chunk = NULL;
        int nc = 0;
        int rc = load_phylip_multi(files[i]->path, &chunk, &nc, idx);
        if (rc != 0 || nc == 0) {
            fprintf(stderr, "Error: failed to load PHYLIP '%s'\n", files[i]->path);
            return 1;
        }
        /* Loci may legitimately carry different numbers of sequences: under
         * the multispecies coalescent each locus is independent, so a study
         * that sequenced different individuals at different loci is normal.
         * BPP's own frogs example has 21/28/28/24/30. Each locus writes its
         * own "<n_seqs> <n_sites>" header, so nothing here needs them equal. */
        all = (LocusAln *)realloc(all, sizeof(LocusAln) * (size_t)(n_all + nc));
        for (int j = 0; j < nc; j++) {
            all[n_all + j] = chunk[j];
            all[n_all + j].source_file = strdup(files[i]->path);
        }
        n_all += nc;
        idx += nc;
        free(chunk);
    }

    return write_alignment_outputs(out_prefix, all, n_all, imap_fi,
                                   opts->min_length, opts->max_missing,
                                   opts->min_snps, opts->keep_invariant, cr);
}
