/*
 * writer.c  –  output file writers
 *
 * BPP sequence file format (per https://bpp.github.io/bpp-manual/)
 * --------------------------------------------------------------
 * Each locus is an independent block with this layout:
 *
 *     <n_seqs> <n_sites>
 *
 *     ^<id1>    ACGT...
 *     ^<id2>    ACGT...
 *
 * The id following the caret (^) is the individual ID; it must be unique
 * within the locus. The Imap file maps those ids to species/populations.
 *
 * There is no file-level header; there is no "^locus_name" delimiter line.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include "bam2bpp.h"

/* Transform an internal sample-name like "ind1^1" into a BPP-safe
 * individual id ("ind1_1") so multiple haplotypes of the same sample
 * don't collide with haplotypes of other samples.  Caller must free. */
static char *to_bpp_id(const char *name)
{
    if (!name) return NULL;
    size_t n = strlen(name);
    char *r = (char *)malloc(n + 1);
    if (!r) return NULL;
    for (size_t i = 0; i < n; i++) {
        r[i] = (name[i] == '^') ? '_' : name[i];
    }
    r[n] = '\0';
    return r;
}

/* -------------------------------------------------------------------------
 * BPP sequence file
 * ---------------------------------------------------------------------- */

void write_bpp_file(const char *prefix,
                    LocusResult *results, int n_results, int n_seqs)
{
    char path[4096];
    snprintf(path, sizeof(path), "%s.txt", prefix);

    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "Error: cannot write '%s': ", path);
        perror(NULL);
        return;
    }

    for (int li = 0; li < n_results; li++) {
        LocusResult *r = &results[li];

        /* Compute id transforms + column width for nicely-aligned output */
        char **ids = (char **)calloc((size_t)r->n_seqs, sizeof(char *));
        int max_id_len = 0;
        for (int i = 0; i < r->n_seqs; i++) {
            ids[i] = to_bpp_id(r->seq_names[i]);
            int len = (int)strlen(ids[i]) + 1;  /* +1 for the caret */
            if (len > max_id_len) max_id_len = len;
        }
        int col_width = max_id_len + 4;   /* 4-space gap before sequence */

        /* Per-locus header */
        if (li > 0) fputc('\n', fp);
        fprintf(fp, "%d %d\n\n", r->n_seqs, r->locus_len);

        for (int i = 0; i < r->n_seqs; i++) {
            /* Build "^<id>" then pad to col_width with spaces */
            int written = fprintf(fp, "^%s", ids[i]);
            int pad = col_width - written;
            if (pad < 1) pad = 1;
            for (int p = 0; p < pad; p++) fputc(' ', fp);
            fprintf(fp, "%s\n", r->seqs[i]);
        }
        for (int i = 0; i < r->n_seqs; i++) free(ids[i]);
        free(ids);
    }

    (void)n_seqs;  /* n_seqs is per-locus in the new format */

    fclose(fp);
    fprintf(stderr, "Wrote BPP sequence file: %s\n", path);
}

/* -------------------------------------------------------------------------
 * BPP Imap file
 *
 * One row per individual ID present in the sequence file.
 * Two whitespace-separated columns: <id>  <population>.
 * ---------------------------------------------------------------------- */

void write_imap_file(const char *prefix,
                     BamFile **bams, int n_bams,
                     const ImapEntry *imap, int n_imap,
                     Phasing phasing)
{
    char path[4096];
    snprintf(path, sizeof(path), "%s.imap", prefix);

    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "Error: cannot write '%s': ", path);
        perror(NULL);
        return;
    }

    for (int i = 0; i < n_bams; i++) {
        const char *pop = imap_lookup(imap, n_imap, bams[i]->sample);
        if (!pop) {
            fprintf(stderr,
                    "Warning: sample '%s' not in Imap; using 'unknown'\n",
                    bams[i]->sample);
            pop = "unknown";
        }

        if (phasing == PHASE_SPLIT) {
            /* Two unique ids per individual: <sample>_1 and <sample>_2. */
            fprintf(fp, "%s_1\t%s\n", bams[i]->sample, pop);
            fprintf(fp, "%s_2\t%s\n", bams[i]->sample, pop);
        } else {
            fprintf(fp, "%s\t%s\n", bams[i]->sample, pop);
        }
    }

    fclose(fp);
    fprintf(stderr, "Wrote Imap file: %s\n", path);
}

/* -------------------------------------------------------------------------
 * Statistics TSV
 * ---------------------------------------------------------------------- */

void write_stats_file(const char *prefix, const LocusStat *stats, int n)
{
    char path[4096];
    snprintf(path, sizeof(path), "%s.stats.tsv", prefix);

    FILE *fp = fopen(path, "w");
    if (!fp) {
        fprintf(stderr, "Error: cannot write '%s': ", path);
        perror(NULL);
        return;
    }

    fprintf(fp, "locus\tlength\tsnps\tmissing_frac\tmean_depth\tstatus\n");

    for (int i = 0; i < n; i++) {
        fprintf(fp, "%s\t%d\t%d\t%.4f\t%.2f\t%s\n",
                stats[i].locus_name,
                stats[i].locus_len,
                stats[i].n_snps,
                stats[i].missing_frac,
                stats[i].mean_depth,
                stats[i].skip_reason ? stats[i].skip_reason : "used");
    }

    fclose(fp);
    fprintf(stderr, "Wrote stats file: %s\n", path);
}
