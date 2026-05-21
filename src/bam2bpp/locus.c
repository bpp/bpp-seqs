/*
 * locus.c  –  BED and Imap file parsers; per-locus sequence statistics
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

#include "bam2bpp.h"

/* -------------------------------------------------------------------------
 * BED parser
 * ---------------------------------------------------------------------- */

Locus *load_bed(const char *path, int *n_out)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Error: cannot open BED file '%s': ", path);
        perror(NULL);
        return NULL;
    }

    int     n   = 0;
    int     cap = 256;
    Locus  *loci = malloc(cap * sizeof(Locus));
    if (!loci) { fclose(f); return NULL; }

    char line[4096];
    int  lineno = 0;

    while (fgets(line, sizeof(line), f)) {
        lineno++;
        /* Skip comments and blank lines */
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;

        if (n == cap) {
            cap *= 2;
            Locus *tmp = realloc(loci, cap * sizeof(Locus));
            if (!tmp) { free(loci); fclose(f); return NULL; }
            loci = tmp;
        }

        char    chrom[512];
        int32_t start, end;
        char    name[512];

        /* BED3 or BED4+ */
        int fields = sscanf(line, "%511s\t%" SCNd32 "\t%" SCNd32 "\t%511s",
                            chrom, &start, &end, name);
        if (fields < 3) {
            fprintf(stderr, "Warning: skipping malformed BED line %d\n", lineno);
            continue;
        }

        if (start >= end) {
            fprintf(stderr, "Warning: degenerate locus at line %d (%s:%d-%d), skipping\n",
                    lineno, chrom, start, end);
            continue;
        }

        loci[n].chrom = strdup(chrom);
        loci[n].start = start;
        loci[n].end   = end;
        loci[n].name  = (fields >= 4) ? strdup(name)
                                       : malloc(64);

        if (fields < 4)
            snprintf(loci[n].name, 64, "locus_%d", n + 1);

        n++;
    }

    fclose(f);
    *n_out = n;
    return loci;
}

void free_loci(Locus *loci, int n)
{
    if (!loci) return;
    for (int i = 0; i < n; i++) {
        free(loci[i].chrom);
        free(loci[i].name);
    }
    free(loci);
}

/* -------------------------------------------------------------------------
 * Imap parser
 * ---------------------------------------------------------------------- */

ImapEntry *load_imap(const char *path, int *n_out)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Error: cannot open Imap file '%s': ", path);
        perror(NULL);
        return NULL;
    }

    int        n   = 0;
    int        cap = 64;
    ImapEntry *imap = malloc(cap * sizeof(ImapEntry));
    if (!imap) { fclose(f); return NULL; }

    char line[1024];

    while (fgets(line, sizeof(line), f)) {
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') continue;

        if (n == cap) {
            cap *= 2;
            ImapEntry *tmp = realloc(imap, cap * sizeof(ImapEntry));
            if (!tmp) { free(imap); fclose(f); return NULL; }
            imap = tmp;
        }

        char sample[256], pop[256];
        if (sscanf(line, "%255s %255s", sample, pop) == 2) {
            imap[n].sample_name = strdup(sample);
            imap[n].population  = strdup(pop);
            n++;
        }
    }

    fclose(f);
    *n_out = n;
    return imap;
}

void free_imap(ImapEntry *imap, int n)
{
    if (!imap) return;
    for (int i = 0; i < n; i++) {
        free(imap[i].sample_name);
        free(imap[i].population);
    }
    free(imap);
}

const char *imap_lookup(const ImapEntry *imap, int n_imap, const char *sample)
{
    for (int i = 0; i < n_imap; i++)
        if (strcmp(imap[i].sample_name, sample) == 0)
            return imap[i].population;
    return NULL;
}

/* -------------------------------------------------------------------------
 * Sequence statistics
 * ---------------------------------------------------------------------- */

/*
 * Fraction of 'N' characters across the entire alignment matrix.
 * A value of 1.0 means every position in every sequence is missing.
 */
double compute_missing_frac(char **seqs, int n_seqs, int locus_len)
{
    if (n_seqs == 0 || locus_len == 0) return 1.0;

    int64_t total   = (int64_t)n_seqs * locus_len;
    int64_t missing = 0;

    for (int i = 0; i < n_seqs; i++)
        for (int j = 0; j < locus_len; j++)
            if (seqs[i][j] == 'N') missing++;

    return (double)missing / (double)total;
}

/*
 * Count segregating sites: columns where at least two distinct non-N bases
 * appear across all sequences.
 */
int count_snps(char **seqs, int n_seqs, int locus_len)
{
    int snps = 0;

    for (int pos = 0; pos < locus_len; pos++) {
        char first = '\0';
        int  is_snp = 0;

        for (int i = 0; i < n_seqs && !is_snp; i++) {
            char b = seqs[i][pos];
            if (b == 'N') continue;
            if (!first) { first = b; continue; }
            if (b != first) is_snp = 1;
        }

        snps += is_snp;
    }

    return snps;
}
