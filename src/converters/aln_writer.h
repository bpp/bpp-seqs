/* aln_writer.h — shared helpers for non-bam2bpp converters. */
#ifndef BPP_SEQS_ALN_WRITER_H
#define BPP_SEQS_ALN_WRITER_H

#include <stdio.h>

#include "inspect.h"
#include "output.h"

double compute_missing_frac(char **seqs, int n_seqs, int locus_len);
int    count_snps          (char **seqs, int n_seqs, int locus_len);

/* Per-locus alignment.  Ownership transfers to write_alignment_outputs
 * (which frees the arrays once written). */
typedef struct {
    char  *name;
    int    length;
    int    n_seqs;
    char **seq_names;     /* size n_seqs */
    char **seqs;          /* size n_seqs, length+1 NUL-terminated */

    /* Provenance (for <prefix>.loci.tsv).  All fields are optional;
     * NULL strings and -1 ints mean "not applicable". */
    char  *source_kind;   /* "BED" | "CHARSET" | "MSA" | "PHYLIP" */
    char  *source_file;
    char  *source_chrom;
    int    source_start;  /* 1-based inclusive, -1 if N/A */
    int    source_end;    /* 1-based inclusive, -1 if N/A */
    int    source_stride; /* default 1 */
} LocusAln;

/* Shared by both converter and bam2bpp paths: write the per-locus
 * provenance to <prefix>.loci.tsv.  Returns 0 on success. */
typedef struct {
    const char *name;
    const char *kind;
    const char *file;
    const char *chrom;
    int start;
    int end;
    int stride;
    int length;
    int n_seqs;
} LocusProv;

int write_loci_tsv(const char *out_prefix,
                   const LocusProv *items, int n);

/* Apply QC filters, write <prefix>.txt and <prefix>.stats.tsv, and emit a
 * BPP-format Imap based on imap_fi. Populates cr with summary information.
 * `loci` and its strings are freed before return. */
int write_alignment_outputs(const char *out_prefix,
                            LocusAln *loci, int n_loci,
                            FileInfo *imap_fi,
                            int min_length, double max_missing,
                            int min_snps, int keep_invariant,
                            ConversionResult *cr);

#endif
