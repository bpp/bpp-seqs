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
} LocusAln;

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
