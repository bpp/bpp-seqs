/* cross_validate.h — cross-file consistency checks. */
#ifndef BPP_SEQS_CROSS_VALIDATE_H
#define BPP_SEQS_CROSS_VALIDATE_H

#include "inspect.h"

typedef struct {
    char *code;        /* e.g. CHROMOSOME_MISMATCH */
    char *severity;    /* info | warning | error | critical */
    char *file;        /* may be NULL */
    char *message;
} ValidationIssue;

typedef struct {
    int  bam_reference_consistent;
    int  bam_reference_matches_fasta;
    int  bed_chromosomes_in_bams;
    int  imap_samples_in_bams;
    int  bam_samples_in_imap;

    char **unmatched_bam_samples;
    int    n_unmatched_bam_samples;
    char **unmatched_imap_samples;
    int    n_unmatched_imap_samples;

    ValidationIssue *issues;
    int              n_issues;
} CrossValidation;

CrossValidation *cross_validate(FileInfo **files, int n);
void              cross_validation_free(CrossValidation *cv);

#endif /* BPP_SEQS_CROSS_VALIDATE_H */
