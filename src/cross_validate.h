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

/* Count issues whose code names an input pairing that cannot be reconciled:
 * reads aligned to a reference other than the one supplied, or a set of
 * BAMs aligned to different references as each other. Converting under either
 * would emit sequences whose coordinates mean nothing, so callers refuse.
 * Other error-severity issues (a missing BED chromosome, say) still only
 * warn -- they narrow the output rather than invalidate it. */
int cross_validation_n_blocking(const CrossValidation *cv);

#endif /* BPP_SEQS_CROSS_VALIDATE_H */
