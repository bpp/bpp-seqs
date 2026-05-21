/* output.h — formatted output (human-readable and JSON) for bpp-seqs. */
#ifndef BPP_SEQS_OUTPUT_H
#define BPP_SEQS_OUTPUT_H

#include <stdio.h>

#include "inspect.h"
#include "cross_validate.h"
#include "workflow.h"

/* Conversion results.  Populated only after a successful conversion. */
typedef struct {
    int has_results;
    char *out_sequences;
    char *out_imap;
    char *out_stats;
    int n_loci_input;
    int n_loci_passed;
    int n_loci_failed;
    int n_sequences;
    int n_too_short;
    int n_high_missing;
    int n_insufficient_snps;
    /* per-locus list (optional, for JSON loci array) */
    struct {
        char *name;
        int   length;
        int   n_snps;
        double missing_fraction;
        double mean_depth;
        char  *status;        /* passed | failed */
        char  *skip_reason;
    } *loci;
    int n_loci;
} ConversionResult;

void conversion_result_init(ConversionResult *r);
void conversion_result_free(ConversionResult *r);
void conversion_result_add_locus(ConversionResult *r,
                                 const char *name, int length, int n_snps,
                                 double missing_fraction, double mean_depth,
                                 const char *status, const char *skip_reason);

/* Status string for the top-level JSON: incomplete | ready | complete | error */
const char *status_string(const WorkflowDecision *d, int conversion_done, int error);

/* Emit human-readable output to fp. */
void print_human(FILE *fp,
                 FileInfo **files, int n_files,
                 const CrossValidation *cv,
                 const WorkflowDecision *d,
                 const ConversionResult *cr,
                 const char *bpp_seqs_version);

/* Emit JSON to fp.  indent>0 for pretty-printing. */
void print_json(FILE *fp, int indent,
                FileInfo **files, int n_files,
                const CrossValidation *cv,
                const WorkflowDecision *d,
                const ConversionResult *cr,
                const char *recommended_command,
                const char *bpp_seqs_version);

#endif /* BPP_SEQS_OUTPUT_H */
