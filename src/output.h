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

        /* Provenance for JSON output (parallels .loci.tsv columns) */
        char *source_kind;
        char *source_file;
        char *source_chrom;
        int   source_start;
        int   source_end;
        int   source_stride;
    } *loci;
    int n_loci;
    char *out_loci;            /* path to <prefix>.loci.tsv */
} ConversionResult;

void conversion_result_init(ConversionResult *r);
void conversion_result_free(ConversionResult *r);
void conversion_result_add_locus(ConversionResult *r,
                                 const char *name, int length, int n_snps,
                                 double missing_fraction, double mean_depth,
                                 const char *status, const char *skip_reason);

/* Set provenance on the most-recently-added locus (call right after
 * conversion_result_add_locus). Any string may be NULL; ints default to
 * start=-1 end=-1 stride=1 if 0 is passed. */
void conversion_result_set_locus_source(ConversionResult *r,
                                        const char *source_kind,
                                        const char *source_file,
                                        const char *source_chrom,
                                        int source_start,
                                        int source_end,
                                        int source_stride);

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
