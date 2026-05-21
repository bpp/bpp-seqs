/* sanity.h — structural sanity checks for BPP sequence data.
 *
 * These are called by every converter right before the BPP writer is
 * invoked.  They catch invariants the writer downstream silently relies on:
 *   - all sequences at a locus must have the same length
 *   - every locus must have at least one sequence
 *   - individual ids (the part after the last '^' in each seq name) must
 *     be unique within a locus — required by BPP
 *
 * Invalid sequence characters are sanitized to 'N' in place and counted
 * as a warning, not a fatal error.
 */
#ifndef BPP_SEQS_SANITY_H
#define BPP_SEQS_SANITY_H

#include <stdio.h>

typedef struct {
    int n_length_mismatches;     /* loci where seqs disagreed in length    */
    int n_duplicate_ids;         /* loci with non-unique individual ids    */
    int n_empty_loci;            /* loci with n_seqs <= 0                  */
    int n_invalid_chars;         /* characters rewritten to 'N'            */
    int n_critical;              /* sum of fatal-class events              */
} SanityReport;

void sanity_report_init(SanityReport *r);

/* Return the individual id portion of a sequence name: text after the
 * last '^', or the whole name if no caret is present. Caller frees. */
char *sanity_individual_id(const char *seq_name);

/* Validate a single locus's alignment. Sanitizes invalid characters in
 * seqs[] to 'N' in place. Critical errors are logged to err.
 *
 * Returns:
 *    0  → OK (proceed)
 *   -1  → critical failure (do not write this locus)
 */
int sanity_check_locus(const char *locus_name,
                       char **seq_names, char **seqs,
                       int n_seqs, int locus_len,
                       SanityReport *report,
                       FILE *err);

/* Imap coverage: every id present in the converted output should have an
 * Imap row, and every Imap row should correspond to some id used in the
 * output. Emits one warning line per orphan to err.
 *
 * Returns the total number of orphans (both directions). 0 = clean. */
int sanity_check_imap_coverage(char **ids, int n_ids,
                               char **imap_samples, int n_imap,
                               FILE *err);

#endif
