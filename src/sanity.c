#include "sanity.h"

#include <ctype.h>
#include <stdlib.h>
#include <string.h>

void sanity_report_init(SanityReport *r) { memset(r, 0, sizeof(*r)); }

const char *sanity_seq_id(const char *seq_name)
{
    if (!seq_name) return NULL;
    const char *caret = strrchr(seq_name, '^');
    return caret ? caret + 1 : seq_name;
}

char *sanity_individual_id(const char *seq_name)
{
    if (!seq_name) return NULL;
    const char *id = sanity_seq_id(seq_name);
    size_t n = strlen(id);
    char *r = (char *)malloc(n + 1);
    if (!r) return NULL;
    memcpy(r, id, n + 1);
    return r;
}

/* Acceptable nucleotide / ambiguity / gap characters. */
static int is_valid_seq_char(char c)
{
    char u = (char)toupper((unsigned char)c);
    switch (u) {
        case 'A': case 'C': case 'G': case 'T': case 'U':
        case 'N': case '-': case '?': case '.':
        case 'M': case 'R': case 'W': case 'S': case 'Y': case 'K':
        case 'B': case 'D': case 'H': case 'V': case 'X':
            return 1;
        default:
            return 0;
    }
}

int sanity_check_locus(const char *locus_name,
                       char **seq_names, char **seqs,
                       int n_seqs, int locus_len,
                       SanityReport *report,
                       FILE *err)
{
    int critical = 0;

    /* Check #3: n_seqs >= 1 */
    if (n_seqs <= 0) {
        report->n_empty_loci++;
        report->n_critical++;
        if (err) fprintf(err,
            "ERROR: locus '%s' has no sequences (n_seqs=%d); skipping.\n",
            locus_name ? locus_name : "?", n_seqs);
        critical = 1;
        return -1;
    }

    /* Check #1: equal lengths */
    for (int i = 0; i < n_seqs; i++) {
        int len = (int)strlen(seqs[i]);
        if (len != locus_len) {
            report->n_length_mismatches++;
            report->n_critical++;
            if (err) fprintf(err,
                "ERROR: locus '%s' sequence #%d ('%s') has length %d but "
                "locus length is %d; skipping locus.\n",
                locus_name ? locus_name : "?", i + 1,
                seq_names[i] ? seq_names[i] : "?", len, locus_len);
            critical = 1;
        }
    }
    if (critical) return -1;

    /* Check #2: unique individual ids (post-caret) */
    char **ids = (char **)calloc((size_t)n_seqs, sizeof(char *));
    for (int i = 0; i < n_seqs; i++) {
        ids[i] = sanity_individual_id(seq_names[i]);
    }
    for (int i = 0; i < n_seqs && !critical; i++) {
        for (int j = i + 1; j < n_seqs; j++) {
            if (ids[i] && ids[j] && strcmp(ids[i], ids[j]) == 0) {
                report->n_duplicate_ids++;
                report->n_critical++;
                if (err) fprintf(err,
                    "ERROR: locus '%s' has duplicate individual id '^%s' "
                    "(sequences #%d and #%d); BPP requires unique ids.\n",
                    locus_name ? locus_name : "?", ids[i], i + 1, j + 1);
                critical = 1;
                break;
            }
        }
    }
    for (int i = 0; i < n_seqs; i++) free(ids[i]);
    free(ids);
    if (critical) return -1;

    /* Check #7: sanitize invalid characters. Non-critical. */
    int invalid_here = 0;
    for (int i = 0; i < n_seqs; i++) {
        for (int k = 0; k < locus_len; k++) {
            if (!is_valid_seq_char(seqs[i][k])) {
                seqs[i][k] = 'N';
                invalid_here++;
            }
        }
    }
    if (invalid_here > 0) {
        report->n_invalid_chars += invalid_here;
        if (err) fprintf(err,
            "WARNING: locus '%s': %d invalid character%s sanitized to 'N'.\n",
            locus_name ? locus_name : "?",
            invalid_here, invalid_here == 1 ? "" : "s");
    }

    return 0;
}

int sanity_check_imap_coverage(char **ids, int n_ids,
                               char **imap_samples, int n_imap,
                               FILE *err)
{
    int orphans = 0;

    /* Each sequence-id should have an Imap row */
    for (int i = 0; i < n_ids; i++) {
        int found = 0;
        for (int j = 0; j < n_imap; j++) {
            if (strcmp(ids[i], imap_samples[j]) == 0) { found = 1; break; }
        }
        if (!found) {
            orphans++;
            if (err) fprintf(err,
                "WARNING: individual id '%s' present in sequence file but "
                "not in Imap. BPP will reject this input.\n", ids[i]);
        }
    }

    /* Each Imap row should be used by some sequence */
    for (int j = 0; j < n_imap; j++) {
        int found = 0;
        for (int i = 0; i < n_ids; i++) {
            if (strcmp(imap_samples[j], ids[i]) == 0) { found = 1; break; }
        }
        if (!found) {
            orphans++;
            if (err) fprintf(err,
                "WARNING: Imap entry '%s' has no matching sequence in "
                "any locus; it will be unused.\n", imap_samples[j]);
        }
    }

    return orphans;
}
