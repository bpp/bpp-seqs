#include "cross_validate.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *xstrdup_local(const char *s)
{
    if (!s) return NULL;
    size_t n = strlen(s);
    char *r = (char *)malloc(n + 1);
    if (!r) return NULL;
    memcpy(r, s, n + 1);
    return r;
}

static void add_issue(CrossValidation *cv,
                      const char *code, const char *severity,
                      const char *file, const char *message)
{
    cv->issues = (ValidationIssue *)realloc(cv->issues,
                    sizeof(ValidationIssue) * (size_t)(cv->n_issues + 1));
    ValidationIssue *iss = &cv->issues[cv->n_issues++];
    iss->code     = xstrdup_local(code);
    iss->severity = xstrdup_local(severity);
    iss->file     = xstrdup_local(file);
    iss->message  = xstrdup_local(message);
}

static int seqref_arrays_equal(const SeqRef *a, int na, const SeqRef *b, int nb)
{
    if (na != nb) return 0;
    for (int i = 0; i < na; i++) {
        if (!a[i].name || !b[i].name) return 0;
        if (strcmp(a[i].name, b[i].name) != 0) return 0;
        if (a[i].length != b[i].length) return 0;
    }
    return 1;
}

static int seqref_has(const SeqRef *a, int na, const char *name)
{
    for (int i = 0; i < na; i++) {
        if (a[i].name && strcmp(a[i].name, name) == 0) return 1;
    }
    return 0;
}

static int sample_in_bams(FileInfo **files, int n, const char *name)
{
    for (int i = 0; i < n; i++) {
        FileInfo *fi = files[i];
        if ((fi->ft == BS_BAM || fi->ft == BS_CRAM) && fi->sample_name &&
            strcmp(fi->sample_name, name) == 0) return 1;
    }
    return 0;
}

/* True if any sequence-bearing file lists this sample/sequence name. */
static int sample_in_seq_files(FileInfo **files, int n, const char *name)
{
    for (int i = 0; i < n; i++) {
        FileInfo *fi = files[i];
        if (fi->ft == BS_VCF || fi->ft == BS_GVCF) {
            for (int j = 0; j < fi->n_sample_names; j++) {
                if (fi->sample_names[j] && strcmp(fi->sample_names[j], name) == 0) return 1;
            }
        } else if (fi->ft == BS_FASTA_MSA || fi->ft == BS_PHYLIP || fi->ft == BS_NEXUS) {
            for (int j = 0; j < fi->n_sequence_names; j++) {
                if (fi->sequence_names[j] && strcmp(fi->sequence_names[j], name) == 0) return 1;
            }
        }
    }
    return 0;
}

static int sample_in_imap(FileInfo **files, int n, const char *name)
{
    for (int i = 0; i < n; i++) {
        FileInfo *fi = files[i];
        if (fi->ft != BS_IMAP) continue;
        for (int j = 0; j < fi->imap_n_entries; j++) {
            if (fi->imap_samples[j] && strcmp(fi->imap_samples[j], name) == 0) return 1;
        }
    }
    return 0;
}

CrossValidation *cross_validate(FileInfo **files, int n)
{
    CrossValidation *cv = (CrossValidation *)calloc(1, sizeof(CrossValidation));
    if (!cv) return NULL;
    cv->bam_reference_consistent   = 1;
    cv->bam_reference_matches_fasta = 1;
    cv->bed_chromosomes_in_bams    = 1;
    cv->imap_samples_in_bams       = 1;
    cv->bam_samples_in_imap        = 1;

    /* 1. BAM @SQ consistency across BAMs */
    FileInfo *first_bam = NULL;
    for (int i = 0; i < n; i++) {
        FileInfo *fi = files[i];
        if (fi->ft == BS_BAM || fi->ft == BS_CRAM) {
            if (!first_bam) first_bam = fi;
            else if (!seqref_arrays_equal(first_bam->seq_refs, first_bam->n_seq_refs,
                                          fi->seq_refs, fi->n_seq_refs)) {
                cv->bam_reference_consistent = 0;
                char msg[512];
                snprintf(msg, sizeof(msg),
                    "BAMs aligned to different references: %s and %s have different @SQ sets.",
                    first_bam->path, fi->path);
                add_issue(cv, "BAM_REF_INCONSISTENT", "error", fi->path, msg);
            }
        }
    }

    /* 2. BAM @SQ vs FASTA reference */
    FileInfo *ref_fa = NULL;
    for (int i = 0; i < n; i++) {
        if (files[i]->ft == BS_FASTA_REFERENCE) { ref_fa = files[i]; break; }
    }
    if (ref_fa && first_bam) {
        /* For each BAM @SQ name, require it to be one of the FASTA sequence
         * names. We don't compare lengths against the FASTA (would need to
         * parse it again); just names. */
        for (int i = 0; i < first_bam->n_seq_refs; i++) {
            const char *sq = first_bam->seq_refs[i].name;
            int found = 0;
            for (int j = 0; j < ref_fa->n_sequence_names; j++) {
                if (ref_fa->sequence_names[j] && strcmp(ref_fa->sequence_names[j], sq) == 0) {
                    found = 1; break;
                }
            }
            if (!found) {
                cv->bam_reference_matches_fasta = 0;
                char msg[512];
                snprintf(msg, sizeof(msg),
                    "BAM @SQ name '%s' not found in FASTA reference '%s'.",
                    sq, ref_fa->path);
                add_issue(cv, "BAM_FASTA_MISMATCH", "error", ref_fa->path, msg);
                break;
            }
        }
    }

    /* 3. BED chromosomes in BAM @SQ */
    FileInfo *bed = NULL;
    for (int i = 0; i < n; i++) {
        if (files[i]->ft == BS_BED) { bed = files[i]; break; }
    }
    if (bed && first_bam) {
        for (int i = 0; i < bed->n_chromosomes; i++) {
            if (!seqref_has(first_bam->seq_refs, first_bam->n_seq_refs, bed->chromosomes[i])) {
                cv->bed_chromosomes_in_bams = 0;
                char msg[512];
                snprintf(msg, sizeof(msg),
                    "BED chromosome '%s' not found in BAM @SQ headers.",
                    bed->chromosomes[i]);
                add_issue(cv, "CHROMOSOME_MISMATCH", "error", bed->path, msg);
            }
        }
    }

    /* 4. Sample-name reconciliation between sequence files and IMAP.
     * For BAM workflows we compare to BAM @RG SM; otherwise we compare
     * to sequence/sample names from MSA/PHYLIP/NEXUS/VCF. */
    int has_imap = 0;
    for (int i = 0; i < n; i++) if (files[i]->ft == BS_IMAP) { has_imap = 1; break; }

    /* Pick which sample-source family is present. */
    int have_seq_source = 0;
    for (int i = 0; i < n; i++) {
        FileType t = files[i]->ft;
        if (t == BS_FASTA_MSA || t == BS_PHYLIP || t == BS_NEXUS ||
            t == BS_VCF || t == BS_GVCF) { have_seq_source = 1; break; }
    }

    if (has_imap && first_bam) {
        /* unmatched_bam_samples: BAM SM tags not in Imap */
        for (int i = 0; i < n; i++) {
            FileInfo *fi = files[i];
            if ((fi->ft == BS_BAM || fi->ft == BS_CRAM) && fi->sample_name) {
                if (!sample_in_imap(files, n, fi->sample_name)) {
                    cv->bam_samples_in_imap = 0;
                    cv->unmatched_bam_samples = (char **)realloc(
                        cv->unmatched_bam_samples,
                        sizeof(char *) * (size_t)(cv->n_unmatched_bam_samples + 1));
                    cv->unmatched_bam_samples[cv->n_unmatched_bam_samples++] =
                        xstrdup_local(fi->sample_name);
                    char msg[256];
                    snprintf(msg, sizeof(msg),
                        "BAM sample '%s' not present in Imap.", fi->sample_name);
                    add_issue(cv, "UNMATCHED_BAM_SAMPLE", "warning", fi->path, msg);
                }
            }
        }
        /* unmatched_imap_samples: Imap names not in any BAM */
        for (int i = 0; i < n; i++) {
            FileInfo *fi = files[i];
            if (fi->ft != BS_IMAP) continue;
            for (int j = 0; j < fi->imap_n_entries; j++) {
                const char *s = fi->imap_samples[j];
                if (!sample_in_bams(files, n, s)) {
                    cv->imap_samples_in_bams = 0;
                    cv->unmatched_imap_samples = (char **)realloc(
                        cv->unmatched_imap_samples,
                        sizeof(char *) * (size_t)(cv->n_unmatched_imap_samples + 1));
                    cv->unmatched_imap_samples[cv->n_unmatched_imap_samples++] =
                        xstrdup_local(s);
                    char msg[256];
                    snprintf(msg, sizeof(msg),
                        "Imap sample '%s' not present in any BAM.", s);
                    add_issue(cv, "UNMATCHED_IMAP_SAMPLE", "warning", fi->path, msg);
                }
            }
        }
    } else if (has_imap && have_seq_source) {
        /* Non-BAM workflow: reconcile against sequence/sample names from
         * MSA/PHYLIP/NEXUS/VCF instead. */
        for (int i = 0; i < n; i++) {
            FileInfo *fi = files[i];
            if (fi->ft != BS_IMAP) continue;
            for (int j = 0; j < fi->imap_n_entries; j++) {
                const char *s = fi->imap_samples[j];
                if (!sample_in_seq_files(files, n, s)) {
                    cv->imap_samples_in_bams = 0;
                    cv->unmatched_imap_samples = (char **)realloc(
                        cv->unmatched_imap_samples,
                        sizeof(char *) * (size_t)(cv->n_unmatched_imap_samples + 1));
                    cv->unmatched_imap_samples[cv->n_unmatched_imap_samples++] =
                        xstrdup_local(s);
                    char msg[256];
                    snprintf(msg, sizeof(msg),
                        "Imap sample '%s' not found in any provided sequence file.", s);
                    add_issue(cv, "UNMATCHED_IMAP_SAMPLE", "warning", fi->path, msg);
                }
            }
        }
    }

    return cv;
}

void cross_validation_free(CrossValidation *cv)
{
    if (!cv) return;
    for (int i = 0; i < cv->n_unmatched_bam_samples; i++) free(cv->unmatched_bam_samples[i]);
    free(cv->unmatched_bam_samples);
    for (int i = 0; i < cv->n_unmatched_imap_samples; i++) free(cv->unmatched_imap_samples[i]);
    free(cv->unmatched_imap_samples);
    for (int i = 0; i < cv->n_issues; i++) {
        free(cv->issues[i].code);
        free(cv->issues[i].severity);
        free(cv->issues[i].file);
        free(cv->issues[i].message);
    }
    free(cv->issues);
    free(cv);
}
