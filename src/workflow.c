#include "workflow.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const char *workflow_name(WorkflowKind w)
{
    switch (w) {
        case WF_BAM2BPP:                return "bam2bpp";
        case WF_BAM2BPP_NEEDS_BED:      return "bam2bpp_needs_bed";
        case WF_BAM2BPP_NEEDS_REF:      return "bam2bpp_needs_ref";
        case WF_NEEDS_ALIGNMENT_FIRST:  return "needs_alignment_first";
        case WF_FASTA2BPP:              return "fasta2bpp";
        case WF_PHYLIP2BPP:             return "phylip2bpp";
        case WF_NEXUS2BPP:              return "nexus2bpp";
        case WF_GVCF2BPP:               return "gvcf2bpp";
        case WF_GVCF2BPP_NEEDS_BED:     return "gvcf2bpp_needs_bed";
        case WF_VCF_NOT_RECOMMENDED:    return "vcf_not_recommended";
        case WF_ASSEMBLY_NOT_SUPPORTED: return "assembly_not_supported";
        case WF_NONE:
        default:                        return "unknown";
    }
}

static char *xstrdup_local(const char *s)
{
    if (!s) return NULL;
    size_t n = strlen(s);
    char *r = (char *)malloc(n + 1);
    if (!r) return NULL;
    memcpy(r, s, n + 1);
    return r;
}

static int count_type(FileInfo **files, int n, FileType t)
{
    int c = 0;
    for (int i = 0; i < n; i++) if (files[i]->ft == t) c++;
    return c;
}

static int has_imap(FileInfo **files, int n)
{
    return count_type(files, n, BS_IMAP) > 0;
}

/* Collect sample names that need an Imap assignment. From BAMs we use
 * @RG SM tag; from MSA/PHYLIP/NEXUS/VCF we use sequence/sample names. */
static void gather_samples_needing_assignment(FileInfo **files, int n,
                                              char ***out, int *n_out)
{
    char **list = NULL;
    int   ln = 0;

    for (int i = 0; i < n; i++) {
        FileInfo *fi = files[i];
        if ((fi->ft == BS_BAM || fi->ft == BS_CRAM) && fi->sample_name) {
            /* dedupe */
            int found = 0;
            for (int k = 0; k < ln; k++) if (strcmp(list[k], fi->sample_name) == 0) { found = 1; break; }
            if (!found) {
                list = (char **)realloc(list, sizeof(char *) * (size_t)(ln + 1));
                list[ln++] = xstrdup_local(fi->sample_name);
            }
        } else if (fi->ft == BS_VCF || fi->ft == BS_GVCF) {
            for (int j = 0; j < fi->n_sample_names; j++) {
                int found = 0;
                for (int k = 0; k < ln; k++) if (strcmp(list[k], fi->sample_names[j]) == 0) { found = 1; break; }
                if (!found) {
                    list = (char **)realloc(list, sizeof(char *) * (size_t)(ln + 1));
                    list[ln++] = xstrdup_local(fi->sample_names[j]);
                }
            }
        } else if (fi->ft == BS_FASTA_MSA || fi->ft == BS_PHYLIP || fi->ft == BS_NEXUS) {
            for (int j = 0; j < fi->n_sequence_names; j++) {
                int found = 0;
                for (int k = 0; k < ln; k++) if (strcmp(list[k], fi->sequence_names[j]) == 0) { found = 1; break; }
                if (!found) {
                    list = (char **)realloc(list, sizeof(char *) * (size_t)(ln + 1));
                    list[ln++] = xstrdup_local(fi->sequence_names[j]);
                }
            }
        }
    }
    *out = list;
    *n_out = ln;
}

static void add_imap_missing(WorkflowDecision *d, FileInfo **files, int n)
{
    char **samps = NULL; int ns = 0;
    gather_samples_needing_assignment(files, n, &samps, &ns);

    d->missing = (MissingItem *)realloc(d->missing,
        sizeof(MissingItem) * (size_t)(d->n_missing + 1));
    MissingItem *mi = &d->missing[d->n_missing++];
    mi->item = xstrdup_local("imap_file");
    mi->required = 1;
    mi->description = xstrdup_local(
        "Tab-separated file mapping each sample name to a population or species. "
        "Two columns, no header. Sample names must match BAM @RG SM tags.");
    mi->samples_needing_assignment = samps;
    mi->n_samples_needing_assignment = ns;
    mi->format_example = xstrdup_local("ind1\tpopA\nind2\tpopA\nind3\tpopB\nind4\tpopB");
}

static void add_missing(WorkflowDecision *d, const char *item, const char *desc)
{
    d->missing = (MissingItem *)realloc(d->missing,
        sizeof(MissingItem) * (size_t)(d->n_missing + 1));
    MissingItem *mi = &d->missing[d->n_missing++];
    mi->item = xstrdup_local(item);
    mi->required = 1;
    mi->description = xstrdup_local(desc);
    mi->samples_needing_assignment = NULL;
    mi->n_samples_needing_assignment = 0;
    mi->format_example = NULL;
}

WorkflowDecision *workflow_decide(FileInfo **files, int n)
{
    WorkflowDecision *d = (WorkflowDecision *)calloc(1, sizeof(WorkflowDecision));
    if (!d) return NULL;
    d->workflow = WF_NONE;

    int n_bam     = count_type(files, n, BS_BAM) + count_type(files, n, BS_CRAM);
    int n_fastq   = count_type(files, n, BS_FASTQ);
    int n_ref     = count_type(files, n, BS_FASTA_REFERENCE);
    int n_bed     = count_type(files, n, BS_BED);
    int n_msa     = count_type(files, n, BS_FASTA_MSA);
    int n_contigs = count_type(files, n, BS_FASTA_CONTIGS);
    int n_phylip  = count_type(files, n, BS_PHYLIP);
    int n_nexus   = count_type(files, n, BS_NEXUS);
    int n_vcf     = count_type(files, n, BS_VCF);
    int n_gvcf    = count_type(files, n, BS_GVCF);

    /* Priority order: most specific actionable workflow first. */
    if (n_bam > 0 && n_ref > 0 && n_bed > 0) {
        d->workflow = WF_BAM2BPP;
    } else if (n_bam > 0 && n_bed == 0) {
        d->workflow = WF_BAM2BPP_NEEDS_BED;
    } else if (n_bam > 0 && n_ref == 0 && n_bed > 0) {
        d->workflow = WF_BAM2BPP_NEEDS_REF;
    } else if (n_fastq > 0) {
        /* Any FASTQ present means raw reads: they must be aligned to a reference
         * before BPP conversion, whether or not a reference/BED is also given.
         * (The old condition required ref+BED, so a bare FASTQ -- the commonest
         * "I have raw reads" case -- fell through to WF_NONE/"unknown".) */
        d->workflow = WF_NEEDS_ALIGNMENT_FIRST;
    } else if (n_msa > 0) {
        d->workflow = WF_FASTA2BPP;
    } else if (n_phylip > 0) {
        d->workflow = WF_PHYLIP2BPP;
    } else if (n_nexus > 0) {
        d->workflow = WF_NEXUS2BPP;
    } else if (n_gvcf > 0 && n_bed > 0) {
        d->workflow = WF_GVCF2BPP;
    } else if (n_gvcf > 0) {
        /* gVCF present but no BED: convertible once loci are defined, exactly
         * like a BAM without a BED (WF_BAM2BPP_NEEDS_BED). Report the gap. */
        d->workflow = WF_GVCF2BPP_NEEDS_BED;
    } else if (n_vcf > 0 && n_bed > 0) {
        d->workflow = WF_VCF_NOT_RECOMMENDED;
    } else if (n_contigs > 0) {
        d->workflow = WF_ASSEMBLY_NOT_SUPPORTED;
    } else {
        d->workflow = WF_NONE;
    }

    /* The only thing that goes in missing[] is an absent Imap, and only
     * for workflows that need it. */
    int needs_imap = 0;
    switch (d->workflow) {
        case WF_BAM2BPP:
        case WF_FASTA2BPP:
        case WF_PHYLIP2BPP:
        case WF_NEXUS2BPP:
        case WF_GVCF2BPP:
        case WF_GVCF2BPP_NEEDS_BED:
            needs_imap = 1;
            break;
        default:
            needs_imap = 0;
    }

    if (needs_imap && !has_imap(files, n)) {
        add_imap_missing(d, files, n);
    }

    /* Partial-BAM workflows: list the still-missing piece. */
    if (d->workflow == WF_BAM2BPP_NEEDS_BED) {
        add_missing(d, "bed_file",
            "BED file defining locus intervals: chrom<TAB>start<TAB>end[<TAB>name]. "
            "Chrom names must match the BAM @SQ headers; start is 0-based, end is exclusive.");
    } else if (d->workflow == WF_BAM2BPP_NEEDS_REF) {
        add_missing(d, "reference_fasta",
            "Reference FASTA the BAMs were aligned to, with a .fai index "
            "(produced by `samtools faidx ref.fa`).");
    } else if (d->workflow == WF_GVCF2BPP_NEEDS_BED) {
        add_missing(d, "bed_file",
            "BED file defining locus intervals: chrom<TAB>start<TAB>end[<TAB>name]. "
            "Chrom names must match the gVCF contigs; start is 0-based, end is "
            "exclusive. Each interval becomes one BPP locus.");
    }

    /* Workflow-level advisories for diagnostic-only states. */
    switch (d->workflow) {
        case WF_VCF_NOT_RECOMMENDED:
            d->advisory = xstrdup_local(
                "Standard (non-gVCF) VCF cannot distinguish 'invariant and covered' "
                "from 'not covered' at non-variant sites. BPP requires that "
                "distinction. Use the BAM files this VCF was called from, or a gVCF "
                "with coverage bands, instead.");
            break;
        case WF_ASSEMBLY_NOT_SUPPORTED:
            d->advisory = xstrdup_local(
                "Assembled contigs cannot be used directly as BPP input. "
                "Loci must first be defined and sequences aligned across individuals. "
                "Either align reads to a reference and provide BAM+BED, or build an "
                "MSA per locus and provide it as FASTA.");
            break;
        case WF_NEEDS_ALIGNMENT_FIRST:
            d->advisory = xstrdup_local(
                "FASTQ input requires alignment before BPP conversion. "
                "Run an aligner (bwa mem for short reads, minimap2 for long reads) "
                "against the reference FASTA, sort and index the BAMs, and re-run "
                "with BAM+reference+BED.");
            break;
        default: break;
    }

    /* Ready to run iff: have a real conversion workflow AND no missing items. */
    int convertable = (d->workflow == WF_BAM2BPP ||
                       d->workflow == WF_FASTA2BPP ||
                       d->workflow == WF_PHYLIP2BPP ||
                       d->workflow == WF_NEXUS2BPP ||
                       d->workflow == WF_GVCF2BPP);
    d->ready_to_run = (convertable && d->n_missing == 0);

    return d;
}

void workflow_decision_free(WorkflowDecision *d)
{
    if (!d) return;
    for (int i = 0; i < d->n_missing; i++) {
        free(d->missing[i].item);
        free(d->missing[i].description);
        free(d->missing[i].format_example);
        for (int j = 0; j < d->missing[i].n_samples_needing_assignment; j++) {
            free(d->missing[i].samples_needing_assignment[j]);
        }
        free(d->missing[i].samples_needing_assignment);
    }
    free(d->missing);
    free(d->advisory);
    free(d);
}
