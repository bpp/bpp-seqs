/* workflow.h — choose a conversion workflow from a set of inspected files,
 * and compute what (if anything) is missing. */
#ifndef BPP_SEQS_WORKFLOW_H
#define BPP_SEQS_WORKFLOW_H

#include "inspect.h"

typedef enum {
    WF_NONE = 0,
    WF_BAM2BPP,
    WF_BAM2BPP_NEEDS_BED,
    WF_BAM2BPP_NEEDS_REF,
    WF_NEEDS_ALIGNMENT_FIRST,
    WF_FASTA2BPP,
    WF_PHYLIP2BPP,
    WF_NEXUS2BPP,
    WF_GVCF2BPP,
    WF_GVCF2BPP_NEEDS_BED,
    WF_VCF_NOT_RECOMMENDED,
    WF_ASSEMBLY_NOT_SUPPORTED
} WorkflowKind;

const char *workflow_name(WorkflowKind w);

typedef struct {
    char *item;          /* "imap_file" */
    int   required;
    char *description;
    char **samples_needing_assignment;
    int    n_samples_needing_assignment;
    char *format_example;
} MissingItem;

typedef struct {
    WorkflowKind  workflow;
    MissingItem  *missing;
    int           n_missing;
    int           ready_to_run;

    /* Workflow-level advisory for diagnostic-only states (e.g.
     * vcf_not_recommended, assembly_not_supported, needs_alignment_first).
     * NULL when the workflow is actionable or unknown. */
    char         *advisory;
} WorkflowDecision;

WorkflowDecision *workflow_decide(FileInfo **files, int n);
void              workflow_decision_free(WorkflowDecision *d);

#endif /* BPP_SEQS_WORKFLOW_H */
