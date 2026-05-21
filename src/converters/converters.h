/* converters.h — entry points for non-bam2bpp workflows. */
#ifndef BPP_SEQS_CONVERTERS_H
#define BPP_SEQS_CONVERTERS_H

#include "inspect.h"
#include "output.h"

typedef struct {
    int    min_length;
    double max_missing;
    int    min_snps;
    int    keep_invariant;
    int    quiet;
} ConvertOpts;

/* Each converter writes <prefix>.txt (sequences), <prefix>.imap, and
 * <prefix>.stats.tsv. Returns 0 on success, non-zero on error. Fills cr
 * with stats. */
int convert_fasta_msa(FileInfo **files, int n_files,
                      FileInfo *imap_fi,
                      const char *out_prefix,
                      const ConvertOpts *opts,
                      ConversionResult *cr);

int convert_phylip   (FileInfo **files, int n_files,
                      FileInfo *imap_fi,
                      const char *out_prefix,
                      const ConvertOpts *opts,
                      ConversionResult *cr);

int convert_nexus    (FileInfo **files, int n_files,
                      FileInfo *imap_fi,
                      const char *out_prefix,
                      const ConvertOpts *opts,
                      ConversionResult *cr);

int convert_gvcf     (FileInfo **files, int n_files,
                      FileInfo *imap_fi,
                      const char *out_prefix,
                      const ConvertOpts *opts,
                      ConversionResult *cr);

#endif /* BPP_SEQS_CONVERTERS_H */
