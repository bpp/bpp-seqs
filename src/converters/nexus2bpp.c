#define _POSIX_C_SOURCE 200809L

#include "aln_writer.h"
#include "converters.h"
#include "nexus.h"       /* shared spec-based NEXUS reader (via -Isrc) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int convert_nexus(FileInfo **files, int n_files,
                  FileInfo *imap_fi,
                  const char *out_prefix,
                  const ConvertOpts *opts,
                  ConversionResult *cr)
{
    /* Single NEXUS file expected (the spec doesn't require multi-file). */
    FileInfo *nx = NULL;
    for (int i = 0; i < n_files; i++) if (files[i]->ft == BS_NEXUS) { nx = files[i]; break; }
    if (!nx || !imap_fi) {
        fprintf(stderr, "Error: nexus2bpp requires a NEXUS file and an Imap.\n");
        return 1;
    }

    NexusDoc d;
    if (nexus_parse(nx->path, &d) != 0) {
        fprintf(stderr, "Error: NEXUS parse failed for '%s': %s\n",
                nx->path, d.err[0] ? d.err : "unknown error");
        nexus_free(&d);
        return 1;
    }
    int ntax = d.ntax, nchar = d.nchar;

    /* Decide loci: one per charset, or one over all sites if none. Each locus
     * gets its own copies of names/seqs so the shared NexusDoc can be freed. */
    int n_loci_out = d.n_charsets ? d.n_charsets : 1;
    LocusAln *loci = (LocusAln *)calloc((size_t)n_loci_out, sizeof(LocusAln));

    if (d.n_charsets == 0) {
        loci[0].name   = strdup("locus1");
        loci[0].length = nchar;
        loci[0].n_seqs = ntax;
        loci[0].seq_names = (char **)calloc((size_t)ntax, sizeof(char *));
        loci[0].seqs      = (char **)calloc((size_t)ntax, sizeof(char *));
        for (int t = 0; t < ntax; t++) {
            loci[0].seq_names[t] = strdup(d.taxa[t]);
            loci[0].seqs[t]      = strdup(d.seq[t]);
        }
        loci[0].source_kind   = strdup("CHARSET");
        loci[0].source_file   = strdup(nx->path);
        loci[0].source_chrom  = NULL;
        loci[0].source_start  = 1;
        loci[0].source_end    = nchar;
        loci[0].source_stride = 1;
    } else {
        for (int c = 0; c < d.n_charsets; c++) {
            int s  = d.cs_start[c];                 /* 1-based inclusive */
            int e  = d.cs_end[c];
            int st = d.cs_stride[c] ? d.cs_stride[c] : 1;
            if (s < 1) s = 1;
            if (e < 1 || e > nchar) e = nchar;      /* end == -1 ("to last") clamps here */
            int len = 0;
            for (int x = s; x <= e; x += st) len++;
            loci[c].name   = strdup(d.cs_name[c]);
            loci[c].length = len;
            loci[c].n_seqs = ntax;
            loci[c].seq_names = (char **)calloc((size_t)ntax, sizeof(char *));
            loci[c].seqs      = (char **)calloc((size_t)ntax, sizeof(char *));
            loci[c].source_kind   = strdup("CHARSET");
            loci[c].source_file   = strdup(nx->path);
            loci[c].source_chrom  = NULL;
            loci[c].source_start  = s;
            loci[c].source_end    = e;
            loci[c].source_stride = st;
            for (int t = 0; t < ntax; t++) {
                loci[c].seq_names[t] = strdup(d.taxa[t]);
                loci[c].seqs[t]      = (char *)malloc((size_t)len + 1);
                int k = 0;
                for (int x = s; x <= e; x += st) loci[c].seqs[t][k++] = d.seq[t][x - 1];
                loci[c].seqs[t][len] = '\0';
            }
        }
    }

    nexus_free(&d);

    return write_alignment_outputs(out_prefix, loci, n_loci_out, imap_fi,
                                   opts->min_length, opts->max_missing,
                                   opts->min_snps, opts->keep_invariant, cr);
}
