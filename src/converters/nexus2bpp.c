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
    if (!imap_fi) {
        fprintf(stderr, "Error: nexus2bpp requires a NEXUS file and an Imap.\n");
        return 1;
    }

    /* Accumulate loci across ALL NEXUS inputs. Two multilocus conventions are
     * both honoured: one file with CHARSETs (each charset -> a locus), and the
     * common one-NEXUS-file-per-locus layout (each file -> a locus). */
    LocusAln *loci = NULL;
    int n_loci_out = 0;
    int n_nexus = 0;

    for (int i = 0; i < n_files; i++) {
        if (files[i]->ft != BS_NEXUS) continue;
        n_nexus++;
        FileInfo *nx = files[i];

        NexusDoc d;
        if (nexus_parse(nx->path, &d) != 0) {
            fprintf(stderr, "Error: NEXUS parse failed for '%s': %s\n",
                    nx->path, d.err[0] ? d.err : "unknown error");
            nexus_free(&d);
            free(loci);
            return 1;
        }
        int ntax = d.ntax, nchar = d.nchar;
        int add = d.n_charsets ? d.n_charsets : 1;
        loci = (LocusAln *)realloc(loci, (size_t)(n_loci_out + add) * sizeof(LocusAln));

        if (d.n_charsets == 0) {
            int c = n_loci_out++;
            char nm[32];
            snprintf(nm, sizeof nm, "locus%d", c + 1);
            loci[c].name   = strdup(nm);
            loci[c].length = nchar;
            loci[c].n_seqs = ntax;
            loci[c].seq_names = (char **)calloc((size_t)ntax, sizeof(char *));
            loci[c].seqs      = (char **)calloc((size_t)ntax, sizeof(char *));
            for (int t = 0; t < ntax; t++) {
                loci[c].seq_names[t] = strdup(d.taxa[t]);
                loci[c].seqs[t]      = strdup(d.seq[t]);
            }
            loci[c].source_kind   = strdup("FILE");
            loci[c].source_file   = strdup(nx->path);
            loci[c].source_chrom  = NULL;
            loci[c].source_start  = 1;
            loci[c].source_end    = nchar;
            loci[c].source_stride = 1;
        } else {
            for (int cs = 0; cs < d.n_charsets; cs++) {
                int c  = n_loci_out++;
                int s  = d.cs_start[cs];             /* 1-based inclusive */
                int e  = d.cs_end[cs];
                int st = d.cs_stride[cs] ? d.cs_stride[cs] : 1;
                if (s < 1) s = 1;
                if (e < 1 || e > nchar) e = nchar;   /* end == -1 ("to last") clamps here */
                int len = 0;
                for (int x = s; x <= e; x += st) len++;
                loci[c].name   = strdup(d.cs_name[cs]);
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
    }

    if (n_nexus == 0) {
        fprintf(stderr, "Error: nexus2bpp requires a NEXUS file and an Imap.\n");
        free(loci);
        return 1;
    }

    return write_alignment_outputs(out_prefix, loci, n_loci_out, imap_fi,
                                   opts->min_length, opts->max_missing,
                                   opts->min_snps, opts->keep_invariant, cr);
}
