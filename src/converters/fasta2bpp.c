#define _POSIX_C_SOURCE 200809L

#include "aln_writer.h"
#include "converters.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

/* Load an MSA FASTA into parallel arrays. All seqs must share a length.
 * Returns 0 on success, non-zero on error. */
static int load_msa(const char *path, char ***names_out, char ***seqs_out,
                    int *n_out, int *len_out)
{
    gzFile gz = gzopen(path, "rb");
    if (!gz) return -1;
    char line[1 << 16];
    char **names = NULL, **seqs = NULL;
    int n = 0;
    int cap = 0;
    char *cur = NULL;
    size_t cur_len = 0, cur_cap = 0;
    int expected_len = -1;

    while (gzgets(gz, line, sizeof(line)) != NULL) {
        if (line[0] == '>') {
            if (cur) {
                if (expected_len < 0) expected_len = (int)cur_len;
                else if ((int)cur_len != expected_len) {
                    free(cur); for (int i = 0; i < n; i++) { free(names[i]); free(seqs[i]); }
                    free(names); free(seqs); gzclose(gz);
                    return -2;
                }
                cur[cur_len] = '\0';
                seqs[n - 1] = cur;
                cur = NULL; cur_len = 0; cur_cap = 0;
            }
            if (n >= cap) { cap = cap ? cap * 2 : 16;
                names = realloc(names, cap * sizeof(char*));
                seqs  = realloc(seqs,  cap * sizeof(char*)); }
            char *p = line + 1;
            char *end = p;
            while (*end && *end != ' ' && *end != '\t' && *end != '\n' && *end != '\r') end++;
            *end = '\0';
            names[n] = strdup(p);
            n++;
        } else {
            for (char *q = line; *q && *q != '\n' && *q != '\r'; q++) {
                if (cur_len + 1 >= cur_cap) {
                    cur_cap = cur_cap ? cur_cap * 2 : 256;
                    cur = realloc(cur, cur_cap);
                }
                cur[cur_len++] = (char)toupper((unsigned char)*q);
            }
        }
    }
    if (cur) {
        if (expected_len < 0) expected_len = (int)cur_len;
        else if ((int)cur_len != expected_len) {
            free(cur); for (int i = 0; i < n; i++) { free(names[i]); free(seqs[i]); }
            free(names); free(seqs); gzclose(gz);
            return -2;
        }
        cur[cur_len] = '\0';
        if (n - 1 >= 0) seqs[n - 1] = cur;
    }
    gzclose(gz);
    *names_out = names;
    *seqs_out  = seqs;
    *n_out     = n;
    *len_out   = expected_len > 0 ? expected_len : 0;
    return 0;
}

int convert_fasta_msa(FileInfo **files, int n_files,
                      FileInfo *imap_fi,
                      const char *out_prefix,
                      const ConvertOpts *opts,
                      ConversionResult *cr)
{
    /* Collect MSA file paths */
    FileInfo **msas = NULL;
    int n_msa = 0;
    for (int i = 0; i < n_files; i++) {
        if (files[i]->ft == BS_FASTA_MSA) {
            msas = realloc(msas, sizeof(*msas) * (size_t)(n_msa + 1));
            msas[n_msa++] = files[i];
        }
    }
    if (n_msa == 0 || !imap_fi) {
        free(msas);
        fprintf(stderr, "Error: fasta2bpp requires at least one MSA FASTA and an Imap.\n");
        return 1;
    }

    LocusAln *loci = (LocusAln *)calloc((size_t)n_msa, sizeof(LocusAln));
    int n_seqs_expected = -1;

    for (int i = 0; i < n_msa; i++) {
        char **names = NULL, **seqs = NULL;
        int n = 0, len = 0;
        if (load_msa(msas[i]->path, &names, &seqs, &n, &len) != 0) {
            fprintf(stderr, "Error: failed to load MSA '%s'\n", msas[i]->path);
            free(msas); free(loci); return 1;
        }
        if (n_seqs_expected < 0) n_seqs_expected = n;
        else if (n != n_seqs_expected) {
            fprintf(stderr, "Error: MSA '%s' has %d sequences but %d expected\n",
                msas[i]->path, n, n_seqs_expected);
            for (int j = 0; j < n; j++) { free(names[j]); free(seqs[j]); }
            free(names); free(seqs); free(msas); free(loci); return 1;
        }

        const char *bn = strrchr(msas[i]->path, '/');
        bn = bn ? bn + 1 : msas[i]->path;
        char nm[256];
        snprintf(nm, sizeof(nm), "%s", bn);
        char *dot = strrchr(nm, '.'); if (dot) *dot = '\0';

        loci[i].name      = strdup(nm);
        loci[i].length    = len;
        loci[i].n_seqs    = n;
        loci[i].seq_names = names;
        loci[i].seqs      = seqs;
        loci[i].source_kind   = strdup("MSA");
        loci[i].source_file   = strdup(msas[i]->path);
        loci[i].source_chrom  = NULL;
        loci[i].source_start  = -1;
        loci[i].source_end    = -1;
        loci[i].source_stride = 1;
    }

    return write_alignment_outputs(out_prefix, loci, n_msa, imap_fi,
                                   opts->min_length, opts->max_missing,
                                   opts->min_snps, opts->keep_invariant, cr);
}
