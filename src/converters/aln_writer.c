#include "aln_writer.h"
#include "sanity.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char *xdup(const char *s)
{
    if (!s) return NULL;
    size_t n = strlen(s);
    char *r = (char *)malloc(n + 1);
    if (!r) return NULL;
    memcpy(r, s, n + 1);
    return r;
}

int write_loci_tsv(const char *out_prefix, const LocusProv *items, int n)
{
    char path[1024];
    snprintf(path, sizeof(path), "%s.loci.tsv", out_prefix);
    FILE *f = fopen(path, "w");
    if (!f) return -1;
    fprintf(f, "locus_name\tsource_kind\tsource_file\tsource_chrom\t"
               "source_start\tsource_end\tsource_stride\tlength\tn_seqs\n");
    for (int i = 0; i < n; i++) {
        const LocusProv *it = &items[i];
        fprintf(f, "%s\t%s\t%s\t%s\t",
                it->name  ? it->name  : ".",
                it->kind  ? it->kind  : ".",
                it->file  ? it->file  : ".",
                it->chrom ? it->chrom : ".");
        if (it->start > 0) fprintf(f, "%d\t", it->start); else fprintf(f, ".\t");
        if (it->end   > 0) fprintf(f, "%d\t", it->end);   else fprintf(f, ".\t");
        fprintf(f, "%d\t%d\t%d\n",
                it->stride > 0 ? it->stride : 1,
                it->length,
                it->n_seqs);
    }
    fclose(f);
    return 0;
}

int write_alignment_outputs(const char *out_prefix,
                            LocusAln *loci, int n_loci,
                            FileInfo *imap_fi,
                            int min_length, double max_missing,
                            int min_snps, int keep_invariant,
                            ConversionResult *cr)
{
    /* Sanity-check every locus's structure and sanitize bad characters. */
    SanityReport sr; sanity_report_init(&sr);
    int *sane = (int *)calloc((size_t)n_loci, sizeof(int));
    for (int i = 0; i < n_loci; i++) {
        sane[i] = (sanity_check_locus(loci[i].name,
                                      loci[i].seq_names, loci[i].seqs,
                                      loci[i].n_seqs, loci[i].length,
                                      &sr, stderr) == 0);
    }

    /* Decide which loci pass QC. */
    int *passes = (int *)calloc((size_t)n_loci, sizeof(int));
    int n_pass = 0;
    int n_too_short = 0, n_high_missing = 0, n_insufficient = 0;
    int n_seqs = (n_loci > 0) ? loci[0].n_seqs : 0;

    for (int i = 0; i < n_loci; i++) {
        const char *skip = NULL;
        if (!sane[i])                                 { skip = "sanity_failed"; }
        if (!skip && loci[i].length < min_length)     { skip = "too_short";          n_too_short++; }
        double miss = compute_missing_frac(loci[i].seqs, loci[i].n_seqs, loci[i].length);
        int    nsnp = count_snps         (loci[i].seqs, loci[i].n_seqs, loci[i].length);
        if (!skip && miss > max_missing)              { skip = "high_missing";       n_high_missing++; }
        if (!skip && !keep_invariant && nsnp < min_snps) { skip = "insufficient_snps"; n_insufficient++; }
        passes[i] = (skip == NULL);
        if (passes[i]) n_pass++;
        conversion_result_add_locus(cr, loci[i].name, loci[i].length, nsnp, miss, 0.0,
                                    passes[i] ? "passed" : "failed", skip);
    }
    free(sane);

    if (n_pass == 0) {
        free(passes);
        fprintf(stderr, "Error: no loci passed QC.\n");
        for (int i = 0; i < n_loci; i++) {
            free(loci[i].name);
            for (int j = 0; j < loci[i].n_seqs; j++) {
                free(loci[i].seq_names[j]);
                free(loci[i].seqs[j]);
            }
            free(loci[i].seq_names);
            free(loci[i].seqs);
        }
        return 1;
    }

    /* Write sequences file in BPP per-locus format:
     *     <n_seqs> <n_sites>
     *
     *     ^<id>    SEQUENCE
     */
    char path[1024];
    snprintf(path, sizeof(path), "%s.txt", out_prefix);
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); free(passes); return 1; }
    int first = 1;
    for (int i = 0; i < n_loci; i++) {
        if (!passes[i]) continue;
        if (!first) fputc('\n', f);
        first = 0;
        fprintf(f, "%d %d\n\n", loci[i].n_seqs, loci[i].length);
        /* compute column width */
        int max_id = 0;
        for (int j = 0; j < loci[i].n_seqs; j++) {
            const char *nm = loci[i].seq_names[j];
            /* Caret already present (pre-formatted BPP name like "label^id")
             * → no extra caret prepended. Otherwise we add one. */
            int already = (strchr(nm, '^') != NULL);
            int l = (int)strlen(nm) + (already ? 0 : 1);
            if (l > max_id) max_id = l;
        }
        int col_width = max_id + 4;
        for (int j = 0; j < loci[i].n_seqs; j++) {
            const char *nm = loci[i].seq_names[j];
            int already = (strchr(nm, '^') != NULL);
            int written = already ? fprintf(f, "%s", nm)
                                  : fprintf(f, "^%s", nm);
            int pad = col_width - written;
            if (pad < 1) pad = 1;
            for (int p = 0; p < pad; p++) fputc(' ', f);
            fprintf(f, "%s\n", loci[i].seqs[j]);
        }
    }
    fclose(f);
    cr->out_sequences = xdup(path);

    /* Write Imap from imap_fi */
    snprintf(path, sizeof(path), "%s.imap", out_prefix);
    f = fopen(path, "w");
    if (!f) { perror(path); free(passes); return 1; }
    if (imap_fi) {
        for (int i = 0; i < imap_fi->imap_n_entries; i++) {
            fprintf(f, "%s\t%s\n", imap_fi->imap_samples[i], imap_fi->imap_pops[i]);
        }
    }
    fclose(f);
    cr->out_imap = xdup(path);

    /* Stats */
    snprintf(path, sizeof(path), "%s.stats.tsv", out_prefix);
    f = fopen(path, "w");
    if (f) {
        fprintf(f, "locus\tlength\tn_snps\tmissing_frac\tmean_depth\tstatus\tskip_reason\n");
        for (int i = 0; i < cr->n_loci; i++) {
            fprintf(f, "%s\t%d\t%d\t%.4f\t%.2f\t%s\t%s\n",
                cr->loci[i].name, cr->loci[i].length, cr->loci[i].n_snps,
                cr->loci[i].missing_fraction, cr->loci[i].mean_depth,
                cr->loci[i].status,
                cr->loci[i].skip_reason ? cr->loci[i].skip_reason : "");
        }
        fclose(f);
        cr->out_stats = xdup(path);
    }

    cr->has_results        = 1;
    cr->n_loci_input       = n_loci;
    cr->n_loci_passed      = n_pass;
    cr->n_loci_failed      = n_loci - n_pass;
    cr->n_sequences        = n_seqs;
    cr->n_too_short        = n_too_short;
    cr->n_high_missing     = n_high_missing;
    cr->n_insufficient_snps = n_insufficient;

    /* Imap coverage check: collect all unique post-caret ids across all
     * passing loci and reconcile against the Imap. Warnings only. */
    if (imap_fi) {
        char **ids = NULL; int n_ids = 0;
        for (int i = 0; i < n_loci; i++) {
            if (!passes[i]) continue;
            for (int j = 0; j < loci[i].n_seqs; j++) {
                char *id = sanity_individual_id(loci[i].seq_names[j]);
                if (!id) continue;
                int dup = 0;
                for (int k = 0; k < n_ids; k++) if (strcmp(ids[k], id) == 0) { dup = 1; break; }
                if (dup) { free(id); continue; }
                ids = (char **)realloc(ids, sizeof(char *) * (size_t)(n_ids + 1));
                ids[n_ids++] = id;
            }
        }
        sanity_check_imap_coverage(ids, n_ids,
                                   imap_fi->imap_samples, imap_fi->imap_n_entries,
                                   stderr);
        for (int i = 0; i < n_ids; i++) free(ids[i]);
        free(ids);
    }

    /* Write per-locus provenance for every passing locus. */
    {
        LocusProv *items = (LocusProv *)calloc((size_t)n_pass, sizeof(LocusProv));
        int k = 0;
        for (int i = 0; i < n_loci; i++) {
            if (!passes[i]) continue;
            items[k].name   = loci[i].name;
            items[k].kind   = loci[i].source_kind   ? loci[i].source_kind   : ".";
            items[k].file   = loci[i].source_file   ? loci[i].source_file   : ".";
            items[k].chrom  = loci[i].source_chrom  ? loci[i].source_chrom  : ".";
            items[k].start  = loci[i].source_start;
            items[k].end    = loci[i].source_end;
            items[k].stride = loci[i].source_stride > 0 ? loci[i].source_stride : 1;
            items[k].length = loci[i].length;
            items[k].n_seqs = loci[i].n_seqs;
            k++;
        }
        if (write_loci_tsv(out_prefix, items, n_pass) == 0) {
            snprintf(path, sizeof(path), "%s.loci.tsv", out_prefix);
            cr->out_loci = xdup(path);
        }
        free(items);
    }

    /* Mirror provenance into the ConversionResult.loci entries (in the
     * same order they were added above), so the JSON output carries the
     * same fields the .loci.tsv has. */
    for (int i = 0; i < cr->n_loci && i < n_loci; i++) {
        free(cr->loci[i].source_kind);
        free(cr->loci[i].source_file);
        free(cr->loci[i].source_chrom);
        cr->loci[i].source_kind   = loci[i].source_kind  ? strdup(loci[i].source_kind)  : NULL;
        cr->loci[i].source_file   = loci[i].source_file  ? strdup(loci[i].source_file)  : NULL;
        cr->loci[i].source_chrom  = loci[i].source_chrom ? strdup(loci[i].source_chrom) : NULL;
        cr->loci[i].source_start  = loci[i].source_start;
        cr->loci[i].source_end    = loci[i].source_end;
        cr->loci[i].source_stride = loci[i].source_stride > 0 ? loci[i].source_stride : 1;
    }

    /* Free input loci */
    for (int i = 0; i < n_loci; i++) {
        free(loci[i].name);
        free(loci[i].source_kind);
        free(loci[i].source_file);
        free(loci[i].source_chrom);
        for (int j = 0; j < loci[i].n_seqs; j++) {
            free(loci[i].seq_names[j]);
            free(loci[i].seqs[j]);
        }
        free(loci[i].seq_names);
        free(loci[i].seqs);
    }
    free(passes);
    return 0;
}
