#include "output.h"
#include "json_writer.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ────────────────────────────────────────────────────────────────────────
 * ConversionResult lifecycle
 * ───────────────────────────────────────────────────────────────────── */

void conversion_result_init(ConversionResult *r)
{
    memset(r, 0, sizeof(*r));
}

void conversion_result_free(ConversionResult *r)
{
    free(r->out_sequences);
    free(r->out_imap);
    free(r->out_stats);
    free(r->out_loci);
    for (int i = 0; i < r->n_loci; i++) {
        free(r->loci[i].name);
        free(r->loci[i].status);
        free(r->loci[i].skip_reason);
        free(r->loci[i].source_kind);
        free(r->loci[i].source_file);
        free(r->loci[i].source_chrom);
    }
    free(r->loci);
    memset(r, 0, sizeof(*r));
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

void conversion_result_add_locus(ConversionResult *r,
                                 const char *name, int length, int n_snps,
                                 double missing_fraction, double mean_depth,
                                 const char *status, const char *skip_reason)
{
    r->loci = realloc(r->loci, sizeof(*r->loci) * (size_t)(r->n_loci + 1));
    r->loci[r->n_loci].name             = xstrdup_local(name);
    r->loci[r->n_loci].length           = length;
    r->loci[r->n_loci].n_snps           = n_snps;
    r->loci[r->n_loci].missing_fraction = missing_fraction;
    r->loci[r->n_loci].mean_depth       = mean_depth;
    r->loci[r->n_loci].status           = xstrdup_local(status);
    r->loci[r->n_loci].skip_reason      = xstrdup_local(skip_reason);
    r->loci[r->n_loci].source_kind      = NULL;
    r->loci[r->n_loci].source_file      = NULL;
    r->loci[r->n_loci].source_chrom     = NULL;
    r->loci[r->n_loci].source_start     = -1;
    r->loci[r->n_loci].source_end       = -1;
    r->loci[r->n_loci].source_stride    = 1;
    r->n_loci++;
}

void conversion_result_set_locus_source(ConversionResult *r,
                                        const char *source_kind,
                                        const char *source_file,
                                        const char *source_chrom,
                                        int source_start,
                                        int source_end,
                                        int source_stride)
{
    if (r->n_loci == 0) return;
    int i = r->n_loci - 1;
    free(r->loci[i].source_kind);
    free(r->loci[i].source_file);
    free(r->loci[i].source_chrom);
    r->loci[i].source_kind   = xstrdup_local(source_kind);
    r->loci[i].source_file   = xstrdup_local(source_file);
    r->loci[i].source_chrom  = xstrdup_local(source_chrom);
    r->loci[i].source_start  = source_start;
    r->loci[i].source_end    = source_end;
    r->loci[i].source_stride = source_stride > 0 ? source_stride : 1;
}

const char *status_string(const WorkflowDecision *d, int conversion_done, int error)
{
    if (error) return "error";
    if (conversion_done) return "complete";
    if (d && d->ready_to_run) return "ready";
    return "incomplete";
}

/* ────────────────────────────────────────────────────────────────────────
 * Human-readable
 * ───────────────────────────────────────────────────────────────────── */

static const char *basename_of(const char *path)
{
    const char *p = strrchr(path, '/');
    return p ? p + 1 : path;
}

static void format_size_bp(int64_t v, char *out, size_t cap)
{
    if (v >= 1000000000LL)      snprintf(out, cap, "%.2f Gb", (double)v / 1e9);
    else if (v >= 1000000LL)    snprintf(out, cap, "%.2f Mb", (double)v / 1e6);
    else if (v >= 1000LL)       snprintf(out, cap, "%.1f kb", (double)v / 1e3);
    else                        snprintf(out, cap, "%lld bp", (long long)v);
}

static void print_file_line_human(FILE *fp, FileInfo *fi)
{
    const char *name = basename_of(fi->path);
    switch (fi->ft) {
        case BS_BAM:
        case BS_CRAM: {
            char buf[256];
            snprintf(buf, sizeof(buf),
                "%s%s%s%s, sample=%s, %dbp %s, %.1f%% mapped",
                fi->ft == BS_CRAM ? "CRAM" : "BAM",
                fi->unaligned ? " (unaligned)" : "",
                fi->aligner ? ", " : "",
                fi->aligner ? fi->aligner : "",
                fi->sample_name ? fi->sample_name : "-",
                fi->read_length_mean,
                fi->read_type ? fi->read_type : "-",
                fi->mapping_rate * 100.0);
            fprintf(fp, "  %-24s %s\n", name, buf);
            break;
        }
        case BS_FASTQ:
            fprintf(fp, "  %-24s FASTQ, %d reads sampled, %d bp mean (range %d-%d), platform guess %s\n",
                name, fi->n_reads_sampled, fi->read_length_mean,
                fi->read_length_min, fi->read_length_max,
                fi->platform_guess ? fi->platform_guess : "?");
            break;
        case BS_FASTA_REFERENCE: {
            char sz[32]; format_size_bp(fi->total_length, sz, sizeof(sz));
            fprintf(fp, "  %-24s FASTA reference, %d sequence%s, %s%s\n",
                name, fi->n_sequences, fi->n_sequences == 1 ? "" : "s",
                sz, fi->indexed ? ", indexed" : ", NOT indexed");
            break;
        }
        case BS_FASTA_MSA:
            fprintf(fp, "  %-24s FASTA MSA, %d sequences, %d bp aligned%s\n",
                name, fi->n_sequences, fi->alignment_length,
                fi->fasta_has_gaps ? " (with gaps)" : "");
            break;
        case BS_FASTA_CONTIGS:
            fprintf(fp, "  %-24s FASTA contigs, %d sequences, mean=%lld, N50=%lld\n",
                name, fi->n_sequences,
                (long long)fi->length_mean, (long long)fi->length_n50);
            break;
        case BS_VCF:
            fprintf(fp, "  %-24s VCF, %d samples, %d records sampled%s\n",
                name, fi->n_samples, fi->n_records_sampled,
                fi->has_phase_info ? ", phased" : "");
            break;
        case BS_GVCF:
            fprintf(fp, "  %-24s gVCF, %d samples, %d records sampled\n",
                name, fi->n_samples, fi->n_records_sampled);
            break;
        case BS_BED:
            fprintf(fp, "  %-24s BED, %d loci, mean length %lld bp\n",
                name, fi->n_loci, (long long)fi->bed_length_mean);
            break;
        case BS_IMAP:
            fprintf(fp, "  %-24s Imap, %d samples → %d populations\n",
                name, fi->imap_n_samples, fi->imap_n_population_names);
            break;
        case BS_PHYLIP:
            fprintf(fp, "  %-24s PHYLIP %s, %d sequences × %d sites\n",
                name, fi->phylip_format ? fi->phylip_format : "?",
                fi->phylip_n_sequences, fi->phylip_n_sites);
            break;
        case BS_NEXUS:
            fprintf(fp, "  %-24s NEXUS, %d sequences × %d sites%s, %d locus block%s\n",
                name, fi->nexus_n_sequences, fi->nexus_n_sites,
                fi->nexus_has_charsets ? " (with charsets)" : "",
                fi->nexus_n_loci, fi->nexus_n_loci == 1 ? "" : "s");
            break;
        default:
            fprintf(fp, "  %-24s UNKNOWN file type\n", name);
            break;
    }
}

void print_human(FILE *fp,
                 FileInfo **files, int n_files,
                 const CrossValidation *cv,
                 const WorkflowDecision *d,
                 const ConversionResult *cr,
                 const char *bpp_seqs_version)
{
    (void)bpp_seqs_version;
    fprintf(fp, "bpp-seqs: inspecting %d file%s\n\n", n_files, n_files == 1 ? "" : "s");

    for (int i = 0; i < n_files; i++) {
        print_file_line_human(fp, files[i]);
    }

    fprintf(fp, "\nWorkflow: %s\n", workflow_name(d->workflow));

    /* Cross-validation */
    int cv_ok = (cv->bam_reference_consistent && cv->bam_reference_matches_fasta &&
                 cv->bed_chromosomes_in_bams && cv->imap_samples_in_bams &&
                 cv->bam_samples_in_imap);
    fprintf(fp, "\nCross-validation: %s\n", cv_ok ? "OK" : "issues found");
    for (int i = 0; i < cv->n_issues; i++) {
        const ValidationIssue *iss = &cv->issues[i];
        fprintf(fp, "  %s [%s]: %s\n",
            iss->severity ? iss->severity : "warning",
            iss->code ? iss->code : "?",
            iss->message ? iss->message : "");
    }

    /* Aggregate warnings from files */
    int any_w = 0;
    for (int i = 0; i < n_files; i++) if (files[i]->n_warnings > 0) { any_w = 1; break; }
    if (any_w) {
        fprintf(fp, "\nWarnings:\n");
        for (int i = 0; i < n_files; i++) {
            for (int j = 0; j < files[i]->n_warnings; j++) {
                Warning *w = &files[i]->warnings[j];
                fprintf(fp, "  %s %s [%s]: %s\n",
                    w->severity ? w->severity : "warning",
                    basename_of(files[i]->path),
                    w->code ? w->code : "?",
                    w->message ? w->message : "");
            }
        }
    }

    /* Missing items */
    if (d->n_missing > 0) {
        fprintf(fp, "\nMissing:\n");
        for (int i = 0; i < d->n_missing; i++) {
            MissingItem *mi = &d->missing[i];
            fprintf(fp, "  %s — %s\n", mi->item, mi->description ? mi->description : "");
            if (mi->n_samples_needing_assignment > 0) {
                fprintf(fp, "  Samples needing assignment: ");
                for (int j = 0; j < mi->n_samples_needing_assignment; j++) {
                    fprintf(fp, "%s%s", mi->samples_needing_assignment[j],
                        j + 1 < mi->n_samples_needing_assignment ? ", " : "");
                }
                fputc('\n', fp);
            }
            if (mi->format_example) {
                fprintf(fp, "  Format: two tab-separated columns, one sample per line, e.g.\n");
                /* indent each line of the example */
                const char *p = mi->format_example;
                while (*p) {
                    const char *nl = strchr(p, '\n');
                    if (nl) { fprintf(fp, "    %.*s\n", (int)(nl - p), p); p = nl + 1; }
                    else    { fprintf(fp, "    %s\n", p); break; }
                }
            }
        }
        fprintf(fp, "\nRe-run with the Imap file to proceed.\n");
    }

    /* Conversion results */
    if (cr && cr->has_results) {
        fprintf(fp, "\nResults:\n");
        fprintf(fp, "  %d / %d loci passed QC\n", cr->n_loci_passed, cr->n_loci_input);
        if (cr->n_loci_failed > 0) {
            fprintf(fp, "  %d loci excluded: %d too short, %d high missing data, %d insufficient SNPs\n",
                cr->n_loci_failed, cr->n_too_short, cr->n_high_missing, cr->n_insufficient_snps);
        }
        fprintf(fp, "\nOutput files:\n");
        if (cr->out_sequences) fprintf(fp, "  %s  BPP sequence file (%d loci, %d sequences)\n",
            cr->out_sequences, cr->n_loci_passed, cr->n_sequences);
        if (cr->out_imap)      fprintf(fp, "  %s  BPP Imap file\n", cr->out_imap);
        if (cr->out_stats)     fprintf(fp, "  %s  Per-locus statistics\n", cr->out_stats);
    }
}

/* ────────────────────────────────────────────────────────────────────────
 * JSON
 * ───────────────────────────────────────────────────────────────────── */

static void jw_str_array(JsonWriter *w, char **arr, int n)
{
    jw_arr_open(w);
    for (int i = 0; i < n; i++) jw_str(w, arr[i]);
    jw_arr_close(w);
}

static void file_info_to_json(JsonWriter *w, FileInfo *fi)
{
    jw_obj_open(w);
    jw_kv_str(w, "path", fi->path);
    jw_kv_str(w, "type", file_type_name(fi->ft));

    switch (fi->ft) {
        case BS_BAM:
        case BS_CRAM:
            jw_kv_str (w, "subtype", fi->unaligned ? "unaligned" : "aligned");
            jw_kv_str (w, "aligner", fi->aligner);
            jw_kv_str (w, "aligner_version", fi->aligner_version);
            jw_kv_str (w, "aligner_command", fi->aligner_command);
            jw_kv_str (w, "sample_name", fi->sample_name);
            jw_kv_str (w, "platform", fi->platform);
            jw_kv_int (w, "read_length_mean", fi->read_length_mean);
            jw_kv_str (w, "read_type", fi->read_type);
            jw_kv_dbl (w, "mapping_rate", fi->mapping_rate);
            jw_kv_dbl (w, "mean_depth_estimate", fi->mean_depth_estimate);
            jw_kv_str (w, "depth_estimate_basis", "reads_x_length_div_reflen");
            jw_kv_int (w, "n_reads_sampled", fi->n_reads_sampled);
            jw_key(w, "reference_sequences"); jw_arr_open(w);
            for (int i = 0; i < fi->n_seq_refs; i++) {
                jw_obj_open(w);
                jw_kv_str(w, "name", fi->seq_refs[i].name);
                jw_kv_int(w, "length", fi->seq_refs[i].length);
                jw_obj_close(w);
            }
            jw_arr_close(w);
            jw_kv_str (w, "reference_path_in_header", fi->reference_path_in_header);
            jw_kv_bool(w, "indexed", fi->indexed);
            break;
        case BS_FASTQ:
            jw_kv_int (w, "read_length_mean", fi->read_length_mean);
            jw_kv_int (w, "read_length_min", fi->read_length_min);
            jw_kv_int (w, "read_length_max", fi->read_length_max);
            jw_kv_str (w, "read_type", fi->read_type);
            jw_kv_int (w, "n_reads_sampled", fi->n_reads_sampled);
            jw_kv_str (w, "platform_guess", fi->platform_guess);
            jw_kv_str (w, "recommended_aligner", fi->recommended_aligner);
            break;
        case BS_FASTA_REFERENCE:
            jw_kv_int (w, "n_sequences", fi->n_sequences);
            jw_kv_int (w, "total_length", fi->total_length);
            jw_key(w, "sequence_names"); {
                int show = fi->n_sequence_names < 10 ? fi->n_sequence_names : 10;
                jw_arr_open(w);
                for (int i = 0; i < show; i++) jw_str(w, fi->sequence_names[i]);
                jw_arr_close(w);
            }
            jw_kv_int (w, "n_sequence_names_total", fi->n_sequence_names);
            jw_kv_bool(w, "indexed", fi->indexed);
            break;
        case BS_FASTA_MSA:
            jw_kv_int (w, "n_sequences", fi->n_sequences);
            jw_kv_int (w, "alignment_length", fi->alignment_length);
            jw_key(w, "sequence_names"); jw_str_array(w, fi->sequence_names, fi->n_sequence_names);
            jw_kv_int (w, "n_loci", 1);
            jw_kv_bool(w, "has_gaps", fi->fasta_has_gaps);
            jw_kv_dbl (w, "missing_fraction", fi->fasta_missing_fraction);
            break;
        case BS_FASTA_CONTIGS:
            jw_kv_int (w, "n_sequences", fi->n_sequences);
            jw_kv_int (w, "length_min", fi->length_min);
            jw_kv_int (w, "length_max", fi->length_max);
            jw_kv_int (w, "length_mean", fi->length_mean);
            jw_kv_int (w, "length_n50", fi->length_n50);
            jw_kv_bool(w, "has_n_runs", fi->has_n_runs);
            break;
        case BS_VCF:
        case BS_GVCF:
            jw_kv_int (w, "n_samples", fi->n_samples);
            jw_key(w, "sample_names"); jw_str_array(w, fi->sample_names, fi->n_sample_names);
            jw_kv_str (w, "reference_in_header", fi->vcf_reference_in_header);
            jw_kv_int (w, "n_records_sampled", fi->n_records_sampled);
            jw_kv_bool(w, "has_phase_info", fi->has_phase_info);
            jw_kv_dbl (w, "phased_fraction", fi->phased_fraction);
            if (fi->ft == BS_GVCF) jw_kv_bool(w, "has_coverage_bands", fi->has_coverage_bands);
            break;
        case BS_BED:
            jw_kv_int (w, "n_loci", fi->n_loci);
            jw_key(w, "chromosomes"); jw_str_array(w, fi->chromosomes, fi->n_chromosomes);
            jw_kv_int (w, "length_min", fi->bed_length_min);
            jw_kv_int (w, "length_max", fi->bed_length_max);
            jw_kv_int (w, "length_mean", fi->bed_length_mean);
            jw_kv_int (w, "length_median", fi->bed_length_median);
            jw_kv_bool(w, "has_names", fi->bed_has_names);
            break;
        case BS_IMAP:
            jw_kv_int (w, "n_samples", fi->imap_n_samples);
            jw_kv_int (w, "n_populations", fi->imap_n_population_names);
            jw_key(w, "sample_names"); jw_str_array(w, fi->imap_sample_names, fi->imap_n_sample_names);
            jw_key(w, "population_names"); jw_str_array(w, fi->imap_population_names, fi->imap_n_population_names);
            jw_key(w, "samples_per_population"); jw_obj_open(w);
            for (int p = 0; p < fi->imap_n_population_names; p++) {
                int c = 0;
                for (int i = 0; i < fi->imap_n_entries; i++) {
                    if (strcmp(fi->imap_pops[i], fi->imap_population_names[p]) == 0) c++;
                }
                jw_kv_int(w, fi->imap_population_names[p], c);
            }
            jw_obj_close(w);
            break;
        case BS_PHYLIP:
            jw_kv_str (w, "format", fi->phylip_format);
            jw_kv_int (w, "n_sequences", fi->phylip_n_sequences);
            jw_kv_int (w, "n_sites", fi->phylip_n_sites);
            jw_kv_dbl (w, "missing_fraction", fi->phylip_missing_fraction);
            break;
        case BS_NEXUS:
            jw_kv_int (w, "n_sequences", fi->nexus_n_sequences);
            jw_kv_int (w, "n_sites", fi->nexus_n_sites);
            jw_kv_int (w, "n_loci", fi->nexus_n_loci);
            jw_kv_bool(w, "has_charsets", fi->nexus_has_charsets);
            jw_key(w, "charset_names"); jw_str_array(w, fi->nexus_charset_names, fi->nexus_n_charset_names);
            jw_kv_dbl (w, "missing_fraction", fi->nexus_missing_fraction);
            break;
        default:
            break;
    }

    /* warnings */
    jw_key(w, "warnings"); jw_arr_open(w);
    for (int i = 0; i < fi->n_warnings; i++) {
        jw_obj_open(w);
        jw_kv_str(w, "code", fi->warnings[i].code);
        jw_kv_str(w, "severity", fi->warnings[i].severity);
        jw_kv_str(w, "message", fi->warnings[i].message);
        jw_obj_close(w);
    }
    jw_arr_close(w);

    jw_obj_close(w);
}

void print_json(FILE *fp, int indent,
                FileInfo **files, int n_files,
                const CrossValidation *cv,
                const WorkflowDecision *d,
                const ConversionResult *cr,
                const char *recommended_command,
                const char *bpp_seqs_version)
{
    JsonWriter w; jw_init(&w, fp, indent);
    int complete = cr && cr->has_results;
    const char *status = status_string(d, complete, 0);

    jw_obj_open(&w);
    jw_kv_str (&w, "bpp_seqs_version", bpp_seqs_version);
    jw_kv_str (&w, "status", status);
    jw_kv_str (&w, "workflow", workflow_name(d->workflow));
    jw_kv_bool(&w, "ready_to_run", d->ready_to_run);

    jw_key(&w, "files_provided"); jw_arr_open(&w);
    for (int i = 0; i < n_files; i++) file_info_to_json(&w, files[i]);
    jw_arr_close(&w);

    /* cross_validation */
    jw_key(&w, "cross_validation"); jw_obj_open(&w);
    jw_kv_bool(&w, "bam_reference_consistent",    cv->bam_reference_consistent);
    jw_kv_bool(&w, "bam_reference_matches_fasta", cv->bam_reference_matches_fasta);
    jw_kv_bool(&w, "bed_chromosomes_in_bams",     cv->bed_chromosomes_in_bams);
    jw_kv_bool(&w, "imap_samples_in_bams",        cv->imap_samples_in_bams);
    jw_kv_bool(&w, "bam_samples_in_imap",         cv->bam_samples_in_imap);
    jw_key(&w, "unmatched_bam_samples");  jw_str_array(&w, cv->unmatched_bam_samples,  cv->n_unmatched_bam_samples);
    jw_key(&w, "unmatched_imap_samples"); jw_str_array(&w, cv->unmatched_imap_samples, cv->n_unmatched_imap_samples);
    jw_key(&w, "issues"); jw_arr_open(&w);
    for (int i = 0; i < cv->n_issues; i++) {
        jw_obj_open(&w);
        jw_kv_str(&w, "code",     cv->issues[i].code);
        jw_kv_str(&w, "severity", cv->issues[i].severity);
        jw_kv_str(&w, "file",     cv->issues[i].file);
        jw_kv_str(&w, "message",  cv->issues[i].message);
        jw_obj_close(&w);
    }
    jw_arr_close(&w);
    jw_obj_close(&w);

    /* missing */
    jw_key(&w, "missing"); jw_arr_open(&w);
    for (int i = 0; i < d->n_missing; i++) {
        MissingItem *mi = &d->missing[i];
        jw_obj_open(&w);
        jw_kv_str (&w, "item", mi->item);
        jw_kv_bool(&w, "required", mi->required);
        jw_kv_str (&w, "description", mi->description);
        jw_key(&w, "samples_needing_assignment");
        jw_str_array(&w, mi->samples_needing_assignment, mi->n_samples_needing_assignment);
        jw_kv_str (&w, "format_example", mi->format_example);
        jw_obj_close(&w);
    }
    jw_arr_close(&w);

    /* top-level warnings (aggregated from files + cross-validation issues
     * at warning+ severity).  Keeping it simple: just per-file warnings. */
    jw_key(&w, "warnings"); jw_arr_open(&w);
    for (int i = 0; i < n_files; i++) {
        for (int j = 0; j < files[i]->n_warnings; j++) {
            jw_obj_open(&w);
            jw_kv_str(&w, "code",     files[i]->warnings[j].code);
            jw_kv_str(&w, "severity", files[i]->warnings[j].severity);
            jw_kv_str(&w, "file",     files[i]->path);
            jw_kv_str(&w, "message",  files[i]->warnings[j].message);
            jw_obj_close(&w);
        }
    }
    jw_arr_close(&w);

    /* recommended_command (only meaningful when ready_to_run but not yet
     * run, i.e. status="ready") */
    if (recommended_command) jw_kv_str(&w, "recommended_command", recommended_command);

    /* output_files + summary + loci, populated only when conversion ran */
    if (cr && cr->has_results) {
        jw_key(&w, "output_files"); jw_obj_open(&w);
        jw_kv_str(&w, "sequences", cr->out_sequences);
        jw_kv_str(&w, "imap",      cr->out_imap);
        jw_kv_str(&w, "stats",     cr->out_stats);
        jw_kv_str(&w, "loci",      cr->out_loci);
        jw_obj_close(&w);

        jw_key(&w, "summary"); jw_obj_open(&w);
        jw_kv_int(&w, "n_loci_input",  cr->n_loci_input);
        jw_kv_int(&w, "n_loci_passed", cr->n_loci_passed);
        jw_kv_int(&w, "n_loci_failed", cr->n_loci_failed);
        jw_kv_int(&w, "n_sequences",   cr->n_sequences);
        jw_key(&w, "failure_reasons"); jw_obj_open(&w);
        jw_kv_int(&w, "too_short",          cr->n_too_short);
        jw_kv_int(&w, "high_missing",       cr->n_high_missing);
        jw_kv_int(&w, "insufficient_snps",  cr->n_insufficient_snps);
        jw_obj_close(&w);
        jw_obj_close(&w);

        jw_key(&w, "loci"); jw_arr_open(&w);
        for (int i = 0; i < cr->n_loci; i++) {
            jw_obj_open(&w);
            jw_kv_str(&w, "name",             cr->loci[i].name);
            jw_kv_int(&w, "length",           cr->loci[i].length);
            jw_kv_int(&w, "n_snps",           cr->loci[i].n_snps);
            jw_kv_dbl(&w, "missing_fraction", cr->loci[i].missing_fraction);
            jw_kv_dbl(&w, "mean_depth",       cr->loci[i].mean_depth);
            jw_kv_str(&w, "status",           cr->loci[i].status);
            jw_kv_str(&w, "skip_reason",      cr->loci[i].skip_reason);
            jw_kv_str(&w, "source_kind",      cr->loci[i].source_kind);
            jw_kv_str(&w, "source_file",      cr->loci[i].source_file);
            jw_kv_str(&w, "source_chrom",     cr->loci[i].source_chrom);
            if (cr->loci[i].source_start > 0) jw_kv_int(&w, "source_start", cr->loci[i].source_start);
            else jw_kv_null(&w, "source_start");
            if (cr->loci[i].source_end > 0)   jw_kv_int(&w, "source_end",   cr->loci[i].source_end);
            else jw_kv_null(&w, "source_end");
            jw_kv_int(&w, "source_stride",    cr->loci[i].source_stride);
            jw_obj_close(&w);
        }
        jw_arr_close(&w);
    } else {
        jw_kv_null(&w, "output_files");
    }

    jw_obj_close(&w);
    jw_finish(&w);
}
