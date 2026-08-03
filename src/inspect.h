/* inspect.h — content-based file type detection and per-type inspection.
 *
 * One entry point: inspect_file(path) → FileInfo*.  Returns NULL on
 * "could not open / read" errors; populates a FileInfo with ft = FT_UNKNOWN
 * (plus a warning) for unrecognised content.
 *
 * All allocations are owned by the returned FileInfo; free with file_info_free().
 */
#ifndef BPP_SEQS_INSPECT_H
#define BPP_SEQS_INSPECT_H

#include <stdint.h>
#include <stddef.h>

typedef enum {
    BS_UNKNOWN = 0,
    BS_BAM,
    BS_CRAM,
    BS_FASTQ,
    BS_FASTA_REFERENCE,
    BS_FASTA_MSA,
    BS_FASTA_CONTIGS,
    BS_VCF,
    BS_GVCF,
    BS_BED,
    BS_IMAP,
    BS_PHYLIP,
    BS_NEXUS,
    BS_CONTROL        /* a BPP control file (not sequence data) */
} FileType;

const char *file_type_name(FileType t);

typedef struct {
    char *code;        /* short uppercase code, e.g. LOW_MAPPING_RATE  */
    char *severity;    /* info | warning | error | critical            */
    char *message;     /* human-readable */
} Warning;

typedef struct {
    char *name;
    int64_t length;
} SeqRef;

typedef struct FileInfo {
    char     *path;
    FileType  ft;

    /* Common */
    Warning  *warnings;
    int       n_warnings;

    /* BAM / CRAM */
    char     *aligner;
    char     *aligner_version;
    char     *aligner_command;
    char     *sample_name;
    char     *platform;
    int       read_length_mean;
    char     *read_type;          /* "paired" | "single" */
    double    mapping_rate;
    double    mean_depth_estimate;
    int       n_reads_sampled;
    SeqRef   *seq_refs;           /* @SQ list */
    int       n_seq_refs;
    char     *reference_path_in_header;
    int       indexed;
    int       unaligned;          /* > 50% reads unmapped */

    /* FASTQ */
    int       read_length_min;
    int       read_length_max;
    char     *platform_guess;
    char     *recommended_aligner;

    /* FASTA_* */
    int       n_sequences;
    int64_t   total_length;
    int       alignment_length;   /* MSA: shared length */
    char    **sequence_names;
    int       n_sequence_names;
    int       fasta_has_gaps;
    double    fasta_missing_fraction;
    /* contigs only */
    int64_t   length_min;
    int64_t   length_max;
    int64_t   length_mean;
    int64_t   length_n50;
    int       has_n_runs;

    /* VCF / GVCF */
    int       n_samples;
    char    **sample_names;
    int       n_sample_names;
    char     *vcf_reference_in_header;
    int       n_records_sampled;
    int       has_phase_info;
    double    phased_fraction;
    int       has_coverage_bands; /* GVCF only */
    SeqRef   *vcf_contigs;        /* ##contig=&lt;ID=…,length=…&gt; */
    int       n_vcf_contigs;

    /* BED */
    int       n_loci;
    char    **chromosomes;
    int       n_chromosomes;
    int64_t   bed_length_min;
    int64_t   bed_length_max;
    int64_t   bed_length_mean;
    int64_t   bed_length_median;
    int       bed_has_names;

    /* IMAP */
    int       imap_n_samples;
    char    **imap_sample_names;
    int       imap_n_sample_names;
    char    **imap_population_names;
    int       imap_n_population_names;
    /* parallel arrays: imap_samples[i] → imap_pops[i] */
    char    **imap_samples;        /* duplicate of sample_names for direct use */
    char    **imap_pops;
    int       imap_n_entries;

    /* PHYLIP */
    char     *phylip_format;       /* "interleaved" | "sequential" */
    int       phylip_n_sequences;
    int       phylip_n_sites;
    int       phylip_n_loci;       /* count of locus headers; >1 => BPP-native multi-locus */
    double    phylip_missing_fraction;

    /* NEXUS */
    int       nexus_n_sequences;
    int       nexus_n_sites;       /* total sites */
    int       nexus_n_loci;        /* number of charsets, or 1 if none */
    int       nexus_has_charsets;
    char    **nexus_charset_names;
    int       nexus_n_charset_names;
    /* charset definitions: [start, end] inclusive, 1-based; stride; total
     * sites covered by this charset; -1 in start to mean a string-range we
     * could not parse and skipped. */
    int      *nexus_charset_starts;
    int      *nexus_charset_ends;
    int      *nexus_charset_strides;
    double    nexus_missing_fraction;
} FileInfo;

FileInfo *inspect_file(const char *path);
void      file_info_free(FileInfo *fi);

/* Append a warning. Severity is "info"|"warning"|"error"|"critical". */
void file_info_add_warning(FileInfo *fi,
                           const char *code,
                           const char *severity,
                           const char *message);

/* Override a FASTA's classification to BS_FASTA_REFERENCE (e.g. in response
 * to a CLI --reference flag). Removes any ASSEMBLY_NOT_ALIGNABLE warning the
 * auto-classifier may have attached, and sets the .fai-indexed flag with the
 * usual MISSING_FAI warning if the index is absent. No-op if `fi` is already
 * BS_FASTA_REFERENCE; returns 0 if the override applied, -1 if `fi` is not a
 * recognised FASTA type. */
int file_info_force_reference(FileInfo *fi);

#endif /* BPP_SEQS_INSPECT_H */
