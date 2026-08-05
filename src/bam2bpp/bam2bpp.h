/*
 * bam2bpp.h  –  BAM to BPP sequence format converter
 *
 * Converts indexed BAM files + reference FASTA directly into the sequence
 * alignment file required by BPP for species-tree / population-genetics
 * analysis.  No VCF or external tools are involved: htslib pileup is used
 * internally to determine per-sample base calls and coverage at every
 * position in each user-defined locus region.
 *
 * Coverage is always known from the pileup.  Positions with depth below
 * --min-dp are emitted as 'N' regardless of the reference base, so the
 * reference sequence is never used as a silent substitute for missing data.
 */

#ifndef BAM2BPP_H
#define BAM2BPP_H

#include <stdint.h>
#include <stdio.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <htslib/hts.h>

/* -------------------------------------------------------------------------
 * Constants
 * ---------------------------------------------------------------------- */

#define BAM2BPP_VERSION   "0.1.0"
#define DEFAULT_MIN_BQ    20
#define DEFAULT_MIN_MQ    20
#define DEFAULT_MIN_DP     5
#define DEFAULT_HET_FREQ   0.20   /* minor-allele freq threshold for het call */
#define DEFAULT_MIN_LEN  100
#define DEFAULT_MAX_MISS   0.5
#define DEFAULT_MIN_SNPS   1
#define DEFAULT_MAX_DEPTH  8000

/* Base indices used in BaseCounts.counts[] */
#define BASE_A 0
#define BASE_C 1
#define BASE_G 2
#define BASE_T 3

/* -------------------------------------------------------------------------
 * Enumerations
 * ---------------------------------------------------------------------- */

typedef enum {
    PHASE_IUPAC   = 0,   /* het → IUPAC ambiguity code; 1 seq per sample  */
    PHASE_SPLIT   = 1,   /* 2 seqs per sample, arbitrary phase             */
    PHASE_HAPLOID = 2,   /* major allele only; 1 seq per sample            */
    PHASE_VCF     = 3    /* 2 seqs per sample using phased VCF (best)      */
} Phasing;

/* Which base caller turns a pileup column into a genotype.
 *
 * CALLER_CONSENSUS is the samtools/Gap5 Bayesian model, vendored under
 * src/vendor/samtools and verified byte-identical to `samtools consensus -A`
 * by the oracle tests. It weighs each observation by its base quality, caps
 * that by mapping quality, and masks positions it cannot call confidently.
 *
 * CALLER_COUNTS is bpp-seqs' original caller: rank raw allele counts, call a
 * heterozygote when the minor allele clears --het-freq. It has no error model,
 * so at low depth a single miscalled base can carry a site. Kept so results
 * from earlier versions can still be reproduced.
 *
 * Only PHASE_IUPAC can use CALLER_CONSENSUS: a consensus caller emits one
 * sequence per individual, with heterozygotes as IUPAC codes. PHASE_SPLIT and
 * PHASE_VCF need two resolved haplotypes, which it does not produce, so they
 * always use CALLER_COUNTS regardless of this setting. */
typedef enum {
    CALLER_CONSENSUS = 0,
    CALLER_COUNTS    = 1
} Caller;

/* What to do with unphased het calls (GT uses '/' not '|') in VCF mode */
typedef enum {
    UNPHASED_MISSING = 0,   /* write N to both haplotypes (safe default) */
    UNPHASED_IUPAC   = 1,   /* write IUPAC code to both (more data, less safe) */
} UnphasedPolicy;

/* -------------------------------------------------------------------------
 * Core data types
 * ---------------------------------------------------------------------- */

/* One record from the loci BED file */
typedef struct {
    char    *chrom;
    int32_t  start;   /* 0-based, inclusive  */
    int32_t  end;     /* 0-based, exclusive  */
    char    *name;
} Locus;

/* One record from the Imap file (sample → population) */
typedef struct {
    char *sample_name;
    char *population;
} ImapEntry;

/* Filtered base counts at one pileup column for one sample */
typedef struct {
    int32_t counts[4];   /* A=0, C=1, G=2, T=3 */
    int32_t depth;       /* total filtered depth (sum of counts) */
} BaseCounts;

/* Sequences produced for one locus (all samples) */
typedef struct {
    char  *name;         /* locus name (strdup of BED name field) */
    char **seqs;         /* n_seqs strings, each locus_len+1 chars */
    char **seq_names;    /* sequence label for each entry          */
    int    n_seqs;
    int    locus_len;
} LocusResult;

/* Statistics recorded for every locus (used or skipped) */
typedef struct {
    char   *locus_name;
    int     locus_len;
    int     n_snps;
    double  missing_frac;
    double  mean_depth;
    char   *skip_reason;  /* NULL = locus was written to output */
} LocusStat;

/* Open BAM file with its header and index */
typedef struct {
    samFile   *fp;
    sam_hdr_t *hdr;
    hts_idx_t *idx;
    char      *sample;   /* SM tag from first @RG; falls back to filename stem */
    char      *path;
    Phasing    sample_phasing;  /* per-sample emission policy. Mirrors args->phasing
                                 * by default; in --phasing vcf mode the VCF
                                 * classification overrides per sample so that
                                 * phased samples emit 2 haplotypes and unphased
                                 * samples emit 1 IUPAC sequence. */
} BamFile;

/* Parsed command-line arguments */
typedef struct {
    char  **bam_paths;
    int     n_bams;
    char   *ref_path;
    char   *bed_path;
    char   *imap_path;
    char   *out_prefix;
    int     min_bq;
    int     min_mq;
    int     min_dp;
    double  het_freq;       /* minor-allele freq ≥ this → call het */
    Phasing phasing;
    Caller  caller;              /* base caller; PHASE_IUPAC only (see above) */
    char   *phased_vcf;          /* path to phased VCF (--phasing vcf only) */
    UnphasedPolicy unphased_policy; /* what to do with unphased het GTs     */
    int     min_length;
    double  max_missing;
    int     keep_invariant;
    int     min_snps;
    int     max_depth;
} Args;

/* -------------------------------------------------------------------------
 * pileup.c  –  BAM access and per-locus pileup
 * ---------------------------------------------------------------------- */

/* Open all BAMs; returns NULL on any failure */
BamFile **open_bams(char **paths, int n, int max_depth);

/* Close and free all BAMs */
void close_bams(BamFile **bams, int n);

/*
 * process_locus()
 *
 * Run htslib multi-sample pileup over [loc->start, loc->end) in all BAMs.
 * For each position × sample:
 *   depth == 0 or depth < min_dp  →  'N'  (missing / low coverage)
 *   depth ≥ min_dp                →  called base (haploid/iupac/split)
 *
 * The reference sequence is NOT used to fill uncovered sites; every
 * character in the output sequences derives from actual reads.
 *
 * Fills *result (caller provides the LocusResult shell).
 * Writes mean per-position depth across all samples to *mean_dp_out.
 * Returns 0 on success, -1 on error.
 */
int process_locus(BamFile **bams, int n_bams,
                  const Locus *loc, const Args *args,
                  LocusResult *result, double *mean_dp_out);

/* Free all heap memory owned by a LocusResult */
void free_locus_result(LocusResult *r);

/* -------------------------------------------------------------------------
 * genotype.c  –  base calling from allele counts
 * ---------------------------------------------------------------------- */

/* Return the single most-frequent base, or 'N' if depth < min_dp */
char call_haploid(const BaseCounts *bc, int min_dp);

/*
 * Return IUPAC code for the position:
 *   hom        →  A / C / G / T
 *   het        →  M R W S Y K  (minor allele freq ≥ het_freq)
 *   low depth  →  N
 */
char call_iupac(const BaseCounts *bc, int min_dp, double het_freq);

/*
 * Write the two alleles of a diploid call to *b1 and *b2.
 * hom: b1 == b2 == major allele
 * het: b1 = major, b2 = minor  (arbitrary phase)
 * low depth: b1 = b2 = 'N'
 */
void call_split(const BaseCounts *bc, int min_dp, double het_freq,
                char *b1, char *b2);

/* -------------------------------------------------------------------------
 * locus.c  –  I/O helpers and sequence statistics
 * ---------------------------------------------------------------------- */

Locus     *load_bed  (const char *path, int *n_out);
void       free_loci (Locus *loci, int n);

ImapEntry *load_imap (const char *path, int *n_out);
void       free_imap (ImapEntry *imap, int n);

/* Return the population string for a sample, or NULL if not found */
const char *imap_lookup(const ImapEntry *imap, int n_imap, const char *sample);

/* Fraction of 'N' characters across all sequences in a locus result */
double compute_missing_frac(char **seqs, int n_seqs, int locus_len);

/* Number of alignment columns where at least 2 distinct non-N bases appear */
int count_snps(char **seqs, int n_seqs, int locus_len);

/* -------------------------------------------------------------------------
 * writer.c  –  output files
 * ---------------------------------------------------------------------- */

/* When non-zero, the three writers suppress their "Wrote ..." stderr
 * progress lines.  Set by bpp-seqs main.c when --quiet is on. */
void bam2bpp_writer_set_quiet(int q);

/*
 * Write the BPP sequence file (<prefix>.txt).
 * Format:
 *     n_seqs  n_loci          ← right-justified 6-wide integers
 *
 *     ^locus_name
 *     sample^1   ACGT...
 *     sample^2   ACGT...
 *     ...
 */
void write_bpp_file(const char *prefix,
                    LocusResult *results, int n_results, int n_seqs);

/*
 * Write the BPP Imap file (<prefix>.imap).
 * For PHASE_SPLIT each sample gets two entries (sample^1, sample^2).
 */
void write_imap_file(const char *prefix,
                     BamFile **bams, int n_bams,
                     const ImapEntry *imap, int n_imap,
                     Phasing phasing);

/* Write per-locus statistics to <prefix>.stats.tsv */
void write_stats_file(const char *prefix, const LocusStat *stats, int n);

/* -------------------------------------------------------------------------
 * main.c
 * ---------------------------------------------------------------------- */

void usage(const char *prog);

#endif /* BAM2BPP_H */
