/*
 * vcf_phase.h  –  phased VCF + BAM coverage → haplotype sequences
 *
 * This module implements --phasing vcf mode.  For each locus it runs two
 * passes in order:
 *
 *   Pass 1 (BAM pileup)
 *     Every position with depth ≥ min_dp gets the reference base written
 *     into both haplotype strings.  Positions below threshold stay 'N'.
 *     This correctly distinguishes "covered, invariant" from "not covered".
 *
 *   Pass 2 (VCF query)
 *     Every variant record in the locus is read.  The phased GT field
 *     determines which allele goes into hap1 (^1) and which into hap2 (^2).
 *     VCF values overwrite whatever pass 1 wrote, including low-DP
 *     masking based on the FORMAT/DP field.
 *
 * The result is two fully phased haplotype sequences per sample with no
 * reference-allele imputation at uncovered sites.
 */

#ifndef VCF_PHASE_H
#define VCF_PHASE_H

#include <htslib/vcf.h>
#include <htslib/hts.h>
#include <htslib/tbx.h>

#include "bam2bpp.h"

/* Policy for unphased heterozygous calls (GT uses '/' rather than '|').
 * Defined in bam2bpp.h — repeated here for documentation only.
 *   UNPHASED_MISSING  write N to both haplotypes (safe default)
 *   UNPHASED_IUPAC    write IUPAC ambiguity code to both
 */

/* Open phased VCF/BCF and load its index.
 * path must point to a bgzipped, tabix- or CSI-indexed file.
 * Returns NULL on failure. */
typedef struct {
    htsFile   *fp;
    bcf_hdr_t *hdr;
    hts_idx_t *idx;     /* CSI/BAI for BCF input; NULL for VCF.gz */
    tbx_t     *tbx;     /* tabix index for VCF.gz; NULL for BCF */
    char      *path;
    int        is_bcf;  /* 1 if BCF, 0 if VCF.gz */
} VcfPhase;

VcfPhase *open_vcf_phase(const char *path);
void      close_vcf_phase(VcfPhase *v);

/*
 * process_locus_vcf()
 *
 * Two-pass locus processor for --phasing vcf mode.
 *
 * bams[]     : one open BAM per sample (provides coverage for non-variant sites)
 * n_bams     : number of BAMs / samples
 * vcf        : phased VCF opened with open_vcf_phase()
 * loc        : locus definition (chrom, start, end, name)
 * ref_seq    : reference bases for [loc->start, loc->end), uppercase
 * args       : global args (min_dp, min_bq, min_mq, unphased_policy, etc.)
 * result     : output — caller provides the LocusResult shell; function fills it
 * mean_dp_out: mean per-position BAM depth across all samples (for stats)
 *
 * Always emits 2 × n_bams sequences (one haplotype pair per sample).
 * Returns 0 on success, -1 on error.
 */
int process_locus_vcf(BamFile       **bams,
                      int             n_bams,
                      VcfPhase       *vcf,
                      const Locus    *loc,
                      const char     *ref_seq,
                      const Args     *args,
                      LocusResult    *result,
                      double         *mean_dp_out);

/* ────────────────────────────────────────────────────────────────────────
 * Whole-file phasing classification.
 *
 * Used before the per-locus loop in --phasing vcf mode to decide whether
 * the VCF actually carries phase information for the requested samples.
 * The policy is whole-file: if any relevant sample has any unphased het
 * GT within the locus set, we degrade to IUPAC.
 *
 * Returns a malloc'd array of per-sample counts indexed against the
 * BAM order (NOT the VCF order).  Samples not found in the VCF have all
 * counts set to 0 and `not_in_vcf = 1`.
 * ───────────────────────────────────────────────────────────────────── */
typedef struct {
    char  *sample;          /* matches bams[i]->sample */
    int    not_in_vcf;      /* 1 if no column with this name in VCF */
    int    n_het;           /* total heterozygous GTs observed */
    int    n_phased;        /* of those, how many were phased ('|') */
    int    n_unphased;      /* of those, how many were unphased ('/') */
} VcfPhaseSampleStat;

/* Scan the VCF over the given BED loci and classify each sample listed
 * in `bams` by its phasing pattern.  Returns NULL on error. */
VcfPhaseSampleStat *vcf_phase_classify(VcfPhase     *vcf,
                                       BamFile     **bams,
                                       int           n_bams,
                                       const Locus  *loci,
                                       int           n_loci);

void vcf_phase_stats_free(VcfPhaseSampleStat *s, int n);

/* ────────────────────────────────────────────────────────────────────────
 * Reference-bias accounting for --phasing vcf.
 *
 * Pass 1 writes the reference base for PHASE_VCF samples and pass 2 overwrites
 * only where the VCF has a record, so a site the reads call non-reference but
 * the VCF omits silently becomes homozygous reference. That is deliberate --
 * the panel filters sequencing error -- but a filtered panel also omits rare
 * and private variants, so the count is surfaced rather than left invisible.
 *
 * reads_nonref: positions where a PHASE_VCF sample's pileup disagreed with the
 *               reference.  overridden: how many of those the VCF did not cover.
 * Reset before the locus loop; read after it.
 * ───────────────────────────────────────────────────────────────────────── */
void vcf_phase_reset_override_stats(void);
void vcf_phase_get_override_stats(int64_t *reads_nonref, int64_t *overridden);


#endif /* VCF_PHASE_H */
