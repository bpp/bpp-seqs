/*
 * vcf_phase.c  –  phased VCF + BAM coverage → haplotype sequences
 */

#define _POSIX_C_SOURCE 200809L

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <inttypes.h>

#include <htslib/vcf.h>
#include <htslib/sam.h>
#include <htslib/hts.h>

#include "bam2bpp.h"
#include "vcf_phase.h"

/* -------------------------------------------------------------------------
 * IUPAC table (re-used from genotype.c logic, kept local here)
 * ---------------------------------------------------------------------- */

static const char iupac_tbl[4][4] = {
    {'A', 'M', 'R', 'W'},
    {'M', 'C', 'S', 'Y'},
    {'R', 'S', 'G', 'K'},
    {'W', 'Y', 'K', 'T'},
};

static int base_to_idx(char b)
{
    switch (b) {
    case 'A': case 'a': return 0;
    case 'C': case 'c': return 1;
    case 'G': case 'g': return 2;
    case 'T': case 't': return 3;
    default: return -1;
    }
}

static char iupac_pair(char a, char b)
{
    int ia = base_to_idx(a);
    int ib = base_to_idx(b);
    if (ia < 0 || ib < 0) return 'N';
    int lo = (ia < ib) ? ia : ib;
    int hi = (ia < ib) ? ib : ia;
    return iupac_tbl[lo][hi];
}

/* -------------------------------------------------------------------------
 * Open / close VCF
 * ---------------------------------------------------------------------- */

VcfPhase *open_vcf_phase(const char *path)
{
    VcfPhase *v = calloc(1, sizeof(VcfPhase));
    if (!v) return NULL;

    v->path = strdup(path);
    v->fp   = hts_open(path, "r");
    if (!v->fp) {
        fprintf(stderr, "Error: cannot open VCF '%s': ", path);
        perror(NULL);
        free(v->path); free(v);
        return NULL;
    }

    v->hdr = bcf_hdr_read(v->fp);
    if (!v->hdr) {
        fprintf(stderr, "Error: cannot read VCF header from '%s'\n", path);
        hts_close(v->fp); free(v->path); free(v);
        return NULL;
    }

    /* Dispatch on file format: BCF wants a CSI/BAI index queried via
     * bcf_itr_*; VCF.gz wants a tabix index queried via tbx_itr_* with
     * a vcf_parse step.  Detect from the htsFile format field. */
    v->is_bcf = (v->fp->format.format == bcf);

    if (v->is_bcf) {
        v->idx = bcf_index_load3(path, NULL, 0);
        if (!v->idx) {
            fprintf(stderr,
                    "Error: no CSI index for BCF '%s'\n"
                    "       Run: bcftools index %s\n", path, path);
            bcf_hdr_destroy(v->hdr);
            hts_close(v->fp); free(v->path); free(v);
            return NULL;
        }
    } else {
        v->tbx = tbx_index_load3(path, NULL, 0);
        if (!v->tbx) {
            fprintf(stderr,
                    "Error: no tabix index for VCF.gz '%s'\n"
                    "       Run: tabix -p vcf %s\n", path, path);
            bcf_hdr_destroy(v->hdr);
            hts_close(v->fp); free(v->path); free(v);
            return NULL;
        }
    }

    fprintf(stderr, "Opened phased %s: %s  (%d samples)\n",
            v->is_bcf ? "BCF" : "VCF.gz",
            path, bcf_hdr_nsamples(v->hdr));
    return v;
}

void close_vcf_phase(VcfPhase *v)
{
    if (!v) return;
    if (v->idx) hts_idx_destroy(v->idx);
    if (v->tbx) tbx_destroy(v->tbx);
    if (v->hdr) bcf_hdr_destroy(v->hdr);
    if (v->fp)  hts_close(v->fp);
    free(v->path);
    free(v);
}

/* -------------------------------------------------------------------------
 * Per-BAM pileup callback  (mirrors the one in pileup.c)
 * ---------------------------------------------------------------------- */

typedef struct {
    samFile   *fp;
    sam_hdr_t *hdr;
    hts_itr_t *itr;
    int        min_mq;
} CovAux;

static int cov_read(void *data, bam1_t *b)
{
    CovAux *aux = (CovAux *)data;
    int ret;
    while ((ret = sam_itr_next(aux->fp, aux->itr, b)) >= 0) {
        if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                            BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY))
            continue;
        if ((int)b->core.qual < aux->min_mq)
            continue;
        break;
    }
    return ret;
}

/* 4-bit htslib base → ACGT index; -1 for ambiguous */
static const int b4_idx[16] = {
    -1, 0, 1, -1, 2, -1, -1, -1,
     3, -1, -1, -1, -1, -1, -1, -1,
};

/* -------------------------------------------------------------------------
 * Pass 1 (unified): BAM pileup for coverage AND IUPAC consensus.
 *
 *   SAMPLE_PHASED samples:
 *     positions with depth >= min_dp get the reference base in BOTH
 *     haplotypes; variant positions will be overwritten in pass 2.
 *   SAMPLE_IUPAC samples:
 *     call_iupac() runs at every covered position and the result is written
 *     to the single output sequence. Pass 2 will not touch these.
 *
 * seq_offset[i] is the starting index in seqs[] for sample i; the number of
 * sequences for sample i is 2 (PHASED) or 1 (IUPAC).
 *
 * Returns total filtered depth across all samples and positions (for stats).
 * ---------------------------------------------------------------------- */

static int64_t pass1_unified(BamFile         **bams,
                              int               n_bams,
                              const SampleMode *modes,
                              const int        *seq_offset,
                              const Locus      *loc,
                              const char       *ref_seq,
                              const Args       *args,
                              char            **seqs)
{
    char region[1024];
    snprintf(region, sizeof(region),
             "%s:%" PRId32 "-%" PRId32,
             loc->chrom, loc->start + 1, loc->end);

    CovAux  *aux  = calloc(n_bams, sizeof(CovAux));
    void   **data = malloc(n_bams * sizeof(void *));
    if (!aux || !data) { free(aux); free(data); return -1; }

    int locus_tid = sam_hdr_name2tid(bams[0]->hdr, loc->chrom);

    for (int i = 0; i < n_bams; i++) {
        aux[i].fp     = bams[i]->fp;
        aux[i].hdr    = bams[i]->hdr;
        aux[i].min_mq = args->min_mq;
        aux[i].itr    = sam_itr_querys(bams[i]->idx, bams[i]->hdr, region);
        data[i]       = &aux[i];
    }

    int     *n_plp = calloc(n_bams, sizeof(int));
    const bam_pileup1_t **plp = malloc(n_bams * sizeof(const bam_pileup1_t *));
    if (!n_plp || !plp) {
        free(n_plp); free(plp);
        for (int i = 0; i < n_bams; i++)
            if (aux[i].itr) hts_itr_destroy(aux[i].itr);
        free(aux); free(data);
        return -1;
    }

    bam_mplp_t mplp = bam_mplp_init(n_bams, cov_read, data);
    bam_mplp_set_maxcnt(mplp, args->max_depth);

    int64_t total_depth = 0;
    int tid, pos;

    while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) {
        if (tid != locus_tid)      continue;
        if (pos < (int)loc->start) continue;
        if (pos >= (int)loc->end)  break;

        int offset = pos - (int)loc->start;

        for (int i = 0; i < n_bams; i++) {
            /* Filter and count bases for this pileup column. We need per-base
             * counts for IUPAC samples; PHASED samples only need total depth,
             * but the per-base loop is identical so we count for both. */
            BaseCounts bc;
            memset(&bc, 0, sizeof(bc));
            for (int j = 0; j < n_plp[i]; j++) {
                const bam_pileup1_t *p = &plp[i][j];
                if (p->is_del || p->is_refskip) continue;
                uint8_t bq = bam_get_qual(p->b)[p->qpos];
                if ((int)bq < args->min_bq) continue;
                int bidx = b4_idx[bam_seqi(bam_get_seq(p->b), p->qpos)];
                if (bidx < 0) continue;
                bc.counts[bidx]++;
                bc.depth++;
            }
            total_depth += bc.depth;

            int off = seq_offset[i];
            if (modes[i] == SAMPLE_PHASED) {
                if (bc.depth >= args->min_dp) {
                    /* Covered invariant position: write reference base to
                     * both haplotypes. Pass 2 will overwrite if a VCF record
                     * exists here. */
                    char rb = ref_seq[offset];
                    seqs[off][offset]     = rb;
                    seqs[off + 1][offset] = rb;
                }
            } else {
                /* SAMPLE_IUPAC: write the IUPAC consensus from BAM pileup.
                 * call_iupac returns 'N' when depth < min_dp. */
                seqs[off][offset] = call_iupac(&bc, args->min_dp, args->het_freq);
            }
        }
    }

    bam_mplp_destroy(mplp);
    free(n_plp);
    free(plp);
    for (int i = 0; i < n_bams; i++)
        if (aux[i].itr) hts_itr_destroy(aux[i].itr);
    free(aux);
    free(data);

    return total_depth;
}

/* -------------------------------------------------------------------------
 * Pass 2: VCF phased genotypes
 *
 * Overwrites variant positions with the correct phased allele.
 * For each record:
 *   phased GT (a|b)  →  hap1 = allele[a], hap2 = allele[b]
 *   unphased GT (a/b) → both haplotypes = N  (or IUPAC, flag-controlled)
 *   missing GT        → both haplotypes = N
 *   DP < min_dp       → both haplotypes = N (mask low-confidence calls)
 *
 * Phase block changes within a locus generate a warning — the phase
 * assignment of ^1 and ^2 may be arbitrarily flipped across the boundary.
 * ---------------------------------------------------------------------- */

static void pass2_vcf(VcfPhase         *vcf,
                      int               n_bams,
                      const SampleMode *modes,
                      const int        *bam_to_vcf,   /* bam[i] → vcf column, -1 if absent */
                      const int        *seq_offset,
                      const Locus      *loc,
                      const Args       *args,
                      char            **seqs)
{
    char region[1024];
    snprintf(region, sizeof(region),
             "%s:%" PRId32 "-%" PRId32,
             loc->chrom, loc->start + 1, loc->end);

    hts_itr_t *itr = vcf->is_bcf
        ? bcf_itr_querys(vcf->idx, vcf->hdr, region)
        : tbx_itr_querys(vcf->tbx, region);
    if (!itr) return;   /* no records in this region */

    bcf1_t  *rec     = bcf_init();
    int32_t *gt_arr  = NULL;
    int      ngt_arr = 0;
    int32_t *dp_arr  = NULL;
    int      ndp_arr = 0;
    int32_t *ps_arr  = NULL;
    int      nps_arr = 0;
    kstring_t kst = {0, 0, NULL};   /* used for VCF.gz text rows */

    /* Track phase-set per sample to warn on within-locus PS changes */
    int32_t *prev_ps = calloc(n_bams, sizeof(int32_t));
    int      ps_warned = 0;
    for (int i = 0; i < n_bams; i++) prev_ps[i] = bcf_int32_missing;

    /* Iteration loop: read one record per pass, dispatching on format. */
    while (1) {
        if (vcf->is_bcf) {
            if (bcf_itr_next(vcf->fp, itr, rec) < 0) break;
        } else {
            if (tbx_itr_next(vcf->fp, vcf->tbx, itr, &kst) < 0) break;
            if (vcf_parse(&kst, vcf->hdr, rec) != 0) continue;
        }
        bcf_unpack(rec, BCF_UN_ALL);

        /* Only biallelic SNPs: skip indels, MNPs, multi-allelic sites */
        if (rec->n_allele != 2) continue;
        if (strlen(rec->d.allele[0]) != 1 ||
            strlen(rec->d.allele[1]) != 1) continue;

        char ref_base = (char)rec->d.allele[0][0];
        char alt_base = (char)rec->d.allele[1][0];
        if (base_to_idx(ref_base) < 0 || base_to_idx(alt_base) < 0) continue;

        int offset = (int)rec->pos - (int)loc->start;
        if (offset < 0 || offset >= (int)(loc->end - loc->start)) continue;

        /* Fetch GT array for all samples */
        int ngt = bcf_get_genotypes(vcf->hdr, rec, &gt_arr, &ngt_arr);
        if (ngt <= 0) continue;   /* no genotype data */

        /* Fetch per-sample DP (optional; -1 if absent) */
        bcf_get_format_int32(vcf->hdr, rec, "DP", &dp_arr, &ndp_arr);

        /* Fetch phase-set tag PS (optional) */
        bcf_get_format_int32(vcf->hdr, rec, "PS", &ps_arr, &nps_arr);

        /* ploidy = ngt / n_vcf_samples */
        int n_vcf = bcf_hdr_nsamples(vcf->hdr);
        int ploidy = (n_vcf > 0) ? ngt / n_vcf : 2;
        if (ploidy < 2) continue;   /* not diploid */

        for (int bi = 0; bi < n_bams; bi++) {
            if (modes[bi] != SAMPLE_PHASED) continue;   /* IUPAC samples filled by pass 1 */
            int vi = bam_to_vcf[bi];
            if (vi < 0) continue;                        /* defensive: shouldn't happen */

            int off = seq_offset[bi];

            /* Per-sample DP check */
            int dp = (dp_arr && vi < ndp_arr &&
                      dp_arr[vi] != bcf_int32_missing)
                     ? dp_arr[vi] : INT32_MAX;
            if (dp < args->min_dp) {
                seqs[off][offset]     = 'N';
                seqs[off + 1][offset] = 'N';
                continue;
            }

            int32_t gt0 = gt_arr[vi * ploidy];
            int32_t gt1 = gt_arr[vi * ploidy + 1];

            /* Missing genotype */
            if (bcf_gt_is_missing(gt0) || gt0 == bcf_int32_vector_end ||
                bcf_gt_is_missing(gt1) || gt1 == bcf_int32_vector_end) {
                seqs[off][offset]     = 'N';
                seqs[off + 1][offset] = 'N';
                continue;
            }

            int a0 = bcf_gt_allele(gt0);   /* 0 = ref, 1 = alt */
            int a1 = bcf_gt_allele(gt1);
            int phased = bcf_gt_is_phased(gt1);   /* gt1 carries phase flag */

            char base0 = (a0 == 0) ? ref_base : alt_base;
            char base1 = (a1 == 0) ? ref_base : alt_base;

            if (!phased && base0 != base1) {
                /* Unphased heterozygote */
                if (args->unphased_policy == UNPHASED_IUPAC) {
                    char code = iupac_pair(base0, base1);
                    seqs[off][offset]     = code;
                    seqs[off + 1][offset] = code;
                } else {
                    seqs[off][offset]     = 'N';
                    seqs[off + 1][offset] = 'N';
                }
                continue;
            }

            /* Phased (or homozygous): assign directly */
            seqs[off][offset]     = base0;
            seqs[off + 1][offset] = base1;

            /* Phase-block continuity check */
            if (!ps_warned && ps_arr && vi < nps_arr &&
                ps_arr[vi] != bcf_int32_missing) {
                int32_t ps = ps_arr[vi];
                if (prev_ps[bi] != bcf_int32_missing &&
                    prev_ps[bi] != ps) {
                    fprintf(stderr,
                        "Warning: phase-set boundary within locus '%s' "
                        "for sample '%s' (PS %d → %d). "
                        "^1/^2 assignment may flip across the boundary.\n",
                        loc->name,
                        /* BAM sample name not stored here; use index */
                        "(see locus stats)", prev_ps[bi], ps);
                    ps_warned = 1;
                }
                prev_ps[bi] = ps;
            }
        }
    }

    free(gt_arr);
    free(dp_arr);
    free(ps_arr);
    free(prev_ps);
    if (kst.s) free(kst.s);
    bcf_destroy(rec);
    hts_itr_destroy(itr);
}

/* -------------------------------------------------------------------------
 * Public: process_locus_vcf
 * ---------------------------------------------------------------------- */

int process_locus_vcf(BamFile         **bams,
                      int               n_bams,
                      const SampleMode *modes,
                      VcfPhase         *vcf,
                      const Locus      *loc,
                      const char       *ref_seq,
                      const Args       *args,
                      LocusResult      *result,
                      double           *mean_dp_out)
{
    int locus_len = (int)(loc->end - loc->start);

    /* ------------------------------------------------------------------
     * Per-sample seq layout: PHASED samples take 2 slots (^1/^2), IUPAC
     * samples take 1.  seq_offset[i] is the starting index in seqs[].
     * ------------------------------------------------------------------*/

    int *seq_offset = malloc(n_bams * sizeof(int));
    if (!seq_offset) return -1;

    int n_seqs = 0;
    for (int i = 0; i < n_bams; i++) {
        seq_offset[i] = n_seqs;
        n_seqs += (modes[i] == SAMPLE_PHASED) ? 2 : 1;
    }

    /* ------------------------------------------------------------------
     * Build BAM sample → VCF column index mapping (PHASED samples only).
     * IUPAC samples leave bam_to_vcf[i] = -1 even if present in the VCF.
     * ------------------------------------------------------------------*/

    int *bam_to_vcf = malloc(n_bams * sizeof(int));
    if (!bam_to_vcf) { free(seq_offset); return -1; }

    int n_vcf_samples = bcf_hdr_nsamples(vcf->hdr);
    for (int i = 0; i < n_bams; i++) {
        bam_to_vcf[i] = -1;
        if (modes[i] != SAMPLE_PHASED) continue;
        for (int j = 0; j < n_vcf_samples; j++) {
            if (strcmp(vcf->hdr->samples[j], bams[i]->sample) == 0) {
                bam_to_vcf[i] = j;
                break;
            }
        }
        /* PHASED implies present in VCF (classified upstream). If absent
         * here it's a caller bug; leave as -1 and pass2 will skip. */
    }

    /* ------------------------------------------------------------------
     * Allocate sequences — all 'N' until a pass writes something.
     * ------------------------------------------------------------------*/

    char **seqs      = calloc(n_seqs, sizeof(char *));
    char **seq_names = calloc(n_seqs, sizeof(char *));
    if (!seqs || !seq_names) goto fail;

    for (int i = 0; i < n_seqs; i++) {
        seqs[i] = malloc(locus_len + 1);
        if (!seqs[i]) goto fail;
        memset(seqs[i], 'N', locus_len);
        seqs[i][locus_len] = '\0';
    }

    for (int i = 0; i < n_bams; i++) {
        size_t slen = strlen(bams[i]->sample);
        int off = seq_offset[i];
        if (modes[i] == SAMPLE_PHASED) {
            seq_names[off]     = malloc(slen + 3);
            seq_names[off + 1] = malloc(slen + 3);
            if (!seq_names[off] || !seq_names[off + 1]) goto fail;
            sprintf(seq_names[off],     "%s^1", bams[i]->sample);
            sprintf(seq_names[off + 1], "%s^2", bams[i]->sample);
        } else {
            seq_names[off] = strdup(bams[i]->sample);
            if (!seq_names[off]) goto fail;
        }
    }

    /* ------------------------------------------------------------------
     * Pass 1 (unified): BAM pileup — REF for PHASED, IUPAC consensus for
     *                                 SAMPLE_IUPAC.
     * Pass 2 (VCF):     phased GTs overwrite variant positions for PHASED.
     * ------------------------------------------------------------------*/

    int64_t total_depth = pass1_unified(bams, n_bams, modes, seq_offset,
                                         loc, ref_seq, args, seqs);
    if (total_depth < 0) goto fail;

    pass2_vcf(vcf, n_bams, modes, bam_to_vcf, seq_offset, loc, args, seqs);

    /* ------------------------------------------------------------------
     * Fill result.
     * ------------------------------------------------------------------*/

    result->name      = strdup(loc->name);
    result->seqs      = seqs;
    result->seq_names = seq_names;
    result->n_seqs    = n_seqs;
    result->locus_len = locus_len;

    *mean_dp_out = (n_bams > 0 && locus_len > 0)
                   ? (double)total_depth / ((double)n_bams * locus_len)
                   : 0.0;

    free(seq_offset);
    free(bam_to_vcf);
    return 0;

fail:
    fprintf(stderr, "Error: out of memory in process_locus_vcf for '%s'\n",
            loc->name);
    if (seqs) {
        for (int i = 0; i < n_seqs; i++) free(seqs[i]);
        free(seqs);
    }
    if (seq_names) {
        for (int i = 0; i < n_seqs; i++) free(seq_names[i]);
        free(seq_names);
    }
    free(seq_offset);
    free(bam_to_vcf);
    return -1;
}

/* ────────────────────────────────────────────────────────────────────────
 * Whole-file phasing classification
 * ───────────────────────────────────────────────────────────────────── */

VcfPhaseSampleStat *vcf_phase_classify(VcfPhase     *vcf,
                                       BamFile     **bams,
                                       int           n_bams,
                                       const Locus  *loci,
                                       int           n_loci)
{
    if (!vcf || !bams || n_bams <= 0) return NULL;

    VcfPhaseSampleStat *stats =
        (VcfPhaseSampleStat *)calloc((size_t)n_bams, sizeof(VcfPhaseSampleStat));
    if (!stats) return NULL;

    int *bam_to_vcf = (int *)malloc((size_t)n_bams * sizeof(int));
    int  n_vcf_samples = bcf_hdr_nsamples(vcf->hdr);
    for (int i = 0; i < n_bams; i++) {
        stats[i].sample = strdup(bams[i]->sample);
        bam_to_vcf[i] = -1;
        for (int j = 0; j < n_vcf_samples; j++) {
            if (strcmp(vcf->hdr->samples[j], bams[i]->sample) == 0) {
                bam_to_vcf[i] = j;
                break;
            }
        }
        if (bam_to_vcf[i] < 0) stats[i].not_in_vcf = 1;
    }

    bcf1_t   *rec     = bcf_init();
    int32_t  *gt_arr  = NULL;
    int       ngt_arr = 0;
    char      region[1024];
    kstring_t kst = {0, 0, NULL};

    for (int li = 0; li < n_loci; li++) {
        snprintf(region, sizeof(region),
                 "%s:%" PRId32 "-%" PRId32,
                 loci[li].chrom, loci[li].start + 1, loci[li].end);
        hts_itr_t *itr = vcf->is_bcf
            ? bcf_itr_querys(vcf->idx, vcf->hdr, region)
            : tbx_itr_querys(vcf->tbx, region);
        if (!itr) continue;

        while (1) {
            if (vcf->is_bcf) {
                if (bcf_itr_next(vcf->fp, itr, rec) < 0) break;
            } else {
                if (tbx_itr_next(vcf->fp, vcf->tbx, itr, &kst) < 0) break;
                if (vcf_parse(&kst, vcf->hdr, rec) != 0) continue;
            }
            bcf_unpack(rec, BCF_UN_ALL);
            if (rec->n_allele < 2) continue;

            int ngt = bcf_get_genotypes(vcf->hdr, rec, &gt_arr, &ngt_arr);
            if (ngt <= 0) continue;
            int ploidy = (n_vcf_samples > 0) ? ngt / n_vcf_samples : 2;
            if (ploidy < 2) continue;

            for (int bi = 0; bi < n_bams; bi++) {
                int vi = bam_to_vcf[bi];
                if (vi < 0) continue;
                int32_t gt0 = gt_arr[vi * ploidy];
                int32_t gt1 = gt_arr[vi * ploidy + 1];
                if (bcf_gt_is_missing(gt0) || bcf_gt_is_missing(gt1) ||
                    gt0 == bcf_int32_vector_end || gt1 == bcf_int32_vector_end) {
                    continue;
                }
                int a0 = bcf_gt_allele(gt0);
                int a1 = bcf_gt_allele(gt1);
                if (a0 == a1) continue;  /* only count heterozygous calls */
                stats[bi].n_het++;
                if (bcf_gt_is_phased(gt1)) stats[bi].n_phased++;
                else                       stats[bi].n_unphased++;
            }
        }
        hts_itr_destroy(itr);
    }

    free(gt_arr);
    if (kst.s) free(kst.s);
    bcf_destroy(rec);
    free(bam_to_vcf);
    return stats;
}

void vcf_phase_stats_free(VcfPhaseSampleStat *s, int n)
{
    if (!s) return;
    for (int i = 0; i < n; i++) free(s[i].sample);
    free(s);
}
