/*
 * pileup.c  –  BAM access and per-locus multi-sample pileup
 *
 * Key design principle
 * --------------------
 * Every sequence in the output is initialised to all-'N'.  Characters are
 * overwritten only when the pileup yields a position with depth ≥ min_dp.
 * Positions with zero coverage (never yielded by bam_mplp_auto) and
 * positions with depth below min_dp both remain 'N'.  The reference
 * sequence is never used to fill uncovered sites.
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <inttypes.h>

#include <htslib/sam.h>
#include <htslib/hts.h>

#include "bam2bpp.h"

/* -------------------------------------------------------------------------
 * Per-BAM data passed to the pileup callback
 * ---------------------------------------------------------------------- */

typedef struct {
    samFile   *fp;
    sam_hdr_t *hdr;
    hts_itr_t *itr;    /* region iterator for the current locus */
    int        min_mq;
} PlpAux;

/*
 * Pileup read callback.  htslib calls this to fetch the next read for one
 * BAM file.  We skip reads that fail standard QC flags or fall below the
 * mapping-quality threshold.
 */
static int plp_read(void *data, bam1_t *b)
{
    PlpAux *aux = (PlpAux *)data;
    int ret;

    while ((ret = sam_itr_next(aux->fp, aux->itr, b)) >= 0) {
        /* Skip unmapped, secondary, QC-fail, duplicate, supplementary */
        if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY |
                            BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY))
            continue;

        if ((int)b->core.qual < aux->min_mq)
            continue;

        break;
    }

    return ret;
}

/* -------------------------------------------------------------------------
 * Extract sample name from BAM header @RG SM tag.
 * Falls back to the filename stem when no @RG line is present.
 * ---------------------------------------------------------------------- */

static char *extract_sample_name(sam_hdr_t *hdr, const char *path)
{
    /* Walk the header text looking for @RG\t...\tSM:<value> */
    const char *text = sam_hdr_str(hdr);
    if (text) {
        const char *rg = strstr(text, "@RG\t");
        while (rg) {
            /* Stay within this @RG line */
            const char *eol = strchr(rg, '\n');
            const char *sm  = strstr(rg, "\tSM:");
            if (sm && (!eol || sm < eol)) {
                sm += 4;   /* skip "\tSM:" */
                const char *end = sm;
                while (*end && *end != '\t' && *end != '\n') end++;
                char *name = malloc(end - sm + 1);
                if (name) {
                    memcpy(name, sm, end - sm);
                    name[end - sm] = '\0';
                    return name;
                }
            }
            rg = strstr(rg + 4, "@RG\t");
        }
    }

    /* Fallback: filename stem (strip path and extension) */
    const char *base = strrchr(path, '/');
    base = base ? base + 1 : path;
    char *name = strdup(base);
    if (name) {
        char *dot = strrchr(name, '.');
        if (dot) *dot = '\0';
    }
    return name;
}

/* -------------------------------------------------------------------------
 * Open / close BAMs
 * ---------------------------------------------------------------------- */

BamFile **open_bams(char **paths, int n, int max_depth)
{
    (void)max_depth;   /* used per-locus in bam_mplp_set_maxcnt */

    BamFile **bams = calloc(n, sizeof(BamFile *));
    if (!bams) return NULL;

    for (int i = 0; i < n; i++) {
        BamFile *b = calloc(1, sizeof(BamFile));
        if (!b) goto fail;
        bams[i] = b;

        b->path = strdup(paths[i]);
        b->fp   = sam_open(paths[i], "rb");
        if (!b->fp) {
            fprintf(stderr, "Error: cannot open BAM '%s': ", paths[i]);
            perror(NULL);
            goto fail;
        }

        b->hdr = sam_hdr_read(b->fp);
        if (!b->hdr) {
            fprintf(stderr, "Error: cannot read header from '%s'\n", paths[i]);
            goto fail;
        }

        b->idx = sam_index_load(b->fp, paths[i]);
        if (!b->idx) {
            fprintf(stderr,
                    "Error: no index for '%s' – run 'samtools index' first\n",
                    paths[i]);
            goto fail;
        }

        b->sample = extract_sample_name(b->hdr, paths[i]);
    }

    return bams;

fail:
    close_bams(bams, n);
    return NULL;
}

void close_bams(BamFile **bams, int n)
{
    if (!bams) return;
    for (int i = 0; i < n; i++) {
        BamFile *b = bams[i];
        if (!b) continue;
        if (b->idx) hts_idx_destroy(b->idx);
        if (b->hdr) sam_hdr_destroy(b->hdr);
        if (b->fp)  sam_close(b->fp);
        free(b->sample);
        free(b->path);
        free(b);
    }
    free(bams);
}

/* -------------------------------------------------------------------------
 * free_locus_result
 * ---------------------------------------------------------------------- */

void free_locus_result(LocusResult *r)
{
    if (!r) return;
    free(r->name);
    if (r->seqs) {
        for (int i = 0; i < r->n_seqs; i++) free(r->seqs[i]);
        free(r->seqs);
    }
    if (r->seq_names) {
        for (int i = 0; i < r->n_seqs; i++) free(r->seq_names[i]);
        free(r->seq_names);
    }
    r->name      = NULL;
    r->seqs      = NULL;
    r->seq_names = NULL;
}

/* -------------------------------------------------------------------------
 * process_locus  –  the main pileup engine
 * ---------------------------------------------------------------------- */

/*
 * Lookup table: 4-bit htslib base encoding → ACGT index (0-3), or -1.
 *   htslib encoding: A=1, C=2, G=4, T=8, N=15, ambiguous=others
 */
static const int b4_to_idx[16] = {
    -1,       /* 0  = padding/unknown  */
     0,       /* 1  = A               */
     1,       /* 2  = C               */
    -1,       /* 3  = M (AC)          */
     2,       /* 4  = G               */
    -1,       /* 5  = R (AG)          */
    -1,       /* 6  = S (CG)          */
    -1,       /* 7  = V (ACG)         */
     3,       /* 8  = T               */
    -1,       /* 9  = W (AT)          */
    -1,       /* 10 = Y (CT)          */
    -1,       /* 11 = H (ACT)         */
    -1,       /* 12 = K (GT)          */
    -1,       /* 13 = D (AGT)         */
    -1,       /* 14 = B (CGT)         */
    -1,       /* 15 = N               */
};

int process_locus(BamFile **bams, int n_bams,
                  const Locus *loc, const Args *args,
                  LocusResult *result, double *mean_dp_out)
{
    int locus_len = (int)(loc->end - loc->start);
    int n_seqs    = (args->phasing == PHASE_SPLIT) ? 2 * n_bams : n_bams;

    /* Region string for htslib iterators: 1-based, inclusive */
    char region[1024];
    snprintf(region, sizeof(region),
             "%s:%" PRId32 "-%" PRId32, loc->chrom,
             loc->start + 1, loc->end);

    /* ------------------------------------------------------------------
     * Allocate output sequences and names; initialise all bases to 'N'.
     * Positions left as 'N' after the pileup loop have zero or low
     * coverage – they are never overwritten with the reference base.
     * ------------------------------------------------------------------ */

    char **seqs      = calloc(n_seqs, sizeof(char *));
    char **seq_names = calloc(n_seqs, sizeof(char *));
    if (!seqs || !seq_names) goto fail_alloc;

    for (int i = 0; i < n_seqs; i++) {
        seqs[i] = malloc(locus_len + 1);
        if (!seqs[i]) goto fail_alloc;
        memset(seqs[i], 'N', locus_len);
        seqs[i][locus_len] = '\0';
    }

    /* Assign sequence names */
    for (int i = 0; i < n_bams; i++) {
        if (args->phasing == PHASE_SPLIT) {
            size_t slen = strlen(bams[i]->sample);
            seq_names[2 * i]     = malloc(slen + 3);
            seq_names[2 * i + 1] = malloc(slen + 3);
            if (!seq_names[2 * i] || !seq_names[2 * i + 1]) goto fail_alloc;
            sprintf(seq_names[2 * i],     "%s^1", bams[i]->sample);
            sprintf(seq_names[2 * i + 1], "%s^2", bams[i]->sample);
        } else {
            seq_names[i] = strdup(bams[i]->sample);
            if (!seq_names[i]) goto fail_alloc;
        }
    }

    /* ------------------------------------------------------------------
     * Set up per-BAM pileup auxiliary structures and region iterators.
     * ------------------------------------------------------------------ */

    PlpAux  *aux  = calloc(n_bams, sizeof(PlpAux));
    void   **data = malloc(n_bams * sizeof(void *));
    if (!aux || !data) { free(aux); free(data); goto fail_alloc; }

    int locus_tid = sam_hdr_name2tid(bams[0]->hdr, loc->chrom);
    if (locus_tid < 0) {
        fprintf(stderr, "Error: contig '%s' not found in BAM header\n",
                loc->chrom);
        free(aux); free(data);
        goto fail_alloc;
    }

    for (int i = 0; i < n_bams; i++) {
        aux[i].fp    = bams[i]->fp;
        aux[i].hdr   = bams[i]->hdr;
        aux[i].min_mq = args->min_mq;
        aux[i].itr   = sam_itr_querys(bams[i]->idx, bams[i]->hdr, region);
        data[i]      = &aux[i];

        if (!aux[i].itr) {
            /* Region absent from this BAM – all positions stay 'N' */
            fprintf(stderr,
                    "Warning: region %s not found in %s; all positions will be N\n",
                    region, bams[i]->path);
        }
    }

    /* ------------------------------------------------------------------
     * Multi-sample pileup loop
     * ------------------------------------------------------------------ */

    int *n_plp = calloc(n_bams, sizeof(int));
    const bam_pileup1_t **plp =
        malloc(n_bams * sizeof(const bam_pileup1_t *));
    if (!n_plp || !plp) {
        free(n_plp); free(plp);
        for (int i = 0; i < n_bams; i++)
            if (aux[i].itr) hts_itr_destroy(aux[i].itr);
        free(aux); free(data);
        goto fail_alloc;
    }

    bam_mplp_t mplp = bam_mplp_init(n_bams, plp_read, data);
    bam_mplp_set_maxcnt(mplp, args->max_depth);

    int64_t total_depth = 0;
    int64_t n_observed  = 0;   /* sample × position pairs yielded by mplp */
    int     tid, pos;

    while (bam_mplp_auto(mplp, &tid, &pos, n_plp, plp) > 0) {

        /* Stay within the locus window */
        if (tid != locus_tid)          continue;
        if (pos < (int)loc->start)     continue;
        if (pos >= (int)loc->end)      break;

        int offset = pos - (int)loc->start;

        for (int i = 0; i < n_bams; i++) {

            BaseCounts bc;
            memset(&bc, 0, sizeof(bc));

            /* Count filtered bases in this pileup column for sample i */
            for (int j = 0; j < n_plp[i]; j++) {
                const bam_pileup1_t *p = &plp[i][j];

                /* Skip deletions and reference skips (introns) */
                if (p->is_del || p->is_refskip) continue;

                /* Base-quality filter */
                uint8_t bq = bam_get_qual(p->b)[p->qpos];
                if ((int)bq < args->min_bq) continue;

                /* Convert 4-bit encoded base to ACGT index */
                int bidx = b4_to_idx[bam_seqi(bam_get_seq(p->b), p->qpos)];
                if (bidx < 0) continue;   /* ambiguous base – skip */

                bc.counts[bidx]++;
                bc.depth++;
            }

            total_depth += bc.depth;
            n_observed++;

            /*
             * Call and write the base(s).
             *
             * If bc.depth < min_dp the call functions return 'N', which
             * leaves the sequence unchanged (already 'N' from memset).
             * We write explicitly for clarity.
             */
            switch (args->phasing) {
                case PHASE_SPLIT: {
                    char b1, b2;
                    call_split(&bc, args->min_dp, args->het_freq, &b1, &b2);
                    seqs[2 * i][offset]     = b1;
                    seqs[2 * i + 1][offset] = b2;
                    break;
                }
                case PHASE_IUPAC:
                    seqs[i][offset] = call_iupac(&bc, args->min_dp,
                                                 args->het_freq);
                    break;
                case PHASE_HAPLOID:
                    seqs[i][offset] = call_haploid(&bc, args->min_dp);
                    break;
                default:
                    /* PHASE_VCF should never reach process_locus —
                     * main() dispatches to process_locus_vcf instead.
                     * Fall back to haploid (coverage pass only). */
                    seqs[i][offset] = call_haploid(&bc, args->min_dp);
                    break;
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

    /* Fill result struct */
    result->name      = strdup(loc->name);
    result->seqs      = seqs;
    result->seq_names = seq_names;
    result->n_seqs    = n_seqs;
    result->locus_len = locus_len;

    *mean_dp_out = (n_observed > 0)
                   ? (double)total_depth / ((double)n_bams * locus_len)
                   : 0.0;

    return 0;

fail_alloc:
    fprintf(stderr, "Error: out of memory processing locus '%s'\n", loc->name);
    if (seqs) {
        for (int i = 0; i < n_seqs; i++) free(seqs[i]);
        free(seqs);
    }
    if (seq_names) {
        for (int i = 0; i < n_seqs; i++) free(seq_names[i]);
        free(seq_names);
    }
    return -1;
}
