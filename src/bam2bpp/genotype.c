/*
 * genotype.c  –  base calling from raw allele counts
 *
 * All three calling modes share the same allele-ranking logic; they differ
 * only in what they do with a heterozygous site.
 *
 * Homozygous call: major allele frequency ≥ (1 - het_freq)
 * Heterozygous call: minor allele frequency ≥ het_freq
 * Low-depth or no-data: 'N'
 */

#include <string.h>
#include "bam2bpp.h"

/* IUPAC ambiguity codes for pairs of ACGT bases.
 * Indexed iupac_tbl[lo][hi] where lo < hi (both in 0-3 = ACGT order).
 * Diagonal entries are canonical bases (homozygous). */
static const char iupac_tbl[4][4] = {
    /* second base: A    C    G    T  */
    /* A */       {'A', 'M', 'R', 'W'},
    /* C */       {'M', 'C', 'S', 'Y'},
    /* G */       {'R', 'S', 'G', 'K'},
    /* T */       {'W', 'Y', 'K', 'T'},
};

static const char base_char[4] = {'A', 'C', 'G', 'T'};

/*
 * top_two() – find the indices of the two most-frequent alleles.
 * *a1 receives the most-frequent index, *a2 the second-most.
 * If the depth is 0, *a1 is set to -1.
 */
static void top_two(const BaseCounts *bc, int *a1, int *a2)
{
    *a1 = -1;
    *a2 = -1;

    for (int i = 0; i < 4; i++) {
        if (bc->counts[i] == 0) continue;
        if (*a1 < 0 || bc->counts[i] > bc->counts[*a1]) {
            *a2 = *a1;
            *a1 = i;
        } else if (*a2 < 0 || bc->counts[i] > bc->counts[*a2]) {
            *a2 = i;
        }
    }
}

/* -------------------------------------------------------------------------
 * Public interface
 * ---------------------------------------------------------------------- */

char call_haploid(const BaseCounts *bc, int min_dp)
{
    if (bc->depth < min_dp) return 'N';
    int a1, a2;
    top_two(bc, &a1, &a2);
    if (a1 < 0) return 'N';
    return base_char[a1];
}

char call_iupac(const BaseCounts *bc, int min_dp, double het_freq)
{
    if (bc->depth < min_dp) return 'N';

    int a1, a2;
    top_two(bc, &a1, &a2);
    if (a1 < 0) return 'N';

    /* Heterozygous? */
    if (a2 >= 0 && bc->counts[a2] > 0) {
        double maf = (double)bc->counts[a2] / bc->depth;
        if (maf >= het_freq) {
            int lo = (a1 < a2) ? a1 : a2;
            int hi = (a1 < a2) ? a2 : a1;
            return iupac_tbl[lo][hi];
        }
    }

    /* Homozygous */
    return base_char[a1];
}

void call_split(const BaseCounts *bc, int min_dp, double het_freq,
                char *b1, char *b2)
{
    if (bc->depth < min_dp) {
        *b1 = *b2 = 'N';
        return;
    }

    int a1, a2;
    top_two(bc, &a1, &a2);
    if (a1 < 0) {
        *b1 = *b2 = 'N';
        return;
    }

    if (a2 >= 0 && bc->counts[a2] > 0) {
        double maf = (double)bc->counts[a2] / bc->depth;
        if (maf >= het_freq) {
            /* Return major as b1, minor as b2 (arbitrary but deterministic) */
            *b1 = base_char[a1];
            *b2 = base_char[a2];
            return;
        }
    }

    *b1 = *b2 = base_char[a1];
}
