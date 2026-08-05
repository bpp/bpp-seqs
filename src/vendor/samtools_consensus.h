/* samtools_consensus.h — bpp-seqs interface to the vendored samtools
 * consensus caller.
 *
 * bpp-seqs' own caller ranks raw allele counts against a fixed minor-allele
 * fraction: no error model, base qualities used only as a cutoff. This wraps
 * the Gap5/samtools Bayesian consensus model instead, which weighs each
 * observation by its base quality, caps that by mapping quality, adjusts for
 * local quality minima and NM density, and returns a phred-scaled confidence
 * so poorly supported positions can become 'N' rather than a confident guess.
 *
 * The model is samtools' `consensus` subcommand, vendored verbatim under
 * vendor/samtools/ -- see that directory's README.md for the pinned version
 * and the sync procedure. The point of vendoring rather than reimplementing is
 * that `samtools consensus -A` can then serve as a test oracle for output
 * equality, which a clean-room reimplementation could only approximate.
 *
 * Defaults below mirror `samtools consensus` so that, given the same options,
 * this produces the same bases. Where bpp-seqs needs different behaviour
 * (fixed reference-coordinate output for per-locus alignments) it is stated.
 */
#ifndef BPP_SEQS_SAMTOOLS_CONSENSUS_H
#define BPP_SEQS_SAMTOOLS_CONSENSUS_H

#include <htslib/sam.h>

typedef struct BppsCons BppsCons;

typedef struct {
    /* Emit IUPAC ambiguity codes at heterozygous sites (samtools -A/--ambig).
     * BPP reads these as unphased heterozygotes under `phase = 1`. */
    int    ambig;
    /* Phred confidence below which a call becomes 'N' (samtools -C/--cutoff). */
    int    cons_cutoff;
    /* Positions covered by fewer than this many reads become 'N'
     * (samtools -d/--min-depth). */
    int    min_depth;
    /* Prior probability that a site is heterozygous (samtools --P-het).
     * A fixed scalar, never estimated from the data, so a call does not
     * depend on how many other individuals happen to be sequenced at the
     * locus -- which matters here because locus sampling varies. */
    double p_het;
    /* Hard filters applied before pileup (samtools --min-MQ / --min-BQ). */
    int    min_mqual;
    int    min_qual;
    /* Character written where a deletion is involved: samtools writes '*' for
     * a homozygous deletion and a lowercase base for a base/deletion
     * heterozygote. BPP can read neither (IUPAC has no base/gap code), so both
     * become this character -- 'N' or '-'. Set it to '*' to pass upstream's
     * bytes through unchanged, which is what the oracle test compares. */
    char   del_char;
} BppsConsOpts;

/* Fill `o` with the same defaults `samtools consensus` uses, except ambig,
 * which bpp-seqs turns on (BPP wants IUPAC heterozygotes). */
void bpps_cons_opts_init(BppsConsOpts *o);

/* Build the model's probability tables. Returns NULL on allocation failure. */
BppsCons *bpps_cons_init(const BppsConsOpts *o);
void      bpps_cons_free(BppsCons *c);

/* Call consensus over [beg, end) of one already-open, indexed BAM/CRAM.
 *
 * beg/end are 0-based half-open, matching BED and bpp-seqs' own Locus. (Note
 * that `samtools consensus -r` takes 1-based inclusive coordinates instead, so
 * the equivalent command line is `-r chrom:beg-(end-1)` -- and position 0
 * cannot be expressed there at all, which is a limit of that CLI, not of this
 * function.)
 *
 * Writes exactly (end - beg) bases into seq_out, in reference coordinates:
 * insertions are skipped and deletions are written as the reference-position
 * gap they occupy, so the result lines up across samples without alignment.
 * That differs from `samtools consensus` run bare, which by default reports
 * insertions and so returns a sequence longer than the span; the oracle test
 * compares against `--show-ins no --show-del yes`, which matches this.
 *
 * Positions with no read coverage are left as 'N'; seq_out is pre-filled.
 *
 * Returns 0 on success, -1 on error. */
int bpps_cons_region(BppsCons *c,
                     samFile *fp, sam_hdr_t *h, hts_idx_t *idx,
                     const char *chrom, hts_pos_t beg, hts_pos_t end,
                     char *seq_out);

/* Version of the vendored samtools source, for --version and provenance. */
const char *bpps_cons_samtools_version(void);

#endif /* BPP_SEQS_SAMTOOLS_CONSENSUS_H */
