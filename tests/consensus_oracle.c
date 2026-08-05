/* consensus_oracle.c — test harness for the vendored samtools consensus model.
 *
 * Prints the consensus of one region as FASTA so it can be diffed against
 * `samtools consensus`. Built on demand by tests/run_tests.sh; not part of the
 * bpp-seqs binary.
 *
 * usage: consensus_oracle BAM CHROM BEG END [cutoff] [min_depth] [del_char]
 *        BEG/END are 0-based half-open, i.e. the same convention as a BED.
 */
#include <stdio.h>
#include <stdlib.h>

#include <htslib/sam.h>

#include "samtools_consensus.h"

int main(int argc, char **argv)
{
    if (argc < 5) {
        fprintf(stderr, "usage: %s BAM CHROM BEG END "
                        "[cutoff] [min_depth] [del_char]\n", argv[0]);
        return 2;
    }
    const char *bam = argv[1], *chrom = argv[2];
    hts_pos_t beg = atoll(argv[3]), end = atoll(argv[4]);
    if (end <= beg) { fprintf(stderr, "empty region\n"); return 2; }

    BppsConsOpts o;
    bpps_cons_opts_init(&o);
    if (argc > 5) o.cons_cutoff = atoi(argv[5]);
    if (argc > 6) o.min_depth   = atoi(argv[6]);
    if (argc > 7) o.del_char    = argv[7][0];

    samFile *fp = sam_open(bam, "r");
    if (!fp) { perror(bam); return 1; }
    sam_hdr_t *h   = sam_hdr_read(fp);
    hts_idx_t *idx = sam_index_load(fp, bam);
    if (!h || !idx) { fprintf(stderr, "cannot read header or index\n"); return 1; }

    char *seq = (char *)malloc((size_t)(end - beg) + 1);
    if (!seq) return 1;
    seq[end - beg] = '\0';

    BppsCons *c = bpps_cons_init(&o);
    if (!c || bpps_cons_region(c, fp, h, idx, chrom, beg, end, seq) != 0) {
        fprintf(stderr, "consensus failed\n");
        return 1;
    }

    printf(">%s\n", chrom);
    for (hts_pos_t i = 0; i < end - beg; i += 70) {
        hts_pos_t n = end - beg - i;
        printf("%.*s\n", (int)(n < 70 ? n : 70), seq + i);
    }

    bpps_cons_free(c);
    free(seq);
    hts_idx_destroy(idx);
    sam_hdr_destroy(h);
    sam_close(fp);
    return 0;
}
