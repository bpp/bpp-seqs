/* cmd_mask.h — `bpp-seqs mask`: apply per-sample mappability masks to an
 * existing BPP sequence file, replacing masked positions with 'N'.
 *
 * Kept out of the conversion path deliberately. Masking is a text rewrite over
 * an already-called alignment, so doing it here means an existing build can be
 * re-masked without re-running the base caller, different masks can be
 * compared without redoing any calling, and the caller itself -- whose output
 * is pinned byte-for-byte against samtools -- is not touched at all.
 */
#ifndef BPP_SEQS_CMD_MASK_H
#define BPP_SEQS_CMD_MASK_H

int cmd_mask(int argc, char **argv);

#endif /* BPP_SEQS_CMD_MASK_H */
