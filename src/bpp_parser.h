/* bpp_parser.h — parse an existing BPP-format sequence file.
 *
 * The file format (per https://bpp.github.io/bpp-manual/) consists of one
 * or more locus blocks, each beginning with `<n_seqs> <n_sites>` on a
 * line by itself, followed by n_seqs sequence lines of the form
 *      <optional_label>^<id>    SEQUENCE
 * or
 *      ^<id>    SEQUENCE
 * Sequence data may span multiple lines (interleaved); the parser
 * accumulates characters until each sample's sequence reaches n_sites.
 *
 * Optionally, a sibling <prefix>.loci.tsv file produced by bpp-seqs
 * carries per-locus provenance (source chrom/start/end/stride etc.); it
 * can be attached to the parsed array with bpp_attach_loci_tsv().
 */
#ifndef BPP_SEQS_BPP_PARSER_H
#define BPP_SEQS_BPP_PARSER_H

#include <stdio.h>

typedef struct {
    int    n_seqs;
    int    n_sites;
    char **seq_names;       /* raw sequence names as they appear in the file */
    char **ids;              /* the part after the last '^' (the BPP id) */
    char **seqs;             /* size n_seqs; each NUL-terminated, length n_sites */

    /* Optional metadata, populated by bpp_attach_loci_tsv() when available. */
    char  *name;             /* locus_name from .loci.tsv (else "locus<N>") */
    char  *source_kind;
    char  *source_file;
    char  *source_chrom;
    int    source_start;     /* -1 if not provided */
    int    source_end;
    int    source_stride;
} BppLocus;

/* Parse the BPP file at `path`. Returns a malloc'd array of n loci and
 * sets *n_out. Returns NULL on error. Caller frees with bpp_loci_free(). */
BppLocus *bpp_parse_file(const char *path, int *n_out);

/* Attach .loci.tsv metadata to a parsed array. Matches by locus_name when
 * possible; otherwise by index. Returns the number of rows attached. */
int bpp_attach_loci_tsv(BppLocus *loci, int n_loci, const char *tsv_path);

void bpp_loci_free(BppLocus *arr, int n);

#endif
