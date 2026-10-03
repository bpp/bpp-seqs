/* lineio.h — read whole lines of any length from a (possibly gzipped) file.
 *
 * A BPP sequence file keeps each sequence on one line, and loci of hundreds
 * of kilobases are ordinary. Reading those through a fixed buffer splits a
 * row, and the remainder is then misread as a new row: an extra taxon, or a
 * malformed locus that ends the parse early. Anything that reads sequence
 * rows goes through this instead. */
#ifndef BPP_SEQS_LINEIO_H
#define BPP_SEQS_LINEIO_H

#include <stddef.h>
#include <zlib.h>

/* Read one whole line, newline included, into *buf (allocated or grown as
 * needed; the caller frees it). Returns the line's length, or -1 at EOF. */
long gz_getline(gzFile gz, char **buf, size_t *cap);

#endif
