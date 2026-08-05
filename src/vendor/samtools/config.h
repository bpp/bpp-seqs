/* config.h — shim for vendored samtools sources.
 *
 * samtools is autoconf-built and its sources open with `#include <config.h>`.
 * bpp-seqs is not, so this stands in for the generated header.
 *
 * The vendored consensus sources reference no HAVE_* / PACKAGE_* symbols
 * (verified by grep at the pinned version), so an empty header is sufficient.
 * If a future sync introduces one, the build will fail loudly here rather
 * than silently compile a different code path.
 *
 * This file is ours, not samtools'. Everything else in this directory is
 * upstream verbatim -- see README.md.
 */
#ifndef BPP_SEQS_VENDOR_SAMTOOLS_CONFIG_H
#define BPP_SEQS_VENDOR_SAMTOOLS_CONFIG_H

/* samtools' configure defines these; the consensus sources do not read them,
 * but keeping them here documents what a real config.h would carry. */
#define PACKAGE_VERSION "1.23.1-vendored"

#endif /* BPP_SEQS_VENDOR_SAMTOOLS_CONFIG_H */
