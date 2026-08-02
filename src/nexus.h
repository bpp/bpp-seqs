#ifndef BPP_SEQS_NEXUS_H
#define BPP_SEQS_NEXUS_H

/* ---------------------------------------------------------------------------
 * NEXUS reader (Maddison, Swofford & Maddison 1997, Syst. Biol. 46:590-621).
 *
 * One spec-based reader shared by the inspector (inspect.c) and the converter
 * (converters/nexus2bpp.c), replacing the two ad-hoc line/substring scanners
 * that disagreed with each other and mishandled comments, quoted labels and
 * whitespace around '='.
 *
 * The reader is whitespace-insensitive by construction: it strips nested
 * `[ ... ]` comments (honouring single-quoted tokens), then parses the
 * DIMENSIONS / FORMAT / MATRIX / CHARSET commands from the cleaned token
 * stream. Single-quoted labels are unquoted (interior spaces -> '_', '' -> ').
 * INTERLEAVE and MATCHCHAR are honoured; {..}/(..) ambiguity groups count as
 * one site ('N').
 * ------------------------------------------------------------------------- */

typedef struct {
    int     ntax;
    int     nchar;
    char  **taxa;        /* [ntax] labels (owned) */
    char  **seq;         /* [ntax] sequences, each nchar long, uppercased,
                            matchchar-expanded (owned) */

    int     n_charsets;
    char  **cs_name;     /* [n_charsets] (owned) */
    int    *cs_start;    /* 1-based inclusive; end == -1 means "to last site" */
    int    *cs_end;
    int    *cs_stride;

    double  missing_fraction;
    char    datatype[16];   /* lowercased, e.g. "dna" ("" if unspecified) */
    int     interleaved;
    char    err[256];       /* "" on success, else a diagnostic message */
} NexusDoc;

/* Parse a (possibly gzipped) NEXUS file. Returns 0 on success. On failure
 * returns non-zero and sets doc->err. Call nexus_free(doc) either way (it is
 * safe on a partially-populated or zeroed doc). */
int  nexus_parse(const char *path, NexusDoc *doc);
void nexus_free(NexusDoc *doc);

#endif /* BPP_SEQS_NEXUS_H */
