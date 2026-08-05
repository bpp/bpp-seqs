/* samtools_consensus.c — bpp-seqs shim around the vendored samtools
 * consensus caller.
 *
 * The model lives in vendor/samtools/bam_consensus.c, which is upstream
 * verbatim. Two of the things we need from it -- consensus_init() and the
 * cons_prob_recall / cons_prob_precise probability tables -- are `static`,
 * so they cannot be reached by linking against a separately compiled object.
 * Including the .c file textually puts this shim in the same translation unit
 * and makes them reachable, which is what lets the upstream file stay
 * byte-for-byte unmodified: syncing is then "drop in the new file and re-run
 * the oracle test", with no patch to rebase.
 *
 * The cost of that choice is that upstream's CLI entry point (main_consensus)
 * is compiled too, so this file must satisfy the handful of samtools-internal
 * symbols it references. Those are stubbed below; none are on any path we
 * call.
 *
 * Everything upstream needs is included before the textual include, so the
 * vendored code sees the environment it expects.
 */

#define _POSIX_C_SOURCE 200809L

#include "samtools_consensus.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ── the vendored model ───────────────────────────────────────────────── */

/* Upstream is warning-clean under its own build flags, not under ours
 * (-Wall -Wextra -Wpedantic). Silence what we cannot fix without editing it;
 * our own code keeps the strict flags. */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#pragma GCC diagnostic ignored "-Wsign-compare"
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#pragma GCC diagnostic ignored "-Wpedantic"
#ifdef __clang__
#pragma GCC diagnostic ignored "-Wunused-but-set-variable"
#endif

#include "samtools/bam_consensus.c"

#pragma GCC diagnostic pop

/* Defined after the include so upstream's attributed declarations in
 * sam_utils.h are seen first. */
void print_error(const char *subcommand, const char *format, ...);
void print_error_errno(const char *subcommand, const char *format, ...);

void print_error(const char *subcommand, const char *format, ...)
{
    va_list ap;
    fprintf(stderr, "bpp-seqs (vendored samtools %s): ", subcommand);
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    fputc('\n', stderr);
}

void print_error_errno(const char *subcommand, const char *format, ...)
{
    va_list ap;
    fprintf(stderr, "bpp-seqs (vendored samtools %s): ", subcommand);
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    fputc('\n', stderr);
}

/* Referenced only by main_consensus()'s option parsing, which we never call.
 * Present so the translation unit links; deliberately inert. */
int parse_sam_global_opt(int c, const char *optarg, const struct option *lopt,
                         sam_global_args *ga)
{
    (void)c; (void)optarg; (void)lopt; (void)ga;
    return -1;
}

void sam_global_args_free(sam_global_args *ga)
{
    (void)ga;
}


/* ── bpp-seqs API ─────────────────────────────────────────────────────── */

struct BppsCons {
    consensus_opts opts;
    char           del_char;
};

void bpps_cons_opts_init(BppsConsOpts *o)
{
    if (!o) return;
    memset(o, 0, sizeof(*o));
    /* Mirrors upstream's own defaults (see the consensus_opts initialiser in
     * bam_consensus.c), except ambig: BPP wants IUPAC heterozygotes. */
    o->ambig       = 1;
    o->cons_cutoff = 10;
    o->min_depth   = 1;
    o->p_het       = P_HET;      /* 1e-3, from bam_consensus.c */
    o->min_mqual   = 0;
    o->min_qual    = 0;
    o->del_char    = 'N';   /* '*' is not a base BPP can use */
}

BppsCons *bpps_cons_init(const BppsConsOpts *in)
{
    BppsCons *c = (BppsCons *)calloc(1, sizeof(*c));
    if (!c) return NULL;

    BppsConsOpts o;
    if (in) o = *in; else bpps_cons_opts_init(&o);

    consensus_opts *p = &c->opts;
    /* Upstream's defaults, field for field. Anything left zero here is a
     * field upstream also defaults to zero. */
    p->mode           = MODE_RECALL;
    p->use_qual       = 0;
    p->min_qual       = o.min_qual;
    p->adj_qual       = 1;
    p->use_mqual      = 1;
    p->scale_mqual    = 1.00;
    p->nm_adjust      = 1;
    p->nm_halo        = 50;
    p->sc_cost        = 60;
    p->low_mqual      = 1;
    p->high_mqual     = 60;
    p->min_depth      = o.min_depth;
    p->call_fract     = 0.75;
    p->het_fract      = 0.5;
    p->het_only       = 0;
    p->fmt            = FASTA;
    p->cons_cutoff    = o.cons_cutoff;
    p->ambig          = o.ambig;
    p->line_len       = 70;
    p->default_qual   = 10;
    p->all_bases      = 1;      /* we want every reference position back */
    p->show_del       = 1;      /* keep reference coordinates ... */
    p->show_ins       = 0;      /* ... so insertions must not shift them */
    p->mark_ins       = 0;
    p->incl_flags     = 0;
    p->excl_flags     = BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP;
    p->min_mqual      = o.min_mqual;
    p->P_het          = o.p_het;
    p->P_indel        = P_INDEL;
    p->het_scale      = P_HET_SCALE;
    p->homopoly_fix   = 0;
    p->homopoly_redux = 0.01;
    p->ref_qual       = 0;
    p->span           = 500000;
    p->fp_out         = NULL;   /* we never write through upstream's output */

    c->del_char = o.del_char ? o.del_char : 'N';

    set_qcal(&p->qcal, QCAL_FLAT);

    /* Build the probability tables exactly as main_consensus() does for
     * MODE_RECALL. These are file-static in the vendored source, reachable
     * only because of the textual include above. */
    consensus_init(p->P_het, p->P_indel, p->het_scale, 0.01,
                   &p->qcal, MODE_RECALL, &cons_prob_recall);

    return c;
}

void bpps_cons_free(BppsCons *c)
{
    free(c);
}

/* Per-region state handed to the pileup callbacks. */
typedef struct {
    ctx        base;        /* upstream's own context; callbacks expect it */
    char      *seq_out;
    hts_pos_t  beg, end;
    char       del_char;
} region_ctx;

/* Column callback: replaces upstream's basic_pileup(), which formats output
 * we do not want. The calling itself is unchanged -- consensus_base() is
 * upstream's, and carries the ambig / cutoff / min-depth logic with it. */
static int bpps_column(void *cd, samFile *fp, sam_hdr_t *h, pileup_t *p,
                       int depth, hts_pos_t pos, int nth, int is_insert)
{
    (void)fp; (void)h; (void)is_insert;
    region_ctx *r = (region_ctx *)cd;

    if (nth) return 0;                       /* inserted base: no ref coord */
    if (pos < r->beg || pos >= r->end) return 0;

    int base = 'N', qual = 0;
    if (consensus_base(r->base.opts, p, pos, depth, &base, &qual) < 0)
        return -1;

    /* A deletion occupies its reference position. samtools writes '*', which
     * is not a base BPP can read, so it is remapped -- but keeping the column
     * is what holds samples in register with each other. */
    if (base == '*') base = r->del_char;

    r->seq_out[pos - r->beg] = (char)base;
    return 0;
}

int bpps_cons_region(BppsCons *c,
                     samFile *fp, sam_hdr_t *h, hts_idx_t *idx,
                     const char *chrom, hts_pos_t beg, hts_pos_t end,
                     char *seq_out)
{
    if (!c || !fp || !h || !idx || !chrom || !seq_out || end <= beg) return -1;

    memset(seq_out, 'N', (size_t)(end - beg));

    char region[1024];
    snprintf(region, sizeof(region), "%s:%" PRIhts_pos "-%" PRIhts_pos,
             chrom, beg + 1, end);

    region_ctx r;
    memset(&r, 0, sizeof(r));
    r.base.opts = &c->opts;
    r.base.h    = h;
    r.seq_out   = seq_out;
    r.beg       = beg;
    r.end       = end;
    r.del_char  = c->del_char;

    r.base.iter = sam_itr_querys(idx, h, region);
    if (!r.base.iter) return -1;

    int rc = pileup_loop(fp, h, readaln2, nm_init, bpps_column, nm_free, &r);

    hts_itr_destroy(r.base.iter);
    return rc < 0 ? -1 : 0;
}

const char *bpps_cons_samtools_version(void)
{
    return "1.23.1";
}
