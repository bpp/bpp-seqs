#include "json_writer.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static void jw_emit_indent(JsonWriter *w)
{
    if (w->indent <= 0) return;
    fputc('\n', w->fp);
    for (int i = 0; i < w->depth; i++) {
        for (int j = 0; j < w->indent; j++) fputc(' ', w->fp);
    }
}

/* Called at the start of every value or key emission. Emits comma between
 * siblings and indentation. Does nothing if the writer is in "after_key"
 * state (the value immediately follows the colon written by jw_key). */
static void jw_prelude(JsonWriter *w)
{
    if (w->after_key) {
        w->after_key = 0;
        return;
    }
    if (w->depth > 0) {
        if (w->first[w->depth - 1]) {
            w->first[w->depth - 1] = 0;
        } else {
            fputc(',', w->fp);
        }
        jw_emit_indent(w);
    }
}

void jw_init(JsonWriter *w, FILE *fp, int indent)
{
    w->fp = fp;
    w->indent = indent;
    w->depth = 0;
    w->after_key = 0;
    memset(w->first, 0, sizeof(w->first));
    memset(w->in_obj, 0, sizeof(w->in_obj));
}

static void jw_scope_open(JsonWriter *w, char ch, int is_obj)
{
    jw_prelude(w);
    fputc(ch, w->fp);
    if (w->depth < JW_MAX_DEPTH) {
        w->first[w->depth] = 1;
        w->in_obj[w->depth] = (unsigned char)is_obj;
    }
    w->depth++;
}

static void jw_scope_close(JsonWriter *w, char ch)
{
    int was_empty = (w->depth > 0 && w->first[w->depth - 1]);
    w->depth--;
    if (!was_empty) jw_emit_indent(w);
    fputc(ch, w->fp);
}

void jw_obj_open (JsonWriter *w) { jw_scope_open (w, '{', 1); }
void jw_obj_close(JsonWriter *w) { jw_scope_close(w, '}'); }
void jw_arr_open (JsonWriter *w) { jw_scope_open (w, '[', 0); }
void jw_arr_close(JsonWriter *w) { jw_scope_close(w, ']'); }

static void jw_write_escaped(JsonWriter *w, const char *s, size_t n)
{
    fputc('"', w->fp);
    for (size_t i = 0; i < n; i++) {
        unsigned char c = (unsigned char)s[i];
        switch (c) {
            case '"':  fputs("\\\"", w->fp); break;
            case '\\': fputs("\\\\", w->fp); break;
            case '\b': fputs("\\b",  w->fp); break;
            case '\f': fputs("\\f",  w->fp); break;
            case '\n': fputs("\\n",  w->fp); break;
            case '\r': fputs("\\r",  w->fp); break;
            case '\t': fputs("\\t",  w->fp); break;
            default:
                if (c < 0x20) fprintf(w->fp, "\\u%04x", c);
                else          fputc((int)c, w->fp);
        }
    }
    fputc('"', w->fp);
}

void jw_key(JsonWriter *w, const char *key)
{
    jw_prelude(w);
    jw_write_escaped(w, key, strlen(key));
    fputc(':', w->fp);
    if (w->indent > 0) fputc(' ', w->fp);
    w->after_key = 1;
}

void jw_str(JsonWriter *w, const char *s)
{
    if (!s) { jw_null(w); return; }
    jw_prelude(w);
    jw_write_escaped(w, s, strlen(s));
}

void jw_str_n(JsonWriter *w, const char *s, size_t n)
{
    if (!s) { jw_null(w); return; }
    jw_prelude(w);
    jw_write_escaped(w, s, n);
}

void jw_int(JsonWriter *w, long long v)
{
    jw_prelude(w);
    fprintf(w->fp, "%lld", v);
}

void jw_uint(JsonWriter *w, unsigned long long v)
{
    jw_prelude(w);
    fprintf(w->fp, "%llu", v);
}

void jw_dbl(JsonWriter *w, double v)
{
    jw_prelude(w);
    if (isnan(v) || isinf(v)) {
        fputs("null", w->fp);
    } else {
        /* %.6g is plenty for the fractions/depths/lengths we emit. */
        fprintf(w->fp, "%.6g", v);
    }
}

void jw_bool(JsonWriter *w, int v)
{
    jw_prelude(w);
    fputs(v ? "true" : "false", w->fp);
}

void jw_null(JsonWriter *w)
{
    jw_prelude(w);
    fputs("null", w->fp);
}

void jw_kv_str (JsonWriter *w, const char *key, const char *val) { jw_key(w, key); jw_str (w, val); }
void jw_kv_int (JsonWriter *w, const char *key, long long val)   { jw_key(w, key); jw_int (w, val); }
void jw_kv_dbl (JsonWriter *w, const char *key, double val)      { jw_key(w, key); jw_dbl (w, val); }
void jw_kv_bool(JsonWriter *w, const char *key, int val)         { jw_key(w, key); jw_bool(w, val); }
void jw_kv_null(JsonWriter *w, const char *key)                  { jw_key(w, key); jw_null(w); }

void jw_finish(JsonWriter *w)
{
    if (w->indent > 0) fputc('\n', w->fp);
    fflush(w->fp);
}
