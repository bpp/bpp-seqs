/* json_writer.h — minimal streaming JSON writer for bpp-seqs.
 *
 * No DOM, no allocation per token. Emits to a FILE*.
 * Caller is responsible for matching obj_open/obj_close and arr_open/arr_close.
 * Comma placement is handled internally via a per-scope first-element flag.
 *
 * NaN/Inf doubles are emitted as JSON null.
 */
#ifndef BPP_SEQS_JSON_WRITER_H
#define BPP_SEQS_JSON_WRITER_H

#include <stdio.h>
#include <stdint.h>

#define JW_MAX_DEPTH 64

typedef struct {
    FILE *fp;
    int   indent;             /* spaces per level; 0 = compact */
    int   depth;
    /* For each open scope, 1 = next item is first (no leading comma). */
    unsigned char first[JW_MAX_DEPTH];
    /* For each open scope, 1 = inside an object (so a key must precede value). */
    unsigned char in_obj[JW_MAX_DEPTH];
    /* Flag set right after jw_key: the next value belongs to that key,
     * so it does not emit indent/comma logic for itself. */
    unsigned char after_key;
} JsonWriter;

void jw_init       (JsonWriter *w, FILE *fp, int indent);

void jw_obj_open   (JsonWriter *w);
void jw_obj_close  (JsonWriter *w);

void jw_arr_open   (JsonWriter *w);
void jw_arr_close  (JsonWriter *w);

/* Inside an object: emit "key": (no comma logic — that ran with the value). */
void jw_key        (JsonWriter *w, const char *key);

void jw_str        (JsonWriter *w, const char *s);          /* NULL → null */
void jw_str_n      (JsonWriter *w, const char *s, size_t n);
void jw_int        (JsonWriter *w, long long v);
void jw_uint       (JsonWriter *w, unsigned long long v);
void jw_dbl        (JsonWriter *w, double v);
void jw_bool       (JsonWriter *w, int v);
void jw_null       (JsonWriter *w);

/* Convenience helpers for object members: emit key then value. */
void jw_kv_str     (JsonWriter *w, const char *key, const char *val);
void jw_kv_int     (JsonWriter *w, const char *key, long long val);
void jw_kv_dbl     (JsonWriter *w, const char *key, double val);
void jw_kv_bool    (JsonWriter *w, const char *key, int val);
void jw_kv_null    (JsonWriter *w, const char *key);

/* Finish output: newline if pretty-printing, flush. */
void jw_finish     (JsonWriter *w);

#endif /* BPP_SEQS_JSON_WRITER_H */
