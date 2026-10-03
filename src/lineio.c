#include "lineio.h"

#include <stdlib.h>
#include <string.h>

long gz_getline(gzFile gz, char **buf, size_t *cap)
{
    if (!*buf || *cap < 2) {
        *cap = 1 << 16;
        *buf = (char *)realloc(*buf, *cap);
    }
    size_t len = 0;
    for (;;) {
        if (gzgets(gz, *buf + len, (int)(*cap - len)) == NULL)
            return len ? (long)len : -1;
        len += strlen(*buf + len);
        if (len > 0 && (*buf)[len - 1] == '\n') return (long)len;
        if (len + 1 < *cap) return (long)len;          /* last line, no newline */
        *cap *= 2;
        *buf = (char *)realloc(*buf, *cap);
    }
}
