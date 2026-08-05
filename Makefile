# Makefile for bpp-seqs
#
# Requires htslib (https://github.com/samtools/htslib).
#
# Detection order for htslib:
#   1. HTSLIB_A=/path/to/libhts.a  -- STATIC link (self-contained binary; used
#                                     by the release CI). Also set HTSLIB_INC.
#   2. pkg-config --libs htslib
#   3. HTSLIB=/path prefix override (make HTSLIB=/opt/htslib)
#   4. Bare -lhts
#
# Static build recipe (produces a binary with no libhts dependency):
#   ( cd htslib && ./configure --disable-libcurl --disable-plugins && make )
#   make HTSLIB_A=$PWD/htslib/libhts.a HTSLIB_INC=$PWD/htslib
#   ldd bpp-seqs   # confirms libhts is absent
#
# EXTRA_CFLAGS / EXTRA_LDFLAGS are appended after the defaults (e.g. for
# cross-compilation: make EXTRA_CFLAGS=-arch\ x86_64 EXTRA_LDFLAGS=-arch\ x86_64).
#
# Targets:
#   make            release build
#   make debug      debug build (-g, no -O2, ASan/UBSan)
#   make clean      remove object files and binary
#   make test       run tests/run_tests.sh
#   make install    install to $(PREFIX)/bin (default PREFIX=/usr/local)

CC      := cc
CSTD    := -std=c11
WARN    := -Wall -Wextra -Wpedantic
CFLAGS  := -O2 $(CSTD) $(WARN) -D_GNU_SOURCE
LDFLAGS :=
PREFIX  ?= /usr/local

# ── htslib detection ────────────────────────────────────────────────────────

# Compression/threading libs the static libhts.a needs pulled in after it.
# Override if your htslib was configured with libdeflate or libcurl.
HTSLIB_STATIC_DEPS ?= -lz -lbz2 -llzma -lpthread -lm

ifdef HTSLIB_A
  HTSLIB_CFLAGS := -I$(HTSLIB_INC)
  HTSLIB_LIBS   := $(HTSLIB_A) $(HTSLIB_STATIC_DEPS)
else ifdef HTSLIB
  HTSLIB_CFLAGS := -I$(HTSLIB)/include
  HTSLIB_LIBS   := -L$(HTSLIB)/lib -lhts -Wl,-rpath,$(HTSLIB)/lib
else
  HTSLIB_CFLAGS := $(shell pkg-config --cflags htslib 2>/dev/null)
  HTSLIB_LIBS   := $(shell pkg-config --libs   htslib 2>/dev/null || echo "-lhts")
endif

INCLUDES := -Isrc -Isrc/bam2bpp -Isrc/vendor -Isrc/vendor/samtools $(HTSLIB_CFLAGS)
CFLAGS   += $(INCLUDES) $(EXTRA_CFLAGS)
LDFLAGS  += $(EXTRA_LDFLAGS)

# Same as CFLAGS but without the strict warning set, for vendored third-party
# code we deliberately keep unmodified.
CFLAGS_VENDOR := $(filter-out $(WARN),$(CFLAGS))
LDLIBS   := $(HTSLIB_LIBS) -lz -lm -lpthread

# ── sources ────────────────────────────────────────────────────────────────

SRCS := \
  src/main.c \
  src/inspect.c \
  src/nexus.c \
  src/cross_validate.c \
  src/workflow.c \
  src/output.c \
  src/json_writer.c \
  src/sanity.c \
  src/bpp_parser.c \
  src/cmd_extract.c \
  src/cmd_windows.c \
  src/bam2bpp/pileup.c \
  src/bam2bpp/genotype.c \
  src/bam2bpp/locus.c \
  src/bam2bpp/writer.c \
  src/bam2bpp/vcf_phase.c \
  src/converters/aln_writer.c \
  src/converters/fasta2bpp.c \
  src/converters/phylip2bpp.c \
  src/converters/nexus2bpp.c \
  src/converters/gvcf2bpp.c \
  src/vendor/samtools_consensus.c \
  src/vendor/samtools/consensus_pileup.c \
  src/vendor/samtools/bam_plbuf.c

OBJS := $(SRCS:.c=.o)
BIN  := bpp-seqs

# ── targets ────────────────────────────────────────────────────────────────

.PHONY: all debug clean test install

all: $(BIN)

debug: CFLAGS := -O0 -g $(CSTD) $(WARN) -D_GNU_SOURCE $(INCLUDES) -fsanitize=address,undefined
debug: LDFLAGS += -fsanitize=address,undefined
debug: $(BIN)

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) $(LDLIBS)

# -MMD -MP writes a .d file of header dependencies beside each object, so a
# header edit rebuilds everything that includes it. Without this a change to
# a struct in a header relinks stale objects compiled against the old layout,
# producing a binary that misbehaves in ways no source inspection explains.
# Vendored samtools sources are upstream verbatim and are not warning-clean
# under our flags (-Wall -Wextra -Wpedantic). Compile them without those so
# a sync never has to touch the code to keep the build quiet; our own sources
# stay strict. Must precede the generic %.o rule: make takes the first
# matching pattern, not the most specific. See src/vendor/samtools/README.md.
src/vendor/samtools/%.o: src/vendor/samtools/%.c
	$(CC) $(CFLAGS_VENDOR) -MMD -MP -c -o $@ $<

%.o: %.c
	$(CC) $(CFLAGS) -MMD -MP -c -o $@ $<

-include $(OBJS:.o=.d)

test: $(BIN)
	./tests/run_tests.sh

install: $(BIN)
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 755 $(BIN) $(DESTDIR)$(PREFIX)/bin/$(BIN)

clean:
	rm -f $(OBJS) $(OBJS:.o=.d) $(BIN)
