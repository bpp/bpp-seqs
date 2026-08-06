# Vendored samtools consensus caller

Everything in this directory **except `config.h`** is upstream samtools source,
byte-for-byte unmodified.

| | |
|---|---|
| **Upstream** | https://github.com/samtools/samtools |
| **Pinned version** | `1.23.1` |
| **Licence** | MIT (Genome Research Ltd.; portions derived from Gap4/5, relicensed by the copyright holder — see the header of `bam_consensus.c`) |
| **Author** | James Bonfield and the samtools contributors |

## Why these files are here

bpp-seqs' own base caller ranks raw allele counts against a fixed minor-allele
fraction: base qualities act only as a cutoff, mapping quality only as a filter,
and there is no error model. That calls a heterozygote from a single miscalled
base at modest depth.

`samtools consensus` implements the Gap5 Bayesian consensus model, which weighs
each observation by its base quality, caps it by mapping quality, adjusts for
local quality minima and NM density, and returns a phred-scaled confidence so
weakly supported positions become `N`. The model is not in htslib — it lives in
the samtools binary — so using it means either shelling out to `samtools`
(which would break the self-contained static binary the release CI builds) or
vendoring it.

Vendoring was chosen over reimplementing because it makes
`samtools consensus -A` a usable **test oracle**: output can be required to
match byte-for-byte, which a clean-room reimplementation could only approximate.
See tests 84+ in `tests/run_tests.sh`.

## Files

| File | Role |
|---|---|
| `bam_consensus.c` | the consensus model (`calculate_consensus_gap5`, `consensus_base`) |
| `bam_consensus_tab.h` | quality-calibration tables |
| `consensus_pileup.c/h` | samtools' pileup engine (depends only on `htslib/sam.h`) |
| `bam_plbuf.c/h` | pileup buffer |
| `samtools.h`, `sam_utils.h`, `sam_opts.h` | headers the above include |
| `config.h` | **ours**, not upstream — a shim for the autoconf header samtools generates |

The wrapper that drives them is `../samtools_consensus.c` / `.h`, which is
bpp-seqs code.

## How to sync to a newer samtools

```sh
TAG=1.24                      # whichever release you are moving to
cd src/vendor/samtools
for f in bam_consensus.c bam_consensus_tab.h consensus_pileup.c consensus_pileup.h \
         bam_plbuf.c bam_plbuf.h samtools.h sam_utils.h sam_opts.h; do
    curl -fsSL -o "$f" "https://raw.githubusercontent.com/samtools/samtools/$TAG/$f"
done
cd ../../..
make clean && make && ./tests/run_tests.sh
```

Then update the pinned version in all four places that record it — a test
(93) fails if they drift apart:

| Where | What |
|---|---|
| this README | the **Pinned version** row above |
| `../samtools_consensus.c` | `bpps_cons_samtools_version()` |
| `tests/run_tests.sh` | `PINNED_SAMTOOLS` |
| `.github/workflows/ci.yml` | `SAMTOOLS_VERSION` |

CI builds that exact release from source rather than taking whatever apt or
brew ships, and runs the suite with `BPP_SEQS_STRICT=1` so a skipped oracle
test fails the build. Without that, an unpinned samtools would make the oracle
tests skip and CI would report a clean run having verified nothing here.

**Do not edit these files.** The wrapper reaches upstream's `static` internals
(`consensus_init`, `cons_prob_recall`) by including `bam_consensus.c`
*textually*, precisely so that no patch has to be maintained across syncs.

### What can break, and how you will know

`calculate_consensus_gap5`, `consensus_base` and `consensus_init` are samtools
*internals* with no API-stability promise, so a future release may rename them
or change their signatures. When that happens you will get a **compile error**
in `../samtools_consensus.c`, or a **failing oracle test** — never a silent
change in results. Fix the wrapper, not these files.

The oracle tests skip themselves when the installed `samtools` version differs
from the pinned one, since only the matching version is guaranteed to agree.
