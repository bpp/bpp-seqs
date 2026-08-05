#!/usr/bin/env bash
# Provenance for examples/02-phylip-ascidian.
#
# The loci/ and imap.txt files here are already committed, so you do NOT need
# to run this to use the example. It records where the data came from and how
# to obtain the full 179-locus set.
#
# Source: Dryad doi:10.5061/dryad.3r2280gf7  (CC0-1.0)
#         Nydam et al. (2021), Scientific Reports 11,
#         doi:10.1038/s41598-021-87255-2
#
# File:   Nydam2021_T500_FinalAlignments.zip   (file id 2705416, ~728 KB)
#
# Dryad's API requires a bearer token for downloads, so either
#
#   (a) download from the dataset page in a browser -- no token needed:
#         https://doi.org/10.5061/dryad.3r2280gf7
#
#   (b) or, with a Dryad API token (see your Dryad account page):
#
#         TOKEN=$(cat ~/.config/dryad_token)
#         curl -L -H "Authorization: Bearer $TOKEN" \
#              https://datadryad.org/api/v2/files/2705416/download \
#              -o Nydam2021_T500_FinalAlignments.zip
#
# The committed loci/ holds the 20 alignments with the most complete taxon
# sampling (all 55 specimens), copied unmodified. To use all 179 instead:
#
#         unzip -q Nydam2021_T500_FinalAlignments.zip
#         cp Nydam2021_T500_FinalAlignments/*.phylip loci/
#         python3 make-imap.py > imap.txt
#
set -euo pipefail
sed -n '2,33p' "$0" | sed 's/^# \{0,1\}//'
