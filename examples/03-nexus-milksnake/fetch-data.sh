#!/usr/bin/env bash
# Provenance for examples/03-nexus-milksnake.
#
# The loci/ and imap.txt files here are already committed, so you do NOT need
# to run this to use the example. It records where the data came from.
#
# Source: Dryad doi:10.5061/dryad.420h7  (CC0-1.0)
#         Ruane et al. (2014), Systematic Biology 63:231-250,
#         doi:10.1093/sysbio/syt099
#
# File:   Lamp_Data.zip   (file id 37134, ~130 KB)
#
# Dryad's API requires a bearer token for downloads, so either
#
#   (a) download from the dataset page in a browser -- no token needed:
#         https://doi.org/10.5061/dryad.420h7
#
#   (b) or, with a Dryad API token (see your Dryad account page):
#
#         TOKEN=$(cat ~/.config/dryad_token)
#         curl -L -H "Authorization: Bearer $TOKEN" \
#              https://datadryad.org/api/v2/files/37134/download -o Lamp_Data.zip
#
# loci/ holds the 11 nuclear alignments plus CYTB_124_Set.nex, unmodified;
# only the filenames were tidied (the "COMPLETE " prefix and a stray ".txt"
# were dropped, and CYTB_124_Set.nex became CYTB.nex). The archive also holds
# Cytb_329_Set.nex, a larger 329-taxon cytb sampling that does not line up with
# the 124 individuals in the nuclear loci, so it is not used here.
#
#         unzip -q Lamp_Data.zip
#         python3 make-imap.py > imap.txt
#
set -euo pipefail
sed -n '2,33p' "$0" | sed 's/^# \{0,1\}//'
