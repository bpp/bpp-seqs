#!/usr/bin/env bash
# Provenance for examples/01-fasta-anastrepha.
#
# The loci/ and imap.txt files in this directory are already committed, so you
# do NOT need to run this to use the example. It is here so the subset can be
# regenerated or widened from the published source.
#
# Source: Dryad doi:10.5061/dryad.7d7wm3864  (CC0-1.0)
#         Norrbom et al. (2025), Systematic Entomology, doi:10.1111/syen.70003
#
# Two files are needed:
#     Individual_Locus_Alignments.zip      (file id 4331360, ~22 MB)
#     AnastrephaNames_LabCodes_*.csv       (file id 4331369, ~47 KB)
#
# Dryad's API requires a bearer token for downloads, so either
#
#   (a) download both from the dataset page in a browser -- no token needed:
#         https://doi.org/10.5061/dryad.7d7wm3864
#
#   (b) or, with a Dryad API token (see
#       https://datadryad.org/api/v2/docs and your Dryad account page):
#
#         TOKEN=$(cat ~/.config/dryad_token)
#         curl -L -H "Authorization: Bearer $TOKEN" \
#              https://datadryad.org/api/v2/files/4331360/download \
#              -o Individual_Locus_Alignments.zip
#         curl -L -H "Authorization: Bearer $TOKEN" \
#              https://datadryad.org/api/v2/files/4331369/download \
#              -o AnastrephaNames_LabCodes.csv
#
# Then unzip and rebuild this example's subset:
#
#         unzip -q Individual_Locus_Alignments.zip
#         python3 make-subset.py DS1_aln AnastrephaNames_LabCodes.csv
#
set -euo pipefail
sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'
