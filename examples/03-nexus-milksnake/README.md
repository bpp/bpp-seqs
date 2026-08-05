# 03 — milksnakes: per-locus NEXUS → BPP (`nexus2bpp`)

Real multilocus data from the milksnake (*Lampropeltis*) species-delimitation
study of Ruane et al. (2014) — the same empirical dataset later reused in the
BPP species-delimitation literature.

| | |
|---|---|
| **Workflow** | `nexus2bpp` |
| **Data** | 12 loci: 11 nuclear + mtDNA cytochrome *b* |
| **Sequences** | 248 (124 individuals × 2 alleles; cytb is haploid, 124) |
| **Groups** | 23 — named species, five geographic *triangulum* lineages, two outgroup genera |
| **Source** | Dryad [doi:10.5061/dryad.420h7](https://doi.org/10.5061/dryad.420h7) |
| **Paper** | Ruane S., Bryson R.W. Jr, Pyron R.A., Burbrink F.T. (2014) Coalescent species delimitation in milksnakes (genus *Lampropeltis*) and impacts on phylogenetic comparative analyses. *Systematic Biology* 63:231–250. [doi:10.1093/sysbio/syt099](https://doi.org/10.1093/sysbio/syt099) |
| **License** | CC0-1.0 (public domain dedication) |

The NEXUS files in `loci/` are the published files, byte-for-byte, renamed only
to drop the `COMPLETE ` prefix and the stray `.txt` in one filename. Nothing in
the sequence data has been altered.

## Run it

```sh
# 1. Inspect: what does conversion still need?
../../bpp-seqs --dry-run loci/*.nex --imap imap.txt

# 2. Convert to BPP format.
../../bpp-seqs --out milksnake loci/*.nex --imap imap.txt
```

All 12 loci pass QC, giving `milksnake.txt` (the BPP sequence file),
`milksnake.imap`, `milksnake.stats.tsv`, and `milksnake.loci.tsv`.

## Where the Imap comes from

`imap.txt` is not hand-written — it is derived from the published taxon labels
by `make-imap.py`, which you can re-run at any time:

```sh
python3 make-imap.py > imap.txt
```

Each label carries the authors' own lineage assignment as its leading token:

```
getula_FTB966_NJ_Cumberland        allele 1 of voucher FTB966  → getula
getula_a_FTB966_NJ_Cumberland      allele 2 of the same snake  → getula
SA_FHGO3090_Ecuador_Pinchincha     a South American lineage    → SA
```

The `_a_` infix marks the second allele of a diploid individual, so both
sequences of an individual map to the same group. `MX`, `CA`, `SA`, `Tam` and
`WEST` are the geographic lineages within the *L. triangulum* complex that the
paper set out to test; `Arizona_elegans` and `Cemophora_coccinea` are outgroups.

## Phasing

These sequences are already resolved alleles, so no `--phasing` option applies —
that flag only affects BAM/CRAM and gVCF input. Use `phase = 0` for every
species in the BPP control file.
