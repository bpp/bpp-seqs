seed = 42

seqfile  = result.txt
Imapfile = result.imap
jobname  = bpp_out

speciesdelimitation = 0
speciestree = 0

# Conventional human topology: AMR sister to EAS (peopled via Beringia).
species&tree = 4  AFR  EUR  EAS  AMR
                  2    2    2    2
                 (AFR, (EUR, (EAS, AMR)));

# Data are phased upstream by the NYGC pipeline.  bpp-seqs --phasing vcf
# emitted two haplotypes per individual, so BPP must NOT re-phase.
phase = 0 0 0 0

usedata = 1
nloci = 17
cleandata = 0

thetaprior = gamma 2 2000
tauprior   = gamma 2 1000

finetune = 1
print = 1 0 0 0

burnin   = 200
sampfreq = 2
nsample  = 200
