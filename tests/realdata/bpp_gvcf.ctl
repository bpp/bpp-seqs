seed = 42
seqfile  = gv_real.txt
Imapfile = gv_real.imap
jobname  = bpp_gvcf_out
speciesdelimitation = 0
speciestree = 0
species&tree = 4  AFR  EUR  EAS  AMR
                  2    2    2    2
                 (AFR, (EUR, (EAS, AMR)));
phase = 1 1 1 1
usedata = 1
nloci = 17
cleandata = 0
thetaprior = gamma 2 2000
tauprior   = gamma 2 1000
finetune = 1
print = 1 0 0 0
burnin   = 2000
sampfreq = 2
nsample  = 2000
