seed = 42
seqfile  = result_medium.txt
Imapfile = result_medium.imap
jobname  = bpp_out

speciesdelimitation = 0
speciestree = 0
species&tree = 4  AFR  EUR  EAS  AMR
                  4    4    4    4
                 (AFR, (EUR, (EAS, AMR)));
phase = 0 0 0 0
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
