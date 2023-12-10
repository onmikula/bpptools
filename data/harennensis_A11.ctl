seed = -1

seqfile = harennensis.phy
Imapfile = harennensis_A11.imap
outfile = harennensis_A11.out
mcmcfile = harennensis_A11.mcmc

speciesdelimitation = 1 0 2
speciestree = 1
speciesmodelprior = 1

species&tree = 6 A B C D E H
                 4 4 2 2 1 5
                 ((((E,(D,C)),B),A),H);
phase = 1 1 1 1 1 1

usedata = 1
nloci = 6
model = HKY
alphaprior = 1 1 4

cleandata = 0

thetaprior = 3 0.05734 e
tauprior = 3 0.02867

locusrate = 1 0 0 1 dir
clock = 1
heredity = 0

finetune = 1: 0.01 0.02 0.03 0.04 0.05 0.01 0.01
print = 1 0 0 0 0
burnin = 100000
sampfreq = 100
nsample = 1000
threads = 4
