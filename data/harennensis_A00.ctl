seed = -1

seqfile = harennensis.phy
Imapfile = harennensis_A00.imap
outfile = harennensis_A00.out
mcmcfile = harennensis_A00.mcmc

speciesdelimitation = 0
speciestree = 0

species&tree = 5 A B CD E H
                 4 4 4 1 5
                 ((A,((E,CD),B)),H);
phase = 1 1 1 1 1

usedata = 1
nloci = 6
model = HKY
alphaprior = 1 1 4

cleandata = 0

thetaprior = 3 0.04278 e
tauprior = 3 0.02983

locusrate = 1 0 0 1 dir
clock = 1
heredity = 0

finetune = 1: 0.01 0.02 0.03 0.04 0.05 0.01 0.01
print = 1 0 0 0 0
burnin = 100000
sampfreq = 100
nsample = 1000
threads = 4
