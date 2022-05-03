### TOOLS
library(ape)
library(phangorn)
library(invgamma)
source("func/bpptools.R")
source("func/readNexus.R")
source("func/concatenate.R")
source("func/find_mcct.R")
source("func/find_spdelim.R")
write.delim <- function(x, file) write.table(x, file, row.names=FALSE, quote=FALSE, sep="\t")



### FILE NAMES
# these names are case-specific
# others (defined in inputs / outputs sections) are generated using the case-specific 'stem' keyword
# 'alignments' are assumed to be a partitioned nexus file
# 'infofile' is assumed to be a tab-delimited file with individual identifiers in the first column and (candidate) species classifiers in the second
# requires data files containing sequences and information about mapping (=classification) of sequences to individuals and individuals to (candidate) species

alignments <- "data/harennensis.nex"
infofile <- "data/harennensis_info.txt"
stem <- "harennensis"



### A11 INPUTS
# file names
seqfile <- paste0("data/", stem, ".phy")
imapfile <- paste0("data/", stem, "_A11.imap")
ctlfile <- paste0("data/", stem, "_A11.ctl")
mcmcfile <- paste0(stem, "_A11.mcmc")
outfile <- paste0(stem, "_A11.out")

# importing information about indiiduals
info <- read.delim(infofile)
# importing list of alignments
seq <- readNexus(alignments)
# input data files, for 
write_bpp_input(seq=seq, info=info, seqfile=seqfile, imapfile=imapfile)

# MSC priors
# a, b = alpha, beta parameterse of inverse gamma priors of tau0 (root time of species tree) and mean theta (population size of the species)
# the purpose of the below code is to obtain reasonable values of these parameters, which should make priors sufficiently diffuse (no too much informative), but centered around values which could approximately fit the data
# alpha = 3 is recommended, beta can be calculated from alpha and a supposed mean of the distribution 
# The mean of tau0 prior can be about half of the maximum observed distance between candidate species (if we believe there is more than one actual species)
# The mean of theta prior can be about mean of genetics distances within pre-defined OTUs (which may or may not correspond to candidate species specified in imap file). It can be also based on previous analyses or on external information about population size and mutation rate.
# Here, we provisionally delimited operational taxonomic units (OTUs) to guide our choice of tau0 prior. The delimitaion was made by function 'estimate_otus' which internally uses an ad hoc created average linkege (a.k.a. UPGMA) tree and looks for a cross section through the tree with the longest average branch length.
# In any case, different approaches are possible and should be tailored to the data set in hand. Sometimes, it may be useful to export the matrix of average distances (intraspecific on diagonal, interspecific off diagonal) and work on it directly (example included).
# All calculations are based on concatenated sequences of all loci (obtained by function 'concatenate')
# Alternatively, gamma prior can be used with shape (alpha) and rate (beta) parameters and the mean m = a / b.

# reading in previously created files with individual IDs in pplace of sequences names
seq <- read_bpp_seq(seqfile, indnames=TRUE)
imap <- read_bpp_imap(imapfile)
# concatenation of multi-locus sequences 
concat <- concatenate(seq)
# defining otus
otus <- estimate_otus(concat, model="K80")
# estimating parameters 
params <- estimate_params(concat, otus, model="K80")
#ave <- avedist(concat, otus, model="K80")
#m_tau <- max(ave[lower.tri(ave)], na.rm=TRUE) / 2; m_tau
#m_theta <- median(diag(ave), na.rm=TRUE); m_theta

# calculating beta parameters of inverse gamma priors of tau0 & theta
a <- 3
m_tau <- params["tau"]
b_tau <- round((a - 1) * m_tau, 5)
m_theta <- params["theta"]
b_theta <- round((a - 1) * m_theta, 5)

# plot of the prior probability densities
legend <- expression(legend("topright",legend=c("Tau","Theta"), col=c(1,8), lwd=2))
plot_pdensity(prior=list(c(a, b_tau),c(a, b_theta)), xlab="Tau / Theta", col=c(1,8), lwd=2, cex=1.5, expr=legend)

# parsing control file
tree <- get_starting_tree("nj", seq=concat, imap=imap)
species <- get_maxsamplesize(seqfile, imap)

parse_bpp_ctlfile(con=ctlfile, spdelim=TRUE, sptree=TRUE, seqfile=sub("^.*\\/","",seqfile), imapfile=sub("^.*\\/","",imapfile), outfile=outfile, mcmcfile=mcmcfile, species=species, tree=tree, nloci=length(seq), thetaprior=c(a, b_theta), tauprior=c(a, b_tau), thetadistr="invgamma", taudistr="invgamma", theta=TRUE, speciesmodelprior=1, phase=TRUE, model="HKY", gamma=TRUE, locusrate=c(1, 0, 0, 1, "dir"), clock=1, cleandata=FALSE, finetune=1, print=c(1,0,0,0,0), burnin=100000, sampfreq=100, nsample=1000, threads=4)



### A11 OUTPUTS
# file names
mcmcfile <- paste0("results/", stem, "_A11.mcmc")
logfile <- paste0("results/", stem, "_A11.log")
treefile <- paste0("results/", stem, "_A11.trees")
mcctfile <- paste0("results/", stem, "_A11_mcct.tre")
delimfile <- paste0("results/", stem, "_A11_specdelim.txt")
ctlfile <- paste0("data/", stem, "_A11.ctl")

# read in posterior sample (from .mcmc file) and export it in format readable by, e.g., FigTree or Tracer
mcmc <- read_bpp_output(mcmcfile)
export_bpp_output(mcmc, trees=treefile, log=logfile, model="A11")

# maximum clade credibility tree (also possible to get it in FigTree-readable format by setting 'figtree=TRUE')
mcct <- find_mcct(mcmc$trees, mcctfile, return=TRUE)

# maximum credibility species delimitation
# speciesDA extracts species delimitations from sampled species trees
# speciesPP calculates posterior probabilities from sampled delimitations
# speciesMC selects maximum credibility species delimitation
# the three coomands are also wrapped into 'find_spdelim'
#delim <- speciesDA(mcmc$trees, collapse=0)
#delim_pp <- speciesPP(delim, sep="-")
#delim_sp <- speciesMC(delim_pp, sep="-")
delim_sp <- find_spdelim(mcmc$trees, collapse=0, sep="-"); delim_sp
write.delim(delim_sp, delimfile)

# mapping species onto MCC tree
mcctfile <- paste0("results/", stem, "_A11_mcct.tre")
delimfile <- paste0("results/", stem, "_A11_specdelim.txt")
tree <- ape::read.tree(mcctfile)
delim <- read.delim(delimfile)

# plot MCC tree annotated with posterior probabilities of clades and species
plot_delim_tree(tree, delim)

# plot of posterior probability density of tau0 and mean theta
posterior <- read.delim(logfile)
plot_pdensity(posterior=posterior[,c("tau.0","mtheta")], xlab="Tau / Theta", col=c(1,8), lty=1, lwd=3, cex=1.5, legend=c("Tau","Theta"))

# plot of prior vs. posterior probability density of tau0 and mean theta, respectively.
priors <- read_bpp_prior(ctlfile)
plot_pdensity(prior=priors["tau",1:2], posterior=posterior[,"tau.0"], xlab="Tau", col=1, lty=c(3,1), lwd=3, cex=1.5, points=FALSE, legend=TRUE)
plot_pdensity(prior=priors["theta",1:2], posterior=posterior[,"mtheta"], xlab="Theta", col=1, lty=c(3,1), lwd=3, cex=1.5, points=FALSE, legend=TRUE)



### A00 INPUTS
# file names
seqfile <- paste0("data/", stem, ".phy")
delimfile <- paste0("results/", stem, "_A11_specdelim.txt")
imapA11 <- paste0("data/", stem, "_A11.imap")
imapA00 <- paste0("data/", stem, "_A00.imap")
mcctfile <- paste0("results/", stem, "_A11_mcct.tre")
ctlfile <- paste0("data/", stem, "_A00.ctl")
mcmcfile <- paste0(stem, "_A00.mcmc")
outfile <- paste0(stem, "_A00.out")


# input data can be either written de novo based on source files 'alignments' & 'infofile' (cf. 'A11 INPUTS') or rewritten based of files from A11 analysis

# input data - written de novo
info <- read.delim(infofile)
seq <- readNexus(alignments)
write_bpp_input(seq=seq, info=info, seqfile=seqfile, imapfile=imapA00)

# input data - rewritten
species <- read.delim(delimfile)$Species
delim <- write_bpp_imap(delim=species, return=TRUE, sep="-")
mcct <- get_collapsed_tree(ape::read.tree(mcctfile), delim)
imap <- write_bpp_imap(imap=imapA11, delim=species, file=imapA00, return=TRUE, sep="-")

# MSC priors
seq <- read_bpp_seq(seqfile, indnames=TRUE)
concat <- concatenate(seq)
params <- estimate_params(concat, imap, model="K80")
a <- 3
m_tau <- params["tau"]
b_tau <- round((a - 1) * m_tau, 5)
m_theta <- params["theta"]
b_theta <- round((a - 1) * m_theta, 5)

# parsing control file
tree <- get_topology(mcct)
species <- get_maxsamplesize(seqfile, imap)

parse_bpp_ctlfile(con=ctlfile, spdelim=FALSE, sptree=FALSE, seqfile=sub("^.*\\/","",seqfile), imapfile=sub("^.*\\/","",imapA00), outfile=outfile, mcmcfile=mcmcfile, species=species, tree=tree, nloci=length(seq), thetaprior=c(a, b_theta), tauprior=c(a, b_tau), thetadistr="invgamma", taudistr="invgamma", theta=TRUE, phase=TRUE, model="HKY", gamma=TRUE, locusrate=c(1, 0, 0, 1, "dir"), clock=1, cleandata=FALSE, finetune=1, print=c(1,0,0,0,0), burnin=100000, sampfreq=100, nsample=1000, threads=4)



### A00 OUTPUTS
# file names
ctlfile <- paste0("data/", stem, "_A00.ctl")
mcmcfile <- paste0("results/", stem, "_A00.mcmc")
logfile <- paste0("results/", stem, "_A00.log")
treefile <- paste0("results/", stem, "_A00.trees")
msctree <- paste0("results/", stem, "_A00_msc.tre")
mscsample <- paste0("results/", stem, "_A00_psample.tre")

# read in posterior sample (from .mcmc file) and export it in format readable by, e.g., FigTree or Tracer
topology <- read_bpp_tree(ctlfile)
mcmc <- read_bpp_output(mcmcfile, tree=topology)
export_bpp_output(mcmc, trees=treefile, log=logfile, model="A00")

# multispecies coalescent species tree
export_msc_tree(target=topology, trees=mcmc$trees, params=mcmc$params, file=msctree, treesfile=mscsample, figtree=TRUE)
msct <- read_msc_tree(msctree)
tiplab <- c(H="harennensis", A="triton_A", B="triton_B", "CD"="triton_CD", E="triton_E")

plot_msc_tree(msct, tip.label=tiplab, labels="gdi", digits=2, lab.font=2)
plot_msc_tree(msct, tip.label=tiplab, labels="theta")

# genealogical divergence index
params <- read.delim(logfile)
topology <- read_bpp_tree(ctlfile)
gdis <- gdi(params, topology)
gdis$means
apply(gdis$distributions, 2, HPD, p=0.95)
leg <- expression(legend("topleft", pch=NA, legend="triton CD", bty="n", cex=1.25, text.font=2))
plot_gdi(gdis, "CD", lwd=3, cex=2, expr=leg)

