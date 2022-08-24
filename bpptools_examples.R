### TOOLS
library(ape)
library(phangorn)
library(invgamma)
library(coda)
source("func/bpptools.R")
source("func/readNexus.R")
source("func/find_mcct.R")
source("func/find_spdelim.R")


### EXAMPLES
### Preparing BPP input files (sequence and imap file)
# from partitioned nexus file and tab-delimited info file
seq <- readNexus("data/harennensis.nex")
info <- read.delim("data/harennensis_info.txt")
write_bpp_input(seq=seq, info=info, seqfile="data/harennensis.phy", imapfile="data/harennensis_A11.imap")


### Converting BPP output (.mcmc) file into formats readable by Tracer and FigTree
mcmc <- read_bpp_output("results/harennensis_A11.mcmc")
export_bpp_output(mcmc, trees="results/harennensis_A11.trees", log="results/harennensis_A11.log", model="A11")


### Finding maximum credibility tree & maximum credibility species delimitation;
# exporting the tree in FigTree-readable format and the delimitation as tab-delimited file;
# plotting the tree annotated with two sets of posterior probabilities (PPs of clades and species)
trees <- ape::read.tree("results/harennensis_A11.trees")
mctree <- find_mcct(trees, file="results/FigTree_harennensis_A11_mcct.tre", return=TRUE, figtree=TRUE)
mcspdelim <- find_spdelim(trees, collapse=0, sep="-", file="results/harennensis_A11_specdelim.txt")
plot_delim_tree(mctree, mcspdelim, edge.width=10, clade.cex=6, species.cex=6, text.cex=1, tip.cex=1.5, offset=0.07)


### Export and plotting of the multispecies coalescent tree from A00 analysis (plot labels - gdi or theta)
topology <- read_bpp_tree("data/harennensis_A00.ctl")
mcmc <- read_bpp_output("results/harennensis_A00.mcmc", tree=topology)
export_msc_tree(target=topology, trees=mcmc$trees, params=mcmc$params, file="results/harennensis_A00_msc.tre", figtree=TRUE)

msct <- read_msc_tree("results/harennensis_A00_msc.tre")
tiplab <- c(H="harennensis", A="triton_A", B="triton_B", "CD"="triton_CD", E="triton_E")
plot_msc_tree(msct, tip.label=tiplab, labels="gdi", digits=2, lab.font=2)
plot_msc_tree(msct, tip.label=tiplab, labels="theta", digits=4)


### Calculation of genealogical divergence index for extant & ancestral species (= species tree branches)
posterior <- read.delim("results/harennensis_A00.log")
topology <- read_bpp_tree("data/harennensis_A00.ctl")
gdis <- gdi(params=posterior, tree=topology)
gdis$means


### Confronting prior and posterior of root divergence time
prior <- read_bpp_prior("data/harennensis_A00.ctl", "tau")
prior <- as.numeric(prior[,1:2])
posterior <- read.delim("results/harennensis_A00.log")
posterior <- posterior[,grep("^tau", colnames(posterior))]
plot_pdensity(prior=prior[1:2], posterior=posterior[,1], xlab="Tau", col=c(8,1), lty=1, lwd=4, cex=1.5, points=FALSE, legend=TRUE)


### Checking for convergence of BPP runs (using established rules of thumb)
# and pooling of the posterior samples
topology <- read_bpp_tree("data/harennensis_A00.ctl")
mcmc1 <- read_bpp_output("results/harennensis_A00_run1.mcmc", tree=topology)
mcmc2 <- read_bpp_output("results/harennensis_A00_run2.mcmc", tree=topology)
runs <- list(mcmc1, mcmc2)
params <- lapply(lapply(runs, "[[", "params"), coda::as.mcmc)
psrf <- coda::gelman.diag(params, autoburnin=FALSE, multivariate=FALSE)$psrf
ess <- coda::effectiveSize(params)
(psrf[,"Point est."] - 1)  < 0.025
ess >= 200

mcmc <- combine_bpp_runs(runs)
export_bpp_output(mcmc, trees="results/harennensis_A00_pooled.trees", log="results/harennensis_A00_pooled.log", model="A00")

