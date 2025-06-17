#' BPP control file.
#' 
#' @description
#' Writes BPP control file.
#' 
#' @param con The file name or `stdout()` if you want to write the output on the screen.
#' @param spdelim Logical, whether to estimate species delimitation.
#' @param sptree Logical, whether to estimate species tree.
#' @param seqfile Name of the file with sequences (`seqfile`).
#' @param imapfile Name of the files with mapping of individuals to species (`imapfile`).
#' @param jobname Character string, prefix all output filenames begin with,
#'   including names of files with the program standard output (`outfile`) and MCMC trace (`mcmcfile`).
#' @param species A named numeric vector, numbers are maximum numbers of sequences for particular species,
#'   names are species labels as they appear in the species tree.
#' @param tree The species tree starting topology in the newick format, as a `phylo` object
#'   or a name of file containing the tree.
#' @param nloci Numeric, how many loci to use (in the same order as in seqfile).
#' @param thetaprior,tauprior Numeric vectors with shape & rate parameters of (inverse) gamma priors for theta and tau.
#' @param thetadistr,taudistr The kind of prior distribution for theta & tau, can be `"gamma"` (default) or `"invgamma"`.
#' @param theta Logical, whether to estimate theta explicitly (default is `TRUE`) or to use analytical integration
#'   of population size (if FALSE, it disables also `threads` argument).
#' @param phiprior A numeric vector with alpha & beta parameters of beta prior for phi in the MSC+I model.
#' @param wprior A numeric vector with shape & rate parameters of gamma prior for mutation-scaled migration rates in the MSC+M model.
#' @param migration A matrix or data frame with columns `"from"`, `"to"`, `"shape"` and `"rate"`, which contain
#'   tip or node labels and shape & rate parameters of gamma prior specified for the given migration rate.
#'   If `"shape"` and `"rate"` are missing, the wprior is used instead.
#' @param spdelimalgorithm The speciesdelimitation algorithm. The algorithm is specified by the first number (0 or 1),
#'   while the others are its finetune parameters (one for algorithm 0, two for algorithm 1).
#' @param speciesmodelprior The species tree prior. The options are: 0 = uniform labelled histories,
#'   1 = uniform rooted trees (default), 2 = uniformSLH, 3 = uniformSRooted.
#' @param phase Logical or binary numeric, whether to phase sequence of particular species, recycled if necessary.
#' @param model The nucleotide substitution model, default is `"JC69"`.
#' @param gamma Logical, whether to model between site variation of substitution rates, default is `FALSE`.
#' @param alphaprior If `gamma == TRUE` this is a prior for the shape of the discretized gamma distribution.
#' @param locusrate Specifies substitution rate variation across loci. The options are: 0 = no varitation (default),
#'   c(1, a_mubar, b_mubar, a_mui, <prior>) = locus rates are estimated or
#'   c(2, LocusRateFileName) = locus rates are specified in a file. 
#' @param clock Specifies substitution rate variation across branches. The options are 1 = strict clock (default),
#'   c(2, a_vbar, b_vbar, a_vi, <prior>, <distribution>) = independent rates or
#'   c(3, a_vbar, b_vbar, a_vi, <prior>, <distribution>) = correlated rates.
#' @param heredity Specifies variation in theta across loci. The options are: 0 = no varitation (default),
#'   c(1, a_gamma b_gamma) = estimate theta multipliers from data under gamma prior or
#'   c(2, HeredityFileName) = theta multipliers are specified in a file.
#' @param thetamodel Constraints on theta parameters.
#' @param cleandata Whether to remove ambiguity nucleotides, default is TRUE when phase is TRUE and FALSE otherwise.
#' @param finetune Logical, whether to perform automatic adjustment of MCMC step lengths. 
#' @param print A vector of five binary values denoting whether to print:
#'   MCMC samples, locusrates(mu_i, nu_i), heredityscalars, genetrees and locusrateparameters, respectively.
#' @param burnin The number of burn-in iterations, default is 20000.
#' @param sampfreq The frequency of MCMC sampling, default is 10.
#' @param nsample The total number of MCMC samples to be obtained, default is 20000.
#' @param seed The seed used to initialize the analysis (normally let to be -1, which means random seed).
#' @param usedata Whether to use data (default) or to sample from prior.
#' @param threads The number of threads (CPUs) used for computations.
#' @export

write_bpp_ctl <- function(con, spdelim, sptree, seqfile, imapfile, jobname, species, tree, nloci, thetaprior, tauprior, thetadistr="gamma", taudistr="gamma", theta=TRUE, phiprior=c(1,1), wprior=c(2,1), migration, spdelimalgorithm=c(0, 2), speciesmodelprior=1, phase=FALSE, model="JC69", gamma=FALSE, alphaprior=c(1,1,4), locusrate=0, clock=1, heredity=0, thetamodel="linked-none", cleandata, finetune=TRUE, print=c(1,0,0,0,0), burnin=50000, sampfreq=10, nsample=5000, seed=-1, usedata=TRUE, threads=NULL) {

	spdelim <- as.numeric(spdelim)
	sptree <- as.numeric(sptree)
	nspecies <- length(species)
	species <- sapply(list(names(species), unname(species)), paste, collapse=" ")
	if (thetadistr == "beta") {
		thetaprior <- c(thetadistr, as.numeric(thetaprior[1:4]))
	} else if (thetadistr == "gamma") {
		thetaprior <- as.numeric(thetaprior[1:2])
		thetaprior <- c("gamma", formatC(thetaprior, format="f", digits=0))
	} else {
		thetaprior <- as.numeric(thetaprior[1:2])
		thetaprior <- c(formatC(thetaprior[1], format="f", digits=0), formatC(thetaprior[2], format="f", digits=5))
	}
	if (isTRUE(theta) & thetadistr == "invgamma") {
		thetaprior <- c(thetaprior, "e")
	}
	if (taudistr == "gamma") {
		tauprior <- c("gamma", formatC(as.numeric(tauprior), format="f", digits=0))
	} else {
		tauprior <- as.numeric(tauprior[1:2])
		tauprior <- c(formatC(tauprior[1], format="f", digits=0), formatC(tauprior[2], format="f", digits=5))
	}
	if (spdelim == 1 & spdelimalgorithm[1] == 0 & length(spdelimalgorithm) != 2) {
		spdelimalgorithm <- c(0, 2)
	} else if (spdelim == 1 & spdelimalgorithm[1] == 1 & length(spdelimalgorithm) != 3) {
		spdelimalgorithm <- c(1, 2, 1)
	}
	if (spdelim == 1) {
		spdelim <- c(spdelim, spdelimalgorithm)
	}
	if (sptree == 1) {
		speciesmodelprior <- paste("speciesmodelprior =", speciesmodelprior)
	} else {
		speciesmodelprior <- NULL
	}
	if (isTRUE(grepl("\\[", tree))) {
		phiprior <- paste0("phiprior = ", paste(phiprior, collapse=" "))
	} else {
		phiprior <- NULL
	}
	if (!missing(migration)) {
		if (length(migration) == 2 & !is.matrix(migration)) {
			migration <- matrix(migration, 1, 2)
		}
		wprior <- paste0("wprior = ", paste(wprior, collapse=" "))
		nmigrations <- paste0("migration = ", as.integer(nrow(migration)))
		migration <- c(wprior, nmigrations, paste("          ", migration[,1], migration[,2]))
	} else {
		migration <- NULL
	}
	phase <- rep_len(as.numeric(phase), nspecies)
	if (isTRUE(gamma)) {
		alphaprior <- paste0("alphaprior = ", paste(alphaprior, collapse=" "))
	} else {
		alphaprior <- NULL
	}
	tree <- get_topology(tree)
	cleandata <- ifelse(missing(cleandata), as.numeric(all(as.logical(phase))), as.numeric(cleandata))
	if (isTRUE(as.logical(finetune))) {
		finetune <- "finetune = 1"
	} else {
		finetune <- "finetune = 0"
	}
	print <- as.numeric(print)
	usedata <- as.numeric(usedata)
	if (isFALSE(theta)) {
		threads <- NULL
	}
	if (!is.null(threads)) {
		threads <- paste0("threads = ", paste(threads, collapse=" "))
	}

	space <- paste(rep(" ", 17), collapse="")
	
	ctl <- c(paste0("seed = ", seed, "\n"), paste("seqfile =", seqfile), paste("Imapfile =", imapfile), paste0("jobname = ", jobname, "\n"),
		paste("speciesdelimitation =", paste(spdelim, collapse=" ")), paste("speciestree =", sptree), speciesmodelprior,  
		paste("\nspecies&tree =", nspecies, species[1]), paste0(space, species[2]), paste0(space, tree), paste0("phase = ", paste(phase, collapse=" "), "\n"),
		paste("usedata =", usedata), paste("nloci =", nloci), paste("model =", model), alphaprior, paste0("\ncleandata = ", cleandata, "\n"), 
		paste0("thetaprior = ", paste(thetaprior, collapse=" ")), paste0("tauprior = ", paste(tauprior, collapse=" ")), phiprior, migration, 
		paste0("\nlocusrate = ", paste(locusrate, collapse=" ")), paste0("clock = ", paste(clock, collapse=" ")), paste0("heredity = ", paste(heredity, collapse=" ")), paste0("thetamodel = ", thetamodel),
		paste0("\n", finetune), paste("print =", paste(print, collapse=" ")), 
		paste("burnin =", formatC(burnin, format="f", digits=0)), paste("sampfreq =", formatC(sampfreq, format="f", digits=0)), paste("nsample =", formatC(nsample, format="f", digits=0)),
		threads)

	writeLines(ctl, con)

}



#' Random starting tree.
#' 
#' @description
#' Creates a starting tree for species tree estimation (A11 or A01 analysis), using [ape::rtree()].
#' 
#' @param species A character vector of species names.
#' @param imap A data frame with a component `"Species"` or just two columns (2nd assumed to contain
#'   species labels) or a name of 'imap' file.
#' @param seed Optional numeric, seed for the random generator, it makes results of [ape::rtree()] repeatable.
#' @returns Starting tree topology in Newick format.
#' @export

get_starting_tree <- function(species, imap, seed=sample(1e8, 1)) {
	if (missing(species)) {
		if (is.character(imap)) {
			species <- read_bpp_imap(imapfile)[,"Species"]
		} else if ("Species" %in% colnames(imap)) {
			species <- imap[,"Species"]
		} else if (ncol(imap) == 2) {
			species <- imap[,2]
		}
	}
	species <- unique(species)
	set.seed(seed)
	tree <- ape::rtree(length(species), tip.label=species)
	tree$edge.length <- NULL
	tree$node.label <- NULL	
	ape::write.tree(tree)
}



#' Maximum sample sizes.
#' 
#' @description
#' Finds the maximum sample size (no. of sequences) per locus for each of the species.
#' 
#' @param seq The sequence file name or a list of alignments.
#' @param imap The imap file name or a data frame with the content of imap file.
#' @param allele An allele identifier, regular expression distinguishing sequences from the same individual.
#' @returns A named numeric vector with maximum sample sizes (the names are species labels).
#' @export

get_maxsamplesize <- function(seq, imap, allele=NULL) {
	if (is.character(seq) & length(seq) == 1) {
		seqlines <- readLines(seqfile)
		seqlines <- seqlines[nchar(seqlines) > 0]
		start <- grep("^[[:digit:][:blank:]]+$", seqlines)
		loci <- rep(seq_along(start), diff(c(start, length(seqlines) + 1)))
		seqlines <- strsplit(seqlines, "\\s")
		nams <- sapply(seqlines, "[[", 1)
		zero <- grepl("^[N[:punct:]\\?]+$", sapply(seqlines, "[[", 2))
		nams <- lapply(split(nams[!zero], loci[!zero]), function(x) x[-1])
		nams <- lapply(nams, function(x) gsub("^.+\\^", "", x))	
	} else {
		nams <- lapply(seq, rownames)
		if (!is.null(allele)) {
			nams <- lapply(nams, function(x) sub(allele, "", x))
		}
	}
	if (is.character(imap)) {
		imap <- read_bpp_imap(imap, col.names=c("Individual", "Species"))
	}
	tab <- lapply(nams, function(x, imap) table(imap[match(x, imap[,1]),2]), imap=imap)
	mat <- matrix(,length(unique(imap[,2])), length(tab), dimnames=list(sort(unique(imap[,2])), names(tab)))
	for (i in seq_along(tab)) {
		mat[names(tab[[i]]),i] <- as.numeric(tab[[i]])
	}
	return(apply(mat, 1, max, na.rm=TRUE))
}



#' Maximum number of loci.
#' 
#' @description
#' Finds the maximum number of loci.
#' 
#' @param seqfile The sequence file name.
#' @returns The number of loci.
#' @export

get_nloci <- function(seqfile) {
	seqlines <- readLines(seqfile)
	seqlines <- seqlines[nchar(seqlines) > 0]
	nloci <- length(grep("^[[:digit:][:blank:]]+$", seqlines))
	return(nloci)
}



#' Tree topology.
#' 
#' @description
#' Writes tree topology in Newick format.
#' 
#' @param tree An object of class `phylo` or a character string with tree in Newick format.
#' @returns A character string with tree topology (no branch lengths) in Newick format.
#' @export

get_topology <- function(tree, node.label=TRUE) {
	tree <- import_tree(tree=tree)
	if (inherits(tree, "phylo")) {
		tree <- list(tree)
	}
	for (i in seq_along(tree)) {
		tree[[i]]$edge.length <- NULL
	}
	if (isFALSE(node.label)) {
		for (i in seq_along(tree)) {
			tree[[i]]$node.label <- NULL
		}	
	}
	tree <- unname(sapply(tree, ape::write.tree))
	return(tree)
}

