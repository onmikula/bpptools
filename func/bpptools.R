# Function:
#	write_bpp_input: writes seq (sequence) and imap (classification) files for BPP analyses
# Arguments:
#	seq: list of alignments (objects of 'matrix' or 'DNAbin' class), assumed to contain row names
#	info: data frame (or matrix) defining classification (=mapping) of sequences to individuals and individuals to (candidate) species. For an overview of acceptable arrangements, see Details.
#	seqfile: sequence file name
#	imapfile: imap file name
# Value:
#	writes seq and imap files as inputs for BPP analyses
# Details:
#	Possible arrangements of 'info' data frame: (1) the first three columns contain names of sequences, individuals and species respectively. (2) these three columns can be differently ordered or included at arbitrary places of a larger data frame, if their names are 'Sequence', 'Individual' and 'Species' (or unambiguous abbreviations). (3) if 'Sequence' column is missing, sequence and individual names are assumed identical. (4) if just two (arbitrarily labelled) columns are present, they are assumed to contain individual and species and sequence names are assumed identical to individual ones.
#	All sequence names should be listed in the info file with an exception of names consisting of individual name followed by allele identifier. If, for instance, FR11_0 and FR11_1 sequences are present, only FR11 can be included in the 'Sequence' column.
#	Individual names should not be each other's partial matches (e.g. if an individual FR11 is present, the name FR1 should be avoided).

write_bpp_input <- function(seq, info, seqfile, imapfile) {

	seqlen <- function(x) {
		rowSums(matrix(!x %in% c("N","-","?"), nrow(x), ncol(x)))
	}

	nams <- c("Sequence","Individual","Species")
	info <- as.data.frame(info, stringsAsFactors=FALSE)
	which <- pmatch(tolower(names(info)), tolower(nams), duplicates.ok=FALSE)
	if (all(1:3 %in% which)) {
		imap <- info[,match(c(1,2,3), which)]
	} else if (all(2:3 %in% which)) {
		imap <- info[,match(c(2,2,3), which)]
	} else if (ncol(info) == 2) {
		imap <- info[,c(1,1,2)]
	} else {
		imap <- info[,c(1,2,3)]
	}
	names(imap) <- c("Sequence","Individual","Species")

	loci <- lapply(lapply(seq, as.matrix), toupper)
	for (i in seq_along(loci)) {
		loci[[i]] <- loci[[i]][seqlen(loci[[i]]) > 0,,drop=FALSE]
	}
	indnam <- vector("list", length(loci))
	for (i in seq_along(loci)) {
		seqnam <- rownames(loci[[i]])
		pmatches <- lapply(paste0("^", imap$Sequence), grep, x=seqnam)
		pmatches <- lapply(seq_along(seqnam), function(j) which(mapply("%in%", j, pmatches)))
		if (any(sapply(pmatches, length) == 0)) {
			stop("all sequence names have to be included in the 'info' data frame")
		}
		for (j in which(sapply(pmatches, length) > 1)) {
			nch <- nchar(imap$Sequence[pmatches[[j]]])
			pmatches[[j]] <- pmatches[[j]][nch == min(nch)]
		}
		indnam[[i]] <- imap$Individual[unlist(pmatches)]
		seqnam <- paste(seqnam, indnam[[i]], sep="^")
		seqstrings <- apply(loci[[i]], 1, paste, collapse="")
		loci[[i]] <- c(paste(dim(loci[[i]]), collapse=" "), paste(seqnam, seqstrings, sep=" "), "\n")
	}
	writeLines(unname(unlist(loci)), con=seqfile)

	indnam <- sort(unique(unlist(indnam)))
	imap <- imap[match(indnam, imap$Individual),c("Individual","Species")]
	imap <- unname(apply(imap, 1, paste, collapse=" "))
	if (missing(imapfile)) {
		imapfile <- paste0(sub("\\.[[:alnum:]]+$", "", seqfile), ".imap")
	}
	writeLines(imap, imapfile)
	
}



# Function:
#	write_bpp_imap: writes imap file, possibly merging species according to a previous species delimitation analysis
# Arguments:
#	imap: matrix or data frame (1st column individuals, 2nd column species) or a name of imap file
#	file: file name of the output
#	return: whether to return data frame
#	delim: result of species delimitation (mergers of candidate species) - either vector with concatenated candidate species names or an 'imap' object (see above) or a name of tab-delimited file where 'imap' object is stored
#	sep: separator of original candidate species names merged together, applied if 'delim' is a vector
# Value:
#	imap data frame (if return == TRUE) and/or .imap file (if its name is specified)

write_bpp_imap <- function(imap=NULL, file, return=TRUE, delim=NULL, sep="-") {
	if (!is.null(imap)) {
		if (is.character(imap) & length(imap) == 1) {
			imap <- read_bpp_imap(imap)
		}
		imap <- split(imap[,1], imap[,2])
	}
	if (!is.null(delim)) {
		if (is.character(delim) & length(delim) == 1) {
			delim <- read.delim(delim)
		}
		if (is.matrix(delim) | is.data.frame(delim)) {
			delim <- split(delim[,1], delim[,2])
		}
		if (!is.list(delim)) {
			delim <- setNames(strsplit(delim, split=sep), delim)
		}
	}

	if (!is.null(delim)) {
		if (!is.null(imap)) {
			delim <- lapply(delim, function(d) unlist(imap[d]))
		}
		imap <- data.frame(Ind=unlist(delim), Pop=rep(names(delim), sapply(delim, length)), row.names=NULL)
	} else if (is.null(imap)) {
		stop("at least one of imap and delim arguments has to be specified")
	}
	
	if (!missing(file)) {
		writeLines(unname(apply(imap, 1, paste, collapse=" ")), file)
	}
	if (isTRUE(return)) {
		return(imap)
	}
}



# Function:
#	read_bpp_tree: reads species tree from BPP control file
# Arguments:
#	ctlfile: the name of control file
# Value:
#	the starting species tree topology in newick format with attribute 'tips' which contains tip labels in the same order as in the control file

read_bpp_tree <- function(ctlfile) {
	lin <- readLines(ctlfile)
	tips <- sub("^\\s*species&tree\\s*=\\s*[[:digit:]]+\\s*", "", lin[grepl("species&tree", lin)])
	tips <- unlist(strsplit(gsub("\\s+", " ", tips), " "))
	tree <- regmatches(lin, regexpr("\\(.*;$", lin))		
	tree <- gsub("[[:blank:]]", "", tree)	
	attr(tree, "tips") <- tips		
	return(tree)
}



# Function:
#	read_bpp_seq: reads sequences from BPP sequence file
# Arguments:
#	seqfile: the name of sequence file
#	loci: names of loci (if not provided, the loci are assigned names Lddd, d is for digits)
#	indnames: whether to use individual rather than sequence names to label sequences (useful for functions like 'estimate_otus' or 'estimate_params')
# Value:
#	the list of matrices with alignments of input sequences

read_bpp_seq <- function(seqfile, loci=NULL, indnames=FALSE) {
	lin <- readLines(seqfile)
	lin <- lin[lin != ""]
	header <- which(!grepl("\\^", lin))
	start <- header + 1
	end <- c(header[-1] - 1, length(lin))
	nloci <- length(start)
	if (is.null(loci)) {
		L <- paste("L", formatC(seq(nloci), format="d", flag=0, width=nchar(nloci)), sep="")
	} else {
		L <- loci
	}
	loci <- setNames(vector("list", nloci), L)
	for (i in seq(nloci)) {
		loci[[i]] <- lin[start[i]:end[i]]
		loci[[i]] <- unlist(strsplit(loci[[i]], "\\s+"))
		ii <- seq_along(loci[[i]]) %% 2
		loci[[i]] <- setNames(loci[[i]][ii == 0], loci[[i]][ii == 1])
		nam <- strsplit(names(loci[[i]]), "\\^")
		if (nchar(nam[[1]][1]) == 0 | isTRUE(indnames)) {
			names(loci[[i]]) <- sapply(nam, "[", 2)
		} else {
			names(loci[[i]]) <- sapply(nam, "[", 1)
		}
		
		loci[[i]] <- do.call(rbind, strsplit(loci[[i]], ""))
	}
	return(loci)
}



# Function:
#	read_bpp_imap: reads BPP imap file
# Arguments:
#	imapfile: the name of imap file
# Value:
#	data frame with columns Individual and Species which maps (=classifies) individuals into species/populations.

read_bpp_imap <- function(imapfile) {
	read.table(imapfile, header=FALSE, col.names=c("Individual", "Species"))
}



# Function:
#	read_bpp_prior: reads BPP imap file
# Arguments:
#	ctlfile: the name of control file
#	which: prior identification ("theta" & "tau" by default)
# Value:
#	data frame with columns 'alpha', 'beta', 'distribution' and possibly 'ncat' (for  discretized gamma)

read_bpp_prior <- function(ctlfile, which=c("theta","tau")) {
	spaces <- "[[:punct:][:space:]]+"
	ctl <- readLines(ctlfile)
	ctl <- gsub(paste0("^", spaces, "|", spaces, "$"), "", ctl)
	ctl <- ctl[sapply(paste0("^", which, "prior"), grep, x=ctl)]
	null <- which(sapply(ctl, length) == 0)
	if (length(null) > 0) {
		sep <- rev(rep_len(rep(c("", "&", ","), c(1, 1, length(null))), length(null)))
		missingpriors <- paste(as.character(rbind(which[null], sep))[-2*length(null)], collapse=" ")
		plural <- ifelse(length(null) == 1, "prior is", "priors are")
		stop(paste(missingpriors, plural, "missing"))
	}
	ctl <- gsub(paste0("^[[:alpha:]]+", spaces, "=", spaces), "", ctl)
	ctl <- strsplit(ctl, "[[:space:]]+")
	gammas <- sapply(ctl, "[[", 1) == "gamma"
	if (sum(gammas) > 0) {
		ctl[gammas] <- lapply(ctl[gammas], "[", -1)
	}
	ab <- do.call(rbind, lapply(lapply(ctl, "[", 1:2), as.numeric))
	priors <- data.frame(alpha=ab[,1], beta=ab[,2], distribution=ifelse(gammas, "gamma", "invgamma"), row.names=which)
	alpha <- which == "alpha"
	if (any(alpha)) {
		priors$distribution[alpha] <- "gamma"
		priors$ncat <- ifelse(alpha, as.numeric(sapply(ctl[alpha], "[[", 3)), NA)
	}
	phi <- which == "phi"
	if (any(phi)) {
		priors$distribution[phi] <- "beta"
	}
	return(priors)
}



# Function:
#	read_bpp_output: reads mcmc trace of BPP into R
# Arguments:
#	mcmc: character string giving name of "mcmc file"
#	starting: whether to retain starting state in the trace
#	thinning: proportion of mcmc steps analyzed
#	tree: a target tree upon which divergence times from A00 analysis are projected (in newick format or as 'phylo' object of ape)
# Value:
#	list with components 'trees' ('multiPhylo' object of ape), 'params' (data frame) and 'nspecies' (numeric vector)
# Details:
#	nspecies can included MCMC trace of number of species (A11, A10) or the fixed number of species
# Depends on: ape

read_bpp_output <- function(mcmc, starting=FALSE, thinning=1, tree=NULL) {
	A00 <- substr(readLines(mcmc, n=1), 1, 1) != "("
	if (A00 == FALSE) {
		steps <- readLines(mcmc)
		if (isFALSE(starting)) {
			steps <- steps[-1]
		}
		steps <- steps[seq(1, length(steps), floor(1/thinning))]
		steps <- strsplit(steps, split=";")
		nspecies <- as.numeric(sapply(steps, "[", 2))
		trees <- paste(sapply(steps, "[", 1), ";", sep="")
		trees <- lapply(trees, function(x) ape::read.tree(text=x))
		trees <- lapply(trees, ape::reorder.phylo, order="postorder")
		ntips <- ape::Ntip(trees[[1]])
		params <- vector("list", length(steps))
		for (i in seq_along(steps)) {
			brtimes <- ape::branching.times(trees[[i]])
			nodeord <- order(match(names(brtimes), trees[[i]]$edge[,1]))
			params[[i]] <- data.frame(tau=c(rep(0, ntips), brtimes[nodeord]))
			thetas <- c(trees[[i]]$tip.label, trees[[i]]$node.label)
			thetas <- suppressWarnings(as.numeric(gsub("[[:alnum:]]*#", "", thetas)))
			tipord <- order(match(seq(ntips), trees[[i]]$edge[,2]))
			params[[i]]$theta <- thetas[c(tipord, nodeord + ntips)]
			trees[[i]]$tip.label <- gsub("#[[:alnum:][:punct:]]*$", "", trees[[i]]$tip.label)
			trees[[i]]$node.label <- NULL
			rownames(params[[i]])[seq(ntips)] <- trees[[i]]$tip.label[tipord]
		}
		class(trees) <- "multiPhylo"
	} else {
		params <- read.delim(mcmc, row.names=1)
		if (isFALSE(starting)) {
			params <- params[-1,]
		}
		params <- params[seq(1, nrow(params), floor(1/thinning)),]
		nspecies <- sum(grepl("^tau", names(params))) + 1
		if (!is.null(tree)) {
			settaus <- function(taus, phy) {
				heights <- ifelse(phy$edge <= ape::Ntip(phy), 0, phy$edge - ape::Ntip(phy))
				heights[heights != 0] <- taus[heights[heights != 0]]
				phy$edge.length <- heights[,1] - heights[,2]
				return(phy)		
			}
			if (is.character(tree)) {
				tree <- ape::read.tree(text=tree)
			}
			tree <- ape::reorder.phylo(tree, order="postorder")			
			taus <- lapply(split(params[,grepl("tau", colnames(params))], seq(nrow(params))), as.numeric)
			trees <- lapply(taus, settaus, phy=tree)		
			class(trees) <- "multiPhylo"
		} else {
			trees <- NULL
			warning("trees are not returned as their topology is not specified")
		}
	}
	result <- list(trees=trees, params=params, nspecies=nspecies)
	class(result) <- "bpp"
	return(result)
}



# Function:
#	combine_bpp_runs: combines posterior samples from independent BPP runs
# Arguments:
#	runs: list of "bpp" objects (outputs from 'read_bpp_output')
# Value:
#	"bpp" object with pooled posterior samples

combine_bpp_runs <- function(runs) {
	trees <- Reduce(c, lapply(runs, "[[", "trees"))
	params <- lapply(runs, "[[", "params")
	if (is.null(dim(params))) {
		params <- Reduce(c, params)
		nparams <- length(params)
	} else {
		params <- Reduce(rbind, params)
		nparams <- nrow(params)
	}
	nspecies <- Reduce(c, lapply(runs, "[[", "nspecies"))
	if (length(nspecies) < nparams) {
		nspecies <- nspecies[1]
	}
	result <- list(trees=trees, params=params, nspecies=nspecies)
	class(result) <- "bpp"
	return(result)
}



# Function:
#	export_bpp_output: exports mcmc trace of BPP into '.trees' file (readable by FigTree, ape::read.tree etc.) and '.log' file readable by Tracer
# Arguments:
#	bpp: object of class "bpp" (from 'read_bpp_output')
#	trees: name of 'trees' file
#	log: name of 'log' file
#	model: "A00", "A01", "A10" or "A11"
# Value:
#	exports '.trees' file (text file with trees in newick format) and '.log' file (tab-delimited text file)
# Depends on: ape

export_bpp_output <- function(bpp, trees=NULL, log=NULL, model="A11") {
	spdelimitation <- substring(model, 2, 2) == "1"
	sptree <- substring(model, 3, 3) == "1"
	if (sptree) {
		tau0 <- sapply(lapply(bpp$params, "[[", "tau"), max)
		mtheta <- sapply(lapply(bpp$params, "[[", "theta"), mean, na.rm=TRUE)
		params <- data.frame(tau.0=tau0, mtheta=mtheta)
		if (spdelimitation) {
			params$nspecies <- bpp$nspecies
		}
	} else {
		if (spdelimitation) {
			ntips <- ape::Ntip(bpp$trees[[1]])
			nnode <- ape::Nnode(bpp$trees[[1]])
			tau <- do.call(rbind, lapply(bpp$params, "[[", "tau"))[,-seq(ntips)]
			colnames(tau) <- paste("tau", seq(nnode) - 1, sep=".")
			params <- data.frame(tau)
			params$mtheta <- sapply(lapply(bpp$params, "[[", "theta"), mean, na.rm=TRUE)
			params$nspecies <- bpp$nspecies
		} else {
			params <- bpp$params
		}
	}
	if(!is.null(log)) {
		write.table(cbind(Sample=seq(nrow(params)) - 1, params), file=log, sep="\t", quote=FALSE, row.names=FALSE)
	}
	if (!is.null(trees) & !is.null(bpp$trees)) {
		ape::write.tree(bpp$trees, file=trees)
	}
}



# Function:
#	export_msc_tree: exports tree with MSC parameters stored in nexus metacomment format
# Arguments:
#	target: tree topology tree in newick format obtained from BPP control file (by 'read_bpp_tree') or a name of BPP control file
#	trees: posterior sample of trees ('multiPhylo' or a list of 'phylo' objects)
#	params: data frame with posterior sample of model parameters
#	file: file name for export of the final annotated tree
#	treesfile: file name for export of the annotated posterior sample of trees
#	fun: summary function, 'mean' (default) or 'median'
#	hpd: proportion of posterior sample defining the highest probability density interval
#	figtree: whether to export the final tree and the posterior sample as .nexus files readable by FigTree 
# Depends on: ape

export_msc_tree <- function(target, trees, params, file, treesfile, fun="mean", hpd=0.95, figtree=FALSE) {
	make_figtree <- function(newick, tips) {
		ntip <- length(tips)
		transl <- paste(paste("\t\t", seq(ntip), sep=""), paste0(tips, ","))
		transl[length(transl)] <- sub("\\,$", "", transl[length(transl)])
		transl <- c("\tTranslate", transl, "\t;")
		for (i in seq(ntip)) {
			pattern <- paste0("[[:punct:]]{1}", tips[i], "[[:punct:]]{1}")
			for (j in seq_along(newick)) {
				m <- regexpr(pattern, newick[j]) + 1
				attr(m, "match.length") <- attr(m, "match.length") - 2
				regmatches(newick[j], m) <- i
			}
		}
		taxa <- c("Begin taxa;", paste("\tDimensions ntax=", ntip, ";", sep=""), "\t\tTaxlabels", paste("\t\t", tips, sep=""), "\t;", "End;")
		no <- formatC(seq_along(newick), format="d", flag=0, width=nchar(length(newick)))
		return(c("#NEXUS", "\n", taxa, "\n", "Begin trees;", transl, paste("tree", no, "=", newick), "End;"))
	}
	
	if (substr(target, 1, 1) != "(") {
		target <- read_bpp_tree(target)
	}
	tips <- attr(target, "tips")
	ntip <- length(tips)
	nodes <- seq(ntip - 1) + ntip
	constraint <- ape::read.tree(text=target)$tip.label
	topologies <- sapply(lapply(trees, ape::rotateConstr, constraint=constraint), get_topology)
	subset <- topologies == target
	if (!any(subset)) {
		stop("none of the trees has topology identical to the target")
	}
	nsample <- sum(subset)
	target <- ape::reorder.phylo(ape::read.tree(text=target), order="postorder")
	target$node.label <- nodes
	params <- params[subset,,drop=FALSE]
	params <- params[,grepl("^theta_*[[:digit:]]+|^tau_*[[:digit:]]+", colnames(params)),drop=FALSE]
	colnames(params) <- regmatches(colnames(params), regexpr("^theta_*[[:digit:]]+|^tau_*[[:digit:]]+", colnames(params)))
	params[,seq(ntip)] <- params[,match(tips, target$tip.label)]
	tau <- cbind(as.data.frame(setNames(rep(list(numeric(nsample)), ntip),paste0("tau",seq(ntip)))), params[,grepl("tau",colnames(params))])
	theta <- params[,grepl("theta",colnames(params))]
	anc <- target$edge[match(seq(ncol(theta)), target$edge[,2]),1]
	gdi <- matrix(NA, nsample, ncol(theta), dimnames=list(rownames(params), paste0("gdi_",seq(ncol(theta)))))
	for (i in setdiff(seq(ncol(theta)), ntip+1)) {
		gdi[,i] <- 1 - exp(-2 * (tau[,anc[i]] - tau[,i]) / theta[,i])
	}
	fun <- match.fun(fun)
	mtau <- round(apply(tau, 2, fun), 6)
	mtheta <- round(apply(theta, 2, fun), 6)
	mgdi <- round(apply(gdi, 2, fun), 6)
	itau <- round(apply(tau, 2, HPD, p=hpd), 6)
	itheta <- round(apply(theta, 2, HPD, p=hpd), 6)
	igdi <- round(apply(gdi, 2, HPD, p=hpd), 6)

	tau <- round(tau, 6)
	theta <- round(theta, 6)
	gdi <- round(gdi, 6)
	node_labels <- matrix(paste0("[&tau=", as.matrix(tau[,-seq(ntip)]), "]"), nsample, ncol(tau) - ntip)
	branch_labels <- matrix(paste0("[&theta=", as.matrix(theta), ",", "gdi=", as.matrix(gdi), "]"), nsample, ncol(theta))
	node_labels[,1] <- paste0(gsub("]", "", node_labels[,1]), ",", gsub("\\[&|,gdi=NA", "", branch_labels[,ntip+1]))
	for (i in seq_along(trees)) {
		trees[[i]]$node.label <- nodes
	}
	trees <- unname(sapply(trees, ape::write.tree))
	branch_patterns <- c(paste0(target$tip.label, ":"), paste0(")", target$node.label, ":"))
	node_patterns <- branch_patterns[-seq(ntip)]
	node_patterns[1] <- sub(":", ";", node_patterns[1])
	for (i in seq_along(trees)) {
		m <- regexpr(node_patterns[1], text=trees[i])
		regmatches(trees[i], m) <- paste0(")", node_labels[i,1], ";")
		for (j in seq_along(nodes)[-1]) {
			m <- regexpr(branch_patterns[j+ntip], text=trees[i])
			regmatches(trees[i], m) <- paste0(regmatches(trees[i], m), branch_labels[i,j+ntip])
			m <- regexpr(node_patterns[j], text=trees[i])
			regmatches(trees[i], m) <- paste0(")", node_labels[i,j], ":")
		}
		for (j in seq(ntip)) {
			m <- regexpr(branch_patterns[j], text=trees[i])
			regmatches(trees[i], m) <- paste0(regmatches(trees[i], m), branch_labels[i,j])
		}
	}
	
	target$edge.length <- mtau[target$edge[,1]] - mtau[target$edge[,2]]
	target <- ape::write.tree(target)
	prob <- paste0(substr(hpd, 3, 4), "%")
	itaulab <- paste0("tau", prob, "_HPD={", apply(itau, 2, paste, collapse=","), "}")
	ithetalab <- paste0("theta", prob, "_HPD={", apply(itheta, 2, paste, collapse=","), "}")
	igdilab <- paste0("gdi", prob, "_HPD={", apply(igdi, 2, paste, collapse=","), "}")
	m <- regexpr(node_patterns[1], text=target)
	regmatches(target, m) <- paste0(")", "[&tau=", mtau[ntip+1], ",", itaulab[ntip+1], ",", "theta=", mtheta[ntip+1], ",", ithetalab[ntip+1], "]", ";")
	for (j in seq_along(nodes)[-1]) {
		m <- regexpr(branch_patterns[j+ntip], text=target)
		regmatches(target, m) <- paste0(regmatches(target, m), "[&theta=", mtheta[ntip+j], ",", ithetalab[ntip+j], ",", "gdi=", mgdi[ntip+j], ",", igdilab[ntip+j], "]")
		m <- regexpr(node_patterns[j], text=target)
		regmatches(target, m) <- paste0(")", "[&tau=", mtau[ntip+j], ",", itaulab[ntip+j], "]:")
	}
	for (j in seq(ntip)) {
		m <- regexpr(branch_patterns[j], text=target)
		regmatches(target, m) <- paste0(regmatches(target, m), "[&theta=", mtheta[j], ",", ithetalab[j], ",", "gdi=", mgdi[j], ",", igdilab[j], "]")
	}	

	if (isTRUE(figtree)) {
		if (!missing(file)) {
			writeLines(make_figtree(target, tips), con=file) 
		}
		if (!missing(treesfile)) {
			writeLines(make_figtree(trees, tips), con=treesfile) 
		}
	} else {
		if (!missing(file)) {
			writeLines(target, con=file)
		}
		if (!missing(treesfile)) {
			writeLines(trees, con=treesfile)
		}		
	}
	
}



# Function:
#	read_msc_tree: reads in multispecies coalescent tree written by export_msc_tree
# Arguments:
#	tree: name of the file with MSC-annotated tree
# Value:
#	tree in 'phylo' format with an additional component 'params' which include MSC annotations associated with particular nodes: divergence time (~distance from tips) and theta (associated with the suporting branch of the node). Root population has theta but not tau.

read_msc_tree <- function(tree) {
	tree <- readLines(tree)
	tree <- gsub("^[[:punct:][:space:]]+", "", tree)
	tree <- tree[nchar(tree) > 0]
	lower <- tolower(tree)
	if (grepl("nexus", lower[1])) {
		transl <- tree[(grep("^translate", lower) + 1):(grep("^tree", lower)[1] - 1)]
		transl <- gsub(",", "", transl)
		transl <- do.call(rbind, strsplit(transl, split="[[:space:]]+"))
		tree <- tree[grep("^tree", lower)[1]]
		tree <- regmatches(tree, regexpr("\\(.+;", tree))
		tree <- gsub("^[[:space:]]+", "", tree)
		for (i in seq(nrow(transl))) {
			tree <- sub(paste0(transl[i,1],":"), paste0(transl[i,2],":"), tree)
		}		
	}
	brlab <- gregexpr("]:\\[[^]]*][^:]", tree)
	attr(brlab[[1]],"match.length") <- attr(brlab[[1]],"match.length") - 1
	regmatches(tree, brlab) <- list(paste0(gsub("^.*&", ",", unlist(regmatches(tree, brlab))), ":"))
	brlab <- gregexpr(":\\[[^]]*]", tree)
	regmatches(tree, brlab) <- list(paste0(gsub("^:", "", unlist(regmatches(tree, brlab))), ":"))
	m <- gregexpr("\\[[^]]+]", tree)
	node_labels <- setNames(unlist(regmatches(tree, m)), paste0("Node", seq_along(unlist(m))))
	regmatches(tree, m) <- list(names(node_labels))
	tree <- ape::read.tree(text=tree)
	ntip <- length(tree$tip.label)
	node_labels <- unname(node_labels[c(gsub("^.+Node", "Node", tree$tip.label), tree$node.label)])
	m <- gregexpr("\\{[^}]+}", node_labels)
	regmatches(node_labels, m) <- lapply(regmatches(node_labels, m), function(x) gsub(",", "-", x))
	node_labels <- strsplit(gsub("\\[|&|]", "", node_labels), split=",")
	node_labels <- lapply(node_labels, function(x) x[!grepl("\\{", x)])
	params <- matrix(, length(node_labels), 3, dimnames=list(NULL, c("tau","theta","gdi")))
	params[seq(ntip),] <- cbind(0, gsub("^[[:alpha:]=]+", "", do.call(rbind, node_labels[seq(ntip)])))
	params[ntip+1,1:2] <- gsub("^[[:alpha:]=]+", "", node_labels[[ntip+1]])
	params[-seq(ntip+1),] <- gsub("^[[:alpha:]=]+", "", do.call(rbind, node_labels[-seq(ntip+1)]))
	mode(params) <- "numeric"
	tree$params <- params
	tree$tip.label <- gsub("Node[[:digit:]]+$", "", tree$tip.label)
	tree$node.label <- NULL
	return(tree)	
}



# Function:
#	parse_bpp_ctlfile: parses BPP control file
# Arguments:
#	con: file name or stdout() if you want to write the output on the screen
#	spdelim: logical, whether to estimate species delimitation
#	sptree: logical, whether to estimate species tree
#	seqfile, imapfile: names of the files with sequences (seqfile) a mapping of individuals to species (imapfile)
#	outfile, mcmcfile: names of the files with program output (outfile) and mcmc trace (mcmcfile)
#	species: named numeric vector, numbers are maxiumum numbers of sequences for particular species, names are species labels as they appear in the species tree
#	tree: species tree topology in the newick format
#	nloci: how many loci to use (in the order of seqfile)
#	thetaprior, tauprior: numeric vectors with alpha & beta parameters of (inverse) gamma priors for theta and tau
#	thetadistr, taudistr: type of prior distribution for theta & tau, can be 'invgamma' (default) or 'gamma'
#	theta: whether to estimate theta explicitly (if TRUE, default) or to use analytical integration of population size (if FALSE< disabling also 'threads' argument)
#	phiprior: numeric vectors with alpha & beta parameters of beta prior for phi in the MSci model 
#	spdelimalgorithm: rjMCMC speciesdelimitation algorithm, the first number (0 or 1) specifies the algorithm, the others specify its finetune parameters (one for algorithm zero, two for algorithm 1), see BPP documentation and Rannala & Yang (2010) for more details
#	speciesmodelprior: species tree prior, options: 0 = uniform labelled histories, 1 = uniform rooted trees (default), 2 = uniformSLH, 3 = uniformSRooted
#	phase: whether to phase sequence of particular species (default is FALSE), single value applies to all species
#	model: nucleotide substitution model "JC69", gamma=FALSE, alphaprior,
#	gamma: FALSE, alphaprior
#	locusrate: substitution rate variation across loci; the options are: 0 = no varitation (default), c(1, a_mubar, b_mubar, a_mui, <prior>) = locus rates are estimated or c(2, LocusRateFileName) = locus rates are specified in a file, see BPP documentation for more details 
#	clock: substitution rate variation across branches; the options are 1 = strict clock (default), c(2, a_vbar, b_vbar, a_vi, <prior>, <distribution>) = independent rates or c(3, a_vbar, b_vbar, a_vi, <prior>, <distribution>) = correlated rates, see BPP documentation for more details
#	heredity: variation in theta across loci; the options are: 0 = no varitation (default), c(1, a_gamma b_gamma) = estimate theta multipliers from data under gamma prior or c(2, HeredityFileName) = theta multipliers are specified in a file, see BPP documentation for more details 
#	cleandata: whether to remove ambiguity nucleotides, default is TRUE when phase is TRUE and FALSE otherwise 
#	finetune: numeric vector of length one or eight, depending on whether the lengths of MCMC steps are supplied
#	print: vector of five binary values denoting whether to print: MCMC samples, locusrates(mu_i, nu_i), heredityscalars, genetrees, locusrateparameters
#	burnin: number of burn-in iterations, default = 20000
#	sampfreq: frequency of MCMC sampling, default = 10
#	nsample: total number of MCMC samples, default = 20000
#	seed: random seed used to initialize the analysis
#	usedata: whether to use data (default) or to sample from prior
#	threads: no. of threads (CPUs) used for computations

parse_bpp_ctlfile <- function(con, spdelim, sptree, seqfile, imapfile, outfile, mcmcfile, species, tree, nloci, thetaprior, tauprior, thetadistr ="invgamma", taudistr ="invgamma", theta=TRUE, phiprior=c(1,1), spdelimalgorithm=c(0, 2), speciesmodelprior=1, phase=FALSE, model="JC69", gamma=FALSE, alphaprior=c(1,1,4), locusrate=0, clock=1, heredity=0, cleandata, finetune=1, print=c(1,0,0,0,0), burnin=20000, sampfreq=10, nsample=20000, seed=-1, usedata=TRUE, threads=NULL) {

	spdelim <- as.numeric(spdelim)
	sptree <- as.numeric(sptree)
	nspecies <- length(species)
	species <- sapply(list(names(species), unname(species)), paste, collapse=" ")
	if (isTRUE(theta)) {
		thetaprior <- c(thetaprior, "e")
	}
	if (thetadistr == "gamma") {
		thetaprior <- c("gamma", thetaprior, "e")
	}
	if (taudistr == "gamma") {
		tauprior <- c("gamma", tauprior)
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

	phase <- rep(as.numeric(phase), nspecies)
	if (isTRUE(gamma)) {
		alphaprior <- paste0("alphaprior = ", paste(alphaprior, collapse=" "))
	} else {
		alphaprior <- NULL
	}
	cleandata <- ifelse(missing(cleandata), as.numeric(all(as.logical(phase))), as.numeric(cleandata))
	if (length(finetune) == 8) {
		finetune <- paste("finetune =", paste0(finetune[1], ":"), paste(finetune[-1], collapse=" "))
	} else {
		finetune <- paste("finetune =", paste0(finetune[1], ":"), paste(c(.01, .02, .03, .04, .05, .01, .01), collapse=" "))
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
	
	ctl <- c(paste0("seed = ", seed, "\n"), paste("seqfile =", seqfile), paste("Imapfile =", imapfile), paste("outfile =", outfile), paste0("mcmcfile = ", mcmcfile, "\n"),
		paste("speciesdelimitation =", paste(spdelim, collapse=" ")), paste("speciestree =", sptree), speciesmodelprior,  
		paste("\nspecies&tree =", nspecies, species[1]), paste0(space, species[2]), paste0(space, tree), paste0("phase = ", paste(phase, collapse=" "), "\n"),
		paste("usedata =", usedata), paste("nloci =", nloci), paste("model =", model), alphaprior, paste0("\ncleandata = ", cleandata, "\n"), 
		paste0("thetaprior = ", paste(thetaprior, collapse=" ")), paste0("tauprior = ", paste(tauprior, collapse=" ")), phiprior, 
		paste0("\nlocusrate = ", paste(locusrate, collapse=" ")), paste0("clock = ", paste(clock, collapse=" ")), paste0("heredity = ", paste(heredity, collapse=" ")),
		paste0("\n", finetune), paste("print =", paste(print, collapse=" ")), 
		paste("burnin =", formatC(burnin, format="f", digits=0)), paste("sampfreq =", as.character(sampfreq)), paste("nsample =", as.character(nsample)),
		threads)

	writeLines(ctl, con)

}



# Function:
#	gdi: calculates genealogical divergence index for every branch of a specified species tree
# Arguments:
#	params: data frame with multispecies coalescent parameters
#	tree: tree in newick format obtained from BPP control file (by 'read_bpp_tree') or the name of BPP control file
# Value:
#	list with components 'means' (vector of posterior means) and 'distributions' (matrix with posterior distributions)
# Depends on: ape

gdi <- function(params, tree) {
	A00 <- !is.null(dim(params))
	if (isTRUE(A00)) {
		if (substr(tree, 1, 1) != "(") {
			tree <- read_bpp_tree(tree)
		}
		tips <- attr(tree, "tips")
		ntip <- length(tips)
		tree <- ape::reorder.phylo(ape::read.tree(text=tree), order="postorder")
		tree$node.label <- seq(ntip - 1) + ntip
		params <- params[,grepl("^theta_*[[:digit:]]+|^tau_*[[:digit:]]+", colnames(params)),drop=FALSE]
		colnames(params) <- regmatches(colnames(params), regexpr("^theta_*[[:digit:]]+|^tau_*[[:digit:]]+", colnames(params)))
		nsample <- nrow(params)
		params[,seq(ntip)] <- params[,match(tips, tree$tip.label)]
		tau <- cbind(as.data.frame(setNames(rep(list(numeric(nsample)), ntip),paste0("tau",seq(ntip)))), params[,grepl("tau",colnames(params))])
		theta <- params[,grepl("theta",colnames(params))]
		anc <- tree$edge[match(seq(ncol(theta)), tree$edge[,2]),1]
		pgdi <- matrix(NA, nsample, ncol(theta), dimnames=list(rownames(params), paste0("sp_",seq(ncol(theta)))))
		for (i in setdiff(seq(ncol(theta)), ntip+1)) {
			pgdi[,i] <- 1 - exp(-2 * (tau[,anc[i]] - tau[,i]) / theta[,i])
		}
		colnames(pgdi)[seq(ntip)] <- tree$tip.label
		pgdi <- pgdi[,-(ntip+1)]	
	} else {
		tips <- sort(tree[[1]]$tip.label)
		ntip <- length(tips)
		nsample <- length(params)
		theta <- do.call(rbind, lapply(params, function(x) x[tips,"theta"]))
		anc <- lapply(tree, function(x) x$edge[,1][match(match(tips, x$tip.label), x$edge[,2])])
		tau <- do.call(rbind, lapply(seq(nsample), function(i) params[[i]][anc[[i]],"tau"]))
		pgdi <- 1 - exp(-2 * tau / theta)
		colnames(pgdi) <- tips
	}
	mgdi <- round(apply(pgdi, 2, mean, na.rm=TRUE), 6)
	return(list(means=mgdi, distributions=pgdi))
}



# Function:
#	HPD: estimates highest posterior density interval
# Arguments:
#	x: posterior sample of a model paramater
#	p: the proportion of posterior density be included in the interval (defining 95% HPD by default, but possibly also 50% HPD and so on)
# Value:
#	numeric vector of posterior probabilities for clades ordered according to node numbers of their ancestors
# Details:
# "The 95% HPD stands for highest posterior density interval and represents the most compact interval on the selected parameter that contains 95% of the posterior probability. It can be loosely thought of as a Bayesian analog to a confidence interval." (Drummond & Bouckaert 2015, p. 89)
# Depends on: ape

HPD <- function(x, p=0.95) {
	x <- na.omit(x)
	if (length(x) > 0) {
		n <- round(length(x) * p)
		w <- seq(1, length(x) - n)
		x1 <- unname(sort(x))
		x2 <- rev(x1)
		int1 <- unname(rbind(x1[w + n], x1[w]))
		hpd1 <- sort(int1[, which.min(diff(int1))])
		int2 <- unname(rbind(x2[w + n], x2[w]))
		hpd2 <- sort(int2[, which.min(diff(int2))])
		hpd <- list(hpd1, hpd2)
		hpd <- hpd[[which.min(sapply(hpd, diff))]]	
	} else {
		hpd <- c(NA, NA)
	}
	return(hpd)
}



# Function:
#	get_collapsed_tree: collapses specified clades to single branches
# Arguments:
#	tree: object of class 'phylo'
#	imap: matrix or data frame (1st column individuals, 2nd column species) or a name of imap file
# Value:
#	tree object of class 'phylo' specifying a tree with a single tip per population and the terminal branch length equal to the mean distance between MRCA of the population and its tips in the original tree
# Details:
#	populations are required to be monophyletic in the tree, but they need not to cover the whole tree
# Depends on: ape

get_collapsed_tree <- function(tree, imap) {
	find_common_branch <- function(tip, phy) ifelse(length(tip) == 1, match(tip, phy$tip.label), ape::getMRCA(phy, tip))
	imap <- split(imap[,1], imap[,2])
	dst <- ape::dist.nodes(tree)
	anc <- lapply(imap, find_common_branch, phy=tree)
	tips <- lapply(imap, match, table=tree$tip.label)
	collapsed <- ape::drop.tip(tree, tip=unlist(lapply(imap, "[", -1)))
	ord <- match(sapply(imap, "[", 1), collapsed$tip.label)
	collapsed$tip.label[ord] <- names(imap)
	return(collapsed)
}


# Function:
#	get_topology: exports topology of the given tree in Newick format
# Arguments:
#	tree: 'phylo' object or a string with tree in Newick format
# Value:
#	character string with tree topology (no branch lengths) in Newick format
# Depends on: ape

get_topology <- function(tree) {
	if (!inherits(tree, "phylo")) {
		if (substr(tree, 1, 1) == "(") {
			tree <- ape::read.tree(text=tree)
		} else {
			tree <- ape::read.tree(file=tree)
		}
	}
	tree$edge.length <- NULL
	tree$node.label <- NULL
	ape::write.tree(tree)
}


# Function:
#	get_starting_tree: creates a starting tree for A11 or A01 analysis
# Arguments:
#	type: type of the starting tree, either 'nj' (neighbor joining) or 'random'
#	seq: a single alignment of sequences with row names corresponding to species (possibly concatenated multi-locus data)
#	imap: matrix or data frame (1st column individuals, 2nd column species) or a name of imap file
# Value:
#	starting tree topology in Newick format
# Depends on: ape

get_starting_tree <- function(type=c("nj","random"), seq, imap) {
	if (is.character(imap)) {
		imap <- read.table(imapfile, header=FALSE, col.names=c("Individual", "Species"))
	}
	if (type[1] == "nj") {
		dst <- avedist(seq, imap, model="K80", diag=FALSE)
		nan <- as.numeric(names(sort(table(which(is.nan(dst), arr.ind=TRUE)[,1]))))
		for (i in nan) {
			jj <- which(is.nan(dst[i,]))
			for (j in jj) {
				dst[i,j] <- dst[i,setdiff(which(dst[j,] == min(dst[j,-jj], na.rm=TRUE)),jj)[1]]
			}
		}
		tree <- phangorn::midpoint(ape::bionj(as.dist(dst)))
	}
	if (type[1] == "random") {
		tree <- ape::rtree(length(unique(tip.label)), tip.label=unique(imap[,2]))
	}
	tree$edge.length <- NULL
	tree$node.label <- NULL
	ape::write.tree(tree)
}


# Function:
#	get_maxsamplesize: maximum sample size (no. of sequences) of species at any locus
# Arguments:
#	seqfile: sequence file name
#	imap: matrix or data frame (1st column individuals, 2nd column species) or imap file name
# Value:
#	named numeric vector with maximum sample sizes

get_maxsamplesize <- function(seqfile, imap) {
	seqlines <- readLines(seqfile)
	seqlines <- seqlines[nchar(seqlines) > 0]
	start <- grep("^[[:digit:][:blank:]]+$", seqlines)
	loci <- rep(seq_along(start), diff(c(start, length(seqlines) + 1)))
	seqlines <- strsplit(seqlines, "\\s")
	nams <- sapply(seqlines, "[[", 1)
	zero <- grepl("^[N[:punct:]\\?]+$", sapply(seqlines, "[[", 2))
	nams <- lapply(split(nams[!zero], loci[!zero]), function(x) x[-1])
	nams <- lapply(nams, function(x) gsub("^.+\\^", "", x))
	if (is.character(imap)) {
		imap <- read.table(imapfile, header=FALSE, col.names=c("Individual", "Species"))
	}
	tab <- lapply(nams, function(x, imap) table(imap[match(x, imap[,1]),2]), imap=imap)
	mat <- matrix(,length(unique(imap[,2])), length(tab), dimnames=list(sort(unique(imap[,2])), names(tab)))
	for (i in seq_along(tab)) {
		mat[names(tab[[i]]),i] <- as.numeric(tab[[i]])
	}
	return(apply(mat, 1, max, na.rm=TRUE))
}



# Function:
#	estimate_otus: estimates operational taxonomic units using and maximum cross-section branch length criterion on ad hoc created average linkage (UPGMA) tree
# Arguments:
#	seq: list of alignments with row names indicating individuals or a single alignment of such sequences (possibly concatenated multi-locus data)
#	model: nucleotide substitution model indicated as in 'ape::dist.dna'
# Value:
#	data frame with columns 'Individuals' (individual names) and 'Species' (OTUs)
# Depends on: ape

estimate_otus <- function(seq, model="raw") {
	sc <- function(k, h, b) {
		mean(b[nh[,2] > k & nh[,1] <= k])
	}
	if (is.list(seq)) {
		bp <- cumsum(sapply(seq, ncol))
		part <- cbind(c(1, bp[-length(bp)] + 1), unname(bp))
		rows <- sort(unique(unlist(lapply(seq, rownames))))
		concatenated <- matrix(missing, length(rows), max(bp), dimnames=list(rows, NULL))
		for (i in seq_along(seq)) {
			concatenated[rownames(seq[[i]]), part[i,1]:part[i,2]] <- seq[[i]]
		}
		seq <- concatenated
	}
	dst <- ape::dist.dna(ape::as.DNAbin(seq), pairwise.deletion=TRUE, model=model, as.matrix=TRUE)
	nan <- as.numeric(names(sort(table(which(is.nan(dst), arr.ind=TRUE)[,1]))))
	for (i in nan) {
		jj <- which(is.nan(dst[i,]))
		for (j in jj) {
			dst[i,j] <- dst[i,setdiff(which(dst[j,] == min(dst[j,-jj], na.rm=TRUE)),jj)[1]]
		}
	}
	tre <- ape::as.phylo(stats::hclust(as.dist(dst), method="average"))
	ntip <- ape::Ntip(tre)
	nh <- matrix(ape::plotPhyloCoor(tre, direction="rightwards")[,1][tre$edge], ape::Nedge(tre), 2)
	knots <- sort(setNames(nh[,2], tre$edge[,2])[!duplicated(nh[,2])])
	knots <- setNames(c(0, knots), c(ntip+1, names(knots)))
	k <- knots[which.max(sapply(knots, sc, h=nh, b=tre$edge.length))]
	anc <- tre$edge[nh[,2] > k & nh[,1] <= k, 2]
	single <- unlist(Filter(function(x) length(x) > 0, as.list(anc[anc <= ntip])))
	otus <- lapply(setdiff(anc, single), function(a) ape::extract.clade(tre, a)$tip.label)
	otus <- c(as.list(tre$tip.label[single]), otus)
	otus <- data.frame(Individual=unlist(otus), Species=paste0("OTU", rep(seq_along(otus), sapply(otus, length))))
	return(otus)
}


# Function:
#	estimate_params: estimate putative means of theta and tau prior distributions
# Arguments:
#	seq: list of alignments with row names indicating individuals or a single alignment of such sequences (possibly concatenated multi-locus data)
#	imap: matrix or data frame (1st column individuals, 2nd column species) or a name of imap file
#	model: nucleotide substitution model indicated as in 'ape::dist.dna'
# Value:
#	vector of length two with putative means of tau and theta priors
# Depends on: ape

estimate_params <- function(seq, imap=NULL, model="raw") {
	if (is.list(seq)) {
		bp <- cumsum(sapply(seq, ncol))
		part <- cbind(c(1, bp[-length(bp)] + 1), unname(bp))
		rows <- sort(unique(unlist(lapply(seq, rownames))))
		concatenated <- matrix(missing, length(rows), max(bp), dimnames=list(rows, NULL))
		for (i in seq_along(seq)) {
			concatenated[rownames(seq[[i]]), part[i,1]:part[i,2]] <- seq[[i]]
		}
		seq <- concatenated
	}
	if (is.null(imap)) {
		otus <- setNames(as.list(rownames(seq)), rownames(seq))
	} else {
		if (is.character(imap)) {
			imap <- read.table(imapfile, header=FALSE, col.names=c("Individual", "Species"))
		}
		otus <- split(imap[,1], imap[,2])	
	}
	dst <- ape::dist.dna(ape::as.DNAbin(seq[unlist(otus),]), pairwise.deletion=TRUE, model=model, as.matrix=TRUE)
	dst[dst == 0] <- 1e-8
	tau <- matrix(,length(otus), length(otus))
	for (i in seq(length(otus) - 1)) {
		for (j in (i+1):length(otus)) {
			tau[i,j] <- mean(dst[otus[[i]], otus[[j]]], na.rm=TRUE)
		}
	}
	tau <- max(tau, na.rm=TRUE) / 2
	otus <- rep(seq_along(otus), sapply(otus, length))[match(rownames(dst), unlist(otus))]
	theta <- mean(dst[outer(otus, otus, "!=") & !diag(length(otus))], na.rm=TRUE)
	params <- c(tau=tau, theta=theta)
	return(params)
}



# Function:
#	avedist: calculates average genetic distance betwee the specified species and within them
# Arguments:
#	seq: list of alignments with row names indicating individuals or a single alignment of such sequences (possibly concatenated multi-locus data)
#	imap: matrix or data frame (1st column individuals, 2nd column species) or a name of imap file
#	model: nucleotide substitution model indicated as in 'ape::dist.dna'
#	diag: whether to include mean intra-specific distances (on a diagonal of the resulting matrix)
# Value:
#	square symmetric matrix giving mean genetic distance between species (off-diagonal) and within them (diagonal)
# Depends on: ape

avedist <- function(seq, imap, model="K80", diag=TRUE) {
	if (is.list(seq)) {
		bp <- cumsum(sapply(seq, ncol))
		part <- cbind(c(1, bp[-length(bp)] + 1), unname(bp))
		rows <- sort(unique(unlist(lapply(seq, rownames))))
		concatenated <- matrix(missing, length(rows), max(bp), dimnames=list(rows, NULL))
		for (i in seq_along(seq)) {
			concatenated[rownames(seq[[i]]), part[i,1]:part[i,2]] <- seq[[i]]
		}
		seq <- concatenated
	}
	if (!inherits(seq, "DNAbin")) {
		seq <- ape::as.DNAbin(seq)
	}
	if (is.character(imap)) {
		imap <- read.table(imapfile, header=FALSE, col.names=c("Individual", "Species"))
	}
	imap <- factor(imap[match(rownames(seq), imap[,1]),2])
	dst <- ape::dist.dna(seq, model=model, pairwise.deletion=TRUE, as.matrix=TRUE)
	ave <- matrix(, nlevels(imap), nlevels(imap), dimnames=list(levels(imap), levels(imap)))
	for (i in seq(nlevels(imap) - 1)) {
		ii <- imap == levels(imap)[i]
		for (j in (i + 1):nlevels(imap)) {
			jj <- imap == levels(imap)[j]
			ave[i,j] <- ave[j,i] <- mean(dst[ii,jj], na.rm=TRUE)
		}
	}
	if (isTRUE(diag)) {
		for (i in seq(nlevels(imap))) {
			ii <- imap == levels(imap)[i]
			if (sum(ii) > 1) {
				ave[i,i] <- mean(dst[ii,ii][lower.tri(diag(sum(ii)))], na.rm=TRUE)
			}
		}
	}
	return(ave)
}



# Function:
#	plot_delim_tree: 
# Arguments:
#	tree:	tree object of class 'phylo'
#	delim: data frame with delimited species (delim$Species) and their posterior probabilities (delim$PP)
#	sep:	separator of candidate species labels in delim$Species
#	clade: vector of node labels; if not tree$node.label (default) it requires node numbers as its names 
#	species: vector of edge labels; if not delim$PP (default) it requires node numbers as its names 
#	digits: the desired number of digits after the decimal point for display of posterior probabilities
#	edge.width:	position of clade names (negative values ~ shift to the left)
#	offset: offset of tip labels as a proportion of root to tip distance
#	tip.cex: size of tip label text
#	clade.pos:	node label position as a proportion of the supporting branches
#	clade.cex:	size of node label point
#	clade.lwd:	width of node label point border
#	clade.pch:	shape of node label
#	clade.col:	color of node label point border
#	clade.bg:	color of node label point background
#	species.pos:	edge label position as a proportion of the branch from its start (= root side)
#	species.cex: size of edge label point
#	species.lwd: width of edge label point border
#	species.pch:	shape of edge label
#	species.col:	color of edge label point border
#	species.bg: color of edge label point background
#	text.pos:	adjustment of node / edge label position (adj of 'text')
#	text.cex:	size of node / edge label
#	text.font:	font of node / edge label
#	text.col:	color of node / edge label
#	direction: tree direction (only 'rightwards' is implemented now)
#	device: "quartz", "x11" (or "X11") or "pdf", if NULL, the objects are plotted into the current device
#	file: name of pdf file if device == "pdf"
#	width:	width of graphical device
#	height:	height of graphical device
#	mai:	size of outer margins in inches, recycled if necessary
#	... other parameters of ape::plot.phylo
# Depends on: ape

plot_delim_tree <- function(tree, delim, sep="-", clade, species, digits=2, edge.width=3, offset=0.04, tip.cex=1, clade.pos=0.5, clade.cex=4, clade.lwd=3, clade.pch=21, clade.col="black", clade.bg="cornsilk", species.pos=1, species.cex=4, species.lwd=3, species.pch=21, species.col="black", species.bg="white", text.pos=c(0.5,0.5), text.cex=0.75, text.font=2, text.col="black", direction, device, file, width=7, height=7, mai=0.2, ...) {

	getpaths <- function(phy, from, to) lapply(to, function(x) ape::nodepath(phy, from=from, to=x))
	find_common_branch <- function(tip, phy) ifelse(length(tip) == 1, match(tip, phy$tip.label), ape::getMRCA(phy, tip))

	ntip <- ape::Ntip(tree)
	if (missing(species) & !missing(delim)) {
		sp <- strsplit(delim$Species, split=sep)
		mrcas <- sapply(sp, find_common_branch, phy=tree)
		species <- setNames(formatC(delim$PP, format="f", digits=digits), mrcas)
	} else if (!missing(species)) {
		mrcas <- as.integer(names(species))
	} else {
		mrcas <- seq(ape::Nnode(tree)) + ntip
	}
	if (missing(clade)) {
		anc <- setdiff(sort(unique(unlist(getpaths(tree, ntip+1, mrcas)))), seq(ntip))[-1]
		clade <- setNames(formatC(as.numeric(tree$node.label), format="f", digits=digits)[anc - ntip], anc)	
	} else {
		anc <- as.integer(names(clade))
	}
	
	if (missing(device)) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
	}
	if (isTRUE(device == "pdf")) {
		pdf(ifelse(missing(file), "bpp.pdf", file), width=width, height=height)
	} else if (!is.null(device)) {
		match.fun(tolower(device))(width=width, height=height)	
	}
	par(mai=rep_len(mai, 4))
	offset <- offset * diff(range(ape::plotPhyloCoor(tree, direction="rightwards")[,"xx"]))
	plot(tree, edge.col=species.col, edge.width=edge.width, type="phylogram", direction="rightwards", show.tip.label=TRUE, cex=tip.cex, label.offset=offset, ...)
	xy <- ape::plotPhyloCoor(tree, direction="rightwards")
	cladexy <- xy[anc,]
	cladexy[,1] <- cladexy[,1] - rep(1 - clade.pos, length(anc)) * tree$edge.length[match(anc, tree$edge[,2])]
	points(cladexy, pch=clade.pch, cex=clade.cex, col=clade.col, bg=clade.bg, lwd=clade.lwd)
	text(cladexy, labels=clade, cex=text.cex, font=text.font, col=text.col, adj=text.pos)
	if (!missing(species)) {
		speciesxy <- xy[mrcas,]
		speciesxy[,1] <- speciesxy[,1] - rep(1 - species.pos, length(mrcas)) * tree$edge.length[match(mrcas, tree$edge[,2])]
		points(speciesxy, cex=species.cex, lwd=species.lwd, pch=species.pch, col=species.col, bg=species.bg)
		text(speciesxy, labels= species, cex=text.cex, font=text.font, col=text.col, adj=text.pos)		
	}
	if (isTRUE(device == "pdf")) {
		dev.off()
	}
	
}



# Function:
#	plot_msc_tree: plots the species tree as composed of boxes whose height ~ branch length and width = half of population size parameter (theta)
# Arguments:
#	tree: tree imported by read_msc_tree
#	asp: the y/x aspect ratio, default is 1 (branches in scale), but a meaningful alternative is "mtheta", which indicates ratio of mean branch width to total tree depth, so the tree as a whole is approximately in scale
#	gap: minimum gap between branches, specified as a proportion of maximum branch width, can be automatically adjusted
#	col: color of branch background
#	border: color of branch outline
#	lwd: width of branch outline
#	cex: size of axis labels (cex.lab), serves also for derivation of other expansion factors (e.g. cex.axis)
#	show.tip.label: whether to show tip labels
#	tip.label: vector of tip labels, used if the desired tip labels differ from tree$tip.label, the vector is named by tree$tip.label, however, to establish mutual correspondence of the labels
#	tip.cex:	size of tip labels
#	tip.font:	font of tip labels
#	timescale: whether to display time scale axis
#	axis.lwd: width of axis line
#	labels: what is displayed at branch labels, options are 'theta' (default), 'tau' or 'gdi'
#	digits: precision (no. of decimal places) for parameters to be displayed at node labels
#	lab.col: color of branch label background
#	lab.border: color of branch label outline
#	lab.lwd: width of branch label outline
#	lab.cex: size of branch label
#	lab.font: font of branch label
#	xlab, ylab: labels of x and y axes
#	xlim, ylim: limits of x and y axes
#	mai:	size of outer margins in inches, recycled if necessary
#	direction: tree direction (only 'rightwards' is implemented now)
#	device: "quartz", "x11" (or "X11") or "pdf", if NULL, the objects are plotted into the current device
#	file: name of pdf file if device == "pdf"
#	width:	width of graphical device
#	height:	height of graphical device

plot_msc_tree <- function(tree, asp=1, gap=0.1, col="grey", border="black", lwd=3, cex=1.5, show.tip.label=TRUE, tip.label=NULL, tip.cex=cex, tip.font=3, timescale=TRUE, axis.lwd=lwd/2, labels="theta", digits=4, lab.col="white", lab.border="black", lab.lwd=lwd/2, lab.cex=1, lab.font=1, xlab="Tau", ylab="Theta", xlim, ylim, mai, direction, device, file, width=7, height=7) {
	branch_coords <- function(tree, params, gap) {
		ntip <- length(tree$tip.label)
		xx <- matrix(params[tree$edge,"tau"], nrow(tree$edge), 2)
		xx <- -1 * (xx - max(xx))
		yy <- cumsum(params[seq(ntip),"theta"]) - 0.5 * params[seq(ntip),"theta"]
		yy <- yy + cumsum(rep(gap, ntip))
		yy <- yy[tree$edge[,2]]
		while (any(is.na(yy))) {
			a <- tree$edge[!is.na(yy),1]
			a <- max(intersect(a[which(duplicated(a))], tree$edge[is.na(yy),2]))
			yy[tree$edge[,2] == a] <- mean(yy[tree$edge[,1] == a])
		}
		yy <- cbind(yy - 0.5 * params[tree$edge[,2],"theta"], yy + 0.5 * params[tree$edge[,2],"theta"])
		return(list(xx=xx, yy=yy))
	}
	drawEllipse <- function(x, y, size, w, h, r, res, col, border, lwd, return=FALSE) {
		rs <- seq(0, 2 * pi, len=res)
		pts <- cbind(0.5 * w * cos(rs), 0.5 * h * sin(rs))
		rot <- matrix(c(cos(r), sin(r), -sin(r), cos(r)), 2, 2)
		pts <- size * pts %*% rot + rep(1, res) %*% t(c(x,y))
		if (isTRUE(return)) {
			return(pts)
		} else {
			polygon(pts, col=col, border=border, lwd=lwd)
		}
	}
	rotxy <- function(xy, r) {
		mat <- matrix(xy, 2, 2)
		cen <- c(1, 1) %*% t(colMeans(mat))
		rot <- matrix(c(cos(r), sin(r), -sin(r), cos(r)), 2, 2)
		return(as.numeric((mat - cen) %*% rot + cen))
	}
	
	if (is.character(tree) & length(tree) == 1) {
		tree <- read_msc_tree(tree)
	}
	ntip <- length(tree$tip.label)
	nodes <- seq(ntip - 1) + ntip
	params <- tree$params
	params[,"theta"] <- 0.5 * params[,"theta"]

	gapsize <- gap
	gap <- gapsize * max(params[,"theta"])
	bc <- branch_coords(tree, params, gap)
	ry <- mean(sort(bc$yy[tree$edge[,1]==(ntip+1),])[2:3])
	ec <- drawEllipse(x=0, y=ry, size=1, w=1.5*gap, h=params[ntip+1,"theta"], r=0, res=200, return=TRUE)
	for (i in nodes) {
		ii <- which(tree$edge[,1] == i)
		igap <- diff(sort(bc$yy[ii,])[2:3])
		if (igap < 0) {
			gap <- gap + 1.25 * abs(igap)
		}
	}
	ii <- which(tree$edge[,1] == (ntip + 1))
	igap <- diff(c(sort(bc$yy[ii,])[2], min(ec[,2])))
	if (igap < 0) {
		gap <- gap + 1.25 * abs(igap)
	}
	if (gap > gapsize * max(params[,"theta"])) {
		bc <- branch_coords(tree, params, gap)
	}

	if (isTRUE(show.tip.label)) {
		terminal <- match(seq(ntip), tree$edge[,2])
		if (is.null(tip.label)) {
			tiplab <- gsub("_", " ", tree$tip.label)
		} else {
			tiplab <- gsub("_", " ", tip.label[tree$tip.label])
		}
	}
	if (missing(xlim)) {
		xlim <- range(c(bc$xx, ec[,1]))
	}
	if (missing(ylim)) {
		ylim <- range(bc$yy)
	}

	if (missing(device)) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
	}
	if (isTRUE(device == "pdf") & missing(file)) {
		file <- "msc_tree.pdf"
	}
	if (isTRUE(device == "pdf")) {
		pdf(file, width=width, height=height)
	} else if (!is.null(device)) {
		match.fun(tolower(device))(width=width, height=height)	
	}

	if (missing(mai)) {
		mai <- c(0, 0, 0.02, 0)
		if (isTRUE(timescale)) {
			cex.axis <- ifelse(cex > 1, 1 + (cex-1) / 2, 1)
			strh <- strheight(s=0, units="inches", cex=cex.axis)
			mai[1] <- 0.02 * min(c(width, height)) + strh
		}
		if (nchar(xlab) > 0) {
			strh <- strheight(s=xlab, units="inches", cex=cex)
			maix <- mai[1] + 1.5 * strh
			mai[1] <- mai[1] + 3.5 * strh
		}
		if (nchar(ylab) > 0) {
			strh <- strheight(s=ylab, units="inches", cex=cex)
			maiy <- 1.5 * strh
			mai[2] <- 4 * strh
		}
		if (isTRUE(show.tip.label)) {
			strw <- max(strwidth(s=tiplab, units="inches", cex=cex, font=tip.font))
			mai[4] <- strw + 0.01 * width
		}
	}
	if (asp == "mtheta") {
		asp <- mean(params[,"theta"]) / max(params[,"tau"])
	}

	par(mai=mai)
	plot(0, xlim=xlim, ylim=ylim, asp=asp, bty="n", xaxt="n", yaxt="n", xlab="", ylab="")
	for (i in nodes) {
		ii <- which(tree$edge[,1] == i)
		lines(bc$xx[ii,1], bc$yy[ii,2], lwd=lwd, col=border)
	}
	for (i in seq(nrow(tree$edge))) {
		rect(bc$xx[i,1], bc$yy[i,1], bc$xx[i,2], bc$yy[i,2], col=col, border=border, lwd=lwd)
	}
	ry <- mean(sort(bc$yy[tree$edge[,1]==(ntip+1),])[2:3])
	drawEllipse(x=0, y=ry, size=1, w=1.5*gap, h=params[ntip+1,"theta"], r=0, res=200, col=col, border=border, lwd=lwd)

	if (isTRUE(timescale)) {
		tau <- pretty(xlim)
		tau <- tau[tau <= xlim[2]]
		at <- tau - (max(tau) - xlim[2])
		tau <- rev(tau - min(tau))
		axis(side=1, at=at, labels=tau, lwd=axis.lwd, cex.axis=cex.axis)
	}
	if (nchar(xlab) > 0) {
		mtext(text=xlab, side=1, at=mean(xlim), cex=cex, line=maix/0.2)
	}
	if (nchar(ylab) > 0) {
		mtext(text=ylab, side=2, at=mean(ylim), cex=cex, line=maiy/0.2)
	}
	if (isTRUE(show.tip.label)) {
		mtext(text=tiplab, side=4, at=rowMeans(bc$yy[terminal,]), cex=tip.cex, font=tip.font, las=2, line=(0.005*width)/0.2)
	}

	if (labels %in% colnames(params)) {
		labtext <- formatC(params[,labels], format="f", digits=digits)
		ord <- match(c(seq(ntip), nodes), tree$edge[,2])
		labxy <- cbind(rowMeans(bc$xx), rowMeans(bc$yy))[ord,]
		labxy[ntip + 1,] <- c(0, ry)
		limits <- t(rbind(diff(t(bc$xx)), diff(t(bc$yy))))[ord,]
		limits[ntip + 1,] <- diff(apply(ec, 2, range))		
	}
	if (labels %in% c("tau","theta")) {
		strw <- strwidth(s=labtext, units="user", cex=1, font=lab.font) * lab.cex
		strh <- strheight(s=labtext, units="user", cex=1, font=lab.font) * lab.cex
		m <- t(c(-1, 1))
		rectxy <- labxy[,c(1,1,2,2)] + cbind(t(t(strw/2)) %*% (0.75 * m), t(t(strh/2)) %*% (1.5 * m))
		srt <- ifelse(diff(t(rectxy[,1:2])) < limits[,1], 0, 90)
		for (i in which(srt == 90)) {
			rectxy[i,] <- rotxy(rectxy[i,], r=pi/2)
		}
		if (labels == "tau") {
			ii <- nodes
		} else if (labels == "theta") {
			ii <- c(seq(ntip), nodes)
		}
		for (i in ii) {
			rect(rectxy[i,1], rectxy[i,3], rectxy[i,2], rectxy[i,4], col=lab.col, border=lab.border, lwd=lab.lwd)
			text(labxy[i,1], labxy[i,2], labels=labtext[i], cex=lab.cex, font=lab.font, srt=srt[i])
		}		
	} else if (labels %in% "gdi") {
		strw <- strwidth(s=labtext, units="user", cex=1, font=lab.font) * lab.cex
		strh <- strheight(s=labtext, units="user", cex=1, font=lab.font) * lab.cex
		ii <- c(seq(ntip), nodes[-1])
		for (i in ii) {
			drawEllipse(labxy[i,1], labxy[i,2], w=0.75*strw[i], h=1.75*strh[i], r=0, size=1, col=lab.col, border=lab.border, lwd=lab.lwd, res=200)
			text(labxy[i,1], labxy[i,2], labels=labtext[i], cex=lab.cex, font=lab.font)
		}
	}

	if (isTRUE(device == "pdf")) {
		dev.off()
	}
}



# Function:
#	plot_pdensity: plots probability density (prior and/or posterior) of MSC parameter(s)
# Arguments:
#	prior: vector two parameters (shape and rate) of (inverse) gamma distribution or a list of these vectors
#	posterior: data frame with posterior samples of multispecies coalescent parameters (e.g., 'params' data frame of 'bpp' object created by 'read_bpp_output')
#	distribution: name of prior distribution(s), recycled if necessary; may be 'invgamma' (default) or 'gamma'
#	qmax: maximum quantile of the distributions to display
#	points: whether to display quantiles of the prior as points
#	q: quantiles of the prior to be indicated by points
#	xlab, ylab: optional, x and y axis labels
#	col: colors for display of posterior and prior, can be a list of vectors
#	lty: line type, recycled if necessary, can be a list of vectors
#	lwd: line width, recycled if necessary, can be a list of vectors
#	cex: point size
#	legend: whether to display legend
#	expr: expression allowing to add other elements
# Depends on: invgamma

plot_pdensity <- function(prior, posterior, distribution="invgamma", qmax=0.99, points=FALSE, q=c(0.025, 0.5, 0.975), xlab="", ylab="density", col=c(1,8), lty=c(3,1), lwd=2, cex=1.5, legend=TRUE, expr) {
	if (!missing(posterior)) {
		if (is.numeric(posterior)) {
			posterior <- data.frame(posterior)
		}
		dposterior <- lapply(seq(length(posterior)), function(i) stats::density(posterior[,i], cut=0))
	} else {
		dposterior <- NULL
	}
	if (!missing(prior)) {
		if (is.data.frame(prior)) {
			prior <- apply(prior[,1:2], 1, as.numeric, simplify=FALSE)
		}
		if (!is.list(prior)) {
			prior <- list(prior)
		}
		distribution <- rep_len(distribution, length(prior))
		qdistr <- list(gamma=stats::qgamma, invgamma=invgamma::qinvgamma)[distribution]
		ddistr <- list(gamma=stats::dgamma, invgamma=invgamma::dinvgamma)[distribution]
		from <- lapply(seq_along(prior), function(i) qdistr[[i]](0.00001, prior[[i]][1], prior[[i]][2]))
		to <- lapply(seq_along(prior), function(i) qdistr[[i]](qmax, prior[[i]][1], prior[[i]][2]))
		xprior <- lapply(seq_along(prior), function(i) seq(from[[i]], to[[i]], length=101))
		dprior <- lapply(seq_along(prior), function(i) ddistr[[i]](xprior[[i]], prior[[i]][1], prior[[i]][2]))
		if (isFALSE(points)) {
			qprior <- NULL
		} else {
			qprior <- lapply(seq_along(prior), function(i) qdistr[[i]](q, prior[[i]][1], prior[[i]][2]))
		}
	} else {
		xprior <- NULL
		dprior <- NULL
		qprior <- NULL
	}
	ymax <- max(c(unlist(dprior), unlist(lapply(dposterior, "[[", "y"))))
	xlim <- range(c(unlist(xprior), unlist(lapply(dposterior, "[[", "x"))))
	pp <- c(length(dprior) > 0, length(dposterior) > 0)
	if (all(pp)) {
		col <- as.list(rep_len(col, 2))
		col[1] <- list(rep_len(col[[1]], length(dprior)))
		col[2] <- list(rep_len(col[[2]], length(dposterior)))	
		lty <- as.list(rep_len(lty, 2))
		lty[1] <- list(rep_len(lty[[1]], length(dprior)))
		lty[2] <- list(rep_len(lty[[2]], length(dposterior)))	
		lwd <- as.list(rep_len(lwd, 2))
		lwd[1] <- list(rep_len(lwd[[1]], length(dprior)))
		lwd[2] <- list(rep_len(lwd[[2]], length(dposterior)))	
	} else if (isTRUE(pp[1])) {
		col <- list(rep_len(col, length(dprior)), NULL)
		lty <- list(rep_len(lty, length(dprior)), NULL)
		lwd <- list(rep_len(lwd, length(dprior)), NULL)
	} else if (isTRUE(pp[2])) {
		col <- list(NULL, rep_len(col, length(dposterior)))		
		lty <- list(NULL, rep_len(lty, length(dposterior)))		
		lwd <- list(NULL, rep_len(lwd, length(dposterior)))		
	} else {
		col <- vector("list", 2)
	}

	par(mai=c(1.02, 1.02, 0.42, 0.42))
	plot(0, 0, xlim=xlim, ylim=c(0, ymax), type="n", xlab=xlab, ylab=ylab, cex.lab=1.5, cex.axis=1.25)
	if (!is.null(dprior)) {
		for (i in seq_along(dprior)) {
			lines(xprior[[i]], dprior[[i]], col=col[[1]][i], lty=lty[[1]][i], lwd=lwd[[1]][i])
		}
	}
	if (!is.null(dposterior)) {
		for (i in seq_along(dposterior)) {
			lines(dposterior[[i]]$x, dposterior[[i]]$y, col=col[[2]][i], lty=lty[[2]][i], lwd=lwd[[2]][i])		
		}
	}
	if (!is.null(dprior) & !is.null(qprior)) {
		for (i in seq_along(qprior)) {
			idprior <- ddistr[[i]](qprior[[i]], prior[[i]][1], prior[[i]][2])
			points(qprior[[i]], idprior, pch=21, cex=cex, col=col[[1]][i], lwd=lwd[1], bg="white")
		}
	}

	if (isTRUE(legend) & all(pp)) {
		leg <- c("Prior","Posterior")
	} else if (is.character(legend)) {
		leg <- legend
	} else {
		leg <- NULL
	}
	if (!is.null(leg)) {
		if (all(pp)) {
			col <- sapply(col, "[[", 1)
			lty <- sapply(lty, "[[", 1)
			lwd <- sapply(lwd, "[[", 1)
		} else if (any(pp)) {
			col <- unlist(col[!sapply(col, is.null)])
			lty <- unlist(lty[!sapply(lty, is.null)])
			lwd <- unlist(lwd[!sapply(lwd, is.null)])
		} else {
			leg <- NULL
		}
	}
	if (!is.null(leg)) {
		legend("topright", legend=leg, col=col, lwd=lwd, lty=lty, bty="n")		
	}

	if (!missing(expr)) {
		if (is.expression(expr)) {
			eval(expr)
		}
	}

}



# Function:
#	plot_gdi: plots gdi posterior distributions
# Arguments:
#	gdis: object of class 'gdi' produced by gdi
#	which: numeric, which species to display
#	hpd: proportion of posterior sample defining the highest probability density interval (suppressed by NA)
#	col: color of line & point outline
#	bg: color of point background
#	lwd: line width
#	cex: point size
#	points: whether to display points
#	xlab, ylab: x & y axis labels
#	mai: size of outer margins in inches (cf. par("mai"))
#	cex.lab: size of axis labels
#	cex.axis: size of axis annotations

plot_gdi <- function(gdis, which=1, hpd=0.95, col=1, bg=8, lwd=4, cex=2, points=TRUE, xlab="gdi", ylab="density", mai, cex.lab=1.5, cex.axis=1.25) {
	d <- stats::density(gdis$distributions[,which], from=0, to=1)
	m <- which.min(abs(d$x - gdis$means[which]))
	if (!is.na(hpd)) {
		limits <- HPD(gdis$distributions[,which], p=hpd)
		low <- which.min(abs(d$x - limits[1]))
		upp <- which.min(abs(d$x - limits[2]))	
	} else {
		low <- NULL
		upp <- NULL
	}
	ymax <- max(d$y)
	if (missing(mai)) {
		mai <- c(1.02, 0.92, 0.32, 0.32)
	}
	par(mai=mai)
	plot(0, 0, xlim=c(0,1), ylim=c(0,ymax), type="n", xlab=xlab, ylab=ylab, cex.lab=cex.lab, cex.axis=cex.axis)
	lines(d$x, d$y, col=col, lwd=lwd)
	points(d$x[c(low,m,upp)], d$y[c(low,m,upp)], pch=21, cex=cex, col=col, bg=bg, lwd=lwd/2)
}


