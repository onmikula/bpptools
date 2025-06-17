#' Reading mcmc file.
#' 
#' @description
#' Reads MCMC trace into R.
#' 
#' @param file character string, the name of .mcmc file.
#' @param thinning proportion of mcmc steps to be retained or a subsampling frequency (if `thin > 1`).
#' @param topology a target tree upon which divergence times from A00 analysis are projected;
#'   in newick format, file name or an object of class `phylo`.
#' @returns The list with components `"trees"` (an object of class `multiPhylo`),
#'   `"params"` (data frame) and `"nspecies"` (numeric vector), containing posterior samples 
#'   of trees, multispecies coalescent parameters and species counts, respectively.
#'   The component `nspecies` can contain posterior sample of species counts (A11, A10 analyses)
#'   or a fixed number of species (A01, A00 analyses).
#' @export

read_bpp_mcmc <- function(file, thinning=1, topology=NULL, mode) {
	first <- gsub("^[[:blank:]]+|[[:blank:]]+$", "", readLines(file, n=5))
	first <- first[nchar(first) > 0][1]
	if (substr(first, 1, 1) != "(") {
		if (grepl("np[[:blank:]]tree", first)) {
			mode <- "A10"
		} else {
			mode <- "A00"
		}
	} else {
		if (grepl("[[:digit:]]$", first)) {
			mode <- "A11"
		} else {
			mode <- "A01"
		}
	}
	if (mode == "A00") {
		params <- read.delim(file, row.names=1)
		if (thinning != 1) {
			thinning <- floor(ifelse(thinning < 1, 1 / thinning, thinning))
			thin <- as.integer(seq(0, nrow(params), by=thinning))[-1]
			params <- params[thin,]
		}
		nspecies <- sum(grepl("^tau", names(params))) + 1
		if (!is.null(topology)) {
			if (is.character(topology)) {
				topology <- import_tree(topology)
			}
			trees <- apply(params[,grepl("tau", colnames(params)), drop=FALSE], 1, settaus, phy=topology)
			class(trees) <- "multiPhylo"
		} else {
			trees <- NULL
			warning("trees are not returned as their topology is not specified")
		}
	} else if (mode == "A10") {
		params <- read.delim(file, row.names=1)
		if (thinning != 1) {
			thinning <- floor(ifelse(thinning < 1, 1 / thinning, thinning))
			thin <- as.integer(seq(0, nrow(params), by=thinning))[-1]
			params <- params[thin,]
		}
		nspecies <- sapply(gregexpr("1", as.character(params[,"tree"])), length) + 1
		params <- params[,c(setdiff(colnames(params), c("np","tree","lnL")), c("lnL","tree","np"))]
		if (is.null(topology)) {
			stop("the tree topology must be supplied")
		}
		topology <- import_tree(topology)
		trees <- apply(params[,grepl("tau", colnames(params)),drop=FALSE], 1, settaus, phy=topology)
		class(trees) <- "multiPhylo"
	} else if (mode == "A01") {
		trees <- gsub("[[:blank:]]", "", readLines(file))
		if (thinning != 1) {
			thinning <- floor(ifelse(thinning < 1, 1 / thinning, thinning))
			thin <- as.integer(seq(0, length(trees), by=thinning))[-1]
			trees <- trees[thin]
		}
		trees <- ape::read.tree(text=trees)
		if (inherits(trees, "phylo")) {
			trees <- list(trees)
			class(trees) <- "multiPhylo"
		}
		nspecies <- ape::Ntip(trees[[1]])
		params <- vector("list", length(trees))
		for (i in seq_along(trees)) {
			brtimes <- ape::branching.times(trees[[i]])
			nodeord <- order(match(names(brtimes), trees[[i]]$edge[,1]))
			params[[i]] <- data.frame(tau=c(rep(0, nspecies), brtimes[nodeord]))
			thetas <- c(trees[[i]]$tip.label, trees[[i]]$node.label)
			thetas <- suppressWarnings(as.numeric(gsub("[[:alnum:]]*#", "", thetas)))
			tipord <- order(match(seq(nspecies), trees[[i]]$edge[,2]))
			params[[i]]$theta <- thetas[c(tipord, nodeord + nspecies)]
			trees[[i]]$tip.label <- gsub("#[[:alnum:][:punct:]]*$", "", trees[[i]]$tip.label)
			trees[[i]]$node.label <- NULL
			rownames(params[[i]])[seq(nspecies)] <- trees[[i]]$tip.label[tipord]
		}
	} else if (mode == "A11") {
		trees <- readLines(file)
		trees <- gsub("[[:blank:]]", "", trees)
		if (thinning != 1) {
			thinning <- floor(ifelse(thinning < 1, 1 / thinning, thinning))
			thin <- as.integer(seq(0, length(trees), by=thinning))[-1]
			trees <- trees[thin]
		}
		nspecies <- as.numeric(regmatches(trees, regexpr("[[:digit:]]+$", trees)))
		trees <- ape::read.tree(text=sub("[[:digit:]]+$", "", trees))
		ntip <- ape::Ntip(trees[[1]])
		params <- vector("list", length(trees))
		for (i in seq_along(trees)) {
			brtimes <- ape::branching.times(trees[[i]])
			nodeord <- order(match(names(brtimes), trees[[i]]$edge[,1]))
			params[[i]] <- data.frame(tau=c(rep(0, ntip), brtimes[nodeord]))
			thetas <- c(trees[[i]]$tip.label, trees[[i]]$node.label)
			thetas <- suppressWarnings(as.numeric(gsub("[[:alnum:]]*#", "", thetas)))
			tipord <- order(match(seq(ntip), trees[[i]]$edge[,2]))
			params[[i]]$theta <- thetas[c(tipord, nodeord + ntip)]
			trees[[i]]$tip.label <- gsub("#[[:alnum:][:punct:]]*$", "", trees[[i]]$tip.label)
			trees[[i]]$node.label <- NULL
			rownames(params[[i]])[seq(ntip)] <- trees[[i]]$tip.label[tipord]
		}
	}
	result <- list(trees=trees, params=params, nspecies=nspecies)
	class(result) <- "bpp"
	return(result)
}



#' Combine BPP runs.
#' 
#' @description
#' Combine MCMC traces from independent replicates of the same analysis in BPP.
#' 
#' @param runs The list of lists returned by `read_bpp_mcmc`.
#' @returns The list with the same components as the input objects,
#'   containing pooled posterior samples.
#' @export

combine_bpp_runs <- function(runs) {
	trees <- Reduce(c, lapply(runs, "[[", "trees"))
	params <- lapply(runs, "[[", "params")
	if (is.null(dim(params[[1]]))) {
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



#' Exporting BPP output.
#' 
#' @description
#' Exports mcmc trace of BPP into into '.trees' file (in newick format,
#' readable by FigTree, ape::read.tree etc.) and '.log' file (tab-delimited, readable by Tracer).
#' 
#' @param bpp The list returned `read_bpp_mcmc`.
#' @param trees The name of 'trees' file.
#' @param log The name of 'log' file.
#' @param model One of these: "A00", "A01", "A10" or "A11".
#' @export

export_bpp_mcmc <- function(bpp, trees=NULL, log=NULL, model="A11") {
	spdelimitation <- substring(model, 2, 2) == "1"
	sptree <- substring(model, 3, 3) == "1"
	if (sptree) {
		tau0 <- sapply(lapply(bpp$params, "[[", "tau"), max)
		mtheta <- sapply(lapply(bpp$params, "[[", "theta"), mean, na.rm=TRUE)
		params <- data.frame(tau0=tau0, mtheta=mtheta)
		if (spdelimitation) {
			params$nspecies <- bpp$nspecies
		}
	} else {
		if (spdelimitation) {
			ntip <- ape::Ntip(bpp$trees[[1]])
			nnode <- ape::Nnode(bpp$trees[[1]])
			tau <- do.call(rbind, lapply(bpp$params, "[[", "tau"))[,-seq(ntip)]
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



#' Importing BPP mcmc trace.
#' 
#' @description
#' Imports mcmc trace of BPP from files where its components are stored.
#' 
#' @param trees The name of "trees" file or an object of class `multiPhylo`.
#' @param log The name of "log" file or a data frame.
#' @export

as.bpp <- function(trees, log) {
	if (!(inherits(trees, "multiPhylo") | inherits(trees[[1]], "phylo"))) {
		trees <- import_tree(trees)
	}
	if (!missing(log)) {
		if (!is.data.frame(log)) {
			params <- read.delim(log)
		} else {
			params <- log
		}
	} else {
		params <- NULL
		warning("'log' argument was not specified, only trees are imported")
	}
	nspecies <- apply(speciesDA(trees), 1, max)
	result <- list(trees=trees, params=params, nspecies=nspecies)
	class(result) <- "bpp"
	return(result)
}




settaus <- function(taus, phy) {
	heights <- ifelse(phy$edge <= ape::Ntip(phy), 0, phy$edge - ape::Ntip(phy))
	heights[heights != 0] <- taus[heights[heights != 0]]
	phy$edge.length <- heights[,1] - heights[,2]
	return(phy)
}
