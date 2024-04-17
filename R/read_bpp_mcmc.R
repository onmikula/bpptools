#' Reading mcmc file.
#' 
#' @description
#' Reads MCMC trace into R.
#' 
#' @param mcmcfile A character string giving the name of .mcmc file.
#' @param thinning A proportion of mcmc steps to be retained.
#' @param tree A target tree upon which divergence times from A00 analysis are projected;
#'   in newick format, file name or an object of class `phylo`.
#' @returns The list with components `"trees"` (an object of class `multiPhylo`),
#'   `"params"` (data frame) and `"nspecies"` (numeric vector), containing posterior samples 
#'   of trees, multispecies coalescent parameters and species counts, respectively.
#'   The component `nspecies` can contain posterior sample of species counts (A11, A10 analyses)
#'   or a fixed number of species (A01, A00 analyses).
#' @export

read_bpp_mcmc <- function(mcmcfile, thinning=1, tree=NULL) {
	A00 <- substr(readLines(mcmcfile, n=1), 1, 1) != "("
	if (isTRUE(A00)) {
		params <- read.delim(mcmcfile, row.names=1)
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
			#tree <- ape::reorder.phylo(tree, order="postorder")			
			taus <- lapply(split(params[,grepl("tau", colnames(params))], seq(nrow(params))), as.numeric)
			trees <- lapply(taus, settaus, phy=tree)		
			class(trees) <- "multiPhylo"
		} else {
			trees <- NULL
			warning("trees are not returned as their topology is not specified")
		}
	} else {
		steps <- readLines(mcmcfile)
		steps <- steps[seq(1, length(steps), floor(1/thinning))]
		steps <- strsplit(steps, split=";")
		nspecies <- as.numeric(sapply(steps, "[", 2))
		trees <- paste(sapply(steps, "[", 1), ";", sep="")
		trees <- lapply(trees, function(x) ape::read.tree(text=x))
#		trees <- lapply(trees, ape::reorder.phylo, order="postorder")
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


