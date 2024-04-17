#' Species delimitation.
#' 
#' @description
#' Finds species delimitation representative for a posterior sample.
#'
#' `find_spdelim` a wrapper function joining the following functions into a single pipeline.
#'
#' `speciesDA` takes a posterior sample of trees and extracts a sample of species delimitations out of it.
#'   It returns a numeric matrix with rows for the sampled delimitation and columns for the tree tips.
#'   The tips with the same number in given row are conspecific in that sampled delimitation.
#'
#' `speciesPP` takes an output of `speciesDA` and converts it to a list with components
#'   `"delims"` (data frame with delimitations and their PPs), `"species"` (data frame with species and PPs)
#'   and `"tips"` (candidate species ordered as in column names of `speciesDA` output). 
#'
#' `speciesM50`, `speciesMAP`, `speciesMC` functions choosing majority consensus, the maximum a posteriori
#'   and the maximum credibility species delimtation, respectively. They all return a data frame with species
#'   in the 1st column and their PPs in the 2nd column and with an attribute `"pp"`, which is PP
#'   of the delimitation as a whole
#'
#' @param trees A posterior sample of rooted trees in a `multiPhylo` object or a list of `phylo` objects
#'   or a name of file with the posterior sample.
#' @param file A name of the output file.
#' @param type The type of delimitation, either `"M50"` (majority consensus, default),
#'   `"MAP"` (maximum aposteriori) or `"MC"` (maximum credibility), see Details for explanantion.
#' @param sep Symbol used to concatenate candidate species names.
#' @param collapse The divergence time (in appropriate units) approximating zero. In analyses based
#'   on the rjMCMC algorithm (used by BPP) it is zero. In analyses based on birth-death-collapse prior
#'   (in BEAST2) it is a small value identical to "collapse height" set in the analysis. 
#' @details The majority consensus delimitation adds the most supported species unless they overlap
#'   (they are in conflict) with any of the species that were already included. The maximum a posteriori
#'   delimitation pick the set of species which was most frequently sampled as a whole. The maximum credibility
#'   delimitation scores every sampled delimitation by the product of PPs of its constituent species.
#' @returns A data frame with names of species in the chosen delimitation and their posterior probabilities
#'   and a tab-delimited file with this output (if `file` is specified).
#' @export

find_spdelim <- function(trees, file, type=c("M50","MAP","MC"), sep="-", collapse=0) {
	if (is.character(trees)) {
		trees <- ape::read.tree(trees)
	}
	delim_mat <- speciesDA(trees, collapse=collapse)
	delim_pp <- speciesPP(delim_mat, sep=sep)
	delim_fun <- match.fun(paste0("species", type[1]))
	delim_sp <- delim_fun(delim_pp, sep=sep)
	if (!missing(file)) {
		write.table(delim_sp, file, row.names=FALSE, quote=FALSE, sep="\t")
	}
	return(delim_sp)
}



speciesDA <- function(trees, collapse=0) {
	endpoints <- function(paths, lengths, collapse) {
		pathlengths <- lapply(lapply(lapply(paths, function(x) lengths[x]), rev), cumsum)
		ends <- sapply(pathlengths, function(x) sum(x <= collapse, na.rm=TRUE)) + 1
		return(sort(unique(mapply("[", lapply(paths, rev), ends))))
	}
	offspring <- function(phy, node) {
		tips <- phy$tip.label
		species <- lapply(node, function(n) if (n <= length(tips)) return(tips[n]) else return(ape::extract.clade(phy, n)$tip.label))
		species <- lapply(species, sort)
		return(species[order(sapply(species, "[", 1))])
	}
	candidates <- sort(trees[[1]]$tip.label)
	n <- 2 * length(candidates) - 1
	brlengths <- lapply(trees, function(phy) phy$edge.length[match(seq(n), phy$edge[,2])])
	nodepaths <- lapply(trees, function(phy) ape::nodepath(phy))
	speciesanc <- mapply(endpoints, nodepaths, brlengths, MoreArgs=list(collapse=collapse), SIMPLIFY=FALSE)
	species <- Map(offspring, trees, speciesanc)
	delim <- matrix(, length(trees), length(candidates), dimnames=list(NULL, candidates))
	for (i in seq_along(species)) {
		for (j in seq_along(speciesanc[[i]])) {			
			delim[i, species[[i]][[j]]] <- j
		}
	}
	return(delim)
}



speciesPP <- function(psample, sep="-") {
	delims <- apply(psample, 1, paste, collapse=sep)
	delims <- sort(table(delims), decreasing=TRUE) / length(delims)
	delims <- data.frame(Delimitation=names(delims), PP=as.numeric(delims), stringsAsFactors=FALSE)
	species <- lapply(seq(nrow(psample)), function(i) split(colnames(psample), psample[i,]))
	species <- unlist(species, recursive=FALSE)
	species <- sapply(species, paste, collapse=sep)
	species <- sort(table(species) / nrow(psample), decreasing=TRUE)
	species <- setNames(as.numeric(species), names(species))
	species <- data.frame(Species=names(species), PP=as.numeric(species), stringsAsFactors=FALSE)
	return(list(species=species, delims=delims, tips=colnames(psample)))
}



speciesM50 <- function(pp, sep="-") {
	no <- length(pp$tips)
	i <- 0
	incl <- NULL
	while (length(incl) < no) {
		add <- unlist(strsplit(pp$species[i+1,1], sep))
		if (length(intersect(add, incl)) == 0) {
			incl <- c(incl, add)
		}
		i <- i + 1
	}
	m50 <- pp$species[1:i,]
	delim <- strsplit(m50[,1], split=sep)
	delim <- setNames(rep(seq_along(delim), sapply(delim, length)), unlist(delim))[pp$tips]
	delim <- paste(match(delim, unique(delim)), collapse=sep)
	attr(m50, "pp") <- pp$delims[match(delim, pp$delims[,1]),2]
	return(m50)
}


speciesMAP <- function(pp, sep="-") {
	map <- split(pp$tips, unlist(strsplit(pp$delims[1,1], split=sep)))
	map <- lapply(map, function(y) intersect(pp$tips, y))
	map <- sapply(map, paste, collapse=sep)
	map <- pp$species[match(map, pp$species[,1]),]
	attr(map, "pp") <- pp$delims[1,2]
	return(map)
}



speciesMC <- function(pp, sep="-") {
	delims <- lapply(strsplit(pp$delims[,1], split=sep), as.numeric)
	delims <- lapply(delims, as.numeric)
	delims <- lapply(delims, function(x) sapply(split(pp$tips, x), paste, collapse=sep))
	score <- sapply(delims, function(x) prod(pp$species[match(x, pp$species[,1]),2]))
	maxsc <- which.max(score)
	mc <- split(pp$tips, unlist(strsplit(pp$delims[maxsc,1], split=sep)))
	mc <- lapply(mc, function(y) intersect(pp$tips, y))
	mc <- sapply(mc, paste, collapse=sep)
	mc <- pp$species[match(mc, pp$species[,1]),]
	attr(mc, "pp") <- pp$delims[maxsc,2]
	return(mc)
} 



#' Pairwise posterior probabilities.
#' 
#' @description
#' Calculates posterior probabilities of being conspecific for any pair of tips.
#' These are summed both over the samples and different species where the tips occurred together.
#' 
#' @param psample The posterior sample of species delimitations in the form returned by `speciesDA`.
#' @details This calculation is intended to capture recurrent features of sampled delimitations.
#'   Note, that even mutually incompatible delimitations can contribute to the same pairwise posterior probability.
#' @returns A symmetric matrix with pairwise posterior probabilities.
#' @export

pairwise_pp <- function(psample) {
	pairwise <- matrix(0, ncol(psample), ncol(psample), dimnames=list(colnames(psample), colnames(psample)))
	for (i in seq(nrow(psample))) {
		pool <- outer(psample[i,], psample[i,], "==")
		diag(pool) <- rowSums(pool) == 1
		pairwise <- pairwise + pool
	}
	pairwise <- pairwise / nrow(psample)
	return(pairwise)
}
