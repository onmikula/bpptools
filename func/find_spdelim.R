# Function:
#	speciesDA: extracts species from species trees 
# Arguments:
#	trees: 'multiPhylo' object or a list of 'phylo' objects assumed to represent posterior sample of rooted trees from Bayesian analysis
#	collapse: divergence time (in appropriate units) approximating zero
# Value:
#	matrix whose rows correspond to delimitations sampled by MCC, columns correspond to individuals and integer entries indicate assignement of individuals to species
# Depends on: ape

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



# Function:
#	speciesPP: enumerates delimited species and calculates their posterior probabilities
# Arguments:
#	psample: posterior sample of species delimitations in the form of matrix produced by 'speciesDA'
#	sp: if not NULL, it specifies a subset of candidate species whose merger into a single species is assessed
#	sep: separator symbol used to concatenate candidate species names
# Value:
#	data frame with names of delimited species and their posterior probabilities (is.null(sp)) or the posterior probability of merger specified by sp

speciesPP <- function(psample, sp=NULL, sep="-") {
	if (!is.null(sp)) {
		no <- as.list(apply(psample[, sp, drop=FALSE], 1, unique))
		pp <- 0
		for (i in which(sapply(no, length) == 1)) {
			pp <- pp + (sum(psample[i,] == no[[i]]) == length(sp))
		}
		pp <- pp / nrow(psample)
	} else {
		species <- lapply(seq(nrow(psample)), function(i) split(colnames(psample), psample[i,]))
		species <- unlist(species, recursive=FALSE)
		species <- sapply(species, paste, collapse=sep)
		pp <- sort(table(species) / nrow(psample), decreasing=TRUE)
		pp <- setNames(as.numeric(pp), names(pp))
		pp <- data.frame(Species=names(pp), PP=as.numeric(pp), stringsAsFactors=FALSE)
	}
	return(pp)
}



# Function:
#	speciesMC: extracts maximum credibility species delimitation
# Arguments:
#	pp: data frame with names of delimited species and their posterior probabilities
#	sep: separator symbol used to concatenate candidate species names
# Value:
#	data frame with names of species in the maximum credibility delimitation and their posterior probabilities

speciesMC <- function(pp, sep="-") {
	if (any(!c("Species", "PP") %in% names(pp))) {
		pp <- setNames(pp[,1:2], c("Species", "PP"))
	}
	species <- strsplit(pp$Species, split=sep)
	basic <- sort(unique(unlist(species)))
	mat <- matrix(0, length(species), length(species), dimnames=list(pp$Species, pp$Species))
	for (i in 1:(length(species) - 1)) {
		for (j in (i+1):length(species)) {
			mat[i,j] <- mat[i,j] <- as.numeric(length(intersect(species[[i]], species[[j]])) == 0)
		}
	}
	mcdelim <- 1
	complete <- length(setdiff(basic, unlist(species[mcdelim]))) == 0
	while (complete == FALSE) {		
		mcdelim <- c(mcdelim, min(which(colSums(mat[mcdelim,,drop=FALSE] == 1) == length(mcdelim))))
		complete <- length(setdiff(basic, unlist(species[mcdelim]))) == 0
	}
	return(pp[mcdelim,])	
} 



# Function:
#	find_spdelim: a wrapper function which extracts maximum credibility species delimitation starting from the posterior sample of species trees. It joins 'speciesDA', 'speciesPP' and 'speciesMC' into a single pipeline.
# Arguments:
#	trees: 'multiPhylo' object or a list of 'phylo' objects assumed to represent posterior sample of rooted trees from Bayesian analysis
#	collapse: divergence time (in appropriate units) approximating zero
#	sep: separator symbol used to concatenate candidate species names
#	file: name of the output file 
# Value:
#	data frame with names of species in the maximum credibility delimitation and their posterior probabilities; tab-delimited file with this output (if 'file' is specified)

find_spdelim <- function(trees, collapse=0, sep="-", file) {
	delim <- speciesDA(trees, collapse= collapse)
	delim_pp <- speciesPP(delim, sep=sep)
	delim_sp <- speciesMC(delim_pp, sep=sep)
	if (!missing(file)) {
		write.table(delim_sp, file, row.names=FALSE, quote=FALSE, sep="\t")
	}
	return(delim_sp)
}



# Function:
#	pairwise_delimitation: quantifies pairwise posterior probabilities the candidate species are conspecific
# Arguments:
#	psample: posterior sample of species delimitations in the form of matrix produced by 'speciesDA'
# Value:
#	symmetric matrix with pairwise posterior probabilities
# Details:
#	This calculation can be though as a kind of consensus solution, as it captures recurrent features of sample delimitations. Note, that even mutually incompatible delimitations can contribute to the same pairwise posterior probability.

pairwise_delimitation <- function(psample) {
	pairwise <- matrix(0, ncol(psample), ncol(psample), dimnames=list(colnames(psample), colnames(psample)))
	for (i in seq(nrow(psample))) {
		pool <- outer(psample[i,], psample[i,], "==")
		diag(pool) <- rowSums(pool) == 1
		pairwise <- pairwise + pool
	}
	pairwise <- pairwise / nrow(psample)
	return(pairwise)
}


