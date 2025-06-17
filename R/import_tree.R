#' Importing trees.
#' 
#' @description
#' Imports a phylogenetic tree(s) from a newick text string or a file written in newick or nexus format 
#'   or from BPP-specific files (control, output or mcmc).
#' 
#' @param tree A character string with a tree written in newick format or a file name,
#'   or an object of class `phylo`.
#' @param thin A proportion of mcmc steps to be retained or a subsampling frequency (if `thin > 1`).
#' @returns An object of class `phylo`.
#' @export

import_tree <- function(tree, thin=1) {
	if (inherits(tree, "phylo") | inherits(tree, "multiPhylo")) {
		phy <- tree
	} else if (inherits(tree[[1]], "phylo")) {
		phy <- tree
	} else {
		newick <- all(sapply(tree, function(x) grepl("^\\(.+;$", x)))
		if (isTRUE(newick)) {
			if (length(tree) > 1) {
				phy <- sapply(tree, function(text) ape::read.tree(text=text))
			} else {
				phy <- ape::read.tree(text=tree)
			}
		} else {
			first <- gsub("^[[:blank:]]+|[[:blank:]]+$", "", base::readLines(tree, n=100))
			first <- first[nchar(first) > 0]
			newick <- grepl("^\\(.+;$", first[1])
			nexus <- any(grepl("nexus", first, ignore.case=TRUE))
			ctl <-  any(grepl("speciesdelimitation[[:blank:]]*=", first, ignore.case=TRUE))
			out <- any(grepl("BPP", first, ignore.case=TRUE))
			if (isTRUE(newick)) {
				mcmc <- grepl("#[[:digit:]]+\\.[[:digit:]]+;|;[[:blank:]]+[[:digit:]]+", first[1])
				if (isTRUE(mcmc)) {
					phy <- read_bpp_mcmc(file=tree)$trees
				} else {
					phy <- ape::read.tree(file=tree)
				}
			} else if (isTRUE(nexus)) {
				phy <- ape::read.nexus(file=tree)
			} else if (isTRUE(ctl) | isTRUE(out)) {
				phy <- read_bpp_tree(tree, as.phylo=TRUE)
				tips <- attr(phy, "tips")
				phy <- ape::read.tree(text=tree)
				attr(phy, "tips") <- tips
			}
		}
	}
	if (all(sapply(phy, inherits, what="phylo"))) {
		class(phy) <- "multiPhylo"
	}
	if (inherits(phy, "multiPhylo") & thin != 1) {
		thin <- floor(ifelse(thin < 1, 1 / thin, thin))
		thin <- as.integer(seq(0, length(phy), by=thin))[-1]
		phy <- phy[thin]
	}
	return(phy)
}

