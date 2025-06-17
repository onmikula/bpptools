#' Clades collapsed to branches.
#' 
#' @description
#' Collapses specified clades to single branches.
#' 
#' @param tree An object of class `phylo` or `multiPhylo` or a character string with tree in newick format
#'   or a name of file containing the tree(s).
#' @param clades A data frame listing tip labels (1st column) and clade labels (2nd column) or a name of mapping file.
#' @details The specified clades must be monophyletic in the tree(s), but may cover only a part of the tree(s).
#' @returns An object of class `phylo` of `multiPhylo` specifying a tree with a single tip per clade
#'   and terminal branch lengths equal to the mean root-to-tip distances in the clades.
#' @export

collapse_clades <- function(tree, clades) {
	find_common_anc <- function(tip, phy) {
		ifelse(length(tip) == 1, phy$edge[match(tip, phy$edge[,2]),1], ape::getMRCA(phy, tip))
	}
	find_clade_edges <- function(phy, from, to) {
		edg <- lapply(to, function(x) ape::nodepath(phy, from=from, to=x))
		edg <- lapply(edg, function(x) cbind(x[-length(x)], x[-1]))
		edg <- rbind(phy$edge, unique(do.call(rbind, edg)))
		which(duplicated(edg, fromLast=TRUE))
	}
	tree <- import_tree(tree)
	if (is.character(clades) & length(clades) == 1) {
		clades <- read_bpp_imap(clades)
	}
	clades <- split(clades[,1], clades[,2])
	tips <- sapply(clades, "[", 1)
	drop <- setdiff(unlist(clades), tips)
	if (inherits(tree, "phylo")) {
		tree <- list(tree)
	}
	for (i in seq_along(tree)) {
		phy <- tree[[i]]	
		off <- lapply(clades, match, table=phy$tip.label)
		anc <- sapply(off, find_common_anc, phy=phy)
		edg <- lapply(seq_along(off), function(j) find_clade_edges(phy=phy, from=anc[j], to=off[[j]]))
		if (!is.null(phy$edge.length)) {
			brlen <- phy$edge.length[match(tips, phy$tip.label)]	
			brlen <- sapply(edg, function(jj) mean(phy$edge.length[jj])) - brlen
		}
		tips <- sapply(clades, "[", 1)
		tree[[i]] <- ape::drop.tip(phy, tip=drop)
		ord <- match(tips, tree[[i]]$tip.label)
		tree[[i]]$tip.label[ord] <- names(clades)
		terminal <- match(ord, tree[[i]]$edge[,2])
		if (!is.null(phy$edge.length)) {
			tree[[i]]$edge.length[terminal] <- tree[[i]]$edge.length[terminal] + brlen
		}
	}
	if (length(tree) == 1) {
		tree <- tree[[1]]
	}
	return(tree)
}

