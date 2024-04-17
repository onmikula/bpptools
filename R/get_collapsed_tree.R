#' Clades collapsed to branches.
#' 
#' @description
#' Collapses specified clades to single branches.
#' 
#' @param tree An object of class `phylo` or a character string with tree in Newick format.
#' @param clades A data frame listing tip labels (1st column) and clade labels (2nd column) or a name of 'imap' file.
#' @details The species are required to be monophyletic in the tree, but they need not to cover the whole tree.
#' @returns An object of class `phylo` specifying a tree with a single tip per clade and the terminal branch length
#'   equal to the mean distance between MRCA of the population and its tips in the original tree.
#' @export

get_collapsed_tree <- function(tree, clades) {
	find_common_branch <- function(tip, phy) ifelse(length(tip) == 1, match(tip, phy$tip.label), ape::getMRCA(phy, tip))
	if (!inherits(tree, "phylo")) {
		if (substr(tree, 1, 1) == "(") {
			tree <- ape::read.tree(text=tree)
		} else {
			tree <- ape::read.tree(file=tree)
		}
	}
	if (is.character(clades)) {
		clades <- read_bpp_imap(clades)
	}
	clades <- split(clades[,1], clades[,2])
	dst <- ape::dist.nodes(tree)
	anc <- lapply(clades, find_common_branch, phy=tree)
	tips <- lapply(clades, match, table=tree$tip.label)
	collapsed <- ape::drop.tip(tree, tip=unlist(lapply(clades, "[", -1)))
	ord <- match(sapply(clades, "[", 1), collapsed$tip.label)
	collapsed$tip.label[ord] <- names(clades)
	return(collapsed)
}

