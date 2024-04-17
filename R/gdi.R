#' Genealogical divergence index.
#' 
#' @description
#' Calculates the genealogical divergence index.
#' 
#' @param params A data frame with multispecies coalescent parameters.
#' @param tree A tree in Newick format obtained by `read_bpp_tree`
#'   or the name of BPP control (.ctl) or output (.out) file.
#' @param interval Either `"HPD"` (default) or `"CPD"`.
#' @param p The proportion of posterior density be included in the interval, default is 0.95.
#' @details Currently, it calculates gdi for every branch in the tree (internal or terminal)
#'   under assumption of no postdivergence gene flow. The calculation for terminal branches 
#'   under models with gene flow is not implemented yet.
#' @returns A list with components `"means"` (vector of posterior means), `"distributions"` (matrix
#'   with posterior distributions) and `"interval"` (matrix with posterior density interval limits).
#' @export

gdi <- function(params, tree, interval="HPD", p=0.95) {
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
	thetanum <- gsub("theta_", "", colnames(params)[seq(ntip)])
	thetanum <- intersect(as.numeric(thetanum), seq(ntip))
	thetamis <- setdiff(seq(ntip), thetanum)
	if (length(thetamis) > 0) {
		thetatip <- data.frame(params[,seq_along(thetanum)], matrix(NA, nrow(params), length(thetamis)))
		names(thetatip) <- c(names(params)[seq_along(thetanum)], paste0("theta_", thetamis))
		params <- cbind(thetatip[,order(c(thetanum, thetamis))], params[,-seq_along(thetanum)])
	}
	tau <- cbind(as.data.frame(setNames(rep(list(numeric(nsample)), ntip),paste0("tau",seq(ntip)))), params[,grepl("tau",colnames(params))])
	theta <- params[,grepl("theta",colnames(params))]
	theta[,seq(ntip)] <- theta[,match(tree$tip.label, tips)]
	anc <- tree$edge[match(seq(ncol(theta)), tree$edge[,2]),1]
	pgdi <- matrix(NA, nsample, ncol(theta), dimnames=list(rownames(params), paste0("sp_",seq(ncol(theta)))))
	for (i in setdiff(seq(ncol(theta)), ntip+1)) {
		pgdi[,i] <- 1 - exp(-2 * (tau[,anc[i]] - tau[,i]) / theta[,i])
	}
	colnames(pgdi)[seq(ntip)] <- tree$tip.label
	clades <- c(as.list(seq(ntip)), ape::prop.part(tree))
	for (i in as.numeric(tree$node.label[-1])) {
		off <- tree$edge[tree$edge[,1] == i,2]
		off <- sort(tree$tip.label[sapply(clades[off], "[", 1)])
		colnames(pgdi)[i] <- paste0("anc:", paste(off, collapse="-"))
	}
	pgdi <- pgdi[,-(ntip+1)]	
	mgdi <- round(apply(pgdi, 2, mean, na.rm=TRUE), 6)
	mgdi[is.nan(mgdi)] <- NA
	cint <- round(t(apply(pgdi, 2, match.fun(interval), p=p)), 6)
	
	return(list(means=mgdi, distributions=pgdi, interval=cint))
}
