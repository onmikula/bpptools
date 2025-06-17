#' Additional statistics.
#' 
#' @description
#' Estimates additional statistics from sampled parameter values.
#' 
#' @param bpp An object of class `bpp` as created by [read_bpp_mcmc] or [as.bpp],
#'   or a data frame with postrior sample of parameters.
#' @param which Character, which statistic to estimate.
#'   Currently limited to the names of functions defined in `bpptools`.
#' @param args A named list of arguments to the estimator function.
#' @export

estimate_statistic <- function(bpp, which, args=NULL) {
	isbpp <- inherits(bpp, "bpp")
	if (isbpp) {
		trees <- bpp$trees
		params <- bpp$params
	} else {
		params <- bpp
	}
	which <- tolower(which)
	if (which == "gdi") {
		addit <- gdi(params, tree=args[["tree"]], interval=args[["interval"]], p=args[["p"]])$distributions
		names(addit) <- paste("gdi", names(addit), sep="_")
	}
	if (which == "tl") {
		addit <- setNames(tl(trees), which)
	}
	
	params <- data.frame(params, addit)
	if (isbpp) {
		bpp$trees <- trees
		bpp$params <- params
	}
	return(bpp)
}



tl <- function(phy) {
	if (inherits(phy, "phylo")) {
		phy <- list(phy)
	}
	treelengths <- data.frame(TL=sapply(lapply(phy, "[", "edge.length"), sum))
	return(treelengths)
}
