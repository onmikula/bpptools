#' Starting species tree.
#' 
#' @description
#' Reads starting species tree from BPP control file (.ctl) or output file (.out).
#' 
#' @param file The name of the control or output file.
#' @param as.phylo Logical, whether to return tree as a `phylo` object.
#' @returns The starting (and possibly fixed) species tree topology in newick or `phylo` format
#'   with attribute `"tips"`, which contains tip labels in the same order as in the control file.
#' @export

read_bpp_tree <- function(file, as.phylo=FALSE) {
	lin <- gsub("^[[:blank:]]*|[[:blank:]]*$", "", readLines(file))
	lin <- sub("[[:blank:]]*\\*.*$|[[:blank:]]*#.*$", "", lin)
	lin <- lin[nchar(lin) > 0]
	type <- ifelse(any(grepl("Command:", head(lin, 50), ignore.case=TRUE)), "out", "ctl")
	if (type == "ctl") {
		tree <- regmatches(lin, regexpr("\\(.*;$", lin))		
		tree <- gsub("[[:blank:]]", "", tree)	
		tips <- sub("^\\s*species&tree\\s*=\\s*[[:digit:]]+\\s*", "", lin[grepl("species&tree", lin)])
		tips <- unlist(strsplit(gsub("\\s+", " ", tips), " "))
		attr(tree, "tips") <- tips	
		sptree <- grepl("=1", gsub(" ", "", grep("speciestree", lin, value=TRUE)))
		if (sptree) {
			warning("species tree inference, starting tree is reported")
		}
	}
	if (type == "out") {
		fin <- grep("spent in MCMC", lin)
		tree <- grep("Initial species tree", lin[seq(fin)], value=TRUE, ignore.case=TRUE)
		tree <- gsub("[[:blank:]]", "", regmatches(tree, regexpr("\\(.+\\);", tree)))
		ntip <- length(unlist(gregexpr("\\(", tree))) + 1
		mode <- list(A00=grep("List of nodes", lin, ignore.case=TRUE),
		A01=grep("Species in order", lin, ignore.case=TRUE),
		A10=grep("Order of ancestral nodes", lin, ignore.case=TRUE))
		if (length(mode[["A00"]]) > 0) {
			tips <- lin[mode[["A00"]] + 2:(ntip+1)]
			tips <- gsub("[[:blank:]]+", " ", tips)
			tips <- sapply(strsplit(tips, " "), "[", 4)
		} else if (length(mode[["A01"]]) > 0) {
			tips <- lin[mode[["A01"]] + 1:ntip]
			tips <- sub("^[[:digit:]]+\\.[[:blank:]]*", "", tips)
			warning("species tree inference, starting tree is reported")
		} else if (length(mode[["A10"]]) > 0) {
#			tips <- lin[mode[["A10"]] + seq(ntip-1)]
			tips <- NULL
		} else {
			tips <- NULL
			warning("species tree inference, starting tree is reported")			
		}
		attr(tree, "tips") <- tips	
	}
	if (isTRUE(as.phylo) ) {
		tips <- attr(tree, "tips")
		tree <- ape::read.tree(text=tree)
		attr(tree, "tips") <- tips
	}
	return(tree)
}
