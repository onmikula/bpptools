#' Starting species tree.
#' 
#' @description
#' Reads starting species tree from BPP control file (.ctl) or output file (.out).
#' 
#' @param file The name of the control or output file.
#' @returns The starting species tree topology in newick format with attribute `"tips"`,
#'   which contains tip labels in the same order as in the control file.
#' @export

read_bpp_tree <- function(file) {
	lin <- readLines(file)
	type <- ifelse(any(grepl("Command:", head(lin, 50))), "out", "ctl")
	if (type == "ctl") {
		tree <- regmatches(lin, regexpr("\\(.*;$", lin))		
		tree <- gsub("[[:blank:]]", "", tree)	
		tips <- sub("^\\s*species&tree\\s*=\\s*[[:digit:]]+\\s*", "", lin[grepl("species&tree", lin)])
		tips <- unlist(strsplit(gsub("\\s+", " ", tips), " "))
		attr(tree, "tips") <- tips	
	}
	return(tree)
}

