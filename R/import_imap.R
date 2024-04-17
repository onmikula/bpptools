#' Importing classification information.
#' 
#' @description
#' The function importing information about classification of individuals into (candidate) species
#' and (possibly) sequences to individuals.
#' 
#' @param file A name of file in table format. It must contain either the columns labelled `"Individuals"`
#'   and `"Species"` or these two and an additional column `"Sequence"` or just two arbitrarily labelled columns,
#'   assumed to be `"Individuals"` and `"Species"` in that order.
#'   The column labels can be contained in the file or supplied by the `col.names` argument. 
#' @param header Logical, whether the first line contains a header.
#' @param col.names An optional vector of column names, useful especially if `header == FALSE`.
#' @param sep A field-separating character, default is `"\t"` (=tabulator).
#' @returns A data frame with components `"Sequence"` (optional), `"Individual"` and "Species".
#' @export

import_imap <- function(file, header=TRUE, col.names, sep="\t") {
	if (missing(col.names)) {
		imap <- read.table(file, header=header, sep=sep)
	} else {
		imap <- read.table(file, header=header, col.names=col.names, sep=sep)
	}
	ncol <- length(imap)
	cols <- match(c("Individual", "Species", "Sequence"), colnames(imap))
	if (!any(is.na(cols))) {
		imap <- imap[,cols]
	} else if (!any(is.na(cols[1:2]))) {
		imap <- imap[,cols[1:2]]
	} else if (length(imap) == 2) {
		names(imap) <- c("Individual", "Species")
	} else {
		stop("The table has more than two columns, but they are not properly labelled")
	}
	return(imap)	
}
