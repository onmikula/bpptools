#' Reading BPP sequence file.
#' 
#' @description
#' Reads sequences from BPP sequence file.
#' 
#' @param seqfile The name of sequence file.
#' @param nloci The number of loci to be imported (in the same order as in seqfile). The default means "all loci".
#' @param locnames A character vector with names of loci, if missing, the loci are given names "Lddd" ("d" is for digits).
#' @returns The list of matrices with locus-specific sequence alignments.
#' @export

read_bpp_seq <- function(seqfile, nloci=-1L, locnames) {
	lin <- readLines(seqfile)
	lin <- lin[lin != ""]
	header <- which(!grepl("\\^", lin))
	start <- header + 1
	end <- c(header[-1] - 1, length(lin))
	if (sign(nloci) < 0) {
		nloci <- length(start)
	}
	nloci <- min(c(nloci, length(start)))
	if (missing(locnames)) {
		locnames <- paste("L", formatC(seq(nloci), format="d", flag=0, width=nchar(nloci)), sep="")
	}
	loci <- setNames(vector("list", nloci), locnames)
	for (i in seq(nloci)) {
		loci[[i]] <- lin[start[i]:end[i]]
		loci[[i]] <- unlist(strsplit(loci[[i]], "\\s+"))
		ii <- seq_along(loci[[i]]) %% 2
		loci[[i]] <- setNames(loci[[i]][ii == 0], loci[[i]][ii == 1])
		names(loci[[i]]) <- sub("\\^.*$", "", names(loci[[i]]))		
		loci[[i]] <- do.call(rbind, strsplit(loci[[i]], ""))
	}
	return(loci)
}



#' Reading BPP mapping file.
#' 
#' @description
#' A wrapper for `read.table`; reads in BPP mapping file with classification of individuals to species.
#' 
#' @param imapfile The name of the mapping file.
#' @param col.names The column names, "Individual" & "Species" by default.
#' @returns The data frame with individual IDs (1st column) and species names (2nd column).
#' @export

read_bpp_imap <- function(imapfile, col.names=c("Individual", "Species")) {
	read.table(imapfile, header=FALSE, col.names=col.names[1:2])
}


