#' Writing BPP input files.
#' 
#' @description
#' Writing BPP input files: the sequence (seq) file in a modified phylip format
#' and the mapping (imap) file with classification of individuals into (candidate) species.
#' 
#' @param seq A list of locus-specific alignments of class `matrix` or `DNAbin`.
#' @param imap A matrix or data frame (1st column individuals, 2nd column species) or a name of imap file.
#' @param seqfile The sequence file name.
#' @param imapfile The imap file name.
#' @param seqmap A data frame mapping names of sequences (1st column) to those of individuals (2nd column).
#' @param allele An allele identifier, regular expression distinguishing sequences from the same individual.
#' @details
#' While `seqmap` contains an explicit list of sequence names, `allele` is a regular expression defining
#' allele-specific extensions of individual names. For instance, if an individual name is "USB22" and
#' `allele` is "_[01]", the sequences with names "USB22_0" and "USB22_1" are recognized as belonging
#' to the same individual.
#'
#' Individual names must not be each other's partial matches. For instance, if an individual "USB22"
#' is present, avoid the name "USB2" in any other individual (use "USB02" instead).
#' @export

write_bpp_input <- function(seq, imap, seqfile, imapfile, seqmap=NULL, allele=NULL) {

	seq <- lapply(lapply(seq, as.matrix), toupper)
	for (i in seq_along(seq)) {
		seq[[i]] <- seq[[i]][seqlen(seq[[i]]) > 0,,drop=FALSE]
	}
	if (is.character(imap)) {
		imap <- read_bpp_imap(imapfile, col.names=c("Individual", "Species"))
	} else {
		imap <- setNames(as.data.frame(imap[,1:2]), c("Individual", "Species"))
	}
	pmatches <- setNames(lapply(imap$Individual, grep, x=imap$Individual), imap$Individual)
	pmatches <- pmatches[sapply(pmatches, length) > 1]
	nmatches <- length(pmatches)
	if (nmatches > 0) {
		if (nmatches == 1) {
			message <- "Consider change of the name "
			namelist <- names(pmatches)
		} else {
			message <- "Consider change of the names "			
			namelist <- paste(paste(names(pmatches)[-nmatches], collapse=", "), "and", names(pmatches)[nmatches])
		}
		message <- paste0(message, namelist, ".")
		stop(paste("Individual names must not be each other's partial matches.\n", message))
	}

	sm <- !is.null(seqmap)
	al <- !is.null(allele)
	if (!sm) {
		if (al) {
			allseqnam <- sort(unique(unlist(lapply(seq, rownames))))
			seqnames <- paste0(imap[,1], allele)
			seqnames <- lapply(seqnames, grep, x=allseqnam, value=TRUE)
			seqmap <- data.frame(Sequence=unlist(seqnames), Individual=rep(imap[,1], sapply(seqnames, length)))
		} else {
			seqmap <- data.frame(Sequence=imap[,1], Individual=imap[,1])
		}
	} else if (al) {
		allseqnam <- sort(unique(unlist(lapply(seq, rownames))))
		seqmap <- setNames(as.data.frame(seqmap[,1:2]), c("Sequence", "Individual"))
		seqmap <- unique(rbind(seqmap, data.frame(Sequence=paste0(imap[,1], allele), Individual=imap[,1])))
		seqnames <- lapply(seqmap$Sequence, grep, x=allseqnam, value=TRUE)
		seqmap <- unique(data.frame(Sequence=unlist(seqnames), Individual=rep(seqmap$Individual, sapply(seqnames, length))))
	} else {
		seqmap <- setNames(as.data.frame(seqmap[,1:2]), c("Sequence", "Individual"))
	}

	indnames <- vector("list", length(seq))
	for (i in seq_along(seq)) {
		indnames[[i]] <- seqmap$Individual[match(rownames(seq[[i]]), seqmap$Sequence)]
		if (any(is.na(indnames[[i]]))) {
			message <- paste0("Check sequence names in locus no. ", i, ".")
			stop(paste("All sequence names must be classified to individuals.\n", message))
		}
		seqnames <- paste(rownames(seq[[i]]), indnames[[i]], sep="^")
		seqstrings <- apply(seq[[i]], 1, paste, collapse="")
		seq[[i]] <- c(paste(dim(seq[[i]]), collapse=" "), paste(seqnames, seqstrings, sep=" "), "\n")
	}
	writeLines(unname(unlist(seq)), con=seqfile)

	indnames <- sort(unique(unlist(indnames)))
	imap <- imap[match(indnames, imap$Individual),c("Individual","Species")]
	imap <- unname(apply(imap, 1, paste, collapse=" "))
	if (missing(imapfile)) {
		imapfile <- sub("\\.[[:alnum:]]+$", ".imap", seqfile)
	}
	writeLines(imap, imapfile)
	
}



#' Writing BPP seq file.
#' 
#' @description
#' It writes BPP seq file.
#' 
#' @param seq A list of locus-specific alignments of class `matrix` or `DNAbin`.
#' @param imap A matrix or data frame (1st column individuals, 2nd column species) or a name of imap file.
#' @param seqfile The sequence file name.
#' @param seqmap A data frame mapping names of sequences (1st column) to those of individuals (2nd column).
#' @param allele An allele identifier, regular expression distinguishing sequences from the same individual.
#' @details
#' While `seqmap` contains an explicit list of sequence names, `allele` is a regular expression defining
#' allele-specific extensions of individual names. For instance, if an individual name is "USB22" and
#' `allele` is "_[01]", the sequences with names "USB22_0" and "USB22_1" are recognized as belonging
#' to the same individual.
#'
#' Individual names must not be each other's partial matches. For instance, if an individual "USB22"
#' is present, avoid the name "USB2" in any other individual (use "USB02" instead).
#' @seealso
#'   [write_bpp_input()] for a function writing both seq & imap files.
#' @export

write_bpp_seq <- function(seq, imap, seqfile, seqmap=NULL, allele=NULL) {

	seq <- lapply(lapply(seq, as.matrix), toupper)
	for (i in seq_along(seq)) {
		seq[[i]] <- seq[[i]][seqlen(seq[[i]]) > 0,,drop=FALSE]
	}
	if (is.character(imap)) {
		imap <- read_bpp_imap(imapfile, col.names=c("Individual", "Species"))
	} else {
		imap <- setNames(as.data.frame(imap[,1:2]), c("Individual", "Species"))
	}
	pmatches <- setNames(lapply(imap$Individual, grep, x=imap$Individual), imap$Individual)
	pmatches <- pmatches[sapply(pmatches, length) > 1]
	nmatches <- length(pmatches)
	if (nmatches > 0) {
		if (nmatches == 1) {
			message <- "Consider change of the name "
			namelist <- names(pmatches)
		} else {
			message <- "Consider change of the names "			
			namelist <- paste(paste(names(pmatches)[-nmatches], collapse=", "), "and", names(pmatches)[nmatches])
		}
		message <- paste0(message, namelist, ".")
		stop(paste("Individual names must not be each other's partial matches.\n", message))
	}

	sm <- !is.null(seqmap)
	al <- !is.null(allele)
	if (!sm) {
		if (al) {
			allseqnam <- sort(unique(unlist(lapply(seq, rownames))))
			seqnames <- paste0(imap[,1], allele)
			seqnames <- lapply(seqnames, grep, x=allseqnam, value=TRUE)
			seqmap <- data.frame(Sequence=unlist(seqnames), Individual=rep(imap[,1], sapply(seqnames, length)))
		} else {
			seqmap <- data.frame(Sequence=imap[,1], Individual=imap[,1])
		}
	} else if (al) {
		allseqnam <- sort(unique(unlist(lapply(seq, rownames))))
		seqmap <- setNames(as.data.frame(seqmap[,1:2]), c("Sequence", "Individual"))
		seqmap <- unique(rbind(seqmap, data.frame(Sequence=paste0(imap[,1], allele), Individual=imap[,1])))
		seqnames <- lapply(seqmap$Sequence, grep, x=allseqnam, value=TRUE)
		seqmap <- unique(data.frame(Sequence=unlist(seqnames), Individual=rep(seqmap$Individual, sapply(seqnames, length))))
	}

	indnames <- vector("list", length(seq))
	for (i in seq_along(seq)) {
		indnames[[i]] <- seqmap$Individual[match(rownames(seq[[i]]), seqmap$Sequence)]
		if (any(is.na(indnames[[i]]))) {
			message <- paste0("Check sequence names in locus no. ", i, ".")
			stop(paste("All sequence names must be classified to individuals.\n", message))
		}
		seqnames <- paste(rownames(seq[[i]]), indnames[[i]], sep="^")
		seqstrings <- apply(seq[[i]], 1, paste, collapse="")
		seq[[i]] <- c(paste(dim(seq[[i]]), collapse=" "), paste(seqnames, seqstrings, sep=" "), "\n")
	}
	writeLines(unname(unlist(seq)), con=seqfile)
	
}



#' Writing BPP imap file.
#' 
#' @description
#' It writes BPP imap file, or rewrites it by merging species according to results of a previous analysis.
#' 
#' @param imap A data frame or a matrix defining classification (=mapping) of individuals (1st column)
#'   to species (2nd column) or a name of .imap file. 
#' @param imapfile The imap file name.
#' @param return Logical, whether to return (possibly rewritten) data frame.
#' @param delim A result of species delimitation (mergers of candidate species);
#'   either a vector with concatenated candidate species names or a data frame with candidate species (1st column)
#'   and their mergers (2nd column) or a name of tab-delimited file with such information.
#' @param sep A separator of the original candidate species names, relevant if the vector is used in `delim` argument.
#' @returns The data frame (if `return == TRUE`) and/or .imap file (if `file` is specified).
#' @export

write_bpp_imap <- function(imap=NULL, imapfile, return=TRUE, delim=NULL, sep="-") {
	if (!is.null(imap)) {
		if (is.character(imap) & length(imap) == 1) {
			imap <- read_bpp_imap(imap)
		}
	}
	if (!is.null(delim)) {
		if (is.character(delim) & length(unlist(delim)) == 1) {
			delim <- read.delim(delim)
		}
		if (is.matrix(delim) | is.data.frame(delim)) {
			delim <- split(delim[,1], delim[,2])
		}
		if (!is.list(delim)) {
			delim <- setNames(strsplit(delim, split=sep), delim)
		}
		if (!is.null(imap)) {
			imap <- split(imap[,1], imap[,2])
			delim <- lapply(delim, function(d) unlist(imap[d]))
		}
		imap <- data.frame(Individual=unlist(delim), Species=rep(names(delim), sapply(delim, length)), row.names=NULL)
		imap$Species <- gsub(sep, "", imap$Species)
	} else if (is.null(imap)) {
		stop("at least one of 'imap' and 'delim' arguments has to be specified")
	}
	
	if (!missing(imapfile)) {
		writeLines(unname(apply(imap, 1, paste, collapse=" ")), imapfile)
	}
	if (isTRUE(return)) {
		return(imap)
	}
}



seqlen <- function(x) {
	rowSums(matrix(!x %in% c("N","-","?"), nrow(x), ncol(x)))
}
