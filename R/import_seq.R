#' Importing sequence data.
#' 
#' @description
#' The functions importing sequence data from various file formats.
#'
#' `import_seq` is a wrapper for the following functions.
#'
#' `import_nexus` from .nexus file, loci defined by an appropriate nexus block or an optional argument.
#'
#' `import_fasta` from .fasta file, loci contained in different files or defined by an optional argument.
#'
#' `import_phylip` from .phy file, loci contained in different files or defined by an optional argument.
#'
#' `import_ipyrad` from .alleles file of ipyrad RADseq assembler.
#'
#' `import_stacks` from .alleles.fas file of stacks RADseq assembler.
#' 
#' @param file character, the name(s) of file(s) to be imported.
#' @param format character, the format of the file(s); it can be `"ipyrad"`, `"stacks"`, `"nexus"`, `"fasta"` or `"phylip"`.
#' @param nloci the number of loci to be imported (the default `-1L` means all of them)
#' @param class character, the output class, either `"matrix"` or `"DNAbin"`.
#' @param toupper logical, whether to represent nucleotides by upper case letters.
#' @param part either a two column matrix indicating starts and ends of loci in the concatenated aligments
#'  (in .fasta or .nexus files) or a name of partition file formatted according to RAxML or MrBayes convention.
#'  In [import_nexus] it is overridden by the nexus block defining character sets, if present within the sequence file.
#' @param rm.empty logical, whether to remove empty (e.g. all "-" or all "N") sequences.
#' @param nlines the number of lines to be imported (the default `-1L` means all of them;
#'  relevent for formats `"ipyrad"` and `"stacks"`. An incomplete locus at the end of the imported lines is removed.
#' @param indiv information about individual names, relevant for the format `"stacks"`. It may be `"pop"`
#'  (individuals correspond to populations in 'Stacks'), a vector of individual names in an order corresponding
#'  to sample numbers in Stacks or `NULL` (default), which makes individuals being numbered as in Stacks.
#' @param interleaved Logical, whether in interleaved (not sequential) format; currently only for [import_phylip].
#' @returns A list of locus-specific alignments. If genomic coordinates are included (in ipyrad or Stacks output
#'  or in a "bedtable" block of .nexus file), they are exported in .bed format as an attribute `bedtable`.
#' @export

import_seq <- function(file, format, nloci = -1L, class = "matrix", toupper = TRUE, part = NULL, rm.empty = FALSE, nlines = -1L, indiv = NULL, interleaved = TRUE) {
	fun <- match.fun(paste("import", format, sep="_"))
	fun(file=file, nloci=nloci, part=part, class=class, toupper=toupper, rm.empty=rm.empty, nlines=nlines, indiv=indiv, interleaved=interleaved)
}


#' Importing nexus files.
#' 
#' @description
#' Importing sequence data from files in .nexus format.
#' 
#' @param file Character, the name of .nexus file.
#' @param nloci The number of loci to be imported (the default `-1L` means all of them).
#' @param class Character, the output class, either 'matrix' or 'DNAbin'.
#' @param toupper Logical, whether to represent nucleotides by upper case letters.
#' @param part It can be a two column matrix indicating starts and ends of loci in the concatenated aligments
#'  or a name of partition file formatted according to RAxML or MrBayes convention. It can be overridden 
#'  by the nexus block defining character sets, if present within the sequence file.
#' @param rm.empty Logical, whether to remove empty (e.g. all "-" or all "N") sequences.
#' @param nlines Not relevant, just for compatibility with [import_seq].
#' @param indiv Not relevant, just for compatibility with [import_seq].
#' @param interleaved Not relevant, currently just for compatibility with [import_seq].
#' @returns A list of locus-specific alignments. If genomic coordinates are included in a "bedtable" block of the file
#'  they are exported in .bed format as an attribute `bedtable`.
#' @export

import_nexus <- function(file, nloci=-1L, class="matrix", toupper=TRUE, part=NULL, rm.empty=FALSE, nlines, indiv, interleaved) {

	loci <- readLines(file)
	loci <- loci[!grepl("^[[:space:]]+$", loci)]
	loci <- loci[nchar(loci) > 0]
	start <- grep("MATRIX|Matrix|matrix", loci) + 1
	end <- base::intersect(which(loci == ";"), start:length(loci))[1] - 1
	
	charsets <- grep("charset", loci[(end + 1):length(loci)], ignore.case=TRUE, value=TRUE)
	if(length(charsets) > 0) {
		part <- sub("^[[:blank:]]*[[:alpha:]]+[[:blank:]]+", "", charsets)
		part <- gsub("[[:blank:]]|;", "", part)
		part <- do.call(rbind, strsplit(part, "="))
		part <- cbind(part[,1], do.call(rbind, strsplit(part[,2], "-")))
		part <- data.frame(Start=as.numeric(part[,2]), End=as.numeric(part[,3]), row.names=part[,1])
	} else if (is.character(part)) {
		part <- read_partition_file(part)		
	}
	if (!is.null(part)) {
		nloci <- ifelse(nloci == -1L, nrow(part), nloci)
		if (nloci > nrow(part)) {
			nloci <- nrow(part)
			single <- (nloci > 1) + 1
			warning(paste("There", c("is","are")[single], "only", nloci, c("locus", "loci")[single], "available."))
		} else if (nloci < nrow(part)) {
			part <- part[1:nloci,]
		}		
	} else {
		nloci <- 1
	}

	bedtable <- grep("begin\\s*bedtable", loci, ignore.case=TRUE)
	if (length(bedtable) > 0) {
		gcend <- grep("end;", loci[(bedtable+1):length(loci)], ignore.case=TRUE)[1] - 1
		bedtable <- loci[bedtable + 1:gcend]
		bedtable <- as.data.frame(do.call(rbind, strsplit(bedtable, "[[:space:]]")))
		if (sign(nloci) == 1) {
			bedtable <- bedtable[1:nloci,,drop=FALSE]
		}
		if (length(bedtable) >= 6) {
			names(bedtable)[1:6] <- c("Scaffold", "Start", "End", "Locus", "Score", "Strand")
		} else {
			names(bedtable) <- c("Scaffold", "Start", "End", "Locus", "Score", "Strand")[1:length(bedtable)]
		}
		bedtable$Start <- as.numeric(bedtable$Start)
		bedtable$End <- as.numeric(bedtable$End)
		if (length(bedtable) >= 5) {
			bedtable$Score <- as.numeric(bedtable$Score)
		}
	}
	snpattr <- grep("begin\\s*snps", loci, ignore.case=TRUE)
	if (length(snpattr) > 0) {
		snpend <- grep("end;", loci[(snpattr+1):length(loci)], ignore.case=TRUE)[1] - 1
		header <- unlist(strsplit(loci[snpattr + 1], "[[:space:]]"))
		snpattr <- loci[snpattr + 2:snpend]
		snpattr <- as.data.frame(do.call(rbind, strsplit(snpattr, "[[:space:]]")))
		snpattr[[2]] <- as.numeric(snpattr[[2]])
		snpattr[[3]] <- as.logical(snpattr[[3]])
		snpattr[[4]] <- as.numeric(snpattr[[4]])
		names(snpattr) <- header
	}

	if (grepl("\\s+", loci[start])) {
		loci <- unlist(strsplit(loci[start:end], "\\s+"))
	} else {
		loci <- loci[start:end]
	}
	odd <- seq_along(loci) %% 2
	if(anyDuplicated(loci[odd == 1]) != 0) {
		loci[odd == 0] <- gsub("[[:blank:]]", "", loci[odd == 0])
		seqnam <- unique(loci[odd == 1])
		loci <- sapply(split(loci[odd == 0], loci[odd == 1]), paste, collapse="")[seqnam]
		loci <- c(names(loci), unname(loci))[order(rep(seq_along(seqnam), 2))]
		odd <- seq_along(loci) %% 2
	}
	if (!is.null(part)) {
		loci[odd == 0] <- substr(loci[odd == 0], 1, part$End[nloci])
	}
	if (nloci == 1) {
		part <- data.frame(Start=1, End=nchar(loci[2]), row.names=NULL)
	}
	substrings <- function(i, j, x, n) setNames(substr(x, i, j), n)
	loci <- setNames(base::mapply(substrings, part$Start, part$End, MoreArgs=list(x=loci[odd == 0], n=loci[odd == 1]), SIMPLIFY=FALSE), rownames(part))
	loci <- lapply(loci, function(x) do.call(rbind, strsplit(x, "")))

	if (isTRUE(rm.empty)) {
		loci <- remove_empty(loci, ambig=2)
	}
	if (isTRUE(toupper)) {
		loci <- lapply(loci, toupper)
	}
	if (class == "DNAbin") {
		loci <- lapply(loci, ape::as.DNAbin)
	}

	if (length(snpattr) > 0) {
		snpattr <- snpattr[snpattr$locus %in% names(loci),,drop=FALSE]
		partloc <- setNames(split(seq(nrow(snpattr)), match(snpattr$locus, names(loci))), names(loci))		
		for (i in seq_along(loci)) {
			attr(loci[[i]], "snp") <- snpattr[partloc[[i]], "snp"]
			attr(loci[[i]], "pis") <- as.logical(snpattr[partloc[[i]], "pis"])
			attr(loci[[i]], "nalleles") <- snpattr[partloc[[i]], "nalleles"]
		}
	}
	if (nloci == 1) {
		loci <- loci[[1]]
	}
	if (length(bedtable) > 0) {
		attr(loci, "bedtable") <- bedtable
	}
	
	return(loci)
}



#' Importing ipyrad output.
#' 
#' @description
#' Importing sequence data from ipyrad .alleles files.
#' 
#' @param file Character, the name of the .alleles file.
#' @param class Character, the output class, either `"matrix"` or `"DNAbin"`.
#' @param nloci The number of loci to be imported (the default `-1L` means all of them).
#' @param toupper Logical, whether to import sequences in upper case letters.
#' @param part Not relevant, just for compatibility with [import_seq].
#' @param rm.empty Not relevant, just for compatibility with [import_seq].
#' @param nlines The number of lines to be imported (the default `-1L` means all of them).
#'  An incomplete locus at the end of the imported lines is removed.
#' @param indiv Not relevant, just for compatibility with [import_seq].
#' @param interleaved Not relevant, just for compatibility with [import_seq].
#' @returns A list of locus-specific sequence alignments. If the assembly of RAD loci was by mapping
#'  to a reference genome, the output has an attribute `bedtable`.
#' @export

import_ipyrad <- function(file, nloci=-1L, class="matrix", toupper=TRUE, part, rm.empty, nlines=-1L, indiv, interleaved) {

	loci <- readLines(file, n=nlines)
	ends <- grep("//", loci)
	if (sign(nloci) == 1) {
		ends <- ends[seq(as.integer(nloci))]
		loci <- loci[1:max(ends)]
	}
	if (max(ends) < length(loci)) {
		loci <- loci[1:max(ends)]
		warning("the last locus was discarded")		
	}
	nams <- loci[ends]
	loci <- strsplit(loci[-ends], "\\s+")
	nloci <- length(ends)
	nseq <- diff(c(0, ends)) - 1
	seqnames <- split(sapply(loci, "[", 1), rep(seq(nloci), nseq))
	loci <- split(sapply(loci, "[", 2), rep(seq(nloci), nseq))
	loci <- lapply(loci, function(x) do.call(rbind, strsplit(x, "")))
	for (i in seq(nloci)) {
		rownames(loci[[i]]) <- seqnames[[i]]
	}
	nams <- gsub("[[:punct:][:space:]]*\\||\\|$", "", nams)

	if (class == "DNAbin") {
		loci <- lapply(loci, ape::as.DNAbin)
	}
	
	gc <- grepl(":", nams[1])
	if (isTRUE(gc)) {
		names(loci) <- paste0("loc", sub(":.*$", "", nams))
		bedtable <- do.call(rbind, strsplit(sub("^[[:digit:]]+:", "", nams), "[:-]"))
		bedtable <- data.frame(Scaffold=bedtable[,1], Start=as.numeric(bedtable[,2]), End=as.numeric(bedtable[,3]), Locus=names(loci), Score=0, Strand="+")
		bedtable$Start <- bedtable$Start - 1
		bedtable$End <- bedtable$End - 1		# only correction of apparent contradiction in ipyrad's .alleles output
		attr(loci, "bedtable") <- bedtable
	} else {
		names(loci) <- paste0("loc", nams)
		attr(loci, "bedtable") <- NULL
	}

	return(loci)

}



#' Importing stacks output.
#' 
#' @description
#' Importing sequence data from stacks . alleles.fas files.
#' 
#' @param file Character, the name of the .alleles.fas file.
#' @param class Character, the output class, either `"matrix"` or `"DNAbin"`.
#' @param nloci The number of loci to be imported (the default `-1L` means all of them).
#' @param toupper Logical, whether to import sequences in upper case letters.
#' @param part Not relevant, just for compatibility with [import_seq].
#' @param rm.empty Not relevant, just for compatibility with [import_seq].
#' @param nlines The number of lines to be imported (the default `-1L` means all of them).
#'  An incomplete locus at the end of the imported lines is removed.
#' @param indiv Information about individual names. It may be `"pop"` (individuals correspond to populations in 'Stacks'),
#'  a vector of individual names in an order corresponding to sample numbers in Stacks or `NULL` (default),
#'  which makes individuals being numbered as in Stacks.
#' @param interleaved Not relevant, just for compatibility with [import_seq].
#' @returns A list of locus-specific sequence alignments. If the assembly of RAD loci was by mapping
#'  to a reference genome, the output has an attribute `bedtable`.
#' @export

import_stacks <- function(file, nloci=-1L, class="matrix", toupper=TRUE, part, rm.empty, nlines=-1L, indiv=NULL, interleaved) {

	loci <- readLines(file, n=nlines)[-1]
	lrows <- grep("^>", loci)
	lmain <- do.call(rbind, strsplit(sub("\\s.*" , "", loci[lrows]), "_"))
	lsupp <- do.call(rbind, strsplit(gsub("^.+\\[|]" , "", loci[lrows]), ";\\s*"))
	if (grepl("^>", loci[length(loci)])) {
		discard <- which(lmain[,2] == lmain[nrow(lmain),2])
		loci <- loci[-seq(lrows[discard[1]], length(loci))]
		lrows <- lrows[-discard]
		lmain <- lmain[-discard,]
		lsupp <- lsupp[-discard,]
		warning("the last locus was discarded")
	}
	ends <- lrows[c(match(unique(lmain[,2])[-1], lmain[,2]) - 1, nrow(lmain))] + 1	
	if (sign(nloci) == 1) {
		ends <- ends[seq(as.integer(nloci))]
		loci <- loci[1:max(ends)]
	} else {
		nloci <- length(ends)
	}
	nseq <- diff(c(0, ends)) / 2
	if (is.null(indiv)) {
		seqnames <- paste(paste0("Sample_", lmain[,4]), lmain[,ncol(lmain)], sep="_")
	} else if (tolower(indiv) == "pop") {
		seqnames <- paste(lsupp[,1], lmain[,ncol(lmain)], sep="_")
	} else {
		seqnames <- paste(indiv[as.numeric(lmain[,4])], lmain[,ncol(lmain)], sep="_")
	}

	seqnames <- split(seqnames, rep(seq(nloci), nseq))
	loci <- split(loci[-lrows], rep(seq(nloci), nseq))
	loci <- lapply(loci, function(x) do.call(rbind, strsplit(x, "")))
	for (i in seq(nloci)) {
		rownames(loci[[i]]) <- seqnames[[i]]
	}
	names(loci) <- paste0("loc", unique(lmain[,2]))

	if (class == "DNAbin") {
		loci <- lapply(loci, ape::as.DNAbin)
	}

	gc <- ncol(lsupp) > 1
	if (isTRUE(gc)) {
		nbp <- sapply(loci, ncol)
		bedtable <- gsub("^;\\s|]$|[[:blank:]]", "", unique(lsupp[,2]))
		bedtable <- do.call(rbind, strsplit(bedtable, ","))
		bedtable <- data.frame(Scaffold=bedtable[,1], Start=as.numeric(bedtable[,2]), End=NA, Locus=names(loci), Score=0, Strand=bedtable[,3])
		plus <- bedtable$Strand == "+"
		minus <- !plus
		bedtable$Start[plus] <- bedtable$Start[plus] - 1
		bedtable$End[plus] <- bedtable$Start[plus] + nbp[plus]
		if (any(minus)) {
			start <- bedtable$Start[minus]
			bedtable$Start[minus] <- start - nbp[minus]
			bedtable$End[minus] <- start
		}	
		attr(loci, "bedtable") <- bedtable
	} else {
		attr(loci, "bedtable") <- NULL
	}

	return(loci)
	
}



#' Importing fasta files.
#' 
#' @description
#' Importing sequence data from files in .fasta format.
#' 
#' @param file Character, the name(s) of file(s) to be imported (each file corresponds to a single locus).
#' @param nloci The number of loci to be imported (the default `-1L` means all of them).
#' @param class Character, the output class, either `"matrix"` or `"DNAbin"`.
#' @param toupper Logical, whether to represent nucleotides by upper case letters.
#' @param part It can be a two column matrix indicating starts and ends of loci in the concatenated aligments
#'  or a name of partition file formatted according to RAxML or MrBayes convention.
#' @param rm.empty Logical, whether to remove empty (e.g. all "-" or all "N") sequences.
#' @param nlines Not relevant, just for compatibility with [import_seq].
#' @param indiv Not relevant, just for compatibility with [import_seq].
#' @param interleaved Not relevant, just for compatibility with [import_seq].
#' @returns A list of locus-specific alignments.
#' @export

import_fasta <- function(file, nloci=-1L, class="matrix", toupper=TRUE, part=NULL, rm.empty=FALSE, nlines, indiv, interleaved) {

	read_fasta <- function(file, fixclass, fixcase, aschar) {
		fixclass(fixcase(ape::read.dna(file, format="fasta", as.matrix=TRUE, as.character=aschar)))
	}

	fixcase <- list(base::identity, base::toupper)[[isTRUE(toupper) + 1]]
	fixclass <- list(base::identity, ape::as.DNAbin)[[(class == "DNAbin" & toupper) + 1]]
	aschar <- class == "matrix" | toupper
	if (is.null(part)) {
		nloci <- ifelse(nloci == -1L, length(file), nloci)
		if (nloci > length(file)) {
			nloci <- length(file)
			single <- (nloci > 1) + 1
			warning(paste("There", c("is","are")[single], "only", nloci, c("locus", "loci")[single], "available."))
		}
		loci <- lapply(file[1:nloci], read_fasta, fixclass=fixclass, fixcase=fixcase, aschar=aschar)
	} else {
		if (length(file) > 1) {
			warning("Only the first element of 'file' is used and partitioned according to the 'part' argument.")
		}
		if (is.character(part)) {
			part <- read_partition_file(part)
		}
		nloci <- ifelse(nloci == -1L, nrow(part), nloci)
		if (nloci > nrow(part)) {
			nloci <- nrow(part)
			single <- (nloci > 1) + 1
			warning(paste("There", c("is","are")[single], "only", nloci, c("locus", "loci")[single], "available."))
		} else if (nloci < nrow(part)) {
			part <- part[1:nloci,]
		}
		loci <- readLines(file[1])
		odd <- seq_along(loci) %% 2
		substrings <- function(i, j, x, n) setNames(substr(x, i, j), n)
		loci <- setNames(base::mapply(substrings, part$Start, part$End, MoreArgs=list(x=loci[odd == 0], n=loci[odd == 1]), SIMPLIFY=FALSE), rownames(part))
		loci <- lapply(loci, function(x) do.call(rbind, strsplit(x, "")))
	}

	if (isTRUE(rm.empty)) {
		loci <- remove_empty(loci, ambig=2)
	}
	if (isTRUE(toupper)) {
		loci <- lapply(loci, toupper)
	}
	if (nloci == 1) {
		loci <- unlist(loci, recursive=FALSE)
	}

	return(loci)
		
}



#' Importing phylip files.
#' 
#' @description
#' Importing sequence data from files in .fasta format.
#' 
#' @param file Character, the name(s) of file(s) to be imported (each file corresponds to a single locus).
#' @param nloci The number of loci to be imported (the default `-1L` means all of them).
#' @param class Character, the output class, either `"matrix"` or `"DNAbin"`.
#' @param toupper Logical, whether to represent nucleotides by upper case letters.
#' @param part It can be a two column matrix indicating starts and ends of loci in the concatenated aligments
#'  or a name of partition file formatted according to RAxML or MrBayes convention.
#' @param rm.empty Logical, whether to remove empty (e.g. all "-" or all "N") sequences.
#' @param nlines Not relevant, just for compatibility with [import_seq].
#' @param indiv Not relevant, just for compatibility with [import_seq].
#' @param interleaved Logical, whether in interleaved (not sequential) format.
#' @returns A list of locus-specific alignments.
#' @export

import_phylip <- function(file, nloci=-1L, class="matrix", toupper=TRUE, part=NULL, rm.empty=FALSE, nlines, indiv, interleaved=TRUE) {

	read_phylip <- function(file, fixclass, fixcase, aschar, format) {
		fixclass(fixcase(ape::read.dna(file, format=format, as.matrix=TRUE, as.character=aschar)))
	}

	fixcase <- list(base::identity, base::toupper)[[isTRUE(toupper) + 1]]
	fixclass <- list(base::identity, ape::as.DNAbin)[[(class == "DNAbin" & toupper) + 1]]
	aschar <- class == "matrix" | toupper
	input <- c("sequential", "interleaved")[interleaved + 1]
	if (is.null(part)) {
		nloci <- ifelse(nloci == -1L, length(file), nloci)
		if (nloci > length(file)) {
			nloci <- length(file)
			single <- (nloci > 1) + 1
			warning(paste("There", c("is","are")[single], "only", nloci, c("locus", "loci")[single], "available."))
		}
		loci <- lapply(file[1:nloci], read_phylip, fixclass=fixclass, fixcase=fixcase, aschar=aschar, format=format)
	} else {
		if (length(file) > 1) {
			warning("Only the first element of 'file' is used and partitioned according to the 'part' argument.")
		}
		if (is.character(part)) {
			part <- read_partition_file(part)
		}
		nloci <- ifelse(nloci == -1L, nrow(part), nloci)
		if (nloci > nrow(part)) {
			nloci <- nrow(part)
			single <- (nloci > 1) + 1
			warning(paste("There", c("is","are")[single], "only", nloci, c("locus", "loci")[single], "available."))
		} else if (nloci < nrow(part)) {
			part <- part[1:nloci,]
		}
		loci <- readLines(file[1])
		odd <- seq_along(loci) %% 2
		substrings <- function(i, j, x, n) setNames(substr(x, i, j), n)
		loci <- setNames(base::mapply(substrings, part$Start, part$End, MoreArgs=list(x=loci[odd == 0], n=loci[odd == 1]), SIMPLIFY=FALSE), rownames(part))
		loci <- lapply(loci, function(x) do.call(rbind, strsplit(x, "")))
	}

	if (isTRUE(rm.empty)) {
		loci <- remove_empty(loci, ambig=2)
	}
	if (isTRUE(toupper)) {
		loci <- lapply(loci, toupper)
	}
	if (nloci == 1) {
		loci <- unlist(loci, recursive=FALSE)
	}

	return(loci)

}



#' Sequence partitions.
#' 
#' @description
#' Reads in partition files defining locus boundaries.
#' 
#' @param file Character, the name of the partition file.
#' @details The file can use either RAxML or MrBayes coding of the partitions. Currently, it does not support
#'  partitioning by codon position.
#' @returns A list of locus-specific alignments. If genomic coordinates are included in a "bedtable" block of the file
#'  they are exported in .bed format as an attribute `bedtable`.
#' @export

read_partition_file <- function(file) {
	part <- readLines(file[1])
	part <- Filter(function(x) nchar(x) > 0, part)
	type <- ifelse(any(grepl("nexus", part, ignore.case=TRUE)), "mb", "raxml")
	if (type == "raxml") {
		part <- gsub("^DNA[[:alpha:]]*,|\\s", "", part)
	} else if (type == "mb") {
		part <- grep("charset", part, ignore.case=TRUE, value=TRUE)
	}		
	if (length(part) > 0) {
		part <- gsub("=", "-", part)
		part <- do.call(rbind, strsplit(part, "-"))
		part <- data.frame(start=as.numeric(part[,2]), end=as.numeric(part[,3]), row.names=part[,1])
	} else {
		part <- NULL
	}
	return(part)
}



#' Remove empty sequences.
#' 
#' @description
#' Removes sequencing 
#' 
#' @param loci Locus-specific sequence alignement or a list of them.
#' @param ambig Numeric, how ambiguous can be valid data, i.e., how many different nucleotides can be coded 
#'  by a single character. The default is 2 (~ heterozygote in a diploid genotype). More ambiguous characters
#'  are considered missing data.
#' @returns A list of locus-specific alignments with empty sequences removed.
#'  If there is no non-empty sequence in a locus, an empty matrix is returned.
#' @export

remove_empty <- function(loci, ambig=2) {
	rmempty <- function(x, n) x[rowSums(matrix(as.character(x) %in% n, nrow(x), ncol(x))) > 0,,drop=FALSE]
	nucleotides <- unlist(list(c("A", "C", "G", "T"), c("R", "Y", "S", "W", "K", "M"), c("B", "D", "H", "V"))[1:ambig])
	nucleotides <- c(nucleotides, tolower(nucleotides))
	loci <- lapply(loci, rmempty, n=nucleotides)
	return(loci)
}
