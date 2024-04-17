#' Importing sequence data.
#' 
#' @description
#' The functions importing sequence data from various file formats, typically splitting them into discrete loci.
#'
#' `import_seq` is a wrapper for the following functions.
#'
#' `import_nexus` from .nexus file, loci defined by an appropriate nexus block or an optional argument.
#'
#' `import_fasta` from .fasta file, loci defined by an optional argument.
#'
#' `import_ipyrad` from .alleles file of ipyrad RADseq assembler.
#'
#' `import_stacks` from .alleles file of stacks RADseq assembler.
#' 
#' @param file The input file name.
#' @param format A character string indicating the file format (`"nexus"`, `"fasta"`, `"ipyrad"` or `"stacks"`).
#' @param nloci The number of loci to be imported (counted from the beginning). The default means "all loci".
#' @param class The class of output object, either 'matrix', 'DNAbin' or 'DNAStringSet';
#'   if `class == "matrix"` and sequence lengths are unequal, the list of character vectors is returned.
#' @param toupper Logical, whether to represent bases by upper case letters.
#' @param part It can be a two column matrix indicating starts and ends of loci in the concatenated aligments
#'   (in .fasta or .nexus files) or a name of partition file in format used by RAxML or MrBayes.
#'   In `import_nexus` it is overridden by the nexus block defining charsets (if present within the file).
#' @param nlines The number of lines to be read by `import_ipyrad`. The default means "all lines".
#'   The last locus is omitted, if it is found incomplete.
#' @param indiv The information about individual names for `import_stacks`.
#'   It can be a vector of individual  names in the order corresponding to sample numbers of Stacks
#'   or the character string `"pop"`, if individuals are equivalent to populations of Stacks.
#'   If NULL (default), the sample numbers are used as individual labels.
#' @returns A list of locus-specific alignments. If genomic coordinates are included (in ipyrad or Stacks output),
#'   they are exported in .bed format as an attribute `bedtable`.
#' @export
import_seq <- function(file, format, nloci=-1L, class="matrix", toupper=FALSE, part=NULL, nlines=-1L, indiv=NULL) {
	format <- tolower(format)
	if (format == "nexus") {
		loci <- import_nexus(file=file, nloci=nloci, class=class, toupper=toupper, part=part)
	}
	if (format == "fasta") {
		loci <- import_fasta(file=file, nloci=nloci, class=class, toupper=toupper, part=part)
	}
	if (format == "ipyrad") {
		loci <- import_ipyrad(file=file, nloci=nloci, class=class, toupper=toupper, nlines=nlines)
	}
	if (format == "stacks") {
		loci <- import_stacks(file=file, nloci=nloci, class=class, toupper=toupper, nlines=nlines, indiv=indiv)
	}
	return(loci)	
}


import_nexus <- function(file, nloci=-1L, class="matrix", toupper=FALSE, part=NULL) {

	loci <- readLines(file)
	loci <- loci[!grepl("^[[:space:]]+$", loci)]
	loci <- loci[nchar(loci) > 0]
	start <- which(loci %in% c("MATRIX", "Matrix", "matrix")) + 1
	end <- base::intersect(which(loci == ";"), start:length(loci))[1] - 1
	charsets <- grep("charset", loci[(end + 1):length(loci)], ignore.case=TRUE, value=TRUE)
	if (length(charsets) > 0) {
		if (sign(nloci) == 1) {
			charsets <- charsets[1:nloci]
		}
		part <- gsub("charset|[[:blank:]]|;", "", charsets, ignore.case=TRUE)
		part <- do.call(rbind, strsplit(part, "="))
		part <- setNames(data.frame(part[,1], do.call(rbind, lapply(strsplit(part[,2], "-"), as.numeric)), stringsAsFactors=FALSE), c("Locus", "Start", "End"))
	} else if (is.character(part)) {
		part <- readLines(part[1])
		part <- Filter(function(x) nchar(x) > 0, part)
		type <- ifelse(any(grepl("nexus", part, ignore.case=TRUE)), "mb", "raxml")
		if (type == "raxml") {
			part <- gsub("^[[:alpha:]]*,|\\s", "", part)
			part <- do.call(rbind, strsplit(part, "="))
		} else if (type == "mb") {
			part <- grep("charset", part, ignore.case=TRUE, value=TRUE)
			if (length(part) > 0) {
				part <- gsub("charset|[[:blank:]]|;", "", part, ignore.case=TRUE)
				part <- do.call(rbind, strsplit(part, "="))
				part <- setNames(data.frame(part[,1], do.call(rbind, lapply(strsplit(part[,2], "-"), as.numeric)), stringsAsFactors=FALSE), c("Locus", "Start", "End"))
			}
		}
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
	if (sign(nloci) == 1) {
		last <- part$End[nloci]
		for (i in which(odd == 0)) {
			loci[i] <- substr(loci[i], 1, last)
		}
	}
	if (length(part) == 0) {
		nloci <- 1
		loci <- list(loci)
	} else {
		nloci <- nrow(part)
		substrings <- function(i, j, x, n) setNames(substr(x, i, j), n)
		loci <- setNames(base::mapply(substrings, part$Start, part$End, MoreArgs=list(x=loci[odd == 0], n=loci[odd == 1]), SIMPLIFY=FALSE), part$Locus)
	}
	
	if (class == "matrix") {
		case <- list(base::identity, base::toupper)[[isTRUE(toupper) + 1]]
		for (i in seq_along(loci)) {
			loci[[i]] <- case(as.matrix(do.call(rbind, strsplit(loci[[i]], ""))))
		}
	} else if (class == "DNAbin") {
		for (i in seq_along(loci)) {
			loci[[i]] <- ape::as.DNAbin(strsplit(loci[[i]], ""))
		}			
	} else if (class == "DNAStringSet") {
		for (i in seq_along(loci)) {
			loci[[i]] <- Biostrings::DNAStringSet(loci[[i]])
		}
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
		loci <- unlist(loci, recursive=FALSE)
	}
	
	if (length(bedtable) > 0) {
		attr(loci, "bedtable") <- bedtable
	}
	
	return(loci)
}


import_fasta <- function(file, nloci=-1L, class="matrix", toupper=FALSE, part=NULL) {
	if (class == "DNAStringSet") {
		seq <- Biostrings::readDNAStringSet(file, format="fasta", nrec=nloci, skip=0L, use.names=TRUE)
	} else {
		seq <- ape::read.dna(file, format="fasta", as.character=TRUE)
		if (isTRUE(toupper)) {
			seq <- lapply(seq, toupper)
		}
		if (sign(nloci) == 1) {
			seq <- seq[1:nloci]
		}
		if (class == "matrix" & length(unique(sapply(seq, length))) == 1) {
			seq <- do.call(rbind, seq)
		} else if (class == "DNAbin") {
			seq <- lapply(seq, ape::as.DNAbin)
		} 
	}
	return(seq)		
}


import_ipyrad <- function(file, nloci=-1L, class="matrix", toupper=FALSE, nlines=-1L) {

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

	if (isTRUE(toupper)) {
		loci <- lapply(loci, toupper)
	}
	if (class == "DNAbin") {
		loci <- lapply(loci, ape::as.DNAbin)
	} else if (class == "DNAStringSet") {
		loci <- lapply(loci, function(x) Biostrings::DNAStringSet(unlist(lapply(split(x, rownames(x)), paste, collapse=""))))
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


import_stacks <- function(file, nloci=-1L, class="matrix", toupper=FALSE, nlines=-1L, indiv=NULL) {

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

	if (isTRUE(toupper)) {
		loci <- lapply(loci, toupper)
	}
	if (class == "DNAbin") {
		loci <- lapply(loci, ape::as.DNAbin)
	} else if (class == "DNAStringSet") {
		loci <- lapply(loci, function(x) Biostrings::DNAStringSet(unlist(lapply(split(x, rownames(x)), paste, collapse=""))))
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
