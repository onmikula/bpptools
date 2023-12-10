### FUNCTION
# read_fasta
### ARGUMENTS
# file: name of '.fasta' file
# class: output clas, either 'matrix', 'DNAbin' or 'DNAStringSet' (if class == "matrix" and sequence lengths are unequal, the output is a list of character vectors)
# nloci: number of loci to be read in (default '-1L' meaning 'all')

read_fasta <- function(file, class="matrix", nloci=-1L) {
	if (class == "DNAStringSet") {
		seq <- Biostrings::readDNAStringSet(file, format="fasta", nrec=nloci, skip=0L, use.names=TRUE)
	} else {
		seq <- lapply(ape::read.dna(file, format="fasta", as.character=TRUE), toupper)
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



### FUNCTION
# read_ipyrad
### ARGUMENTS
# file: name of '.alleles' ipyrad output file with phased sequences of discrete ddRAD loci
# class: output clas, either 'matrix', 'DNAbin' or 'DNAStringSet'
# nloci: number of loci to be read in (default '-1L' meaning 'all')
# nlines: number of lines to be read in (default '-1L' meaning 'all')

read_ipyrad <- function(file, class="matrix", nloci=-1L, nlines=-1L) {

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



### FUNCTION
# read_stacks
### ARGUMENTS
# file: name of '.alleles.fas' STACKS output file with phased sequences of discrete ddRAD loci
# class: output clas, either 'matrix', 'DNAbin' or 'DNAStringSet'
# nloci: number of loci to be read in (default '-1L' meaning 'all')
# nlines: number of lines to be read in (default '-1L' meaning 'all')
# indiv: information about individual names, may be "pop" (individuals correspond to populationsin 'Stacks'), a vector of individual  names in the order corresponding to sample numbers in 'Stacks' or NULL (default), which makes individuals being numbered as in 'Stacks'

read_stacks <- function(file, class="matrix", nloci=-1L, nlines=-1L, indiv=NULL) {

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



### FUNCTION
#	read_nexus: reads sequence from (partitioned) .nexus file
### ARGUMENTS
#	file: name of the partitioned .nexus file
#	toupper: whether to export sequences in upper case letters
#	as.DNAbin: whether to export sequences as 'DNAbin' objects
# Value:
#	(a list of) alignment(s) in matrix or DNAbin format (sigle alignment if no partitions are defined in the .nexus file)
# Depends on: ape, Biostrings

read_nexus <- function(file, class="matrix", nloci=-1L, toupper=TRUE) {
	loci <- readLines(file)
	loci <- loci[!grepl("^[[:space:]]+$", loci)]
	loci <- loci[nchar(loci) > 0]
	start <- which(loci %in% c("MATRIX", "Matrix", "matrix")) + 1
	end <- base::intersect(which(loci == ";"), start:length(loci))[1] - 1
	part <- grep("charset", loci[(end + 1):length(loci)], ignore.case=TRUE, value=TRUE)
	if (length(part) > 0) {
		if (sign(nloci) == 1) {
			part <- part[1:nloci]
		}
		part <- gsub("charset|[[:blank:]]|;", "", part, ignore.case=TRUE)
		part <- do.call(rbind, strsplit(part, "="))
		part <- setNames(data.frame(part[,1], do.call(rbind, lapply(strsplit(part[,2], "-"), as.numeric)), stringsAsFactors=FALSE), c("Locus", "Start", "End"))
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
