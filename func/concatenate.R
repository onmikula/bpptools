# Function:
#	concatenate: concatenates multi-locus sequences reads sequence from (partitioned) .nexus file
# Arguments:
#	dna: a list of alignements in matrix or DNAbin format
#	toupper: whether to export sequences in upper case letters
#	as.DNAbin: whether to export the alignment in 'DNAbin' format
#	missing: symbol representing missing nucleotides
# Value:
#	a concatenated alignment in matrix or DNAbin format with 'part' attribute defining locus boundaries
# Depends on: ape

concatenate <- function(dna, toupper=TRUE, as.DNAbin=FALSE, missing="-") {
	bp <- cumsum(sapply(dna, ncol))
	part <- attr(dna, "part")
	if (is.null(part)) {
		part <- cbind(c(1, bp[-length(bp)] + 1), unname(bp))
		rownames(part) <- names(bp)
	}
	dna <- lapply(dna, as.matrix)
	rows <- sort(unique(unlist(lapply(dna, rownames))))
	concatenated <- matrix(missing, length(rows), max(bp), dimnames=list(rows, NULL))
	for (i in seq_along(dna)) {
		concatenated[rownames(dna[[i]]), part[i,1]:part[i,2]] <- dna[[i]]
	}
	if (isTRUE(toupper)) {
		concatenated <- toupper(concatenated)
	}
	if (isTRUE(as.DNAbin)) {
		concatenated <- ape::as.DNAbin(concatenated)
	}
	attr(concatenated, "part") <- part
	return(concatenated)
}
