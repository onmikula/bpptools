#' Sequence concatenation.
#' 
#' @description
#' Concatenates multi-locus sequences into a single alignement, while preserving their atributes.
#' 
#' @param seq The list of alignments with row names indicating individuals
#'   or a single such alignment (possibly concatenated multi-locus data).
#' @param nloci The number of loci to be imported (counted from the beginning). The default means "all loci".
#' @param class The class of output object, either 'matrix', 'DNAbin' or 'DNAStringSet'.
#' @param toupper Logical, whether to represent bases by upper case letters. Always `TRUE` if `allele` is specified.
#' @param missing Character, how to represent missing data.
#' @param allele The allele identifier distinguishing haplotypes of the same individual.
#'   If supplied, the individual diploid genotypes are collapsed into a single sequence.
#' @returns A single concatenated alignment with attributes `part` (a two-column matrix, which defines
#'   paritioning into loci) and possibly `bedtable` and `snps` (if these were included in the input data).
#' @export

concatenate <- function(seq, nloci=-1L, class="matrix", toupper=TRUE, missing="-", allele) {

	bedtable <- attr(seq, "bedtable")
	if (sign(nloci) == 1) {
		nloci <- min(c(nloci, length(loci)))
		seq <- seq[1:nloci]
		if (!is.null(bedtable)) {
			bedtable <- bedtable[1:nloci,,drop=FALSE]
		}
	}
	bp <- cumsum(sapply(seq, ncol))
	part <- attr(seq, "part")
	if (is.null(part)) {
		part <- cbind(c(1, bp[-length(bp)] + 1), unname(bp))
		rownames(part) <- names(bp)
	}
	if (!is.null(attr(seq[[1]], "snp"))) {
		snps <- lapply(seq, function(x) cbind(snp=attr(x, "snp"), pis=as.character(attr(x, "pis")), nalleles=attr(x, "nalleles")))
		snps <- cbind(locus=rep(names(seq), sapply(snps, nrow)), do.call(rbind, snps))
	} else {
		snps <- NULL
	}
	if (inherits(seq[[1]], "DNAbin")) {
		seq <- lapply(seq, as.character)
	}
	rows <- sort(unique(unlist(lapply(seq, rownames))))
	concatenated <- matrix(missing, length(rows), max(bp), dimnames=list(rows, NULL))
	for (i in seq_along(seq)) {
		concatenated[rownames(seq[[i]]), part[i,1]:part[i,2]] <- seq[[i]]
	}

	if(!missing(allele)) {
		glue <- function(x) apply(x, 2, paste, collapse="")
		ambig <- c(AA="A", CC="C", GG="G", TT="T", "A-"="A", "C-"="C", "G-"="G", "T-"="T", "-A"="A", "-C"="C", "-G"="G", "-T"="T", AC="M", CA="M", AG="R", GA="R", AT="W", TA="W", CG="S", GC="S", CT="Y", TC="Y", GT="K", TG="K", "NN"="N", "N-"="N", "-N"="N", "--"="-")
		phased <- grepl(allele, rownames(concatenated))
		nind <- sum(phased) / 2
		if (nind >= 1) {
			unphased <- concatenated[!phased,,drop=FALSE]
			phased <- toupper(concatenated[phased,,drop=FALSE])
			phased <- phased[order(rownames(phased)),,drop=FALSE]
			phased[!phased %in% c("A","C","G","T")] <- "-"
			indnames <- sub(allele, "", rownames(phased)[2*seq(nind)])
			phased <- lapply(seq(nind), function(i) unname(ambig[glue(phased[(2 * i) + c(-1, 0),,drop=FALSE])]))
			phased <- setNames(phased, indnames)
			phased <- do.call(rbind, phased)
			ord <- match(unique(sub(allele, "", rownames(concatenated))), c(rownames(unphased), rownames(phased)))
			concatenated <- rbind(unphased, phased)[ord,,drop=FALSE]
		} else {
			toupper <- TRUE
		}
	}

	if (isTRUE(toupper)) {
		concatenated <- toupper(concatenated)
	}
	if (class == "DNAbin") {
		concatenated <- ape::as.DNAbin(concatenated)
	} else if (class == "DNAStringSet") {
		concatenated <- Biostrings::DNAStringSet(unlist(lapply(split(concatenated, rownames(concatenated)), paste, collapse="")))
	}
	attr(concatenated, "snps") <- snps
	attr(concatenated, "bedtable") <- bedtable
	attr(concatenated, "part") <- part
	return(concatenated)	
}
