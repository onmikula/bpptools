#' Average genetic distances.
#' 
#' @description
#' Calculates average genetic distance between the specified species and within them.
#' 
#' @param seq The list of locus-specific alignments or a single such alignment (possibly 
#'   with concatenated multi-locus data). The row names must be interpretable as individual labels
#'   (possibly with help of `seqmap` or `allele` arguments).
#' @param imap A matrix or data frame (1st column individuals, 2nd column species) or a name of imap file.
#' @param seqmap A data frame mapping names of sequences (1st column) to those of individuals (2nd column).
#' @param allele An allele identifier, regular expression distinguishing sequences from the same individual.
#' @param model The nucleotide substitution model indicated as in [ape::dist.dna()].
#' @param diag Logical, whether to calculate average intraspecific distances, default is `TRUE`.
#' @returns A square symmetric matrix giving mean genetic distance between species (off-diagonal)
#'   and within them (diagonal).
#' @export

avedist <- function(seq, imap=NULL, seqmap=NULL, allele=NULL, model="JC69", diag=TRUE) {

	if (is.list(seq)) {
		seq <- concatenate(seq, class="DNAbin")
	} else if (!inherits(seq, "DNAbin")) {
		seq <- ape::as.DNAbin(seq)
	}

	im <- !is.null(imap)
	sm <- !is.null(seqmap)
	al <- !is.null(allele)
	if (im) {
		if (is.character(imap)) {
			imap <- read_bpp_imap(imapfile, col.names=c("Individual", "Species"))
		} else {
			imap <- setNames(as.data.frame(imap[,1:2]), c("Individual", "Species"))
		}
		if (!sm) {
			if (al) {
				seqmap <- data.frame(Sequence=paste0(imap[,1], allele), Individual=imap[,1])
			} else {
				seqmap <- data.frame(Sequence=imap[,1], Individual=imap[,1])
			}
		} else if (al) {
			seqmap <- setNames(as.data.frame(seqmap[,1:2]), c("Sequence", "Individual"))
			seqmap <- unique(rbind(seqmap, data.frame(Sequence=paste0(imap[,1], allele), Individual=imap[,1])))
			seqnames <- lapply(seqmap$Sequence, grep, x=rownames(seq), value=TRUE)
			seqmap <- unique(data.frame(Sequence=unlist(seqnames), Individual=rep(seqmap$Individual, sapply(seqnames, length))))
		}
	} else {
		if (!sm) {
			if (al) {
				indnames <- sub(allele, "", rownames(seq))
				imap <- data.frame(Individual=unique(indnames), Species=NA)
				seqmap <- data.frame(Sequence=rownames(seq), Individual=indnames)
			} else {
				imap <- data.frame(Individual=rownames(seq), Species=NA)
				seqmap <- data.frame(Sequence=rownames(seq), Individual=rownames(seq))
			}
		} else {
			seqmap <- setNames(as.data.frame(seqmap[,1:2]), c("Sequence", "Individual"))
			imap <- data.frame(Individual=seqmap$Individual, Species=NA)
		}
	}

	indnames <- seqmap$Individual[match(rownames(seq), seqmap$Sequence)]
	sp <- factor(imap$Species[match(indnames, imap$Individual)])
	dst <- ape::dist.dna(seq, model=model, pairwise.deletion=TRUE, as.matrix=TRUE)
#	dst[which(dst == 0)] <- 1e-8
	ave <- matrix(, nlevels(sp), nlevels(sp), dimnames=list(levels(sp), levels(sp)))
	for (i in seq(nlevels(sp) - 1)) {
		ii <- sp == levels(sp)[i]
		for (j in (i + 1):nlevels(sp)) {
			jj <- sp == levels(sp)[j]
			ave[i,j] <- ave[j,i] <- mean(dst[ii,jj], na.rm=TRUE)
		}
	}
	
	if (isTRUE(diag)) {
		for (i in seq(nlevels(sp))) {
			ii <- sp == levels(sp)[i]
			if (sum(ii) > 1) {
				ave[i,i] <- mean(dst[ii,ii][lower.tri(diag(sum(ii)))], na.rm=TRUE)
			}
		}
	}

	return(ave)

}
