#' Prior distribution parameters.
#' 
#' @description
#' Estimate means of theta and tau prior distributions which can be considered reasonable.
#' 
#' @param seq The list of locus-specific alignments or a single such alignment (possibly 
#'   with concatenated multi-locus data). The row names must be interpretable as individual labels
#'   (possibly with help of `seqmap` or `allele` arguments).
#' @param imap A matrix or data frame (1st column individuals, 2nd column species) or a name of imap file.
#' @param d_tau,d_theta The type of prior distributions for tau and theta, either `"gamma"` or `"invgamma"`.
#' @param a_tau,a_theta Alpha parameters of tau and theta priors.
#' @param seqmap A data frame mapping names of sequences (1st column) to those of individuals (2nd column).
#' @param allele An allele identifier, regular expression distinguishing sequences from the same individual.
#' @param model The nucleotide substitution model indicated as in [ape::dist.dna()].
#' @details Any analysis using priors informed by the analyzed data contains some circularity of reasoning.
#'   Still, the strategy is sometimes used to set prior parameters to resonable values, but in that case
#'   the priors should be diffuse (not too informative).
#' @returns A data frame with columns `"alpha"`, `"beta"`, `"mean"` & `"distribution"`,
#'   the latter giving the type of distribution.
#' @export

estimate_priors <- function(seq, imap=NULL, d_tau, d_theta, a_tau, a_theta, seqmap=NULL, allele=NULL, model="JC69") {

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
				seqnames <- paste0(imap[,1], allele)
				seqnames <- lapply(seqnames, grep, x=rownames(seq), value=TRUE)
				seqmap <- data.frame(Sequence=unlist(seqnames), Individual=rep(imap[,1], sapply(seqnames, length)))
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

	if (!im) {
		otus <- estimate_otus(seq, seqmap=seqmap, allele=allele, model=model)
		imap$Species <- otus$Species[match(imap$Individual, otus$Individual)]
	}

	dst <- ape::dist.dna(seq[seqmap$Sequence,], pairwise.deletion=TRUE, model=model, as.matrix=TRUE)
	dst[dst == 0] <- 1e-8
	otus <- split(seqmap$Sequence, imap$Species[match(seqmap$Individual, imap$Individual)])
	otus <- lapply(otus, function(x) sort(unique(x)))
	tau <- matrix(,length(otus), length(otus), dimnames=list(names(otus), names(otus)))
	for (i in seq(length(otus) - 1)) {
		for (j in (i+1):length(otus)) {
			tau[i,j] <- mean(dst[otus[[i]], otus[[j]]], na.rm=TRUE)
		}
	}
	m_tau <- max(tau, na.rm=TRUE) / 2
	otus <- rep(seq_along(otus), sapply(otus, length))[match(rownames(dst), unlist(otus))]
	m_theta <- mean(dst[outer(otus, otus, "!=") & !diag(length(otus))], na.rm=TRUE)
	means <- c(m_tau, m_theta)
	b_tau <- ifelse(d_tau == "invgamma", (a_tau - 1) * m_tau, a_tau / m_tau)
	b_theta <- ifelse(d_theta == "invgamma", (a_theta - 1) * m_theta, a_theta / m_theta)

	priors <- data.frame(alpha=round(c(a_tau, a_theta), 5), beta=round(c(b_tau, b_theta), 5), mean=round(c(m_tau, m_theta), 5), distribution=c(d_tau, d_theta), row.names=c("tau", "theta"))

	return(priors)

}



#' Operational taxonomic units.
#' 
#' @description
#' Provides delimitation of operational taxonomic units (OTUs, provisional species) using
#' the maximum cross-section branch length criterion on ad hoc created average linkage (UPGMA) dendrogram.
#' 
#' @param seq The list of locus-specific alignments or a single such alignment (possibly 
#'   with concatenated multi-locus data). The row names must be interpretable as individual labels
#'   (possibly with helpof `allele` or `seqmap arguments`).
#' @param seqmap A data frame mapping names of sequences (1st column) to those of individuals (2nd column).
#' @param allele An allele identifier, regular expression distinguishing sequences from the same individual.
#' @param model The nucleotide substitution model indicated as in [ape::dist.dna()].
#' @details The function provides provisional delimitation of operational taxonomic units
#'   intended to be used as a rough guide to inform pre-estimation of prior parameters.
#'   It uses a simple "maximum cross-section length" method, which bulids ultrametric UPGMA dendrogram
#'   from the matrix of pairwise inter-individual genetic distances and then cuts it into clusters
#'   at some distance from the tips. The distance is chosen to maximize the average length of branches
#'   that are cut by the cross-section line. 
#' @returns A data frame with columns `"Individuals"` and `"Species"`.
#' @export

estimate_otus <- function(seq, seqmap=NULL, allele=NULL, model="JC69") {
	sc <- function(k, h, b) {
		mean(b[nh[,2] > k & nh[,1] <= k])
	}

	if (is.list(seq)) {
		seq <- concatenate(seq, class="DNAbin")
	} else if (!inherits(seq, "DNAbin")) {
		seq <- ape::as.DNAbin(seq)
	}
	
	sm <- !is.null(seqmap)
	al <- !is.null(allele)
	if (!sm) {
		if (al) {
			indnames <- sub(allele, "", rownames(seq))
			seqmap <- data.frame(Sequence=rownames(seq), Individual=indnames)
		} else {
			seqmap <- data.frame(Sequence=rownames(seq), Individual=rownames(seq))
		}
	} else {
		seqmap <- setNames(as.data.frame(seqmap[,1:2]), c("Sequence", "Individual"))
	}

	dst <- avedist(seq, imap=seqmap, seqmap=NULL, allele=NULL, model=model, diag=FALSE)
	dst[dst == 0] <- 1e-8
	diag(dst) <- 0
	
	tre <- ape::as.phylo(stats::hclust(as.dist(dst), method="average"))
	ntip <- ape::Ntip(tre)
	nh <- matrix(ape::plotPhyloCoor(tre, direction="rightwards")[,1][tre$edge], ape::Nedge(tre), 2)
	knots <- sort(setNames(nh[,2], tre$edge[,2])[!duplicated(nh[,2])])
	knots <- setNames(c(0, knots), c(ntip+1, names(knots)))
	k <- knots[which.max(sapply(knots, sc, h=nh, b=tre$edge.length))]
	anc <- tre$edge[nh[,2] > k & nh[,1] <= k, 2]
	single <- unlist(Filter(function(x) length(x) > 0, as.list(anc[anc <= ntip])))
	otus <- lapply(setdiff(anc, single), function(a) ape::extract.clade(tre, a)$tip.label)
	otus <- c(as.list(tre$tip.label[single]), otus)
	otus <- data.frame(Individual=unlist(otus), Species=paste0("OTU", rep(seq_along(otus), sapply(otus, length))))

	return(otus)
}



#' Starting tree for species tree estimation.
#' 
#' @description
#' Estimates starting tree for species tree estimation (analyses A11 or A01)
#' from average genetic distances between (candidate) species.
#' 
#' @param seq The list of locus-specific alignments or a single such alignment (possibly 
#'   with concatenated multi-locus data). The row names must be interpretable as individual labels
#'   (possibly with helpof `allele` or `seqmap arguments`). Alternatively, it can be a matrix of
#'   average genetic distances between (candidate) species, supplied as an object of class `dist`
#'   or a symmetric matrix of mode `numeric`.
#' @param method The tree estimation method, either `"nj"` (neighbor joining) or `"upgma"`
#'   (unweighted pair group method with arithmetic mean); corresponds to options `"single"`
#'   and `"average"`, respectively, in [stats::hclust()].
#' @param imap A matrix or data frame (1st column individuals, 2nd column species) or a name of imap file.
#' @param seqmap A data frame mapping names of sequences (1st column) to those of individuals (2nd column).
#' @param allele An allele identifier, regular expression distinguishing sequences from the same individual.
#' @param model The nucleotide substitution model indicated as in [ape::dist.dna()].
#' @returns The starting tree topology in Newick format.
#' @export

estimate_starting_tree <- function(seq, method=c("nj","upgma"), imap=NULL, seqmap=NULL, allele=NULL, model="JC69") {

	is_symmetric <- function(m, r) {
		if (!is.matrix(m)) {
			return(FALSE)
		} else if (nrow(m) != ncol(m)) {
			return(FALSE)
		} else {
			m <- round(m, r)
			return(all(m == t(m)))
		}
	}

	if (inherits(seq, "dist")) {
		tau <- seq
	} else if (is_symmetric(m, 8)) {
		tau <- as.dist(seq)
	} else {
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
					seqnames <- paste0(imap[,1], allele)
					seqnames <- lapply(seqnames, grep, x=rownames(seq), value=TRUE)
					seqmap <- data.frame(Sequence=unlist(seqnames), Individual=rep(imap[,1], sapply(seqnames, length)))
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
		ave <- matrix(0, nlevels(sp), nlevels(sp), dimnames=list(levels(sp), levels(sp)))
		for (i in seq(nlevels(sp) - 1)) {
			ii <- sp == levels(sp)[i]
			for (j in (i + 1):nlevels(sp)) {
				jj <- sp == levels(sp)[j]
				ave[i,j] <- ave[j,i] <- mean(dst[ii,jj], na.rm=TRUE)
			}
		}
		ave <- as.dist(ave)
	}

	method <- c(nj="single", upgma="average", single="single", average="average")[method]
	tree <- ape::as.phylo(stats::hclust(ave, method=method))
	tree$edge.length <- NULL
	tree$node.label <- NULL
	ape::write.tree(tree)

}
