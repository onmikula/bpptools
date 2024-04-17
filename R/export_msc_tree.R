#' Multispecies coalescent tree.
#' 
#' @description
#' Exports species tree annotated by the multispecies coalescent parameters in nexus metacomment format.
#' 
#' @param target The tree topology in Newick format obtained by `read_bpp_tree`
#'   or a name of BPP control (.ctl) or output (.out) file.
#' @param trees A posterior sample of rooted trees in a `multiPhylo` object or a list of `phylo` objects
#'   or a name of file with the posterior sample.
#' @param params A data frame with posterior sample of model parameters.
#' @param file A name of the file with the final MSC-annotated tree.
#' @param treesfile A file name for export of the annotated posterior sample of trees.
#' @param fun The summary function, `"mean"` (default) or `"median"`.
#' @param interval Either `"HPD"` (default) or `"CPD"`.
#' @param p The proportion of posterior density be included in the interval, default is 0.95.
#' @param figtree Logical, whether to export the final tree and the posterior sample.
#'   as .nexus files readable by FigTree.
#' @param constraint a character vector specifying bottom-up succession of tips when plotting the tree
#' @returns `read_msc_tree` returns an object of class `phylo` with an additional component `params`, which is
#'   a data frame including  MSC annotations: divergence times (tau) and mutation-scaled population sizes (theta).
#'   Every row corresponds to a single node and the parameters relate to the branch leading to the node.
#'   The root species has its theta but not tau.
#' @export

export_msc_tree <- function(target, trees, params, file, treesfile, fun="mean", interval="HPD", p=0.95, gdi=FALSE, figtree=FALSE, constraint) {
	make_figtree <- function(newick, tips) {
		ntip <- length(tips)
		transl <- paste(paste("\t\t", seq(ntip), sep=""), paste0(tips, ","))
		transl[length(transl)] <- sub("\\,$", "", transl[length(transl)])
		transl <- c("\tTranslate", transl, "\t;")
		for (i in seq(ntip)) {
			pattern <- paste0("[[:punct:]]{1}", tips[i], "[[:punct:]]{1}")
			for (j in seq_along(newick)) {
				m <- regexpr(pattern, newick[j]) + 1
				attr(m, "match.length") <- attr(m, "match.length") - 2
				regmatches(newick[j], m) <- i
			}
		}
		taxa <- c("Begin taxa;", paste("\tDimensions ntax=", ntip, ";", sep=""), "\t\tTaxlabels", paste("\t\t", tips, sep=""), "\t;", "End;")
		no <- formatC(seq_along(newick), format="d", flag=0, width=nchar(length(newick)))
		return(c("#NEXUS", "\n", taxa, "\n", "Begin trees;", transl, paste("tree", no, "=", newick), "End;"))
	}
	fun <- match.fun(fun)
	int <- match.fun(interval)
	
	if (substr(target, 1, 1) != "(") {
		target <- read_bpp_tree(target)
	}
	tips <- attr(target, "tips")
	target <- get_topology(target)
	constraint <- ape::read.tree(text=target)$tip.label
	ntip <- length(tips)
	nodes <- seq(ntip - 1) + ntip
	if (is.character(trees)) {
		trees <- ape::read.tree(trees)
	}
	topologies <- sapply(lapply(trees, ape::rotateConstr, constraint=constraint), get_topology)
	subset <- topologies == target
	if (!any(subset)) {
		stop("none of the trees has topology identical to the target")
	}
	nsample <- sum(subset)
	target <- ape::read.tree(text=target)
	target$node.label <- nodes
	params <- params[subset,,drop=FALSE]
	params <- params[,grepl("^theta_*[[:digit:]]+|^tau_*[[:digit:]]+", colnames(params)),drop=FALSE]
	colnames(params) <- regmatches(colnames(params), regexpr("^theta_*[[:digit:]]+|^tau_*[[:digit:]]+", colnames(params)))
	thetanum <- gsub("theta_", "", colnames(params)[seq(ntip)])
	thetanum <- intersect(as.numeric(thetanum), seq(ntip))
	thetamis <- setdiff(seq(ntip), thetanum)
	if (length(thetamis) > 0) {
		thetatip <- data.frame(params[,seq_along(thetanum)], matrix(NA, nrow(params), length(thetamis)))
		names(thetatip) <- c(names(params)[seq_along(thetanum)], paste0("theta_", thetamis))
		params <- cbind(thetatip[,order(c(thetanum, thetamis))], params[,-seq_along(thetanum)])
	}
	tau <- cbind(as.data.frame(setNames(rep(list(numeric(nsample)), ntip),paste0("tau",seq(ntip)))), params[,grepl("tau",colnames(params))])
	theta <- params[,grepl("theta",colnames(params))]
	theta[,seq(ntip)] <- theta[,match(target$tip.label, tips)]
	estimgdi <- gdi
	if (isTRUE(estimgdi)) {
		anc <- target$edge[match(seq(ncol(theta)), target$edge[,2]),1]
		gdi <- matrix(NA, nsample, ncol(theta), dimnames=list(rownames(params), paste0("gdi_",seq(ncol(theta)))))
		for (i in setdiff(seq(ncol(theta)), ntip+1)) {
			gdi[,i] <- 1 - exp(-2 * (tau[,anc[i]] - tau[,i]) / theta[,i])
		}
		mgdi <- round(apply(gdi, 2, fun), 6)
		igdi <- round(apply(gdi, 2, int, p=p), 6)
		gdi <- round(gdi, 6)
	} else {
		mgdi <- NULL
		igdi <- NULL
		gdi <- NULL
	}
	mtau <- round(apply(tau, 2, fun), 6)
	mtheta <- round(apply(theta, 2, fun), 6)
	itau <- round(apply(tau, 2, int, p=p), 6)
	itheta <- round(apply(theta, 2, int, p=p), 6)

	tau <- round(tau, 6)
	theta <- round(theta, 6)

	node_labels <- matrix(paste0("[&tau=", as.matrix(tau[,-seq(ntip)]), "]"), nsample, ncol(tau) - ntip)
	if (isTRUE(estimgdi)) {
		branch_labels <- matrix(paste0("[&theta=", as.matrix(theta), ",", "gdi=", as.matrix(gdi), "]"), nsample, ncol(theta))
		node_labels[,1] <- paste0(gsub("]", "", node_labels[,1]), ",", gsub("\\[&|,gdi=NA", "", branch_labels[,ntip+1]))
	} else {
		branch_labels <- matrix(paste0("[&theta=", as.matrix(theta), "]"), nsample, ncol(theta))
	}
	for (i in seq_along(trees)) {
		trees[[i]]$node.label <- nodes
	}
	trees <- unname(sapply(trees, ape::write.tree))
	branch_patterns <- c(paste0(target$tip.label, ":"), paste0(")", target$node.label, ":"))
	node_patterns <- branch_patterns[-seq(ntip)]
	node_patterns[1] <- sub(":", ";", node_patterns[1])
	for (i in seq_along(trees)) {
		m <- regexpr(node_patterns[1], text=trees[i])
		regmatches(trees[i], m) <- paste0(")", node_labels[i,1], ";")
		for (j in seq_along(nodes)[-1]) {
			m <- regexpr(branch_patterns[j+ntip], text=trees[i])
			regmatches(trees[i], m) <- paste0(regmatches(trees[i], m), branch_labels[i,j+ntip])
			m <- regexpr(node_patterns[j], text=trees[i])
			regmatches(trees[i], m) <- paste0(")", node_labels[i,j], ":")
		}
		for (j in which(!is.na(theta[1,seq(ntip)]))) {
			m <- regexpr(branch_patterns[j], text=trees[i])
			regmatches(trees[i], m) <- paste0(regmatches(trees[i], m), branch_labels[i,j])
		}
	}
	
	target$edge.length <- mtau[target$edge[,1]] - mtau[target$edge[,2]]
	target <- ape::write.tree(target)
	prob <- paste0(substr(p, 3, 4), "%")
	itaulab <- paste0("tau", prob, "_", interval, "={", apply(itau, 2, paste, collapse=","), "}")
	ithetalab <- paste0("theta", prob, "_", interval, "={", apply(itheta, 2, paste, collapse=","), "}")
	if (isTRUE(estimgdi)) {
		igdilab <- paste0("gdi", prob, "_", interval, "={", apply(igdi, 2, paste, collapse=","), "}")
	}
	m <- regexpr(node_patterns[1], text=target)
	regmatches(target, m) <- paste0(")", "[&tau=", mtau[ntip+1], ",", itaulab[ntip+1], ",", "theta=", mtheta[ntip+1], ",", ithetalab[ntip+1], "]", ";")
	for (j in seq_along(nodes)[-1]) {
		m <- regexpr(branch_patterns[j+ntip], text=target)
		if (isTRUE(estimgdi)) {
			regmatches(target, m) <- paste0(regmatches(target, m), "[&theta=", mtheta[ntip+j], ",", ithetalab[ntip+j], ",", "gdi=", mgdi[ntip+j], ",", igdilab[ntip+j], "]")
		} else {
			regmatches(target, m) <- paste0(regmatches(target, m), "[&theta=", mtheta[ntip+j], ",", ithetalab[ntip+j], "]")
		}
		m <- regexpr(node_patterns[j], text=target)
		regmatches(target, m) <- paste0(")", "[&tau=", mtau[ntip+j], ",", itaulab[ntip+j], "]:")
	}
	for (j in which(!is.na(theta[1,seq(ntip)]))) {
		m <- regexpr(branch_patterns[j], text=target)
		if (isTRUE(estimgdi)) {
			regmatches(target, m) <- paste0(regmatches(target, m), "[&theta=", mtheta[j], ",", ithetalab[j], ",", "gdi=", mgdi[j], ",", igdilab[j], "]")
		} else {
			regmatches(target, m) <- paste0(regmatches(target, m), "[&theta=", mtheta[j], ",", ithetalab[j], "]")
		}
	}

	if (isTRUE(figtree)) {
		if (!missing(file)) {
			writeLines(make_figtree(target, tips), con=file) 
		}
		if (!missing(treesfile)) {
			writeLines(make_figtree(trees, tips), con=treesfile) 
		}
	} else {
		if (!missing(file)) {
			writeLines(target, con=file)
		}
		if (!missing(treesfile)) {
			writeLines(trees, con=treesfile)
		}		
	}
}



#' @export
read_msc_tree <- function(file) {
	tree <- readLines(file)
	tree <- gsub("^[[:blank:]]+", "", tree)
	tree <- tree[nchar(tree) > 0]
	lower <- tolower(tree)
	if (grepl("nexus", lower[1])) {
		transl <- tree[(grep("^translate", lower) + 1):(grep("^tree", lower)[1] - 1)]
		transl <- gsub(",", "", transl)
		transl <- do.call(rbind, strsplit(transl, split="[[:space:]]+"))
		tree <- tree[grep("^tree", lower)[1]]
		tree <- regmatches(tree, regexpr("\\(.+;", tree))
		tree <- gsub("^[[:space:]]+", "", tree)
		for (i in seq(nrow(transl))) {
			tree <- sub(paste0(",",transl[i,1],":"), paste0(",",transl[i,2],":"), tree)
			tree <- sub(paste0("\\(",transl[i,1],":"), paste0("(",transl[i,2],":"), tree)
		}		
	}
	brlab <- gregexpr("]:\\[[^]]*][^:]", tree)
	attr(brlab[[1]],"match.length") <- attr(brlab[[1]],"match.length") - 1
	regmatches(tree, brlab) <- list(paste0(gsub("^.*&", ",", unlist(regmatches(tree, brlab))), ":"))
	brlab <- gregexpr(":\\[[^]]*]", tree)
	regmatches(tree, brlab) <- list(paste0(gsub("^:", "", unlist(regmatches(tree, brlab))), ":"))
	m <- gregexpr("\\[[^]]+]", tree)
	node_labels <- setNames(unlist(regmatches(tree, m)), paste0("Node", seq_along(unlist(m))))
	regmatches(tree, m) <- list(names(node_labels))
	tree <- ape::read.tree(text=tree)
	ntip <- length(tree$tip.label)
	node_labels <- unname(node_labels[c(gsub("^.+Node", "Node", tree$tip.label), tree$node.label)])
	m <- gregexpr("\\{[^}]+}", node_labels)
	regmatches(node_labels, m) <- lapply(regmatches(node_labels, m), function(x) gsub(",", "-", x))
	node_labels <- strsplit(gsub("\\[|&|]", "", node_labels), split=",")
	node_labels <- lapply(node_labels, function(x) x[!grepl("\\{", x)])
	isgdi <- any(grepl("gdi",node_labels[[1]]))
	if (isgdi) {
		params <- matrix(, length(node_labels), 3, dimnames=list(NULL, c("tau","theta","gdi")))
	} else {
		params <- matrix(, length(node_labels), 2, dimnames=list(NULL, c("tau","theta")))
	}
	params[seq(ntip),] <- cbind(0, gsub("^[[:alpha:]=]+", "", do.call(rbind, node_labels[seq(ntip)])))
	params[ntip+1,1:2] <- gsub("^[[:alpha:]=]+", "", node_labels[[ntip+1]])
	params[-seq(ntip+1),] <- gsub("^[[:alpha:]=]+", "", do.call(rbind, node_labels[-seq(ntip+1)]))
	mode(params) <- "numeric"
	tree$params <- as.data.frame(params)
	tree$tip.label <- gsub("Node[[:digit:]]+$", "", tree$tip.label)
	tree$node.label <- NULL
	return(tree)	
}



#' @export
rotate_msc_tree <- function(tree, constraint) {
	ntip <- ape::Ntip(tree)
	rotated <- ape::rotateConstr(tree, constraint)
	rewritten <- ape::read.tree(text=ape::write.tree(rotated))
	tips <- rotated$edge[match(seq(ntip), rewritten$edge[,2]),2]
	root <- ntip + 1
	nodes <- rotated$edge[match(ntip + 2:(ntip - 1), rewritten$edge[,2]),2]
	tree$edge <- rewritten$edge
	tree$edge.length <- rewritten$edge.length
	tree$tip.label <- tree$tip.label[tips]
	tree$params <- tree$params[c(tips, root, nodes),]
	return(tree)
}
