#' Species tree.
#' 
#' @description
#' Finds species tree representative for a posterior sample.
#'
#' `find_sptree` a wrapper function joining the following functions into a single pipeline.
#'
#' `treePP` takes a posterior sample of trees and converts it into a list with components
#'   `"trees"` (data frame with species trees, their posterior probabilities and clade credibilities) 
#'   and `"clades"` (data frame with clades and their posterior probabilities).
#'
#' `treeMAP`, `treeMC`, `treeCON` functions choosing the maximum a posteriori
#'   the maximum credibility or the majority consensus species tree, respectively. They all return a `phylo` object 
#'   with tree whose node labels indicate PPs of the clades.
#'
#' @param trees A posterior sample of rooted trees in a `multiPhylo` object or a list of `phylo` objects
#'   or an object of class `bpp` or a name of file with the posterior sample.
#' @param type The type of summary tree, either `"MAP"` (maximum aposteriori, default),
#'   `"MCC"` (maximum clade credibility) or `"CON"` (majority consensus), see Details for explanantion.
#' @param file A name of the output file.
#' @param thinning A proportion of mcmc steps to be retained or a subsampling frequency (if `thin > 1`).
#' @param digits The number of digits in PP values.
#' @details The maximum a posteriori tree is the most frequently sampled tree.
#'   The maximum clade credibility tree scores every sampled tree by the product of PPs of its constituent clades.
#'   The majority consensus tree adds the most supported clades unless they overlap (they are in conflict)
#'   with any of the clades that were already included. 
#' @returns A `phylo` object with tree whose node labels indicate posterior probabilities of the clades.
#' @export 
                                                                                                                                                                                                                                                                                                                          
find_sptree <- function(trees, type=c("MAP","MCC","CON"), file, thinning=1, digits=2) {
	if (is.character(trees) & length(trees) == 1) {
		trees <- import_tree(trees, thinning=thinning)
	} else if (inherits(trees, "bpp")) {
		trees <- trees$trees
	}
	tree_pp <- treePP(trees)
	tree_fun <- match.fun(paste0("tree", toupper(type[1])))
	tree_sp <- tree_fun(tree_pp, digits=digits)
	if (!missing(file)) {
		ape::write.tree(tree_sp, file)
	}
	return(tree_sp)
}


#' @export 
treePP <- function(trees) {
	nametips <- function(part) lapply(lapply(part, function(x) attr(part, "labels")[x]), sort)

	trees <- lapply(trees, ape::rotateConstr, constraint=trees[[1]]$tip.label)
	topol <- get_topology(trees)
	treepp <- sort(table(topol), decreasing=TRUE) / length(trees)
	treedf <- data.frame(tree=names(treepp), pp=as.numeric(treepp))

	clades <- ape::prop.part(trees)
	cladepp <- attr(clades, "number") / length(trees)
	clades <- nametips(clades)
	if (length(treedf$tree) == 1) {
		credib <- list(nametips(ape::prop.part(import_tree(treedf$tree))))
	} else {
		credib <- lapply(lapply(import_tree(treedf$tree), ape::prop.part), nametips)
	}
	credib <- lapply(credib, match, table=clades)
	treedf$credib <- sapply(seq_along(credib), function(i) sum(log(cladepp[credib[[i]]])))
	clades <- sapply(clades, function(x) paste0("(", paste(x, collapse=","), ")"))
	cladedf <- data.frame(clade=clades, pp=cladepp)
	cladedf <- cladedf[order(cladedf$pp, decreasing=TRUE),]

	result <- list(trees=treedf, clades=cladedf, tips=labels)
	return(result)
}



#' @export 
treeMAP <- function(pp, digits=2) {
	map <- ape::read.tree(text=pp$trees$tree[which.max(pp$trees$pp)])
	mapclades <- ape::prop.part(map)
	mapclades <- lapply(mapclades, function(x, labels) labels[x], labels=attr(mapclades, "labels"))
	mapclades <- lapply(mapclades, sort)
	mapclades <- sapply(mapclades, function(x) paste0("(", paste(x, collapse=","), ")"))
	map$node.label <- formatC(pp$clades$pp[match(mapclades, pp$clades$clade)], format="f", digits=digits)
	return(map)
}


#' @export 
treeMCC <- function(pp, digits=2) {
	mcc <- ape::read.tree(text=pp$trees$tree[which.max(pp$trees$cred)])
	mccclades <- ape::prop.part(mcc)
	mccclades <- lapply(mccclades, function(x, labels) labels[x], labels=attr(mccclades, "labels"))
	mccclades <- lapply(mccclades, sort)
	mccclades <- sapply(mccclades, function(x) paste0("(", paste(x, collapse=","), ")"))
	mcc$node.label <- formatC(pp$clades$pp[match(mccclades, pp$clades$clade)], format="f", digits=digits)
	return(mcc)
}


#' @export 
treeCON <- function(pp, digits=2) {
	clades <- strsplit(gsub("\\(|\\)", "", pp$clades$clade), ",")
	i <- 1
	con <- clades[i]
	while (i < length(clades)) {
		i <- i + 1
		included <- sapply(lapply(con, function(old, new) new %in% old, new=clades[[i]]), all)
		includes <- sapply(lapply(con, function(old, new) old %in% new, new=clades[[i]]), all)
		if (all(included | includes)) {
			con <- append(con, clades[i])
		}
	}
	con <- con[order(sapply(con, length))]
	node.label <- formatC(pp$clades$pp[match(con, clades)], format="f", digits=digits)
	while (length(con) > 1) {
		nwk <- paste0("(", con[[1]][1], ",", con[[1]][2], ")", node.label[1])
		for (i in 2:length(con)) {
			ii <- match(con[[1]][1], con[[i]])
			if (!is.na(ii)) {
				con[[i]][ii] <- nwk
				con[[i]] <- setdiff(con[[i]], con[[1]])
			}
		}
		node.label <- node.label[-1]
		con <- con[-1]
	}
	nwk <- paste0("(", con[[1]][1], ",", con[[1]][2], ")", node.label[1], ";")
	con <- ape::read.tree(text=nwk)
	return(con)	
}

