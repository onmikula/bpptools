#' Reading in BPP output.
#' 
#' @description
#' Reads BPP output report into R.
#' 
#' @param outfile A character string giving the name of .out file.
#' @returns The list with components depending on the mode of the analysis.
#'   It always includes components `"model"` ("MSC", "MSC-M", "MSC-I") and `"mode"` (A00 to A11).
#'   Among others, it reports the maximum aposteriori (MAP) species tree (for A01 and A11),
#'   majority consensus species tree (for A01), MAP species delimitation (for A10 and A11)
#'   and majority consensus species delimitation (A11). 
#' @export

read_bpp_output <- function(outfile) {

	lin <- readLines(outfile)
	lin <- sub("^[[:blank:]]+", "", lin)
	lin <- lin[nchar(lin) > 0]
	fin <- grep("spent in MCMC", lin)

	spdelim <- any(grepl("delim", lin))
	sptree <- any(grepl("tree", gsub("guide tree", "", lin[-seq(fin)], ignore.case=TRUE), ignore.case=TRUE))

	if (spdelim & sptree) {
		lin <- lin[-seq(fin)]
		brackets <- grep("^\\([[:alpha:]]\\)", lin)
		modelprob <- lin[(brackets[1]+1):(brackets[2]-1)]
		delimstr <- regexpr("\\(", modelprob)
		attr(delimstr,"match.length") <- regexpr("\\)", modelprob) - delimstr + 1
		delims <- regmatches(modelprob, delimstr)
		modelprob <- regmatches(modelprob, delimstr, invert=TRUE)
		newick <- gsub("[[:blank:]]+", "", sapply(modelprob, "[", 2))
		modelprob <- sub("[[:blank:]]+$", "", sapply(modelprob, "[", 1))
		modelprob <- gsub("[[:blank:]]+", " ", modelprob)
		modelprob <- setNames(as.data.frame(do.call(rbind, strsplit(modelprob, " "))), c("count", "pp", "cumpp", "nspecies"))
		modelprob$count <- as.integer(modelprob$count)
		modelprob$pp <- as.numeric(modelprob$pp)
		modelprob$cumpp <- as.numeric(modelprob$cumpp)
		modelprob$nspecies <- as.integer(modelprob$nspecies)
		modelprob$spdelim <- delims
		modelprob$sptree <- newick
		maptree <- sort(tapply(modelprob$pp, modelprob$sptree, sum), decreasing=TRUE)[1]
		maptree <- ape::read.tree(text=names(maptree))
		delimprob <- lin[(brackets[2]+1):(brackets[3]-1)]
		delimstr <- regexpr("\\(.*\\)", delimprob)
		delims <- regmatches(delimprob, delimstr)
		delimprob <- sub("[[:blank:]]+$", "", substr(delimprob, 1, delimstr - 1))
		delimprob <- gsub("[[:blank:]]+", " ", delimprob)
		delimprob <- setNames(as.data.frame(do.call(rbind, strsplit(delimprob, " "))), c("count", "pp", "nspecies"))
		delimprob$count <- as.integer(delimprob$count)
		delimprob$pp <- as.numeric(delimprob$pp)
		delimprob$nspecies <- as.integer(delimprob$nspecies)
		delimprob$spdelim <- delims
		speciesprob <- lin[(brackets[3]+1):(brackets[4]-1)]
		speciesprob <- gsub("[[:blank:]]+", " ", speciesprob)
		speciesprob <- setNames(as.data.frame(do.call(rbind, strsplit(speciesprob, " "))), c("count", "pp", "species"))
		speciesprob$count <- as.integer(speciesprob$count)
		speciesprob$pp <- as.numeric(speciesprob$pp)
		nspeciesprob <- lin[(brackets[4]+1):length(lin)]
		nspeciesprob <- gsub(" =\\s+", "=", nspeciesprob)
		nspeciesprob <- gsub("[[:blank:]]+", " ", nspeciesprob)
		nspeciesprob <- setNames(as.data.frame(do.call(rbind, strsplit(nspeciesprob, " "))), c("posterior", "prior"))
		nspeciesprob <- data.frame(nspecies=1:nrow(nspeciesprob), prior=as.numeric(sub("^.*=","", nspeciesprob$prior)), posterior=as.numeric(sub("^.*=","", nspeciesprob$posterior)))
		result <- list(maptree=maptree, modelprob=modelprob, delimprob=delimprob, speciesprob=speciesprob, nspeciesprob=nspeciesprob, model="MSC", mode="A11")
	} else if (spdelim) {
		lin <- lin[-seq(fin)]
		brackets <- c("Number of species-delimitation models", "Order of ancestral nodes", "Guide tree")
		brackets <- sapply(brackets, grep, x=lin, ignore.case=TRUE)
		tree <- lin[grep("Guide tree", lin, ignore.case=TRUE) + 1]
		tree <- sub(";;", ";", tree)
		tree <- ape::read.tree(text=tree)
		tree$node.label <- sub("#", "", tree$node.label)
		ntip <- length(tree$tip.label)
		ancord <- lin[(brackets[2] + 1):(brackets[3] - 1)]
		labels <- data.frame(Node=seq_along(ancord)+ntip, Label=ancord)
		delimprob <- lin[(brackets[1] + 1):(brackets[2] - 1)]
		delimprob <- gsub("[[:blank:]]+", " ", delimprob)
		delimprob <- setNames(as.data.frame(do.call(rbind, strsplit(delimprob[-1], " "))[,-1]), c("model", "prior", "posterior"))
		result <- list(guidetree=tree, delimprob=delimprob, labels=labels, model="MSC", mode="A10")
	} else if (sptree) {
		lin <- lin[-seq(fin)]
		brackets <- grep("^\\([[:alpha:]]\\)", lin)
		spord <- lin[(grep("Species in order", lin, ignore.case=TRUE) + 1):(brackets[1] - 1)]
		spord <- sub("^[[:digit:]]+\\.*[[:blank:]]*", "", spord)
		labels <- data.frame(Node=seq_along(spord), Label=spord)
		treeprob <- lin[(brackets[1]+1):(brackets[2]-1)]
		treestr <- regexpr("\\(.*\\);", treeprob)
		newick <- gsub("[[:blank:]]+", "", regmatches(treeprob, treestr))
		treeprob <- sub("[[:blank:]]+$", "", substr(treeprob, 1, treestr-1))
		treeprob <- gsub("[[:blank:]]+", " ", treeprob)
		treeprob <- setNames(as.data.frame(do.call(rbind, strsplit(treeprob, " "))), c("count", "pp", "cumpp"))
		treeprob$count <- as.integer(treeprob$count)
		treeprob$pp <- as.numeric(treeprob$pp)
		treeprob$cumpp <- as.numeric(treeprob$cumpp)
		treeprob$sptree <- newick
		splitprob <- lin[(brackets[2]+1):(brackets[3]-1)]
		splitprob <- gsub("[[:blank:]]+", " ", splitprob)
		splitprob <- setNames(as.data.frame(do.call(rbind, strsplit(splitprob, " "))), c("count", "pp", "split"))
		splitprob$count <- as.integer(splitprob$count)
		splitprob$pp <- as.numeric(splitprob$pp)
		contree <- gsub("[[:blank:]]", "", lin[brackets[3]+1])
		contree <- ape::read.tree(text=contree)
		contree$node.label <- sub("#", "", contree$node.label)
		maptree <- gsub("[[:blank:]]", "", lin[(brackets[4]+1):length(lin)])
		pp <- as.numeric(gsub("\\[P=|]$", "", regmatches(maptree, regexpr("\\[.*]", maptree))))
		maptree <- gsub("\\[.+$", "", maptree)
		maptree <- ape::read.tree(text=maptree)
		if (inherits(maptree, "multiPhylo")) {
			target <- lapply(maptree, function(x) {x$node.label <- sub("#", "", x$node.label); return(x)})
			maptree <- NULL 			
		} else {
			target <- NULL
			maptree$node.label <- sub("#", "", maptree$node.label)
		}
		result <- list(maptree=maptree, contree=contree, treeprob=treeprob, splitprob=splitprob, labels=labels, tips=spord, target=target, model="MSC", mode="A01")
	} else {
		tree <- grep("Initial species tree", lin[seq(fin)], value=TRUE, ignore.case=TRUE)
		tree <- gsub("[[:blank:]]", "", regmatches(tree, regexpr("\\(.+\\);", tree)))
		if (length(tree) == 0) tree <- NULL else tree <- ape::read.tree(text=tree)
		lin <- lin[-seq(fin)]
		params <- lin[(grep("list", lin, ignore.case=TRUE) + 1):length(lin)]
		params <- gsub("\\(.*\\)", "", params)
		squarebr <- regexpr("\\[.+]$", params[-1])
		species <- gsub("\\[\\s*|\\s*]", "", regmatches(params[-1], squarebr))
		params <- sub("[[:blank:]]+$", "", gsub("[[:blank:]]+", " ", substr(params, 1, squarebr-1)))
		params <- do.call(rbind, strsplit(params, " "))
		params <- setNames(data.frame(as.numeric(params[-1,1])+1, as.numeric(params[-1,2]), as.numeric(params[-1,3]), params[-1,4]), params[1,])
		stats <- lin[1:(grep("list", lin, ignore.case=TRUE) - 1)]
		stats <- gsub("[[:blank:]]+", " ", stats)
		stats <- strsplit(stats, " ")
		stats <- setNames(data.frame(do.call(rbind, stats[-1]), row.names=1), stats[[1]])
		if (!is.null(tree)) {
			ntip <- length(tree$tip.label)
			taus <- params$Tau[-seq(ntip)]
			heights <- ifelse(tree$edge <= ntip, 0, tree$edge - ntip)
			heights[heights != 0] <- taus[heights[heights != 0]]
			tree$edge.length <- heights[,1] - heights[,2]
			reord <- match(tree$tip.label, species[seq(ntip)])
			params[seq(ntip),] <- params[reord,]
			species[seq(ntip)] <- species[reord]
			tree$params <- params[,c("Tau","Theta")]
		}
		result <- list(tree=tree, paramstat=stats, labels=params[,c("Node","Label")], species=species, model="MSC", mode="A00")
	}

	return(result)
}

