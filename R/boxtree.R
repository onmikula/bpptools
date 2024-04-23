#' Box depiction of species tree.
#' 
#' @description
#' Plots species tree with branches represented by rectangles (boxes) whose sides are proportional to multispecies coalescent parameters ('tau' and 'theta'), expressed in comparable units (expected number of mutations, see Details).
#' 
#' @param tree Tree imported by `read_msc_tree`, i.e., `phylo` object with the additional component `params`,
#'   or name of a file, which can be imported by `read_msc_tree`.
#' @param gap The minimum gap between branches, specified as a proportion of the maximum branch width,
#'   but then automatically adjusted.
#' @param col Color of branch background, can be a vector ordered as branches in `tree$edge`.
#' @param border Color of branch outline.
#' @param lwd Width of branch outline.
#' @param cex Size of axis labels, serves also for derivation of size of axis annotation. 
#' @param round How much rounded are branch corners. Specified as a proportion of the minimum of branch
#'    lengths and widths, which is then used as a diameter of the circle defining the curvature.
#' @param show.tip.label Logical, whether to show tip labels
#' @param tip.label Vector of tip labels, used if the desired tip labels differ from `tree$tip.label`,
#'   the vector is named by `tree$tip.label`, however, to establish mutual correspondence of the labels.
#' @param tip.cex Size of tip labels.
#' @param tip.font Font of tip labels.
#' @param tip.offset A multiplayer of space between tree tips and species labels.
#' @param root.lwd A multiplayer of the width of the root line.
#' @param timescale Logical, whether to display time scale axis.
#' @param axis.lwd Width of the axis line.
#' @param labels Character, what is displayed at branch labels, options are `"theta"`, `"tau"` or `"gdi"`.
#' @param digits Numeric, the number of decimal places for parameters to be displayed at node labels.
#' @param lab.col Color of branch label background.
#' @param lab.border Color of branch label outline.
#' @param lab.lwd Width of branch label outline
#' @param lab.cex Size of branch label
#' @param lab.font Font of branch label
#' @param xlab,ylab Labels of x and y axes
#' @param xlim,ylim Limits of x and y axes
#' @param mai Size of outer margins in inches, recycled if necessary
#' @param direction Tree direction (only 'rightwards' is implemented now)
#' @param device Character, type of the graphical device. If not specified (default), either `"quartz"` or `"x11"`
#'   is chosen (depending on the platform). If `NULL`, the objects are plotted into the current device. If it is
#'   `"pdf"`, the pdf file is created, whose name is specified by the `pdf` argument.
#' @param pdf Character, the name of pdf file, if specified, it overrides `device` argument by setting it to `"pdf"`.
#' @param width Width of graphical device.
#' @param height Height of graphical device.
#' @details Both divergence times ('taus') and population size parameters ('thetas') are
#'   in units of the expected number of substitutions per site. However, taus correspond to 
#'   tip-node distances, whereas thetas to tip-node-tip distances, so the thetas are divided
#'   by two to be in scale in the figure.
#' @export

boxtree <- function(tree, gap = 0.1, col = "grey", border="black", lwd=3, cex=1.5, round=0.75, show.tip.label=TRUE, tip.label=NULL, tip.cex=cex, tip.font=3, tip.offset=1.1, root.lwd=1, timescale=TRUE, axis.lwd=lwd/2, labels=NULL, digits=4, lab.col="white", lab.border="black", lab.lwd=lwd/2, lab.cex=1, lab.font=1, xlab="Tau", ylab="Theta", xlim, ylim, mai, direction, device, pdf, width=NULL, height=NULL) {

# tree data
	if (is.character(tree) & length(tree) == 1) {
		tree <- read_msc_tree(tree)
	}
	ntip <- length(tree$tip.label)
	root <- ntip + 1
	nodes <- seq(ntip - 1) + ntip
	params <- tree$params
	params[,"theta"] <- 0.5 * params[,"theta"]

# provisional branch coordinates
	bgap <- gap * max(params[,"theta"], na.rm=TRUE)
	ord <- order(match(seq(ntip), tree$edge[,2]))
	thetas <- ifelse(is.na(params[seq(ntip),"theta"]), 0, params[seq(ntip),"theta"])
	xx <- matrix(max(params[,"tau"]) - params[tree$edge,"tau"], nrow(tree$edge), 2)
	yy <- cumsum(thetas) - 0.5 * thetas
	yy <- yy + cumsum(rep(bgap, ntip))
	yy <- yy[tree$edge[,2]]
	while (any(is.na(yy))) {
		a <- tree$edge[!is.na(yy),1]
		a <- max(intersect(a[which(duplicated(a))], tree$edge[is.na(yy),2]))
		yy[tree$edge[,2] == a] <- mean(yy[tree$edge[,1] == a])
	}
	thetas <- ifelse(is.na(params[tree$edge[,2],"theta"]), 0, params[tree$edge[,2],"theta"])
	yy <- cbind(yy - 0.5 * thetas, yy + 0.5 * thetas)

# root branch coordinates & possibly shift of other branch coordinates
	edg <- lapply(nodes, get_clade_edges, phy=tree)
	edg <- c(as.list(match(seq(ntip), tree$edge[,2])), edg)
	off <- which(tree$edge[,1] == root)
	off <- off[order(rowMeans(yy[off,]))]
	bound <- as.numeric(t(yy[off,]))
	ry <- mean(bound) + c(-1/2, 1/2) * params[root,"theta"]
	rgap <- gap * params[root,"theta"]
	rshift <- rgap + c(bound[2] - ry[1], ry[2] - bound[3])
	if (rshift[1] > 0) {
		branches <- edg[[tree$edge[off[1],2]]]
		yy[branches,] <- yy[branches,] - rshift[1]
	}
	if (rshift[2] > 0) {
		branches <- edg[[tree$edge[off[2],2]]]
		yy[branches,] <- yy[branches,] + rshift[2]
	}

# shifts of species branches
	for (anc in which(tree$edge[,2] %in% nodes)) {
		off <- which(tree$edge[,1] == tree$edge[anc,2])
		off <- off[order(rowMeans(yy[off,]))]
		ancgap <- gap * params[tree$edge[anc,2],"theta"]
		shifts <- c(yy[off[1],1] - yy[anc,1], yy[anc,2] - yy[off[2],2])
		shifts <- ifelse(shifts > 0, shifts + ancgap, shifts)
		branches <- edg[[tree$edge[anc,2]]]
		off1 <- ifelse(tree$edge[off[1],2] > ntip, tree$edge[off[1],2], NA)
		off2 <- ifelse(tree$edge[off[2],2] > ntip, tree$edge[off[2],2], NA)
		if (shifts[1] > 0) {
			branches1 <- setdiff(branches, c(anc, off[2], unlist(edg[off2]))) 
			yy[branches1,] <- yy[branches1,] - shifts[1]
		}
		if (shifts[2] > 0) {
			branches2 <- setdiff(branches, c(anc, off[1], unlist(edg[off1]))) 
			yy[branches2,] <- yy[branches2,] + shifts[2]
		}
		offgap <- gap * max(params[tree$edge[off,2],"theta"])
		shifts <- diff(yy[off,][2:3])
		shifts <- 0.5 * ifelse(shifts > 0 & shifts < offgap, offgap, shifts)
		if (shifts > 0) {
			branches1 <- unique(c(off[1], unlist(edg[off1])))
			yy[branches1,] <- yy[branches1,] - shifts[1]
			branches2 <- unique(c(off[2], unlist(edg[off2])))
			yy[branches2,] <- yy[branches2,] + shifts[2]
		}
	}

# overlaps of both sister & non-sister branches & shifts of the subclades they belong to
	# sectimes = 'cross-section times' - tau - max length of the daughters
	# branchpairs = indices of branches within this timeframe -> their non-ancestral pairs 
	# overlaps = overlaps of branches along y axis
	sectimes <- unique(sapply(nodes, function(i) min(xx[which(tree$edge[,1] == i),2])))
	for (sec in sectimes) {
		branchpairs <- sort(which(xx[,1] <= sec & xx[,2] >= sec))
		branchpairs <- utils::combn(branchpairs, 2, simplify=FALSE)
		branchpairs <- branchpairs[sapply(branchpairs, function(ii) !any(tree$edge[ii,1] %in% tree$edge[ii,2]))]
		overlaps <- sapply(branchpairs, function(ii) diff(yy[ii,][order(rowMeans(yy[ii,])),][2:3]))
		while (any(overlaps > 0)) {
			mrca <- ape::getMRCA(tree, tree$edge[branchpairs[[which.max(overlaps)]],2])
			basal <- edg[tree$edge[tree$edge[,1] == mrca,2] - ntip]
			basal <- basal[order(rowMeans(yy[tree$edge[,1] == mrca,]))]
			shifts <- 0.5 * max(overlaps) + gap * max(params[tree$edge[ii,2],"theta"])
			yy[basal[[1]],] <- yy[basal[[1]],] - shifts
			yy[basal[[2]],] <- yy[basal[[2]],] + shifts
			overlaps <- sapply(branchpairs, function(ii) diff(yy[ii,][order(rowMeans(yy[ii,])),][2:3]))
		}
	}

# tip-to-root check to shift ancestral branches to gaps between their offspring
	sisters <- which(tree$edge[,2] <= ntip)
	sisters <- split(sisters, tree$edge[sisters,1])
	sisters <- sisters[sapply(sisters, length) == 2]
	anc <- match(tree$edge[sapply(sisters, "[", 1),1], tree$edge[,2])
	while (length(sisters) > 0) {
		ygap <- yy[sisters[[1]],]
		yanc <- mean(yy[anc[1],])
		ydif <- mean(ygap[2:3]) - yanc
		if (ydif > 0) {
			yspace <- max(ygap) - yy[anc[1],2]
			yshift <- ifelse(yspace > 0, min(c(ydif, yspace)), 0)
		} else {
			yspace <- min(ygap) - yy[anc[1],1]
			yshift <- ifelse(yspace < 0, max(c(ydif, yspace)), 0)
		}
		yy[anc[1],] <- yy[anc[1],] + yshift
		newsisters <- which(tree$edge[,1] == tree$edge[anc[1],1])
		newanc <- match(tree$edge[newsisters[1],1], tree$edge[,2])
		if (length(newsisters) == 2 & !is.na(newanc)) {
			sisters <- c(sisters, list(newsisters))
			anc <- c(anc, newanc)
		}
		sisters <- sisters[-1]
		anc <- anc[-1]
	}

# preparing plot functions & arguments
	if (isTRUE(show.tip.label)) {
		terminal <- match(seq(ntip), tree$edge[,2])
		if (is.null(tip.label)) {
			tiplab <- gsub("_", " ", tree$tip.label)
		} else {
			tiplab <- gsub("_", " ", tip.label[tree$tip.label])
		}
	}
	if (missing(xlim)) {
		xlim <- range(xx)
	}
	if (missing(ylim)) {
		ylim <- range(yy)
	}
	if (!missing(pdf)) {
		device <- "pdf"
	}
	if (missing(device)) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
	}
	if (isTRUE(device == "pdf") & missing(pdf)) {
		pdf <- "msc_tree.pdf"
	}
	col <- rep(col, nrow(tree$edge))


# device dimensions
	#strh <- strheight(s=0, units="inches", cex=cex.axis)
	#strw <- max(strwidth(s=tiplab, units="inches", cex=tip.cex, font=tip.font))
	strh <- 0.13888890
	strw <- 0.09269206

	nulldim <- c(w=is.null(width), h=is.null(height))
	if (nulldim["w"]) {
		width <- 7
	}
	if (nulldim["h"]) {
		height <- 7
	}

# margins
	missmai <- missing(mai)
	if (missmai) {
		mai <- c(0, 0, 0.02, 0)
		if (isTRUE(timescale)) {
			cex.axis <- ifelse(cex > 1, 1 + (cex-1) / 2, 1)
			mai[1] <- 0.02 * min(c(width, height)) + strh * cex.axis
		}
		if (nchar(xlab) > 0) {
			maix <- mai[1] + 1.5 * strh * cex
#			mai[1] <- mai[1] + 3.5 * strh * cex
			mai[1] <- mai[1] + 3.0 * strh * cex
		}
		if (nchar(ylab) > 0) {
			maiy <- mai[2] + 1.5 * strh * cex
#			mai[2] <- mai[2] + 4.0 * strh * cex
			mai[2] <- mai[2] + 3.5 * strh * cex
		}
		if (isTRUE(show.tip.label)) {
			mai[4] <- strw * tip.cex * max(nchar(tiplab)) + 0.01 * width
		}
	} else {
		maix <- mai[1]
		maiy <- mai[2]
	}

	if (any(nulldim)) {
		oldwh <- c(w=width, h=height)
		innerw <- width - sum(mai[c(2,4)])
		innerh <- height - sum(mai[c(1,3)])
		lrat <- diff(ylim) / diff(xlim)
		irat <- innerh / innerw
		rat <- lrat / irat
		if (nulldim["w"] & nulldim["h"]) {
			if (rat > 1) {
				width <- sum(mai[c(2,4)]) + innerw / rat
			} else {
				height <- sum(mai[c(1,3)]) + rat * innerh
			}
		} else if (nulldim["h"]) {
			height <- sum(mai[c(1,3)]) + rat * innerh
		} else {
			width <- sum(mai[c(2,4)]) + innerw / rat
		}
	}

# graphical device
	if (isTRUE(device == "pdf")) {
		pdf(file, width=width, height=height)
	} else if (!is.null(device)) {
		match.fun(tolower(device))(width=width, height=height)	
	}
	par(mai=mai)

# plotting
	plot(0, xlim=xlim, ylim=ylim, asp=1, bty="n", axes=FALSE, ann=FALSE)
	for (i in nodes) {
		ii <- which(tree$edge[,1] == i)
		lines(xx[ii,1], rowMeans(yy[ii,]), lwd=lwd, col=border)
	}
	d <- round * min(c(diff(t(xx)), diff(t(yy))))
	zero <- which(diff(t(yy)) == 0)
	for (i in setdiff(seq(nrow(tree$edge)), zero)) {
		drawStadium(xx[i,], yy[i,], d, col=col[i], border=border, lwd=lwd)
	}
	for (i in zero) {
		lines(xx[i,], yy[i,], lwd=lwd, col=border, lty=1)
	}
	basal <- sort(yy[tree$edge[,1]==(ntip+1),])[2:3]
	ry <- mean(basal) + c(-1/2, 1/2) * params[ntip + 1,"theta"]
	lines(cbind(0,ry), col=border, lty=1, lwd=root.lwd*2.5*lwd)

# time scale
	if (isTRUE(timescale)) {
		tau <- pretty(xlim)
		tau <- tau[tau <= xlim[2]]
		at <- tau - (max(tau) - xlim[2])
		tau <- rev(tau - min(tau))
		tau <- formatC(tau, format="f", digits=max(nchar(tau)))
		tau <- formatC(as.numeric(tau), format="f", digits=max(nchar(gsub("^.*\\.|0+$", "", tau))))
		axis(side=1, at=at, labels=tau, lwd=axis.lwd, cex.axis=cex.axis)
	}
	if (nchar(xlab) > 0) {
		mtext(text=xlab, side=1, at=mean(xlim), cex=cex, line=maix/0.2)
	}
	if (nchar(ylab) > 0) {
		mtext(text=ylab, side=2, at=mean(ylim), cex=cex, line=maiy/0.2)
	}
	if (isTRUE(show.tip.label)) {
		mtext(text=tiplab, side=4, at=rowMeans(yy[terminal,]), cex=tip.cex, font=tip.font, las=2, line=tip.offset*maix)
	}

# labels
	if (is.null(labels)) labels <- NA
	if (labels %in% colnames(params)) {
		labtext <- params[,labels]
		if (labels == "theta") {
			labtext <- 2 * labtext
		}
		labtext <- formatC(labtext, format="f", digits=digits)
		ord <- match(c(seq(ntip), nodes), tree$edge[,2])
		labxy <- cbind(rowMeans(xx), rowMeans(yy))[ord,]
		labxy[ntip + 1,] <- c(0, mean(ry))
		limits <- t(rbind(diff(t(xx)), diff(t(yy))))[ord,]
		limits[ntip + 1,] <- c(0, diff(ry))		
	}

	if (labels %in% c("tau","theta")) {
		labw <- strw * nchar(labtext) * lab.cex
		labh <- strh * (nchar(labtext) > 0) * lab.cex
		m <- t(c(-1, 1))
		rectxy <- labxy[,c(1,1,2,2)] + cbind(t(t(labw/2)) %*% (0.75 * m), t(t(labh/2)) %*% (1.5 * m))
		srt <- ifelse(diff(t(rectxy[,1:2])) < limits[,1], 0, 90)
		for (i in which(srt == 90)) {
			rectxy[i,] <- rotxy(rectxy[i,], r=pi/2)
		}
		if (labels == "tau") {
			ii <- nodes
		} else if (labels == "theta") {
			ii <- c(seq(ntip), nodes)
		}
		ii <- setdiff(ii, which(is.na(params[,labels])))
		for (i in ii) {
			rect(rectxy[i,1], rectxy[i,3], rectxy[i,2], rectxy[i,4], col=lab.col, border=lab.border, lwd=lab.lwd)
			text(labxy[i,1], labxy[i,2], labels=labtext[i], cex=lab.cex, font=lab.font, srt=srt[i])
		}		
	} else if (labels %in% "gdi") {
		labw <- strw * nchar(labtext) * lab.cex
		labh <- strh * (nchar(labtext) > 0) * lab.cex
		ii <- c(seq(ntip), nodes[-1])
		ii <- setdiff(ii, which(is.na(params[,labels])))
		for (i in ii) {
			drawEllipse(labxy[i,1], labxy[i,2], w=0.75*labw[i], h=1.75*strh[i], r=0, size=1, col=lab.col, border=lab.border, lwd=lab.lwd, res=200)
			text(labxy[i,1], labxy[i,2], labels=labtext[i], cex=lab.cex, font=lab.font)
		}
	}

# graphical device
	if (isTRUE(device == "pdf")) {
		invisible(dev.off())
	}
	
}



drawEllipse <- function(x, y, size, w, h, r, res, col, border, lwd, return=FALSE) {
	rs <- seq(0, 2 * pi, len=res)
	pts <- cbind(0.5 * w * cos(rs), 0.5 * h * sin(rs))
	rot <- matrix(c(cos(r), sin(r), -sin(r), cos(r)), 2, 2)
	pts <- size * pts %*% rot + rep(1, res) %*% t(c(x,y))
	if (isTRUE(return)) {
		return(pts)
	} else {
		polygon(pts, col=col, border=border, lwd=lwd)
	}
}


drawStadium <- function(xx, yy, d, col, border, lwd, return=FALSE) {
	len <- 100
	rmat <- function(r, pts) pts %*% matrix(c(cos(r), sin(r), -sin(r), cos(r)), 2, 2)
	rs <- seq(pi / 2, 0, len=len)
	pts <- d * cbind(0.5 * cos(rs), 0.5 * sin(rs))
	pts <- c(list(pts), lapply(c(pi/2, pi, 3*pi/2), rmat, pts=pts))
	xx <- xx + c(0.5,-0.5) * d
	yy <- yy + c(0.5,-0.5) * d
	pts[[1]] <- pts[[1]] + rep(1,len) %*% t(c(xx[2], yy[2]))
	pts[[2]] <- pts[[2]] + rep(1,len) %*% t(c(xx[2], yy[1]))
	pts[[3]] <- pts[[3]] + rep(1,len) %*% t(c(xx[1], yy[1]))
	pts[[4]] <- pts[[4]] + rep(1,len) %*% t(c(xx[1], yy[2]))
	pts <- do.call(rbind, pts)
	if (isTRUE(return)) {
		return(pts)
	} else {
		polygon(pts, col=col, border=border, lwd=lwd)
	}
}


rotxy <- function(xy, r) {
	mat <- matrix(xy, 2, 2)
	cen <- c(1, 1) %*% t(colMeans(mat))
	rot <- matrix(c(cos(r), sin(r), -sin(r), cos(r)), 2, 2)
	return(as.numeric((mat - cen) %*% rot + cen))
}


get_clade_edges <- function(anc, phy) {
	tips <- match(ape::extract.clade(phy,anc)$tip.label, phy$tip.label)
	last <- phy$edge[,2] %in% tips
	edges <- last
	while(any(last)) {
		edges <- edges | last
		last <- phy$edge[,2] %in% phy$edge[last,1]
	}
	return(which(edges))
}

