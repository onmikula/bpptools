#' Correlations of parameter samples.
#' 
#' @description
#' Plot of pairwise correlations between model parameters as they were sampled by a Markov Chain Monte Carlo.
#' 
#' @param mcmc Either names of two .mcmc files or a list with objects of class `bpp` created by `read_bpp_mcmc`.
#' @param param Name of parameters to be examined (columns in `mcmc$params`), can be specified as a regular expression.
#' @param type Either `"scatter"` (default) or `"ellipse"` (not implemented yet). 
#' @param device Character, type of the graphical device. If not specified (default), either `"quartz"` or `"x11"`
#'   is chosen (depending on the platform). It can be also a name of function creating graphical device
#'   for vector graphic or bitmap file. If `NULL` or `NA`, the objects are plotted into the current device (if exists). 
#' @param file Character, the name of file (for graphical devices like `pdf`, `svd` or `png`).
#' @param size Size of the graphical device (width and height in inches).
#' @param asp Aspect ratio of pairwise scatters.
#' @param thinning Thinning interval.
#' @param mincor The minimum correlation to be highlighted by a colored box.
#' @param labels Axis labels; pretty names of parameters used for display.
#' @param cex.lab Size of axis labels.
#' @param cex.axis Size of tick labels.
#' @param pt.col Point color.
#' @param lin.col Line color.
#' @param lin.lwd Line width.
#' @param box.col Box color.
#' @param box.lwd Box width.
#' @export

plot_correlations <- function(mcmc, param, type=c("scatter", "ellipse"), device, file, size=7, asp=NA, thinning=1, mincor=0.5, labels=NULL, cex.lab=1.5, cex.axis=1, pt.col=4, lin.col=2, lin.lwd=2, box.col="red", box.lwd=3, ...) {
	if (inherits(mcmc, "bpp")) {
		mcmc <- mcmc$params
	}
	if (missing(param)) {
		param <- c("tau", "theta", "phi", "^M")
	}
	if (any(!param %in% colnames(mcmc))) {
		pattern <- paste(setdiff(param, colnames(mcmc)), collapse="|")
		param <- c(intersect(colnames(mcmc), param), grep(pattern, colnames(mcmc), ignore.case=TRUE, value=TRUE))
	}
	if (length(param) == 0) {
		stop("specified parameters not included or the specification was invalid")
	} else if (length(param) == 1) {
		stop("just a single parameter specified")
	}
	mcmc <- mcmc[,param]
	np <- ncol(mcmc)
	if (is.null(labels)) {
		labels <- param
	}
	
	corcoef <- matrix(NA, np, np, dimnames=list(labels, labels))
	for (i in 1:(np-1)) {
		for (j in (i+1):np) {
			corcoef[j,i] <- stats::cor(mcmc[,i], mcmc[,j], method="pearson")
		}
	}
	means <- apply(mcmc, 2, mean)
	sds <- apply(mcmc, 2, sd)

	size <- rep_len(size, 2)
	thinning <- floor(ifelse(thinning < 1, 1 / thinning, thinning))
	thin <- as.integer(seq(0, nrow(mcmc), by=thinning))[-1]
	
	if (missing(device)) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
		devname <- NULL		
	}
	currdevice <- is.null(device) | isTRUE(is.na(device))
	nulldevice <- .Device == "null device"
	if (currdevice & nulldevice) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
		devname <- NULL		
	} else if (currdevice) {
		isfile <- FALSE
		isbitmap <- FALSE
	}
	if (!is.null(device) & !isTRUE(is.na(device))) {
		if (!device %in% c("quartz", "x11")) {
			isfile <- TRUE
			if (missing(file)) {
				file <- paste0("convergence.", sub("^.*_", "", device))
			}
			devname <- file
			isbitmap <- device %in% c("bmp","jpeg","png","tiff")
		} else {
			isfile <- FALSE
			isbitmap <- FALSE
		}
		device <- match.fun(tolower(device))
		if (isTRUE(isbitmap)) {
			device(devname, width=size[1], height=size[2], units="in")	
		} else {
			device(devname, width=size[1], height=size[2])	
		}
	}
	ww <- c(0.2, rep(1, np - 1))
	hh <- c(0.2, rep(1, np - 1))
	din <- par("din")
	win <- din * c(1 / sum(ww), 1 / sum(hh))
	if (win[1] < 2.1875) {
		size <- size * 2.1875 / win[1]
		invisible(dev.off())
		if (currdevice) {
			device <- match.fun(ifelse(.Platform$OS.type == "unix", "quartz", "x11"))
			devname <- NULL		
		}
		if (isTRUE(isbitmap)) {
			device(devname, width=size[1], height=size[2], units="in")	
		} else {
			device(devname, width=size[1], height=size[2])	
		}
		din <- par("din")
		win <- din * c(1 / sum(ww), 1 / sum(hh))
	}
	ww <- din[1] * ww / sum(ww)
	hh <- din[2] * hh / sum(hh)
	mai <- win[c(2,1,2,1)] * c(0.82, 0.82, 0.42, 0.42) / 7
	topmargin <- mai[1] - mai[3]
	hh <- hh + c(0, ifelse(np > 2, -1, 0), rep(0, ifelse(np > 2, np - 3, 0)), rep(1, ifelse(np > 2, 1, 0))) * topmargin
	laymat <- matrix(0, np, np)
	laymat[lower.tri(laymat)] <- seq(choose(np, 2)) + (2 * np - 2)
	laymat <- rbind(c(0, np:(2 * np - 2)), cbind(seq(np - 1), laymat[-1, -np]))
	graphics::layout(laymat, widths=ww, heights=hh)
	par(mai=c(0,0,0,0))
	for (i in 2:np) {
		plot(1, 1, axes=FALSE, ann=FALSE, bty="n", type="n")
		text(1, 1, labels=labels[i], cex=cex.lab, srt=90)
	}
	for (i in 1:(np-1)) {
		plot(1, 1, axes=FALSE, ann=FALSE, bty="n", type="n")
		text(1, 1, labels=labels[i], cex=cex.lab)
	}
	for (j in 1:(np-1)) {
		for (i in (j+1):np) {
			if (np > 2) {
				if (i == 2) par(mai=mai+c(0,0,-topmargin,0))
				if (i > 2 & i < np) par(mai=mai)
				if (i == np) par(mai=mai+c(topmargin,0,0,0))
			} else {
				par(mai=mai)
			}
			b <- corcoef[i,j] * sds[j] / sds[i]
			a <- means[j] - b * means[i]
			plot(mcmc[thin,c(i,j)], xlab="", ylab="", asp=asp, cex.axis=cex.axis, col=pt.col, ...)
			abline(a, b, col=lin.col, lwd=lin.lwd, ...)
			if (abs(corcoef[i,j]) >= mincor) {
				box(which="plot", lty="solid", col=box.col, lwd=box.lwd, ...)
			}
		}
	}
	
	if (isfile) {
		invisible(dev.off())
	}

}
