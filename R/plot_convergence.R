#' Convergence of two MCMC runs.
#' 
#' @description
#' Graphical depiction of convergence between posterior samples obtained from two independent Markov CHains Monte Carlo.
#' 
#' @param mcmc Either names of two .mcmc files or a list with objects of class `bpp` created by `read_bpp_mcmc`.
#' @param param Name of parameter to be examined (a name of column in `mcmc$params`).
#' @param estimate What to estimate from posterior sample (default is `mean`).
#' @param interval Which kind of interval to use, either `"HPD"` (default) or `"CPD"`
#'   meaning the highest or the central posterior density.
#' @param p Numeric, proportion of posterior density encompassed by the interval.
#' @param device Character, type of the graphical device. If not specified (default), either `"quartz"` or `"x11"`
#'   is chosen (depending on the platform). If it is `"pdf"`, the pdf file is created, whose name is specified 
#'   by the `pdf` argument. If `NULL`, the objects are plotted into the current device. 
#' @param pdf Character, the name of pdf file, if specified, it overrides `device` argument by setting it to `"pdf"`.
#' @param size Size of the graphical device (width and height in inches).
#' @param pch Point symbols.
#' @param cex Point size.
#' @param pt.col Point color.
#' @param lin.col Line color.
#' @param pt.lwd Width of point boundary.
#' @param lin.lwd Line width.
#' @param mai Size of outer margins in inches, recycled if necessary
#' @param title Logical, if `TRUE` (default) it uses either `param` or `main` argument (if the latter is specified).
#' @details Both divergence times ('taus') and population size parameters ('thetas') are
#'   in units of the expected number of substitutions per site. However, taus correspond to 
#'   tip-node distances, whereas thetas to tip-node-tip distances, so the thetas are divided
#'   by two to be in scale in the figure.
#' @export

plot_convergence <- function(mcmc, param, estimate="mean", interval="HPD", p=0.95, device, pdf, size=7, pch=21, cex=3, pt.col=4, lin.col=8, pt.lwd=2, lin.lwd=5, mai, title=TRUE, ...) {
	colstat <- function(x, fun) apply(x, 2, fun)
	estimate <- match.fun(estimate)
	interval <- match.fun(interval)
	if (is.character(mcmc)) {
		mcmc <- lapply(mcmc, read_bpp_mcmc) 
	}
	if (inherits(mcmc[[1]], "bpp")) {
		mcmc <- lapply(mcmc, "[[", "params")
	}
	mcmc <- lapply(mcmc, function(x) x[,grep(param, colnames(x), ignore.case=TRUE),drop=FALSE])
	est <- as.matrix(do.call(cbind, lapply(mcmc, colstat, fun=estimate)))
	int <- lapply(lapply(mcmc, colstat, fun=interval), t)
	if (is.null(colnames(est))) {
		colnames(est) <- names(int) <- paste("Run", seq(ncol(est)))
	}
	ellipsis <- list(...)
	if ("main" %in% names(ellipsis)) {
		title <- TRUE
	} else if (isTRUE(title)) {
		main <- paste0(toupper(substr(param, 1, 1)), tolower(substr(param, 2, nchar(param))))
	} else {
		main <- NULL
	}
	if (length(main) > 0 & !"main" %in% names(ellipsis)) {
		cex.main <- 1.5
	} else if (length(main) == 0) {
		cex.main <- 0
	}
	cex.lab <- ifelse(!"cex.lab" %in% names(ellipsis), 1.5, ellipsis[["cex.lab"]])
	cex.axis <- ifelse(!"cex.axis" %in% names(ellipsis), 1.25, ellipsis[["cex.axis"]])
	if (missing(mai) & isTRUE(title)) {
		mai <- c(0.92, 0.92, 0.82, 0.42)	
	} else if (missing(mai)) {
		mai <- c(0.92, 0.92, 0.42, 0.42)
	}
	pt.bg <- ifelse(pch >= 21, pt.col, NA)
	pt.col <- ifelse(pch >= 21, 1, pt.col)
	size <- rep_len(size, 2)
	if (ncol(est) == 2) {
		if (missing(device)) {
			device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
			devname <- NULL
		}
		ispdf <- !missing(pdf) | isTRUE(device == "pdf")
		if (isTRUE(ispdf)) {
			device <- "pdf"
			devname <- ifelse(missing(pdf), "mcmc.pdf", pdf)
		}
		if (!is.null(device)) {
			device <- match.fun(tolower(device))
			device(devname, width=size[1], height=size[2])	
		}
		xlim <- ylim <- range(unlist(int))
		par(mai=mai)
		plot(est, type="n", xlim=xlim, ylim=ylim, asp=1, cex.lab=cex.lab, cex.axis=cex.axis, main=main, cex.main=cex.main, ...)
		abline(0, 1, lwd=lin.lwd, col=lin.col)
		for (i in seq(nrow(int[[1]]))) {
			lines(cbind(int[[1]][i,], est[i,2]), col=1, lwd=pt.lwd)
			lines(cbind(est[i,1], int[[2]][i,]), col=1, lwd=pt.lwd)
		}
		points(est, pch=pch, col=pt.col, bg=pt.bg, lwd=pt.lwd, cex=cex)
	} else {
		warning("For more than two runs is not yet implemented")
	}
	if (ispdf) {
		invisible(dev.off())
	}
}
