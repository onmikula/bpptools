#' Prior and posterior probability densities.
#' 
#' @description
#' Plots prior and/or posterior probability density for chosen parameters of the multispecies coalescent model.
#' 
#' @param prior Numeric vector giving parameters of prior distribution or a list of these vectors.
#'   See Details for notes on parametrization of particular distributions.
#' @param posterior A data frame with posterior samples of multispecies coalescent parameters,
#'   e.g., `params` data frame of `bpp` object created by `read_bpp_mcmc`.
#' @param distribution Character vector giving name(s) of prior distribution(s).
#' @param qmax The maximum quantile of the distributions to display.
#' @param points Logical, whether to display quantiles of the prior as points.
#' @param q Quantiles of the prior to be indicated by points.
#' @param device Character, type of the graphical device. If not specified (default), either `"quartz"` or `"x11"`
#'   is chosen (depending on the platform). If it is `"pdf"`, the pdf file is created, whose name is specified 
#'   by the `pdf` argument. If `NULL`, the objects are plotted into the current device. 
#' @param pdf Character, the name of pdf file, if specified, it overrides `device` argument by setting it to `"pdf"`.
#' @param width Width of graphical device.
#' @param height Height of graphical device.
#' @param xlim Limits of x axis.
#' @param xlab,ylab Labels of x and y axes
#' @param col Line color for prior and posterior (in this order), recycled if necessary.
#'   It can be a list of vectors with different colors for different priors and posteriors, respectively.
#' @param lty Line type for prior and posterior, recycled if necessary, can be a list of vectors.
#' @param lwd Line width for prior and posterior, recycled if necessary, can be a list of vectors.
#' @param cex Point size.
#' @param legend Logical, whether to display legend (symbols of prior and posterior) or character, giving
#'  the labels of displayed distributions (e.g. gamma, invgamma for priors or pamateer names for posterior). 
#' @param expr An expression allowing to add other elements to the plot.
#' @details Currently acceptable prior distributions are `gamma`, `invgamma` and `beta`, which all accept two parameters.
#'   For the gamma distribution, we used shape and rate (1 / scale) parameters, which gives the mean = shape / rate.
#'   For the inverse gamma distribution, shape and scale parameters are specified with the mean = scale / (shape - 1).
#' @export



dinvgamma <- function(x, shape, scale, log=FALSE) {
	d <- scale^shape / gamma(shape) * ((1 / x)^(shape + 1)) * exp(-scale / x)
	if (isTRUE(log)) d <- log(d)
	return(d)
}
qinvgamma <- function(p, shape, scale, lower.tail=TRUE, log.p=FALSE) {
	stats::qgamma(1 - p, shape, 1/scale, lower.tail=lower.tail, log.p=log.p)^(-1)
}



plot_priorposterior <- function(prior, posterior, distribution, qmax=0.99, points=FALSE, q=c(0.025, 0.5, 0.975), device, pdf, width=7, height, xlim=NULL, xlab="", ylab="density", col=c(1,8), lty=c(3,1), lwd=2, cex=1.5, legend=TRUE, expr) {
	if (!missing(posterior)) {
		if (is.numeric(posterior)) {
			posterior <- data.frame(posterior)
		}
		dposterior <- lapply(seq(length(posterior)), function(i) stats::density(posterior[,i], cut=0))
	} else {
		dposterior <- NULL
	}
	if (!missing(prior)) {
		if (is.data.frame(prior)) {
			prior <- apply(prior[,1:2], 1, as.numeric, simplify=FALSE)
		}
		if (!is.list(prior)) {
			prior <- list(prior)
		}
		if (missing(distribution)) {
			stop("name of prior distribution must be specified")
		}
		distribution <- rep_len(distribution, length(prior))
		qdistr <- list(gamma=stats::qgamma, invgamma=qinvgamma, beta=stats::qbeta)[distribution]
		ddistr <- list(gamma=stats::dgamma, invgamma=dinvgamma, beta=stats::dbeta)[distribution]
		from <- lapply(seq_along(prior), function(i) qdistr[[i]](0.00001, prior[[i]][1], prior[[i]][2]))
		to <- lapply(seq_along(prior), function(i) qdistr[[i]](qmax, prior[[i]][1], prior[[i]][2]))
		xprior <- lapply(seq_along(prior), function(i) seq(from[[i]], to[[i]], length=101))
		dprior <- lapply(seq_along(prior), function(i) ddistr[[i]](xprior[[i]], prior[[i]][1], prior[[i]][2]))
		if (isFALSE(points)) {
			qprior <- NULL
		} else {
			qprior <- lapply(seq_along(prior), function(i) qdistr[[i]](q, prior[[i]][1], prior[[i]][2]))
		}
	} else {
		xprior <- NULL
		dprior <- NULL
		qprior <- NULL
	}
	ymax <- max(c(unlist(dprior), unlist(lapply(dposterior, "[[", "y"))))
	if (is.null(xlim)) {
		xlim <- range(c(unlist(xprior), unlist(lapply(dposterior, "[[", "x"))))
	}
	pp <- c(length(dprior) > 0, length(dposterior) > 0)
	if (all(pp)) {
		col <- as.list(rep_len(col, 2))
		col[1] <- list(rep_len(col[[1]], length(dprior)))
		col[2] <- list(rep_len(col[[2]], length(dposterior)))	
		lty <- as.list(rep_len(lty, 2))
		lty[1] <- list(rep_len(lty[[1]], length(dprior)))
		lty[2] <- list(rep_len(lty[[2]], length(dposterior)))	
		lwd <- as.list(rep_len(lwd, 2))
		lwd[1] <- list(rep_len(lwd[[1]], length(dprior)))
		lwd[2] <- list(rep_len(lwd[[2]], length(dposterior)))	
	} else if (isTRUE(pp[1])) {
		col <- list(rep_len(col, length(dprior)), NULL)
		lty <- list(rep_len(lty, length(dprior)), NULL)
		lwd <- list(rep_len(lwd, length(dprior)), NULL)
	} else if (isTRUE(pp[2])) {
		col <- list(NULL, rep_len(col, length(dposterior)))		
		lty <- list(NULL, rep_len(lty, length(dposterior)))		
		lwd <- list(NULL, rep_len(lwd, length(dposterior)))		
	} else {
		col <- vector("list", 2)
	}

	if (missing(device)) {
		device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
		devname <- NULL
	}
	ispdf <- ifelse(missing(pdf), FALSE, !is.na(pdf)) | isTRUE(device == "pdf")
	if (isTRUE(ispdf)) {
		device <- "pdf"
		devname <- ifelse(missing(pdf), "mcmc.pdf", pdf)
	}
	if (!is.null(device)) {
		device <- match.fun(tolower(device))
		device(devname, width=width, height=height)	
	}

	par(mai=c(1.02, 1.02, 0.42, 0.42))
	plot(0, 0, xlim=xlim, ylim=c(0, ymax), type="n", xlab=xlab, ylab=ylab, cex.lab=1.5, cex.axis=1.25)
	if (!is.null(dprior)) {
		for (i in seq_along(dprior)) {
			lines(xprior[[i]], dprior[[i]], col=col[[1]][i], lty=lty[[1]][i], lwd=lwd[[1]][i])
		}
	}
	if (!is.null(dposterior)) {
		for (i in seq_along(dposterior)) {
			lines(dposterior[[i]]$x, dposterior[[i]]$y, col=col[[2]][i], lty=lty[[2]][i], lwd=lwd[[2]][i])		
		}
	}
	if (!is.null(dprior) & !is.null(qprior)) {
		for (i in seq_along(qprior)) {
			idprior <- ddistr[[i]](qprior[[i]], prior[[i]][1], prior[[i]][2])
			points(qprior[[i]], idprior, pch=21, cex=cex, col=col[[1]][i], lwd=lwd[1], bg="white")
		}
	}

	if (isTRUE(legend) & all(pp)) {
		leg <- c("Prior","Posterior")
	} else if (is.character(legend)) {
		leg <- legend
	} else {
		leg <- NULL
	}
	if (!is.null(leg)) {
		if (all(pp)) {
			col <- sapply(col, "[[", 1)
			lty <- sapply(lty, "[[", 1)
			lwd <- sapply(lwd, "[[", 1)
		} else if (any(pp)) {
			col <- unlist(col[!sapply(col, is.null)])
			lty <- unlist(lty[!sapply(lty, is.null)])
			lwd <- unlist(lwd[!sapply(lwd, is.null)])
		} else {
			leg <- NULL
		}
	}
	if (!is.null(leg)) {
		legend("topright", legend=leg, col=col, lwd=lwd, lty=lty, bty="n")		
	}

	if (!missing(expr)) {
		if (is.expression(expr)) {
			eval(expr)
		}
	}

	if (ispdf) {
		invisible(dev.off())
	}

}

