#' Traces of independent MCMC runs.
#' 
#' @description
#' Plotting of independent MCMC traces (series of sampled paramter values) to check their convergence.
#' 
#' @param mcmc Either names of two .mcmc files or a list with objects of class `bpp` created by `read_bpp_mcmc`.
#' @param par Name of parameter to be examined (a name of column in `mcmc$params`) 
#'   or a regular expression defining more of them. 
#' @param device Character, type of the graphical device. If not specified (default), either `"quartz"` or `"x11"`
#'   is chosen (depending on the platform). It can be also a name of function creating graphical device
#'   for vector graphic or bitmap file. If `NULL` or `NA`, the objects are plotted into the current device (if exists). 
#' @param file Character, the name of file (for graphical devices like `pdf`, `svd` or `png`).
#' @param width Width of graphical device (in inches).
#' @param height Height of graphical device (in inches).
#' @param palette Colors used to distinguish MCMC traces.
#' @param mai Size of outer margins in inches, recycled if necessary
#' @param title Logical, if `TRUE` (default) it uses either `par` or `main` argument (if the latter is specified).
#' @param lwd Numeric, line width.
#' @param ... Graphical parameters, accepted by `plot.default`.
#' @export

plot_traces <- function(mcmc, par, device, file, width=7, height, palette, mai, title=TRUE, lwd=1, ...) {

	if (is.character(mcmc)) {
		mcmc <- lapply(mcmc, read_bpp_mcmc) 
	}
	if (inherits(mcmc[[1]], "bpp")) {
		mcmc <- lapply(mcmc, "[[", "params")
	}
	mcmc <- lapply(mcmc, function(x) x[,grep(par, colnames(x), ignore.case=TRUE),drop=FALSE])
	pooled <- do.call(rbind, mcmc)
	nrun <- length(mcmc)
	ngen <- nrow(mcmc[[1]])
	npar <- ncol(mcmc[[1]])
	
	ellipsis <- list(...)
	if ("main" %in% names(ellipsis)) {
		title <- TRUE
	} else if (isTRUE(title)) {
		main <- paste0(toupper(substr(par, 1, 1)), tolower(substr(par, 2, nchar(par))))
	} else {
		main <- NULL
	}
	if (length(main) > 0 & !"cex.main" %in% names(ellipsis)) {
		cex.main <- 2
	} else if (length(main) == 0) {
		cex.main <- 0
	}
	cex.lab <- ifelse(!"cex.lab" %in% names(ellipsis), 1.5, ellipsis[["cex.lab"]])
	cex.axis <- ifelse(!"cex.axis" %in% names(ellipsis), 1.5, ellipsis[["cex.axis"]])

	if (missing(mai)) {
		mai <- c(0.42, 0.62, 0.22, 0.22)
		mai[3] <- ifelse(isTRUE(title), 0.62, mai[3])
	} else {
		mai <- rep_len(mai, 4)
	}
	if (missing(height)) {
		height <- max(c(width, npar * 1.5 + 1))
	}
	defpalette <- ifelse(missing(palette), TRUE, all(is.na(palette)))
	if (!defpalette) {
		if (length(palette) < nrun) {
			defpalette <- TRUE
			warning("there is more traces than colors in the specified palette, the default palette will be used")
		}
	}
	if (defpalette) {
		if (npar <= 11) {
			#RColorBrewer::brewer.pal(11, "Spectral")	
			palette <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")	
		} else if (npar <= 20) {
			#https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors
			palette <- c(Red="#e6194b", Green="#3cb44b", Yellow="#ffe119", Blue="#0082c8", Orange="#f58231", Purple="#911eb4", Cyan="#46f0f0", Magenta="#f032e6", Lime="#d2f53c", Pink="#fabebe", Teal="#008080", Lavender="#e6beff", Brown="#aa6e28", Beige="#fffac8", Maroon="#800000", Mint="#aaffc3", Olive="#808000", Coral="#ffd8b1", Navy="#000080", Grey="#808080")[c(1,3,5,7,9,11,13,15,17,19,20,18,16,14,12,10,8,6,4,2)]
		} 
	}
	palette <- palette[as.integer(seq(1, length(palette), length=nrun))]

# graphical device
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
				file <- paste0("mcmc_traces.", sub("^.*_", "", device))
			}
			devname <- file
			isbitmap <- device %in% c("bmp","jpeg","png","tiff")
		} else {
			isfile <- FALSE
			isbitmap <- FALSE
		}
		device <- match.fun(tolower(device))
		if (isTRUE(isbitmap)) {
			device(devname, width=width, height=height, units="in")	
			
		} else {
			device(devname, width=width, height=height)	
		}
	}

	xlim <- rep(list(c(-10, ngen+10)), npar)	
	ylim <- apply(pooled, 2, range, simplify=FALSE)
	dmai <- max(c(0, mai[3] - mai[1]))
	if (npar == 1) {
		mais <- mai[1]+mai[3]
	} else {
		mais <- c(mai[4]/2+mai[3], rep(mai[4], npar-2), mai[1]+mai[4]/2)
	}
	heights <- rep((height - sum(mais)) / npar, npar) + mais
	layout(matrix(seq(npar), npar, 1), heights=heights)

	for (i in seq(npar)) {
		if (i == 1) {
			par(mai=c(mai[4]/2, mai[2:4]))
			plot(1, pooled[1,i], type="n", xlim=xlim[[i]], ylim=ylim[[i]], xlab="", ylab=names(pooled)[i], cex.lab=cex.lab, cex.axis=cex.axis, main=main, cex.main=cex.main, xaxt="n", xaxs="i", ...)
		} else if (i < npar) {
			par(mai=c(mai[4]/2, mai[2], mai[4]/2, mai[4]))
			plot(1, pooled[1,i], type="n", xlim=xlim[[i]], ylim=ylim[[i]], xlab="", ylab=names(pooled)[i], cex.lab=cex.lab, cex.axis=cex.axis, main=NULL, xaxt="n", xaxs="i", ...)
		} else {
			par(mai=c(mai[1], mai[2], mai[4]/2, mai[4]))
			plot(1, pooled[1,i], type="n", xlim=xlim[[i]], ylim=ylim[[i]], xlab="", ylab=names(pooled)[i], cex.lab=cex.lab, cex.axis=cex.axis, main=NULL, xaxs="i", ...)
		}
		for (j in seq(nrun)) {
			lines(seq(ngen), mcmc[[j]][,i], col=palette[j], lwd=lwd, ...)
		}
	}
	
	if (isfile) {
		invisible(dev.off())
	}
}


#' Samples from independent MCMC runs.
#' 
#' @description
#' Box plot showing posterior samples from independent MCMC runs to check their convergence.
#' 
#' @param mcmc Either names of two .mcmc files or a list with objects of class `bpp` created by `read_bpp_mcmc`.
#' @param par Name of parameter to be examined (a name of column in `mcmc$params`) 
#'   or a regular expression defining more of them. 
#' @param interval Which kind of interval to use, either `"HPD"` (default) or `"CPD"`
#'   meaning the highest or the central posterior density.
#' @param p Numeric, proportion of posterior density encompassed by the interval.
#' @param device Character, type of the graphical device. If not specified (default), either `"quartz"` or `"x11"`
#'   is chosen (depending on the platform). It can be also a name of function creating graphical device
#'   for vector graphic or bitmap file. If `NULL` or `NA`, the objects are plotted into the current device (if exists). 
#' @param file Character, the name of file (for graphical devices like `pdf`, `svd` or `png`).
#' @param width Width of graphical device (in inches).
#' @param height Height of graphical device (in inches).
#' @param palette Colors used to distinguish MCMC traces.
#' @param mai Size of outer margins in inches, recycled if necessary
#' @param title Logical, if `TRUE` (default) it uses either `par` or `main` argument (if the latter is specified).
#' @param lwd Numeric, width of the box and whisker lines (and possibly mean point borders).
#' @param bg Background color of the box (and possibly mean point).
#' @param pch Symbol used for the mean point.
#' @param cex Size of the mean point.
#' @param ... Graphical parameters, accepted by `plot.default`.
#' @export

plot_samples <- function(mcmc, par, interval="HPD", p=0.95, device, file, width, height=7, palette, mai, title=TRUE, lwd=1, bg="white", pch=21, cex=1.5, ...) {

	box <- function(sample, no, col, lwd, bg, pch, cex) {
		ran <- range(sample)
		int <- interval(sample, p=p)
		lines(c(no, no), ran, col=col, lwd=lwd)
		lines(c(no-0.10, no+0.10), ran[c(1,1)], col=col, lwd=lwd)
		lines(c(no-0.10, no+0.10), ran[c(2,2)], col=col, lwd=lwd)
		polygon(c(no-0.10, no+0.10, no+0.10, no-0.10), int[c(1,1,2,2)], density=NULL, angle=45, border=col, col=bg, lty=1, lwd=lwd)
		points(no, mean(sample), pch=pch, cex=cex, bg=bg, lwd=lwd, col=col)
	}		

	if (is.character(mcmc)) {
		mcmc <- lapply(mcmc, read_bpp_mcmc) 
	}
	if (inherits(mcmc[[1]], "bpp")) {
		mcmc <- lapply(mcmc, "[[", "params")
	}
	mcmc <- lapply(mcmc, function(x) x[,grep(par, colnames(x), ignore.case=TRUE),drop=FALSE])
	pooled <- do.call(rbind, mcmc)
	nrun <- length(mcmc)
	ngen <- nrow(mcmc[[1]])
	npar <- ncol(mcmc[[1]])
	interval <- match.fun(interval)
	
	ellipsis <- list(...)
	if ("main" %in% names(ellipsis)) {
		title <- TRUE
	} else if (isTRUE(title)) {
		main <- paste0(toupper(substr(par, 1, 1)), tolower(substr(par, 2, nchar(par))))
	} else {
		main <- NULL
	}
	if (length(main) > 0 & !"cex.main" %in% names(ellipsis)) {
		cex.main <- 2
	} else if (length(main) == 0) {
		cex.main <- 0
	}
	cex.lab <- ifelse(!"cex.lab" %in% names(ellipsis), 1.5, ellipsis[["cex.lab"]])
	cex.axis <- ifelse(!"cex.axis" %in% names(ellipsis), 1.5, ellipsis[["cex.axis"]])

	if (missing(mai)) {
		mai <- c(0.22, 0.62, 0.22, 0.22)
		mai[3] <- ifelse(isTRUE(title), 0.62, mai[3])
	} else {
		mai <- rep_len(mai, 4)
	}
	height <- max(c(height, npar * 1.5 + 1))
	if (missing(width)) {
		width <- nrun * 1 + 0.5
	}

	defpalette <- ifelse(missing(palette), TRUE, all(is.na(palette)))
	if (!defpalette) {
		if (length(palette) < nrun) {
			defpalette <- TRUE
			warning("there is more traces than colors in the specified palette, the default palette will be used")
		}
	}
	if (defpalette) {
		if (npar <= 11) {
			#RColorBrewer::brewer.pal(11, "Spectral")	
			palette <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")	
		} else if (npar <= 20) {
			#https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors
			palette <- c(Red="#e6194b", Green="#3cb44b", Yellow="#ffe119", Blue="#0082c8", Orange="#f58231", Purple="#911eb4", Cyan="#46f0f0", Magenta="#f032e6", Lime="#d2f53c", Pink="#fabebe", Teal="#008080", Lavender="#e6beff", Brown="#aa6e28", Beige="#fffac8", Maroon="#800000", Mint="#aaffc3", Olive="#808000", Coral="#ffd8b1", Navy="#000080", Grey="#808080")[c(1,3,5,7,9,11,13,15,17,19,20,18,16,14,12,10,8,6,4,2)]
		} 
	}
	palette <- palette[as.integer(seq(1, length(palette), length=nrun))]

# graphical device
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
				file <- paste0("mcmc_samples.", sub("^.*_", "", device))
			}
			devname <- file
			isbitmap <- device %in% c("bmp","jpeg","png","tiff")
		} else {
			isfile <- FALSE
			isbitmap <- FALSE
		}
		device <- match.fun(tolower(device))
		if (isTRUE(isbitmap)) {
			device(devname, width=width, height=height, units="in")	
			
		} else {
			device(devname, width=width, height=height)	
		}
	}

	xlim <- rep(list(c(0.85, 1 + (nrun - 1) * 0.30 + 0.15)), npar)	
	ylim <- apply(pooled, 2, range, simplify=FALSE)
	dmai <- max(c(0, mai[3] - mai[1]))
	if (npar == 1) {
		mais <- mai[1]+mai[3]
	} else {
		mais <- c(mai[4]/2+mai[3], rep(mai[4], npar-2), mai[1]+mai[4]/2)
	}
	heights <- rep((height - sum(mais)) / npar, npar) + mais
	layout(matrix(seq(npar), npar, 1), heights=heights)

	for (i in seq(npar)) {
		if (i == 1) {
			par(mai=c(mai[4]/2, mai[2:4]))
			plot(1, pooled[1,i], type="n", xlim=xlim[[i]], ylim=ylim[[i]], xlab="", ylab=names(pooled)[i], cex.lab=cex.lab, cex.axis=cex.axis, main=main, cex.main=cex.main, xaxt="n", xaxs="i", ...)
		} else if (i < npar) {
			par(mai=c(mai[4]/2, mai[2], mai[4]/2, mai[4]))
			plot(1, pooled[1,i], type="n", xlim=xlim[[i]], ylim=ylim[[i]], xlab="", ylab=names(pooled)[i], cex.lab=cex.lab, cex.axis=cex.axis, main=NULL, xaxt="n", xaxs="i", ...)
		} else {
			par(mai=c(mai[1], mai[2], mai[4]/2, mai[4]))
			plot(1, pooled[1,i], type="n", xlim=xlim[[i]], ylim=ylim[[i]], xlab="", ylab=names(pooled)[i], cex.lab=cex.lab, cex.axis=cex.axis, main=NULL, xaxt="n", xaxs="i", ...)
		}		
		for (j in seq(nrun)) {
			box(sample=mcmc[[j]][,i], no=1 + (j - 1) * 0.30, col=palette[j], lwd=lwd, bg=bg, pch=pch, cex=cex)			
		}
	}
	
	if (isfile) {
		invisible(dev.off())
	}
}


#' Comparing independent MCMC runs.
#' 
#' @description
#' Displays traces and/or posterior samples from independent MCMC runs to check their convergence.
#' 
#' @param mcmc Either names of two .mcmc files or a list with objects of class `bpp` created by `read_bpp_mcmc`.
#' @param par Name of parameter to be examined (a name of column in `mcmc$params`) 
#'   or a regular expression defining more of them. 
#' @param type Character, type of graph, either `"traces"` or `"samples"` or both.
#' @param interval Which kind of interval to use, either `"HPD"` (default) or `"CPD"`
#'   meaning the highest or the central posterior density.
#' @param p Numeric, proportion of posterior density encompassed by the interval.
#' @param device Character, type of the graphical device. If not specified (default), either `"quartz"` or `"x11"`
#'   is chosen (depending on the platform). It can be also a name of function creating graphical device
#'   for vector graphic or bitmap file. If `NULL` or `NA`, the objects are plotted into the current device (if exists). 
#' @param file Character, the name of file (for graphical devices like `pdf`, `svd` or `png`).
#' @param width Width of graphical device (in inches).
#' @param height Height of graphical device (in inches).
#' @param palette Colors used to distinguish MCMC traces.
#' @param mai Size of outer margins in inches, recycled if necessary
#' @param title Logical, if `TRUE` (default) it uses either `par` or `main` argument (if the latter is specified).
#' @param lwd Numeric, width of the box and whisker lines (and possibly mean point borders).
#' @param bg Background color of the box (and possibly mean point).
#' @param pch Symbol used for the mean point.
#' @param cex Size of the mean point.
#' @param ... Graphical parameters, accepted by `plot.default`.
#' @export

plot_runs <- function(mcmc, par, type, interval="HPD", p=0.95, device, file, width, height, palette, mai, title=TRUE, lwd=1, bg="white", pch=21, cex=1.5, ...) {

	box <- function(sample, no, col, lwd, bg, pch, cex) {
		ran <- range(sample)
		int <- interval(sample, p=p)
		lines(c(no, no), ran, col=col, lwd=lwd)
		lines(c(no-0.10, no+0.10), ran[c(1,1)], col=col, lwd=lwd)
		lines(c(no-0.10, no+0.10), ran[c(2,2)], col=col, lwd=lwd)
		polygon(c(no-0.10, no+0.10, no+0.10, no-0.10), int[c(1,1,2,2)], density=NULL, angle=45, border=col, col=bg, lty=1, lwd=lwd)
		points(no, mean(sample), pch=pch, cex=cex, bg=bg, lwd=lwd, col=col)
	}		

	if (missing(type)) {
		type <- c("traces", "samples")
	}

	ellipsis <- list(...)
	if ("main" %in% names(ellipsis)) {
		title <- TRUE
	}
	defpalette <- ifelse(missing(palette), TRUE, all(is.na(palette)))

	if (length(type) == 1) {
		if (isTRUE(defpalette)) {
			palette <- NA
		}
		if (missing(device)) {
			device <- ifelse(.Platform$OS.type == "unix", "quartz", "x11")
			pdf <- NA
		}
		isfile <- ifelse(missing(pdf), FALSE, !is.na(pdf)) | isTRUE(device == "pdf")
		if (device == "pdf") {
			device <- "pdf"
			pdf <- ifelse(missing(pdf), "mcmc.pdf", pdf)
			pdf <- ifelse(is.na(pdf), "mcmc.pdf", pdf)
		}
		if (type == "traces") {
			if (missing(width)) {
				width <- 7
			}
			if (missing(height)) {
				height <- max(c(width, npar * 1.5 + 1))
			}
			if (missing(mai)) {
				mai <- c(0.42, 0.62, 0.22, 0.22)
				mai[3] <- ifelse(isTRUE(title), 0.62, mai[3])
			} else {
				mai <- rep_len(mai, 4)
			}	
			plot_traces(mcmc=mcmc, par=par, device=device, pdf=pdf, width=width, height=height, palette=palette, mai=mai, title=title, lwd=lwd, ...)
		}
		if (type == "samples") {
			if (missing(width)) {
				width <- nrun * 1 + 0.5
			}
			if (missing(height)) {
				height <- 7
			}
			if (missing(mai)) {
				mai <- c(0.22, 0.62, 0.22, 0.22)
				mai[3] <- ifelse(isTRUE(title), 0.62, mai[3])
			} else {
				mai <- rep_len(mai, 4)
			}	
			plot_samples(mcmc=mcmc, par=par, interval=interval, p=p, device=device, pdf=pdf, width=width, height=height, palette=palette, mai=mai, title=title, lwd=lwd, bg=bg, pch=pch, cex=cex, ...)
		}
	} else {
		if (is.character(mcmc)) {
			mcmc <- lapply(mcmc, read_bpp_mcmc) 
		}
		if (inherits(mcmc[[1]], "bpp")) {
			mcmc <- lapply(mcmc, "[[", "params")
		}
		mcmc <- lapply(mcmc, function(x) x[,grep(par, colnames(x), ignore.case=TRUE),drop=FALSE])
		pooled <- do.call(rbind, mcmc)
		nrun <- length(mcmc)
		ngen <- nrow(mcmc[[1]])
		npar <- ncol(mcmc[[1]])
		interval <- match.fun(interval)
	
		if (!"main" %in% names(ellipsis) & isTRUE(title)) {
			main <- paste0(toupper(substr(par, 1, 1)), tolower(substr(par, 2, nchar(par))))
		} else {
			main <- NULL
		}
		if (length(main) > 0) {
			main <- paste0(main, ": ", c("traces", "samples"))
			if (!"cex.main" %in% names(ellipsis)) {
				cex.main <- 2
			}
		} else {
			cex.main <- 0
		}
		cex.lab <- ifelse(!"cex.lab" %in% names(ellipsis), 1.5, ellipsis[["cex.lab"]])
		cex.axis <- ifelse(!"cex.axis" %in% names(ellipsis), 1.5, ellipsis[["cex.axis"]])
		if (isFALSE(defpalette)) {
			if (length(palette) < nrun) {
				defpalette <- TRUE
				warning("there is more traces than colors in the specified palette, the default palette will be used")
			}
		}
		if (isTRUE(defpalette)) {
			if (npar <= 11) {
				#RColorBrewer::brewer.pal(11, "Spectral")	
				palette <- c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")	
			} else if (npar <= 20) {
				#https://sashat.me/2017/01/11/list-of-20-simple-distinct-colors
				palette <- c(Red="#e6194b", Green="#3cb44b", Yellow="#ffe119", Blue="#0082c8", Orange="#f58231", Purple="#911eb4", Cyan="#46f0f0", Magenta="#f032e6", Lime="#d2f53c", Pink="#fabebe", Teal="#008080", Lavender="#e6beff", Brown="#aa6e28", Beige="#fffac8", Maroon="#800000", Mint="#aaffc3", Olive="#808000", Coral="#ffd8b1", Navy="#000080", Grey="#808080")[c(1,3,5,7,9,11,13,15,17,19,20,18,16,14,12,10,8,6,4,2)]
			} 
		}
		palette <- palette[as.integer(seq(1, length(palette), length=nrun))]

		if (missing(mai)) {
			mai <- c(0.42, 0.62, 0.22, 0.22)
			mai[3] <- ifelse(isTRUE(title), 0.62, mai[3])
		} else {
			mai <- rep_len(mai, 4)
		}	
		right <- nrun * 1 + 0.5 - (mai[2] - mai[4])
		left <- ifelse(missing(width), 7, width - right)
		width <- c(left, right)		
		if (missing(height)) {
			height <- max(c(left, npar * 1.5 + 1))
		}
		lwd <- rep_len(lwd, 2)

# graphical device
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
					file <- paste0("mcmc_runs.", sub("^.*_", "", device))
				}
				devname <- file
				isbitmap <- device %in% c("bmp","jpeg","png","tiff")
			} else {
				isfile <- FALSE
				isbitmap <- FALSE
			}
			device <- match.fun(tolower(device))
			if (isTRUE(isbitmap)) {
				device(devname, width=width, height=height, units="in")	
			} else {
				device(devname, width=width, height=height)	
			}
		}

		dmai <- max(c(0, mai[3] - mai[1]))
		if (npar == 1) {
			mais <- mai[1]+mai[3]
		} else {
			mais <- c(mai[4]/2+mai[3], rep(mai[4], npar-2), mai[1]+mai[4]/2)
		}
		heights <- rep((height - sum(mais)) / npar, npar) + mais

		layout(matrix(seq(2*npar), npar, 2), widths=width/sum(width), heights=heights)

# plot_traces
		xlim <- rep(list(c(-10, ngen+10)), npar)	
		ylim <- apply(pooled, 2, range, simplify=FALSE)
		if (is.null(main)) main_traces <- NULL else main_traces <- main[1]
		for (i in seq(npar)) {
			if (i == 1) {
				par(mai=c(mai[4]/2, mai[2:4]))
				plot(1, pooled[1,i], type="n", xlim=xlim[[i]], ylim=ylim[[i]], xlab="", ylab=names(pooled)[i], cex.lab=cex.lab, cex.axis=cex.axis, main=main_traces, cex.main=cex.main, xaxt="n", xaxs="i", ...)
			} else if (i < npar) {
				par(mai=c(mai[4]/2, mai[2], mai[4]/2, mai[4]))
				plot(1, pooled[1,i], type="n", xlim=xlim[[i]], ylim=ylim[[i]], xlab="", ylab=names(pooled)[i], cex.lab=cex.lab, cex.axis=cex.axis, main=NULL, xaxt="n", xaxs="i", ...)
			} else {
				par(mai=c(mai[1], mai[2], mai[4]/2, mai[4]))
				plot(1, pooled[1,i], type="n", xlim=xlim[[i]], ylim=ylim[[i]], xlab="", ylab=names(pooled)[i], cex.lab=cex.lab, cex.axis=cex.axis, main=NULL, xaxs="i", ...)
			}
			for (j in seq(nrun)) {
				lines(seq(ngen), mcmc[[j]][,i], col=palette[j], lwd=lwd, ...)
			}
		}

# plot_samples
		xlim <- rep(list(c(0.85, 1 + (nrun - 1) * 0.30 + 0.15)), npar)	
		ylim <- apply(pooled, 2, range, simplify=FALSE)
		if (is.null(main)) main_samples <- NULL else main_samples <- main[2]
		for (i in seq(npar)) {
			if (i == 1) {
				par(mai=c(mai[4]/2, mai[4], mai[3], mai[4]))
				plot(1, pooled[1,i], type="n", xlim=xlim[[i]], ylim=ylim[[i]], xlab="", ylab=names(pooled)[i], cex.lab=cex.lab, cex.axis=cex.axis, main=main_samples, cex.main=cex.main, xaxt="n", xaxs="i", yaxt="n", ...)
			} else if (i < npar) {
				par(mai=c(mai[4]/2, mai[4], mai[4]/2, mai[4]))
				plot(1, pooled[1,i], type="n", xlim=xlim[[i]], ylim=ylim[[i]], xlab="", ylab=names(pooled)[i], cex.lab=cex.lab, cex.axis=cex.axis, main=NULL, xaxt="n", xaxs="i", yaxt="n", ...)
			} else {
				par(mai=c(mai[1], mai[4], mai[4]/2, mai[4]))
				plot(1, pooled[1,i], type="n", xlim=xlim[[i]], ylim=ylim[[i]], xlab="", ylab=names(pooled)[i], cex.lab=cex.lab, cex.axis=cex.axis, main=NULL, xaxt="n", xaxs="i", yaxt="n", ...)
			}		
			for (j in seq(nrun)) {
				box(sample=mcmc[[j]][,i], no=1 + (j - 1) * 0.30, col=palette[j], lwd=lwd, bg=bg, pch=pch, cex=cex)
			}
		}
	}
	
	if (isfile) {
		invisible(dev.off())
	}
	
}

