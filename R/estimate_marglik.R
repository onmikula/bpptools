#' Marginal likelihood by thermodynamic integration.
#' 
#' @description
#' Extracts the estimate of marginal likelihood of a model from the outputs of
#' a series of MCMC runs conducted in the thermodynamic integration procedure.
#' 
#' @param output Names of the .out files.
#' @param weights Gauss-Legendre weights of a name of .csv file containing them.
#' @returns The logarithm of the marginal likelihood.
#' @export

estimate_marglik <- function(output, weights) {
	bflines <- sapply(lapply(output, readLines), function(x) grep("BFbeta", x, ignore.case=TRUE, value=TRUE))
	bflines <- do.call(rbind, strsplit(bflines, split="[[:blank:]=]+"))
	bflines <- data.frame(beta=as.numeric(bflines[,2]), lnL=as.numeric(bflines[,4]))
#	bflines <- data.frame(beta=as.numeric(bflines[,2]), weight=as.numeric(bflines[,4]), lnL=as.numeric(bflines[,6]))
	if (is.character(weights)) {
		bflines$weight <- read.table(weights, header=TRUE, sep=",")$weight
	} else {
		bflines$weight <- weights
	}
	margL <- sum(bflines$weight * bflines$lnL / 2)
	return(margL)
}


#' Bayes factor estimated by Savage-Dickey ratio.
#' 
#' @description
#' Calculates Savage-Dickey ratio as an estimate of Bayes factor comparing the analyzed 
#' model to a nested model, where the parameter in question has a specified null value.
#' 
#' @param prior The type of prior distribution, must be a character string
#'   for which quantile function named `paste0("p", prior)` is available.
#' @param pars A numeric vector with parameters of the prior distribution.
#'   They are used as arguments of the quantile function.
#' @param posterior A posterior sample of the parameter in question.
#' @param null A numeric vector giving the null region, i.e. the interval close to the null value.
#'   If just a single number is supplied, the other limit is supposed to be zero.
#' @returns The value of Savage-Dickey ratio.
#' @export

estimate_sdratio <- function(prior, pars, posterior, null) {

	if (length(null) == 1) {
		null <- c(0, null)
	} else {
		null <- sort(null)
	}
	
	cdf <- match.fun(paste0("p", prior))
	if (length(pars) == 1) {
		num <- cdf(null[2], pars, lower.tail=TRUE) - cdf(null[1], pars, lower.tail=TRUE)
	} else {
		num <- cdf(null[2], pars[1], pars[2], lower.tail=TRUE) - cdf(null[1], pars[1], pars[2], lower.tail=TRUE)		
	}
	
	posterior <- as.numeric(posterior)
	den <- sum(posterior <= null[2] & posterior > null[1]) / length(posterior)
	if (den == 0) {
		dens <- stats::density(posterior, from=0, to=max(posterior), cut=0, bw="nrd0", kernel="gaussian")
		int <- dens$x <= null[2] & dens$x > null[1]
		dx <- diff(dens$x[int])
		y <- dens$y[int]
		den <- sum(dx * (y[-length(y)] - min(y))) + sum(0.5 * dx * diff(y))
	}
	
	return(num / den)
	
}
