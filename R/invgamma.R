#' Inverse gamma distribution.
#' 
#' @description
#' Density function, distribution function, quantile function and random generation
#'   for the Inverse gamma distribution with parameters `shape` and `scale`.
#' 
#' @param x,q Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param shape Shape parameter, must be greater than one.
#' @param scale Scale parameter, must be positive.
#' @param log,log.p Logical, whether logarithms of densities/probabilities are returned.
#' @param lower.tail Logical, if TRUE (default), probabilities are ---, otherwise, ---.
#' @details The `qinvgamma` function is currently a copy of implementation from `invgamma` package.
#' @export

dinvgamma <- function(x, shape, scale, log=FALSE) {
	d <- scale^shape / gamma(shape) * ((1 / x)^(shape + 1)) * exp(-scale / x)
	if (isTRUE(log)) d <- log(d)
	return(d)
}


#' @export

pinvgamma <- function(q, shape, scale, lower.tail=TRUE, log.p=FALSE) {
	if (isTRUE(lower.tail)) {
		p <- integrate(dinvgamma(x, shape, scale), 0, q)
	} else {
		p <- integrate(dinvgamma(x, shape, scale), q, Inf)
	}
	if (isTRUE(log.p)) p <- log(p)
	return(p)
}


#' @export

qinvgamma <- function(p, shape, scale, lower.tail=TRUE, log.p=FALSE) {
	1 / stats::qgamma(1 - p, shape, scale, lower.tail=lower.tail, log.p=log.p)
}


#' @export

rinvgamma <- function(n, shape, scale) {
	1 / (stats::rgamma(n, shape, scale))
}

