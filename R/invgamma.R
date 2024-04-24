#' Inverse gamma distribution.
#' 
#' @description
#' Density and quantile function for the Inverse gamma distribution with parameters `shape` and `scale`.
#' 
#' @param x Vector of quantiles.
#' @param p Vector of probabilities.
#' @param shape Shape parameter, must be greater than one.
#' @param scale Scale parameter, must be positive.
#' @param log,log.p Logical, whether logarithms of densities/probabilities are returned.
#' @details The `qinvgamma` function is currently a copy of implementation from `invgamma` package.
#' @export

dinvgamma <- function(x, shape, scale, log=FALSE) {
	d <- scale^shape / gamma(shape) * ((1 / x)^(shape + 1)) * exp(-scale / x)
	if (isTRUE(log)) d <- log(d)
	return(d)
}


#' @export

qinvgamma <- function(p, shape, scale, lower.tail=TRUE, log.p=FALSE) {
	stats::qgamma(1 - p, shape, 1/scale, lower.tail=lower.tail, log.p=log.p)^(-1)
}

