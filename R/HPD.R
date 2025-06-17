#' The highest posterior density interval.
#' 
#' @description
#' Estimates the highest posterior density interval.
#' 
#' @param x A numeric vector, the posterior sample of a model paramater.
#' @param p The proportion of posterior density be included in the interval, default is 0.95.
#' @details "The 95% HPD stands for highest posterior density interval and represents the most compact interval
#'   on the selected parameter that contains 95% of the posterior probability." (Drummond & Bouckaert 2015, p. 89)
#' @returns A numeric vector with lower and upper limit of HPD interval.
#' @export

HPD <- function(x, p=0.95) {
	x <- na.omit(x)
	if (length(x) > 0) {
		n <- round(length(x) * p)
		w <- seq(1, length(x) - n)
		x1 <- unname(sort(x))
		x2 <- rev(x1)
		int1 <- unname(rbind(x1[w + n], x1[w]))
		hpd1 <- sort(int1[, which.min(diff(int1))])
		int2 <- unname(rbind(x2[w + n], x2[w]))
		hpd2 <- sort(int2[, which.min(diff(int2))])
		hpd <- list(hpd1, hpd2)
		hpd <- hpd[[which.min(sapply(hpd, diff))]]	
	} else {
		hpd <- c(NA, NA)
	}
	return(hpd)
}


#' The central posterior density interval.
#' 
#' @description
#' Estimates the central posterior density interval.
#' 
#' @param x A numeric vector representing posterior sample of a model parameter.
#' @param p The proportion of posterior density be included in the interval, default is 0.95.
#'   It implies the maignitude of exterme quantiles of posterior density, `q = (1 - p) / 2`.
#' @details The central posterior density interval corresponds to the range of values after trimming
#'   of tails corresponding to the extreme quantiles of posterior density.
#' @returns A numeric vector with lower and upper limit of CPD interval.
#' @export

CPD <- function(x, p=0.95) {
	x <- na.omit(x)
	if (length(x) > 0) {
		p <- (1 - p) / 2
		q <- quantile(x, probs=c(p, 1 - p))
		cpd <- range(x[x > q[1] & x <= q[2]])
	} else {
		cpd <- c(NA, NA)
	}
	return(cpd)
}
