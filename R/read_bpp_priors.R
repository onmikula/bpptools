#' Priors from BPP control file.
#' 
#' @description
#' Reads information on priors from BPP control file.
#' 
#' @param ctl The name of the control file.
#' @param which A keyword identifying priors to be read; default is `c("tau", "theta")`.
#'   Other options are `"mig"`, `"phi"` or `"alpha"`.
#' @returns A data frame with columns `"alpha"`, `"beta"`, `"mean"` & `"distribution"`
#'   (the latter giving the type of distribution) and possibly `"ncat"` (for the discretized gamma).
#' @export

read_bpp_priors <- function(ctl, which=c("tau","theta")) {
	spaces <- "[[:punct:][:space:]]+"
	ctl <- readLines(ctl)
	ctl <- gsub(paste0("^", spaces, "|", spaces, "$"), "", ctl)
	ctl <- lapply(paste0("^", which, "prior"), grep, x=ctl, value=TRUE)
	null <- which(sapply(ctl, length) == 0)
	if (length(null) > 0) {
		sep <- rev(rep_len(rep(c("", "&", ","), c(1, 1, length(null))), length(null)))
		missingpriors <- paste(as.character(rbind(which[null], sep))[-2*length(null)], collapse=" ")
		plural <- ifelse(length(null) == 1, "prior is", "priors are")
		if (length(null) < length(ctl)) {
			which <- which[-null]
			warning(paste(missingpriors, plural, "missing"))
		} else {
			stop(paste(missingpriors, plural, "missing"))
		}
	}
	ctl <- unlist(ctl)
	ctl <- gsub(paste0("^[[:alpha:]]+", spaces, "=", spaces), "", ctl)
	ctl <- strsplit(ctl, "[[:space:]]+")
	gammas <- sapply(ctl, "[[", 1) == "gamma" | which == "mig"
	if (any(gammas)) {
		ctl[gammas] <- lapply(ctl[gammas], setdiff, y="gamma")
	}
	ab <- do.call(rbind, lapply(lapply(ctl, "[", 1:2), as.numeric))
	priors <- data.frame(alpha=ab[,1], beta=ab[,2], mean=NA, distribution=ifelse(gammas, "gamma", "invgamma"), row.names=which)
	priors$mean[gammas] <- priors$alpha[gammas] / priors$beta[gammas]
	priors$mean[!gammas] <- priors$beta[!gammas] / (priors$alpha[!gammas] - 1)
	alpha <- which == "alpha"
	if (any(alpha)) {
		priors$distribution[alpha] <- "gamma"
		priors$mean[alpha] <- priors$alpha[alpha] / priors$beta[alpha]
		priors$ncat <- ifelse(alpha, as.numeric(sapply(ctl[alpha], "[[", 3)), NA)
	}
	phi <- which == "phi"
	if (any(phi)) {
		priors$distribution[phi] <- "beta"
		priors$mean[phi] <- priors$alpha[phi] / (priors$alpha[phi] + priors$beta[phi])
	}
	return(priors)
}
