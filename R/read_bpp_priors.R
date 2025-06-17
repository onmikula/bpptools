#' Priors from BPP control file.
#' 
#' @description
#' Reads information on priors from BPP control file.
#' 
#' @param ctl The name of the control file.
#' @param which Optional keyword identifying priors to be read; e.g. `"tau"`, `"theta"`, `"mig"` or `"phi"`.
#' @returns A data frame with columns `"alpha"`, `"beta"`, `"mean"` & `"distribution"`
#'   (the latter giving the type of distribution) and possibly `"ncat"` (for the discretized gamma).
#' @export

read_bpp_priors <- function(ctl, which, digits=5) {
	spaces <- "[[:punct:][:space:]]*"
	ctl <- readLines(ctl)
	ctl <- ctl[!grepl("\\*", ctl)]
	ctl <- gsub(paste0("^", spaces, "|", spaces, "$"), "", ctl)
	ctl <- ctl[nchar(ctl) > 0]
	mig <- grep("migration", ctl)
	if (length(mig) > 0) {
		nmig <- as.numeric(sub("^.+=[[:blank:]]*", "", ctl[mig]))
		mig <- ctl[mig + seq(nmig)]
	} else {
		mig <- character(0)
	}
	ctl <- grep("prior", ctl, value=TRUE)
	ctl <- grep("speciesmodel", ctl, value=TRUE, invert=TRUE)
	if (missing(which)) {
		which <- sub("prior", "", sub(paste0(spaces, "=.*$"), "", ctl))
	}
	priortypes <- c("tau", "theta", "mig", "phi", "alpha")
	which <- c(intersect(priortypes, which), setdiff(which, priortypes))
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
	ctl <- sub(paste0("^[[:alpha:]]+", spaces, "=", spaces), "", ctl)
	ctl <- strsplit(ctl, "[[:space:]]+")
	for (i in which(which %in% c("mig", "alpha"))) {
		ctl[[i]] <- c("gamma", ctl[[i]])
	}
	gammas <- sapply(ctl, "[[", 1) == "gamma"
	if (any(gammas)) {
		ctl[gammas] <- lapply(ctl[gammas], function(x) x[-1])
	}
	ab <- do.call(rbind, lapply(lapply(ctl, "[", 1:2), as.numeric))
	priors <- data.frame(alpha=ab[,1], beta=ab[,2], mean=NA, distribution=ifelse(gammas, "gamma", NA), row.names=which)
	priors[which %in% c("tau", "theta") & !gammas, "distribution"] <- "invgamma"
	priors[which == "phi", "distribution"] <- "beta"
	priors$mean[gammas] <- priors$alpha[gammas] / priors$beta[gammas]
	invgammas <- priors$distribution == "invgamma"
	priors$mean[invgammas] <- priors$beta[invgammas] / (priors$alpha[invgammas] - 1)
	betas <- priors$distribution == "beta"
	priors$mean[betas] <- priors$alpha[betas] / (priors$alpha[betas] + priors$beta[betas])
	mig <- strsplit(mig, "[[:space:]]+")
	if (length(mig) > 0) {		
		ownprior <- sapply(mig, length) == 4
		if (any(ownprior)) {
			general <- as.numeric(priors["mig",1:2])
			mig <- lapply(mig, function(x) rep_len(c(x, general), 4))
			mig <- do.call(rbind, mig)
			wpriors <- data.frame(alpha=as.numeric(mig[,3]), beta=as.numeric(mig[,4]), mean=NA, distribution="gamma", from=mig[,1], to=mig[,2], row.names=paste0("mig", formatC(seq(nmig), format="d", width=nchar(nmig), flag=0)))
			wpriors$mean <- wpriors$alpha / wpriors$beta
			priors$from <- NA
			priors$to <- NA
			priors <- rbind(priors, wpriors)			
		}
	}
	alpha <- rownames(priors) == "alpha"
	if (any(alpha)) {
		priors$ncat <- ifelse(alpha, as.numeric(sapply(ctl[alpha], "[[", 3)), NA)
	}
	priors$mean <- round(priors$mean, digits)
	return(priors)
}
