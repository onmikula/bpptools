% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_priors.R
\name{estimate_priors}
\alias{estimate_priors}
\title{Prior distribution parameters.}
\usage{
estimate_priors(
  dat,
  imap = NULL,
  d_tau,
  d_theta,
  a_tau,
  a_theta,
  seqmap = NULL,
  allele = NULL,
  model = "JC69",
  p = 0.95,
  round = TRUE
)
}
\arguments{
\item{dat}{A data frame or a matrix containing posterior samples of 'tau' (1st column) or 'theta' (2nd column)
or sequence data in the form of a list of locus-specific alignments or a single such alignment (possibly
with concatenated multi-locus data). The row names must be interpretable as individual labels
(possibly with help of \code{seqmap} or \code{allele} arguments).}

\item{imap}{A matrix or data frame (1st column individuals, 2nd column species) or a name of imap file.}

\item{d_tau, d_theta}{The type of prior distributions for tau and theta, either \code{"gamma"} or \code{"invgamma"}.}

\item{a_tau, a_theta}{Alpha parameters of tau and theta priors.}

\item{seqmap}{A data frame mapping names of sequences (1st column) to those of individuals (2nd column).}

\item{allele}{An allele identifier, regular expression distinguishing sequences from the same individual.}

\item{model}{The nucleotide substitution model indicated as in \code{\link[ape:dist.dna]{ape::dist.dna()}}.}

\item{p}{Numeric, the quantile of the highest posterior density used for estimation from posterior samples.}

\item{round}{Logical, whether to round beta parameters of gamma priors to integer values.}
}
\value{
A data frame with columns \code{"alpha"}, \code{"beta"}, \code{"mean"} & \code{"distribution"},
the latter giving the type of distribution.
}
\description{
Estimate means of theta and tau prior distributions which can be considered reasonable.
}
\details{
Any analysis using priors informed by the analyzed data contains some circularity of reasoning.
Still, the strategy is sometimes used to set prior parameters to resonable values, but in that case
the priors should be diffuse (not too informative).
}
