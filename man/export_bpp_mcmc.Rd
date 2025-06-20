% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_bpp_mcmc.R
\name{export_bpp_mcmc}
\alias{export_bpp_mcmc}
\title{Exporting BPP output.}
\usage{
export_bpp_mcmc(bpp, trees = NULL, log = NULL, model = "A11")
}
\arguments{
\item{bpp}{The list returned \code{read_bpp_mcmc}.}

\item{trees}{The name of 'trees' file.}

\item{log}{The name of 'log' file.}

\item{model}{One of these: "A00", "A01", "A10" or "A11".}
}
\description{
Exports mcmc trace of BPP into into '.trees' file (in newick format,
readable by FigTree, ape::read.tree etc.) and '.log' file (tab-delimited, readable by Tracer).
}
