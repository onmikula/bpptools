% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_bpp_tree.R
\name{read_bpp_tree}
\alias{read_bpp_tree}
\title{Starting species tree.}
\usage{
read_bpp_tree(file, as.phylo = FALSE)
}
\arguments{
\item{file}{The name of the control or output file.}

\item{as.phylo}{Logical, whether to return tree as a \code{phylo} object.}
}
\value{
The starting (and possibly fixed) species tree topology in newick or \code{phylo} format
with attribute \code{"tips"}, which contains tip labels in the same order as in the control file.
}
\description{
Reads starting species tree from BPP control file (.ctl) or output file (.out).
}
