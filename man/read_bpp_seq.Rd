% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_bpp_input.R
\name{read_bpp_seq}
\alias{read_bpp_seq}
\title{Reading BPP sequence file.}
\usage{
read_bpp_seq(seqfile, nloci = -1L, locnames)
}
\arguments{
\item{seqfile}{The name of sequence file.}

\item{nloci}{The number of loci to be imported (in the same order as in seqfile). The default means "all loci".}

\item{locnames}{A character vector with names of loci, if missing, the loci are given names "Lddd" ("d" is for digits).}
}
\value{
The list of matrices with locus-specific sequence alignments.
}
\description{
Reads sequences from BPP sequence file.
}
