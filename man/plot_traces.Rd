% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_runs.R
\name{plot_traces}
\alias{plot_traces}
\title{Traces of independent MCMC runs.}
\usage{
plot_traces(
  mcmc,
  par,
  device,
  file,
  width = 7,
  height,
  palette,
  mai,
  title = TRUE,
  lwd = 1,
  ...
)
}
\arguments{
\item{mcmc}{Either names of two .mcmc files or a list with objects of class \code{bpp} created by \code{read_bpp_mcmc}.}

\item{par}{Name of parameter to be examined (a name of column in \code{mcmc$params})
or a regular expression defining more of them.}

\item{device}{Character, type of the graphical device. If not specified (default), either \code{"quartz"} or \code{"x11"}
is chosen (depending on the platform). It can be also a name of function creating graphical device
for vector graphic or bitmap file. If \code{NULL} or \code{NA}, the objects are plotted into the current device (if exists).}

\item{file}{Character, the name of file (for graphical devices like \code{pdf}, \code{svd} or \code{png}).}

\item{width}{Width of graphical device (in inches).}

\item{height}{Height of graphical device (in inches).}

\item{palette}{Colors used to distinguish MCMC traces.}

\item{mai}{Size of outer margins in inches, recycled if necessary}

\item{title}{Logical, if \code{TRUE} (default) it uses either \code{par} or \code{main} argument (if the latter is specified).}

\item{lwd}{Numeric, line width.}

\item{...}{Graphical parameters, accepted by \code{plot.default}.}
}
\description{
Plotting of independent MCMC traces (series of sampled paramter values) to check their convergence.
}
