\name{plot.slices}
\alias{plot.slices}
\title{plot.slices}
\usage{
  plot.slices(file = NULL, outpath = "", type = "png",
    moments = NULL)
}
\arguments{
  \item{file}{path/to/est.slices.RData if not in getwd()}

  \item{outpath}{path to directory where to save plots}

  \item{type}{string indicating file type to save plot as.
  currently png and pdf only.}

  \item{moments}{vector of strings for which moments to
  show on the plot. By default will take all moments in
  mcf$data.moments.}
}
\description{
  generates plots for each moments/parameter combinations
  depends on output of function compute.slices(), which
  saves a dataframe with plotting data and the used mopt
  config into a a file \code{est.slices.RData}
}

