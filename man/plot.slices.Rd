% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fun.slices.r
\name{plot.slices}
\alias{plot.slices}
\title{plot.slices}
\usage{
\method{plot}{slices}(file = NULL, outpath = "", type = "png")
}
\arguments{
\item{file}{path/to/est.slices.RData if not in getwd()}

\item{outpath}{path to directory where to save plots}

\item{type}{string indicating file type to save plot as. currently png and pdf only.}
}
\description{
generates plots for each moments/parameter combinations
depends on output of function compute.slices(), which saves
a dataframe with plotting data and the used mopt config into a 
a file \code{est.slices.RData}
}

