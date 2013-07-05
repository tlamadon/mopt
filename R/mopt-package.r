#' Parallel optimizer for R 
#'
#' implements MCMC in parallel on OpenMP or MPI
#' 
#' The main goal of this package is to provide a library for simulated moment
#' estimation done in parallel on either MPI or openMP. The code is writtern 
#' specificaly for moment estiamtion in the sense that the objective function
#' should return, not only an overall value, but also each moments. This is 
#' then used in the analysis part. For instance the researcher can analyze
#' how each moment is affected by each parameter.
#'
#' The main function is \code{\link{runMOpt}}, and the associated \code{\link{mopt_config}}. This gives
#' full control to the researchers of the methods used to look for the minimum distance estimator. Several
#' algorithms are implemented and we encourage other people to write more of them and leverage the parallelization
#' structure built into this library.
#' 
#' A very useful function is \code{\link{compute.slices}} and the associated \code{\link{plot.slices}} which will generate plots of how
#' the different moments vary as you vary each parameters separately.
#' 
#' @import data.table MASS locfit digest
#' @docType package
#' @name mopt
#' @author Ran Gu \email{lionup@@gmail.com}, Thibaut Lamadon \email{thibaut.lamadon@@gmail.com} , Florian Oswald \email{florian.oswald@@gmail.com}
#' @seealso runMOpt
#' @references
#' \url{https://github.com/tlamadon/mopt}
#' @aliases mopt mopt-package
NULL


