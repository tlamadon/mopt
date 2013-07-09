
#' create an objective function which returns custom
#' output from a bivariate normal density
#' @param norm_mean 2d vector with the true mean
#' @param sample_size number of random draws to compute distance
#' @export
#' @examples
#' myfunc <- objfc.norm2(c(0,0))
#' myfunc(list(m1=0,m2=0))
objfc.norm2.custom <- function(mu=c(0,0),sigma=diag(2),ns=50) {

  obj <- function(p,ms=c()) {

	  t1 <- proc.time()[3]
	  stopifnot(is.list(p))	# we require p be a list

	  # simulate moments 
	  x <- unlist(p[c('x1','x2')]) 	# however rmultnorm doesn't take lists.
	  dat <- data.frame(mvrnorm(n=ns, x, Sigma=sigma))
	  names(dat) <- c("x","y")
	  
	  # say we want a bivariate kernel estimate of this data as an output
	  out <- locfit(~ x+y, data=dat)
	  
	  res = list(
	    p      = p,
	    status = 1,
		output = out)

	   return(res)
   }

   return(obj)
}
