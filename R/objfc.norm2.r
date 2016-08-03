require(MASS)

#' create an objective function which computes
#' the distance between simulated means and true means
#' @param norm_mean 2d vector with the true mean
#' @param sample_size number of random draws to compute distance
#' @export
#' @examples
#' myfunc <- objfc.norm2(c(0,0))
#' myfunc(list(m1=0,m2=0))
objfc.norm2 <- function(mu=c(0,0),sigma=diag(2),ns=50,failp=0) {

  obj <- function(p,ms=c()) {

	  t1 <- proc.time()[3]
	  stopifnot(is.list(p))	# we require p be a list

	  # through an expection sometimes
	  if (runif(1)<failp) {
	    return(NA)
	  }
	  
	  # simulate moments 
	  x <- unlist(p[c('x1','x2')]) 	# however rmultnorm doesn't take lists.
	  moments = apply( mvrnorm(n=ns, x, Sigma=sigma) ,2,mean)
	  value = mean((mu - moments)^2)

	  sm  =data.frame(names=c('m1','m2'),value=moments)
	  sm$names = paste(sm$names)

	  # get system info
	  si <- Sys.info()
	  nname <- si["nodename"]
	  names(nname) <- NULL

	  # get memory usage
	  mem <- paste(gc()[1,2],"MB used")

	  res = list(
				 p      = p,
				 chain  = p$chain,
				 value  = value,
				 status = 1,
				 sm     = sm,
				 infos  = list(node=nname,mem=mem),
				 time   = as.numeric(proc.time()[3] - t1))

   }

   return(obj)
}
