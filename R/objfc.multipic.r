#' create an objective function n dimensions and
#' k local maximums, returns maximum and function
#' @export
#' @examples	
#' x.max = runif(10)
#' myfunc <- objfc.multipic(10,25,x.max=x.max)
#' myfunc(x.max)
objfc.multipic <- function(n,k,x.max) {

	# creating the objective function
	M = array(runif(k*n),dim=c(n,k))
	A = array(runif(k*n),dim=c(1,k)) 
	A[1] = 1.5

	obj <-  function(p) {

	  x <- unlist(p)	
	  # we take the minimum distance from any of the mode
	  r = apply(M,2, function(v) sum((v-x)^2))
	  i = which.min(r)

	  res = list(
	    p=p,
	    value= -A[i]/ (1+r[i]) ,
	    sim.moments=data.frame(names=c('m1'),value=-A[i]/ (1+r[i])),
	    infos = list(),
	    time  = 0)

	  return( res )
	}

	return(obj)
}

