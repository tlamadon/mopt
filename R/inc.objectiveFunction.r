# creating the objective function
n = 5 # number of dimensions
k = 10 # number of pics
# we draw 8 minimum points
M = array(runif(k*n),dim=c(n,k))
A = array(runif(k*n),dim=c(1,k)) 
A[1] = 1.5

obj <-  function(x) {
  # we take the minimum distance from any of the mode
  r = apply(M,2, function(v) sum((v-x)^2))
  i = which.min(r)
  return(  -A[i]/ (1+r[i]) )
}

MOPT_OBJ_FUNC <- function(p) {
  r = list()
  r$time  = 0
  r$value = obj(c(p))
  return(r)
}

