require(mopt)

setwd("~/git/mopt/")

## Create a 5 dimentional function with 10 local maxima.
# creating the objective function   
n = 5 # number of dimensions
k = 10 # number of pics
# we draw 8 minimum points
M = array(runif(k*n),dim=c(n,k))
A = array(runif(k*n),dim=c(1,k)) 
A[1] = 1.5                          # the golabal maxima is at point M[,1]

obj <-  function(x) {
  # we take the minimum distance from any of the mode
  r = apply(M,2, function(v) sum((v-x)^2))
  i = which.min(r)
  return(  -A[i]/ (1+r[i]) )    # the optimizer is looking for 'global' minimum
}

# the solution should be M[,1]


MOPT_OBJ_FUNC <- function(p) {

  if (runif(1) > 0.99) return(NA)   # return NA 1% of the time to check the optimizer deals with it properly

  r = list()
  r$time  = 0
  x = unlist(p[paste('p',1:n,sep='')])  # transforms a list to a vector
  r$value = obj(x)  
  return(r)
}

p = as.list(seq(0.5,0.5,l=n))
names(p) <- paste('p',1:n,sep='')

source('R/inc.mopt.r')
source('R/algo.bgp.r')

# configure optimizer
mcf                  = mopt_config(p)
mcf$wd               = '~/git/mopt/'                  # for objective function to be minimized
mcf$params_all       = paste('p',1:n,sep = '')
mcf$params_to_sample = paste('p',1:n,sep = '')
mcf$iter             = 10000
mcf$mode             = 'serial'
mcf$save_freq        = 200
mcf$shock_var        = 10      # initial value of r shock to parameters   

# adding bounds to the parameters
mcf <- mcf + 
  samplep('p1',0,1) + 
  samplep('p2',0,1) + 
  samplep('p3',0,1) + 
  samplep('p4',0,1) + 
  samplep('p5',0,1)  

# prepare the estimator
mcf = prepare.mopt_config(mcf)  

# starting the optimizer
runMOpt(mcf,FALSE)
  
