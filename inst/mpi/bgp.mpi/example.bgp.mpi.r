# estimate parameters of bivariate normal with MCMC algorithm BGP on cluster.

require(mopt)

# let's take a dummy objective function
MOPT_OBJ_FUNC <- objfc.norm2(c(0,0),ns=2000)

# starting parameters
p <- list(x1=0.5,x2=0.5)
#MOPT_OBJ_FUNC(p)

# then we want to setup the mopt
mcf                  = mopt_config(p)
mcf$wd               = getwd()
mcf$params_to_sample = c('x1','x2')
mcf$moments_to_use   = c('m1','m2')
mcf$mode             = 'mpi'
mcf$source_on_slaves = 'slaves.r'
mcf$iter             = 300
mcf$algo             = algo.bgp
mcf$objfunc          = MOPT_OBJ_FUNC

# set the parameter bounds
mcf <- mcf + 
  samplep('x1',-1,1) +
  samplep('x2',-1,1)

# adding data moment values
mcf <- mcf + datamoments(c('m1','m2'),
                           c(0,0),
                           c(0.1,0.1))

# prepare to run with OpenMP
# require(parallel)
# options(mc.cores = detectCores())

# finalize the preparation
mcf <- prepare.mopt_config(mcf)

# run optimization
res = runMOpt(mcf,FALSE)

print(res)

cat('DONE.')

stopCluster(mcf$cl)


