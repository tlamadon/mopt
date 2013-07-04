# generate custom output for a model
require(mopt)

# let's take a dummy objective function
MOPT_OBJ_FUNC <- objfc.norm2.custom(c(0,0),ns=2000)

# starting parameters
p <- list(x1=0.5,x2=0.5)
MOPT_OBJ_FUNC(p)

# then we want to setup the mopt
mcf                  = mopt_config(p)
mcf$wd               = getwd()
mcf$params_to_sample = c('x1','x2')
mcf$mode             = 'multicore'
mcf$algo             = algo.bgp


# set the parameter bounds
mcf <- mcf + 
  samplep('x1',-1,1) +
  samplep('x2',-1,1)


# prepare to run with OpenMP
require(parallel)
options(mc.cores = detectCores())

# finalize the preparation
mcf <- prepare.mopt_config(mcf)

# compute slices and generate plots
res <- compute.slices2(mcf,ns=30,pad=0.1)


