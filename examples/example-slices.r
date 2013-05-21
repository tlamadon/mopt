# generate slices for a model
require(mopt)

# let's take a dummy objective function
MOPT_OBJ_FUNC <- objfc.norm2(c(0,0),ns=2000)

# starting parameters
p <- list(x1=0.5,x2=0.5)
MOPT_OBJ_FUNC(p)
# then we want to setup the mopt

mcf                  = mopt_config(p)
mcf$wd               = getwd()
mcf$params_to_sample = c('x1','x2')
mcf$moments_to_use   = c('m1','m2')
mcf$iter             = 5000
mcf$mode             = 'multicore'
mcf$np_shock         = 2
mcf$shock_var        = 1
mcf$save_freq        = 10

mcf <- mcf + 
  samplep('x1',-1,1) +
  samplep('x2',-1,1)

# adding data moment values
mcf <- mcf + datamoments(c('m1','m2'),
                           c(0,0),
                           c(0.1,0.1))
# prepare to run with OpenMP
require(parallel)
options(mc.cores = detectCores())

# finalize the preparation
mcf <- prepare.mopt_config(mcf)

compute.slices(mcf)
plot.slices(p,mcf)

