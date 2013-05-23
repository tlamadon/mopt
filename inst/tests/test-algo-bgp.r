require(mopt)
require(testthat)
context("algo.bgp")

test_that("simple checks", {

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
  mcf$mode             = 'serial'
  mcf$iter             = 100
  mcf$algo             = algo.bgp
  mcf$N                = 4

  # set the parameter bounds
  mcf <- mcf + 
    samplep('x1',-1,1) +
    samplep('x2',-1,1)

  # adding data moment values
  mcf <- mcf + datamoments(c('m1','m2'),
                             c(0,0),
                             c(0.1,0.1))

  # get a set of parameters
  ps = mopt:::computeInitialCandidates(mcf$N,mcf)
  mcf <- prepare.mopt_config(mcf)

  # evaluate them
  rd = mopt:::evaluateParameters(ps,mcf)

  # call algo.bgp
  priv=list()
  algo.bgp(rd, 0, mcf, mcf$pdesc, priv)

})