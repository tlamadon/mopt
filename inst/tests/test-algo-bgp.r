require(mopt)
require(testthat)
context("algo.bgp")

test_that("simple checks", {

  # let's take a dummy objective function
  MOPT_OBJ_FUNC <- objfc.norm2(c(0,0),ns=2000)

  # starting parameters
  p <- list(x1=0.5,x2=0.5)
  # MOPT_OBJ_FUNC(p)

  # then we want to setup the mopt
  mcf                  = mopt_config(p)
  mcf$wd               = getwd()
  mcf$params_to_sample = c('x1','x2')
  mcf$moments_to_use   = c('m1','m2')
  mcf$mode             = 'serial'
  mcf$iter             = 100
  mcf$algo             = algo.bgp
  mcf$N                = 4
  mcf$objfunc          = MOPT_OBJ_FUNC

  # set the parameter bounds
  mcf <- mcf + 
    samplep('x1',-1,1) +
    samplep('x2',-1,1)

  # adding data moment values
  mcf <- mcf + datamoments(c('m1','m2'),
                             c(0,0),
                             c(0.1,0.1))

  # get a set of parameters
  mcf <- prepare.mopt_config(mcf)
  ps  <- computeInitialCandidates(mcf$N,mcf)

  rd1  <- evaluateParameters(ps,mcf)
  rd2  <- evaluateParameters(ps,mcf)

  # call algo.bgp
  priv=list()
  res = algo.bgp(rd1,rd2, 0, mcf, mcf$pdesc, priv)

  expect_true(length(res$ps)==mcf$N, label='number of parameter set is correct')
  expect_true(length(res$ps)==mcf$N, label='number of parameter set is correct')
  #p1 = res$ps[[1]]
  #expect_true()

})