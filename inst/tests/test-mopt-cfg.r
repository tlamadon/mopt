require(mopt)
require(testthat)
context("mopt.cfg")

test_that("checking structure creation", {

  p <- list(x1=0.5,x2=0.5)

  # then we want to setup the mopt
  mcf                  = mopt_config(p)
  mcf$params_to_sample = c('x1','x2')
  mcf$moments_to_use   = c('m1','m2')
  mcf$mode             = 'serial'
  mcf$iter             = 100
  mcf$algo             = algo.bgp
  mcf$objfunc          = function(x) x

  # set the parameter bounds
  mcf <- mcf + 
    samplep('x1',-1,1) +
    samplep('x2',-1,1)

  expect_true(nrow(mcf$pdesc)==2,label="pdesc has correct number of elements")
  expect_true(is.numeric(mcf$pdesc$ub),label="pdesc bounds are numerics")

  # adding data moment values
  mcf <- mcf + datamoments(c('m1','m2'),
                             c(0,0),
                             c(0.1,0.1))

  expect_true(nrow(mcf$pdesc)==2,label="pdesc has correct number of elements")
  expect_true(is.numeric(mcf$pdesc$ub),label="pdesc bounds are numerics")
  
})