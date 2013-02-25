# This demonstrate a simple example. It shows how to write an objective 
# function and also how to write a sampler.

# In this example we are going to estimate the mean of a simple joint normal distribution. 
# We are given 2 moments which are the means in the population, and 
# we want to recover those using sampling. This is a silly example since 
# the mean is already the best estimator, however it will allow us to understand 
# the structure of the program.

# Objective function The signature for the objective function is describe [here](Example).

# what libraries are needed to run this example?
library( ? )



true_mean = c(2,-1)
obj <- function(p) {

  stopifnot(is.list(p))	# we require p be a list

  # simulate moments 
  p <- unlist(p) 	# however rmultnorm doesn't take lists.
  moments = apply(rmultnorm(50,mu=p,diag(2)),2,mean)
  value = mean((true_mean - moments)^2)

  res = list(
    p=p,
    value=value,
    sim.moments=data.frame(names=c('m1','m2'),value=moments),
    infos = list(),
    time  = 0)

   return(res)
}

# the minimum an algo has to do is generate new 
# candidates... here we will just randomly sample
# with a normal and give back N parameters
algo.example <- function(chains, last, opts, params, priv) {

  # first we need to compute the acceptance probability
  # by taking the exponential of the difference between 
  # new value and old value within each chain
  for (i in 1:N) {

    # compute the acceptance probability and the draw
    chains[J(i,last)]$acc.pr = exp( chains[J(i,last)]$value - chains[J(i,last-1)]$value)
    chains[J(i,last)]$acc    = (runif() > chains[J(i,last)]$acc)

    # overwrite chain if not accepted
    if (acc==FALSE) {
      chains[J(i,last)]   <- chains[J(i,last-1)]
      chains[J(i,last)]$n = last
    }

    # second we jump the parameters within their bounds
    # here I am just uniformly picking them using the param description
    vals = rnorm(nrow(params))
    p = jumpParams.normalAndMirrored()
  }

  return(p)
}

run.example() <- function() {

  # create the mopt object with its config
  mcf <- mopt(init=p) + 
        opts(N =4,)


}
