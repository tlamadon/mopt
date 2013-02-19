# This demonstrate a simple example. It shows how to write an objective 
# function and also how to write a sampler.

# In this example we are going to estimate the mean of a simple joint normal distribution. 
# We are given 2 moments which are the means in the population, and 
# we want to recover those using sampling. This is a silly example since 
# the mean is already the best estimator, however it will allow us to understand 
# the structure of the program.

# Objective function The signature for the objective function is describe [here](Example).

true_mean = c(2,-1)
obj <- function(p) {

  # simulate moments 
  moments = apply(rmultnorm(50,mu=true_mean,diag(2)),2,mean)
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
algo.example <- function(chains, last, opts, priv) {

}

run.example() <- function() {

# create the mopt object with its config



}