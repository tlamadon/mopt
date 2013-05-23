#' Baragatti Grimaud and Pommeret MCMC chain
#' @export
#' @family algos
algo.bgp <- function(evals, chains, last, cfg, pdesc, priv) {

  # multichain as in Baragatti Grimaud and Pommeret
  # we are going to use N chains with each a different temperiing
  # and we are also going to use energy rings for cross chain jumps

  rr = data.frame()
  ps=list()

  params_to_sample  = cfg$params_to_sample
  params_to_sample2 = paste('p',cfg$params_to_sample,sep='.')
  N                 = cfg$N
  revals = data.frame()

  # compute overall variance/covariance matrix

  # if chains is empty then we accept all
  if (is.null(priv$chain.states)) {
    chain.states = list()
    chain.states$tempering = seq(100,1,l=N)
    chain.states$acc       = seq(0.5,0.5,l=N)
    chain.states$shock_var = seq(cfg$shock_var,cfg$shock_var,l=N)
    priv = list(chain.states = chain.states)
  }
  if (nrow(chains)>0)  {
    # select the last 30 * params^2 observations
    chain.states     = priv$chain.states
    lower_bound_index = pmax(1, nrow(chains) - 30 * length(params_to_sample))
    VV = cov(chains[lower_bound_index:nrow(chains),params_to_sample2]) + 0.001 * diag(mcf$pdesc$ub - mcf$pdesc$lb)
    colnames(VV) <- params_to_sample
    rownames(VV) <- params_to_sample
  }    

  # for each chain, we need the last 2 values
  # and decides whether we accept or reject
  # th new one
  for (c.current in 1:N) {
    # for given current chain, find the latest value
    i.latest = max(subset(chains,chain==c.current)$i)
    # select latest row in the history
    val_old = tail(subset(chains, chain==c.current & i==i.latest),1)
    # select the evaluations
    val_new = subset(evals,chain==c.current)[1,]

    # check for an Energy Ring jump
    if ( 0.05 > runif(1)) {
      # find the list of past values within 10% an energy ring in other chains
      im2 = which( (chains$chain != c.current) &  
                   (abs(chains$value - val_old$value)/
                    abs(chains$value) <0.1))

      if (length(im2)>0) {
       im2 = sample(im2,1)
       val_new = chains[im2,]
       val_new$chain = c.current
      }
    }

    # if there 
    if (all(chains$chain!=c.current)) {
      val_new  = val_old
      next_val = val_new
      ACC = 0
      prob= 0
      status = -1
    # check that old value is NA, in which case we accept new one for sure
    } else if (!is.finite(val_old$value)) {
      next_val = val_new
      ACC  =1
      prob =1
      status = -2
    } else {

      # compute accept reject -- classic Metropolis Hasting
      prob = pmin(1, exp( chain.states$tempering[c.current] * (val_old$value - val_new$value)))
      if (is.na(prob)) prob = 0; 

      # decide on accept reject
      if (prob > runif(1)) {
        next_val = val_new
        ACC = 1
      } else {
        next_val = val_old
        ACC = 0
      }
      status = 1
    }

    cat(sprintf(' value and ratio: %f/%f A=%d prob=%f rate=%f var=%f status=%d\n', val_old$value, val_new$value,  ACC, prob, chain.states$acc[c.current] ,chain.states$shock_var[c.current],status))
    
    # updating sampling variance
    chain.states$acc[c.current]       = 0.9*chain.states$acc[c.current] + 0.1*ACC
    # increase/decrease by 5%
    chain.states$shock_var[c.current] = chain.states$shock_var[c.current] * (1+ 0.05*( 2*(chain.states$acc[c.current]>0.234) -1))

    # change tempering
    # if (abs(cf$acc - 0.234)<0.1)  cf$tempering  =  cf$tempering  * 1.1
    # if (abs(cf$acc - 0.234)>0.22)  cf$tempering =  cf$tempering * 0.9 
    # if (prob > 0.95) cf$theta[i] = cf$theta[i]*1.1
    # if (prob < 0.05) cf$theta[i] = cf$theta[i]*0.9 

    # append to the chain
    rownames(next_val) <- NULL
    revals = rbind(revals,next_val)

    # then we compute a guess for the chain  
    # pick some parameters to update 
    p.new = cfg$initial_value
    p.new$chain = c.current
    p.new[params_to_sample] = next_val[1,params_to_sample2]
    p.new = shockallp(p.new, chain.states$shock_var[c.current], VV, cfg)

    ps[[c.current]] = p.new
  }
 
  priv = list(chain.states = chain.states)
  return(list(ps = ps , priv=priv, evals=revals ))
}


