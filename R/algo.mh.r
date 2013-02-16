#' Metropolis Hasting MCMC chain updating
#' @export
#' @family algos
algo.mh <- function(param_data,rd,cf,N,niter) {

  # just a mulitchain MH

  # then we can compute the Pr of the last evaluations at those
  # values
  rr = data.frame()

  # if param_data is empty then we accept all
  if (nrow(param_data)==0) {
    cf$tempering = 1

    return(list(nps = computeInitialCandidates(N,cf) , ndt = rd,cf=cf ))
  }

  # for each chain, we need the last value, and the new value
  for (i in 1:N) {
    # get last realisation for that chain
    im = which( (param_data$chain == i) & param_data$i == max(param_data$i[param_data$chain==i]))[[1]]
    val_old = param_data[im,]
    val_new = rd[rd$chain==i,]

    # compute accept reject
    prob = pmin(1, exp( cf$tempering * (val_old$value - val_new$value)))
    if (is.na(prob)) prob = 0; 

    # decide on accept reject
    if (prob > runif(1)) {
      next_val = val_new
      ACC = 1
    } else {
      next_val = val_old
      ACC = 0
    }
    #cat(' value and ratio:', epdf(val_old$value),' / ',epdf(val_new$value),'  ',val_old$value, '/' ,val_new$value, ' A=', ACC, '-',prob,' rate=',cf$acc ,'\n')
    cat(' value and ratio:', val_old$value, '/' ,val_new$value, ' A=', ACC, '-',prob,' rate=',cf$acc ,'\n')
    
    # updating sampling variance
    cf$acc = 0.99*cf$acc + 0.01*ACC
    cf$shock_var = cf$shock_var + 1/niter*( 2*(cf$acc>0.234) -1)

    # change tempering
    if (abs(cf$acc - 0.234)<0.1)  cf$tempering  =  cf$tempering  * 1.1
    if (abs(cf$acc - 0.234)>0.22)  cf$tempering =  cf$tempering * 0.9 
    # if (prob > 0.95) cf$theta[i] = cf$theta[i]*1.1
    # if (prob < 0.05) cf$theta[i] = cf$theta[i]*0.9 

    # append to the chain
    rr = rbind(rr,next_val)

    # then we compute a guess for the chain  
    # pick some parameters to update 
    param <- sample(cf$params_to_sample, cf$np_shock)

    val_new = next_val
    for (pp in param) {
       # update value
       val_new[[pp]] = shockp(pp, val_new[[pp]] , cf$shock_var * cf$theta[i] , cf)
    }

    ps[[i]] = val_new
  }
 
   return(list(nps = ps , ndt = rr, cf=cf ))
}