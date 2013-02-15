#' @export
#' @family algos
computeCandidatesWangLandau <- function(param_data,rd,cf,N,niter) {

  # updating sampling distribution
  # ------------------------------

  # we are doing a Wang-Landau
  # so we need to compute the histogram over 'Energies'
  # actually why not do this non parametrically
  # we take all past value and cmpute pdf

  # for each value we can compute how many
  # evaluations have happen close to them
  # and use this to compute the accept reject

  # then we can compute the Pr of the last evaluations at those
  # values
  rr = data.frame()

  # if param_data is empty then we accept all
  if (nrow(param_data)==0) {
    return(list(nps = computeInitialCandidates(N,cf) , ndt = rd,cf=cf ))
  }

  # for each chain, we need the last value, and the new value
  for (i in 1:N) {
    # get last realisation for that chain
    im = which( (param_data$chain == i) & param_data$i == max(param_data$i[param_data$chain==i]))[[1]]
    val_old = param_data[im,]
    val_new = rd[rd$chain==i,]

    epdf <- function(x) {
      hx =  findInterval(x,cf$breaks)
      if (hx==0) return(10)
      return(cf$theta[hx])
    }

    # compute accept reject
    prob = pmin(1, epdf(val_old$value) / epdf(val_new$value))
    # prob = pmin(1, exp((val_old$value - val_new$value)))
    if (is.na(prob)) prob = 0; 

    # decide on accept reject
    if (prob > runif(1)) {
      next_val = val_new
      ACC = 1
    } else {
      next_val = val_old
      ACC = 0
    }
    cat(' value and ratio:', epdf(val_old$value),' / ',epdf(val_new$value),'  ',val_old$value, '/' ,val_new$value, ' A=', ACC, '-',prob,'\n')
    
    # update the distribution
    # first we get the histogram for this value
    hx = findInterval(next_val$value,cf$breaks)
    cf$theta[hx] = cf$theta[hx] * exp(1/cf$kk)
    cf$theta     = cf$theta * exp( -1/(cf$kk * (1:length(cf$theta))))
    cf$theta     = cf$theta / sum(cf$theta)
    cf$nu[hx]    = cf$nu[hx] + 1/niter
    cf$nu        = cf$nu - 1/niter

    # updating sampling variance
    cf$acc = 0.9*cf$acc + 0.1*ACC
    cf$shock_var = cf$shock_var + 1/niter*( 2*(cf$acc>0.234) -1)

    # check flat histogram
    if (max( abs( cf$nu - 1/(1:50)) < 0.01)) {
        cat('updating kk')
        cf$kk = cf$kk + 1
        cf$nu = 0 * cf$nu
    }

    # change tempering
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