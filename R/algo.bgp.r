#' @export
#' @family algos
computeCandidatesBGP <- function(param_data,rd,cf,N,niter) {

  # multichain as in Baragatti Grimaud and Pommeret
  # we are going to use N chains with each a different temperiing
  # and we are also going to use energy rings for cross chain jumps

  rr = data.frame()
  ps=list()

  # compute overall variance/covariance matrix

  # if param_data is empty then we accept all
  if (nrow(param_data)==0) {
    cf$tempering = seq(100,0.001,l=N)
    cf$acc       = seq(0.5,0.5,l=N)
    cf$shock_var = seq(cf$shock_var,cf$shock_var,l=N)

    return(list(nps = computeInitialCandidates(N,cf) , ndt = rd,cf=cf ))
  } else {
    # select the last 30 * params^2 observations
    lower_bound_index = pmax(1, nrow(param_data) - 30 * length(cf$params_to_sample))
    VV = cov(param_data[lower_bound_index:nrow(param_data),cf$params_to_sample])
  }    

  # for each chain, we need the last value, and the new value
  for (i in 1:N) {
    # get last realisation for that chain
    # what is this doing?! looks awful! :-)

    im = which( (param_data$chain == i) & param_data$i == max(param_data$i[param_data$chain==i]))[[1]]

    val_old = param_data[im,]
    val_new = rd[rd$chain==i,]

    # check for an Energy Ring jump
    if ( 0.05 > runif(1)) {
      # find the list of past values within 10% an energy ring in other chains
      im2 = which( (param_data$chain != i) &  
                   (abs(param_data$value - val_old$value)/
                    abs(param_data$value) <0.1))

      if (length(im2)>0) {
       im2 = sample(im2,1)
       val_new = param_data[im2,]
       val_new$chain = i
      }
    }

    # check is new value is NA
    if (all(rd$chain!=i)) {
      val_new = val_old
      next_val = val_new
      ACC = 0
      prob= 0
    # check that old value is NA!
    } else if (!is.finite(val_old$value)) {
      next_val = val_new
      ACC = 1
      prob =1
    } else {

      # compute accept reject -- classic Metropolis Hasting
      prob = pmin(1, exp( cf$tempering[i] * (val_old$value - val_new$value)))
      if (is.na(prob)) prob = 0; 

      # decide on accept reject
      if (prob > runif(1)) {
        next_val = val_new
        ACC = 1
      } else {
        next_val = val_old
        ACC = 0
      }
    }

    cat(' value and ratio:', val_old$value, '/' ,val_new$value, ' A=', ACC, '-',prob,' rate=',cf$acc[i] ,' var=',cf$shock_var[i], '\n')
    
    # updating sampling variance
    cf$acc[i] = 0.9*cf$acc[i] + 0.1*ACC
    # increase/decrease by 5%
    cf$shock_var[i] = cf$shock_var[i] * (1+ 0.05*( 2*(cf$acc[i]>0.234) -1))


    # change tempering
    # if (abs(cf$acc - 0.234)<0.1)  cf$tempering  =  cf$tempering  * 1.1
    # if (abs(cf$acc - 0.234)>0.22)  cf$tempering =  cf$tempering * 0.9 
    # if (prob > 0.95) cf$theta[i] = cf$theta[i]*1.1
    # if (prob < 0.05) cf$theta[i] = cf$theta[i]*0.9 

    # append to the chain
    rownames(next_val) <- NULL
    rr = rbind(rr,next_val)

    # then we compute a guess for the chain  
    # pick some parameters to update 

    val_new = shockallp(next_val, cf$shock_var[i], VV, cf)

    ps[[i]] = val_new
  }
 
   return(list(nps = ps , ndt = rr, cf=cf ))
}


