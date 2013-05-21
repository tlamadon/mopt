#' Runs a serial simplex
#' 
#' @export
#' @example examples/example-serial-optim.r
run.simplex <- function(mcf) {

  # reading configuration
  # =====================
  pdesc = mcf$pdesc
  p     = mcf$initial_value

  infos <- list(count=0,best=10000)

  # create the objective function
  objfunc <- function(x) {
    # get the parameter to sample
    p[mcf$params_to_sample] = x
    rr = MOPT_OBJ_FUNC(p)

    if (rr$value<infos$best) {
      infos$best <<- rr$value
    }
    infos$count <<- infos$count+1

    cat(sprintf("[%i] value=%f best=%f \n",infos$count,rr$value,infos$best))
    return(rr$value)
  }

  # get bounds
  ub = rep(0,length(mcf$params_to_sample))
  lb = rep(0,length(mcf$params_to_sample))
  for (i in 1:length(mcf$params_to_sample)) {
    pp = mcf$params_to_sample[i]
    ub[i] = mcf$pdesc$ub[mcf$pdesc$param==pp]
    lb[i] = mcf$pdesc$lb[mcf$pdesc$param==pp]
  }

  par0 = as.numeric(p[mcf$params_to_sample])
  control = list(maxfeval = mcf$iter,info=TRUE)
  # start the simplex
  res = hjkb(par0, objfunc, lower=lb,upper=ub)

  return(res)
}