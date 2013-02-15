# we want to have a smart optimizer:
# - the optimizer should be able to deal with mulitgoal
# this is because estimation is difficult and usually
# the programmer knows about which moment should
# be controlling which parameter

# first we need the use of an extended objective
# function. In our case, this function should return
# not only the total value, but a decomposition of that value
# for example that might be the value on a set of moments

# Then the user should able to supply a design matrix, that
# defines the order in which both moments and parameters should
# be considered for estimation
# in this mulitspte estimation, the optimizer would start
# by maximizing a subset of moments by shocking only 
# a given set of parameters. It would then move to changing
# deeper moments
# In some sense this says that the optimizer should
# first focus on subspaces of the parameter/objective value state

# As it estimates, the optimizer should learn aobut the influence matrix.
# from evaluations, it can learn which parameters supposively affect
# the moments, and use that infomration to compute steepest descent
# this is a sort of surrogate method

# the objective function takes in relaxation parameters, tolerances
# and max iteration, this can also be controled by the estimator
# the objective function should return a list of achieved tolerance
# and number of iterations
# possiby the funciton should be able to use a cache and return
# whether it was able to use it

# also the optimizer shuold be able to apply transforms to the 
# parameters

# Feedback
# the optimizer can get stuck an a low value (maybe a local min)
# I need to have a better search theory with a smoother acceptance
# drawing a probability, or picking the claue to be shocked from past
# values with higher probabilities on smaller values.... -- keep several bests
# or do multiple chains......
# or ability to change several at the same time 

# TODO
#  - add 'param.' in front of the parameters
#  - store value function distances / tolerances / loop counts
#  - save which node computed which value
#  - try to do some load balancing (clusterApplyLB)

require(plyr)
require(MSBVAR)
source('~/git/Utils/R/inc.utils.r')

# functions that allows to mirror the parameters
# to within the range of values
fitMirror <- function (x,LB=NA,UB=NA) {
  test1 = TRUE
  test2 = TRUE
  while( test1 | test2) {
    if ( !is.na(UB) & (x>UB)) {
      x = 2*UB - x 
    } else {
      test1 = FALSE
    } 
    if ( !is.na(LB) & (x<LB)) {
      x = 2*LB - x
    } else {
      test2 = FALSE
    } 
  }
 return(x)
}

# setting up the config
# should be given 
# list of model parameters including, 
#   * priority, 
#   * moment a priori connection
#   * transform (optional)
#   * bounds (optional)
# list of computational parameters
#   * tolerance
#   * iteration limit

#' create the configuration for the MCMC 
#' chain. 
# @seelalso samplep
mopt_config <- function(p) {
	stopifnot(is.list(p))
 cf = list()

 cf$mode           = 'mpi'
 cf$iter           = 10
 cf$use_last_run   = TRUE
 cf$file_chain     = 'evaluations.dat'
 cf$file_lastparam = 'param_submit.dat'
 cf$wd             = '~/git/ssp/R/'
 cf$source_on_nodes = 'run.modelest.r'
 cf$params_to_sample = c()	# is that a vector of names or indices?
 cf$run            = 0
 cf$shock_var      = 0.1
 cf$moments_to_use = c()
 cf$moments.data   = c()
 cf$moments.sd     = c()
 cf$np_shock       = 1
 cf$save_freq      = 25
 cf$initial_value = p	

 pd = data.frame()
 for (n in names(p)) {
   pd = rbind(pd,data.frame(param=n, lb=NA , ub=NA))
 }
 rownames(pd) <- names(p)
 cf$pdesc = pd

 class(cf) <- 'mopt_config'

 return(cf)
}


samplep <- function(pp,lb,ub) {
  res = list()
  class(res) <- 'mopt_smaplep'
  res$pp = pp
  res$lb = lb
  res$ub = ub
  return(res)
}

datamoments <- function(names,values,sds) {
  res = list()
  rr = data.frame(moment=names,value = values, sd = sds)
  class(res) <- 'mopt_dmoms'
  res$dd = rr
  return(res)
}


"+.mopt_config" <- function(cf,argb) {

 if (class(argb)=='mopt_smaplep') {
    cf$pdesc[argb$pp,'lb'] = argb$lb
    cf$pdesc[argb$pp,'ub'] = argb$ub
 }

 if (class(argb)=='mopt_dmoms') {
    cf$data.moments = argb$dd
 }

 return(cf)
}

mopt_obj_wrapper <- function(p) {
  m = tryCatch( {

    # get result
    r = MOPT_OBJ_FUNC(p)  

    if (!is.list(r)) {
      r = list(value=r)
    }

    # check that there is a status 
    # and that values is not NA
    if (!('status' %in% names(r))) { 
      r$status=1
    }

    if (is.nan(r$value) | is.na(r$value)) {
      r$status=-1
    }

    r
  },error = function(e) {
    list(status=-1,error=e$message)
  } ) 
  return(m) 
  # returns NA if there is any sort of crash
  # later it would be good to get the error message
}

evaluateParameters <- function(ps,cf) {
    
    #save evaluations to file
    #save(ps,file='lasteval.dat')

    cat('Sending parameter evaluations...')
    vals = cf$mylapply(ps,mopt_obj_wrapper)

	# what is ICNOV doing?
    ICONV = rep(FALSE,length(vals)) 
    for (i in 1:length(vals)) {
      if (!('status' %in% names(vals[[i]])))  {
        ICONV[i]=FALSE;
      } else { 
        ICONV[i] = vals[[i]]$status==1
      }
    }
    # we start by removing any NA value
    #ICONV = unlist(Map(function(e) {e$status==1},vals))
    cat('done (vals=',sum(ICONV),'/',length(ps),')\n')    
    vals = vals[ICONV] # finish treating non convergence correctly!!!!!
    
    # transform to array structure
    rvals = list2df(vals) 
    vals = unlist(rvals$value)  # store the objective values

    # saving values to data.frame for later use
    rd       = list2df(ps[])
    rd$i     = cf$i
    rd$run   = cf$run
    rd = merge(rvals , rd[ICONV,c('chain',setdiff(names(rd),names(rvals)))], by='chain') # save all returned values
    return(rd)
}

getParamStructure <- function() {
  p = getStartingValues()
  pd = data.frame()
  for (n in names(p)) {
    pd = rbind(pd,data.frame(param=n, lb=NA , ub=NA))
  }
  rownames(pd) <- names(p)

  pd[pd$param=='rho','lb'] = 0
  pd[pd$param=='rho','ub'] = 1

  return(pd)
}

prepare.mopt_config <- function(cf) {
  if (cf$mode=='mpi') {
    cat('[mode=mpi] USING MPI !!!!! \n')
    require(snow)  
    # creating the cluster
    cl <- makeCluster()
    # setting up the slaves
    eval(parse(text = paste("clusterEvalQ(cl,setwd('",cf$wd,"'))",sep='',collapse=''))) 
    eval(parse(text = paste("clusterEvalQ(cl,source('",cf$source_on_nodes,"'))",sep='',collapse=''))) 

    # adding the normal lapply
    cf$mylapply = function(a,b) { return(parLapply(cl,a,b))}

    # adding the load balanced lapply
    cf$mylbapply = function(a,b) { return(clusterApplyLB(cl,a,b))}

    cf$N = length(cl)
  } else if (cf$mode=='multicore') {
    cat('[mode=mulicore] YEAH !!!!! \n')    
    require(parallel)
    cf$mylapply  = mclapply;
    cf$mylbapply = mclapply;
    cf$N=detectCores()
  } else {
    cf$mode = 'serial'
    cat('[mode=serial] NOT USING MPI !!!!! \n')    
    cf$mylapply  = lapply;
    cf$mylbapply = lapply;
    cf$N=3
  }
  return(cf)
}


# this is the main function, it will run the
# optimizer in parallel if called with MPI=TRUE
# it will search for the minimum
runMOpt <- function(cf,autoload=TRUE) {

  # reading configuration
  # =====================
  pdesc = cf$pdesc
  p     = cf$initial_value

  # setting up the cluster if MPI
  # =============================
  last_time = as.numeric(proc.time()[3])

  # saving all evaluations with param values
  # start from best last value
  if (file.exists(cf$file_chain) & autoload ) {
    load(cf$file_chain)
    cf$run = max(param_data$run) +1

  } else {
    param_data = data.frame()
    cf$run = 0
  }

  cat('Number of nodes: ',cf$N,'\n')

  # what is the meaning of these paramters?
  cf$theta  = seq(1,l=50)
  cf$breaks = seq(-2,0,l=50) 
  cf$acc    = 0.5
  cf$nu     = seq(0,l=50)
  cf$kk     = 2

  # we want to keep a matrix for the chain
  # it should be a data.frame in long format
  # with a chain and iteration indicator
  # each row will also have the all parameters and the 
  # objective value
  # we are going to carry 2 arrays
  #   -- param_data: will store all evaluations
  #   -- chain_data: will store the chain results

  #  ========================   MCMC loop =======================

  # get initial candidates
  cat('Computing intial candidates\n')
  ps = computeInitialCandidates(cf$N,cf)

  cat('Starting main MCMC loop\n')  
  for (i in 1:cf$iter) {

    cf$i=i
    # step 1, evaluate candidates 
    # ----------------------------------------------------------------------
    rd = evaluateParameters(ps,cf)

    if ( (i %% cf$save_freq)==1 & i>10 ) {
      save(param_data,file='evaluations.dat')
      #hist(param_data$value,50)
      #plot(cf$theta)
    }

    # step 3, updating chain and computing guesses
    rr = computeCandidatesBGP(param_data,rd, cf, cf$N, i)
    rr$ndt$i=i
    param_data = rbind(param_data,rr$ndt) # append to previous values
    ps = rr$nps
    cf = rr$cf
    cat('[',i,'/', cf$iter ,'][', as.numeric(proc.time()[3]), '][', as.numeric(proc.time()[3]) - last_time  , '] best value is ' , min(param_data$value,na.rm=TRUE) ,'\n') 
    last_time = as.numeric(proc.time()[3])
    for (pp in cf$params_to_sample) {
      cat(' range for ',pp,' ',range(param_data[,pp]),'\n')
    }
  }

  #saving the data set
  save(param_data,file='evaluations.dat')  

  # stopping the cluster using R snow command
  if (cf$USE_MPI) {
    stopCluster(cl)
  }
}


computeInitialCandidates <- function(N,cf) {
  ps = list()
  # generate N guesses
  for (j in 1:N) {
    np = cf$initial_value 

    # pick some parameters to update 
    param <- sample(cf$params_to_sample, cf$np_shock)

    for (pp in param) {
       # update value
       np[[pp]] = shockp(pp, np[[pp]] , cf$shock_var , cf)
    }

    np$chain = j
    ps[[j]] <- np
  }  
  return(ps)
}  


computeCandidates <- function(param_data,rd,cf,N,niter) {

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


shockp <- function(name,value,shocksd,cf) {
  sh = rnorm(1,0,shocksd)
  # update value
  return(fitMirror( value * (1 + sh/100) ,
                              LB = cf$pdesc[name,'lb'],
                              UB = cf$pdesc[name,'ub']))
}

shockallp <- function(p,shocksd,VV,cf) {
  sh = rmultnorm(1,rep(0,nrow(VV)),VV) * shocksd

  # update value for each param
  for (pp in colnames(sh)) {
    p[[pp]] = fitMirror( p[[pp]] + sh[,pp] ,
                              LB = cf$pdesc[pp,'lb'],
                              UB = cf$pdesc[pp,'ub'])  
  }

  return(p)
}

computeCandidatesMH <- function(param_data,rd,cf,N,niter) {

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
    cat(' value and ratio:', epdf(val_old$value),' / ',epdf(val_new$value),'  ',val_old$value, '/' ,val_new$value, ' A=', ACC, '-',prob,' rate=',cf$acc ,'\n')
    
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


#' will vary each parameter independently from the 
#' starting value and store all values. By default
#' it will run for all parameters, otherwise it uses
#' the argument list
compute.slices <- function(mcf,ns=30,pad=0.2) {

  # reading configuration
  # =====================
  pdesc = mcf$pdesc
  p     = mcf$initial_value

  # storing the time
  # =============================
  last_time = as.numeric(proc.time()[3])

  p2 = mcf$initial_value #predict(mcf,'p.all')
  maxres =  MOPT_OBJ_FUNC(p2)

  cat('Number of nodes: ',mcf$N,'\n')


  rr = data.frame() 
  pp = mcf$params_to_sample[1]
  nn=1

  for (pp in mcf$params_to_sample) {

    # we create a range, let's span the entire domain
    lb = mcf$pdesc[pp,'lb']
    ub = mcf$pdesc[pp,'ub']
    prange = seq( lb+(ub-lb)*pad/2,  lb+(ub-lb)*(1-pad/2) ,l=ns)

    # we put the values in a list of parameters
    ptmp =p2
    ps = list()
    j=0
    for (val in prange) {    
      j=j+1
      ptmp[[pp]] = val
      ptmp$param_value = val
      ptmp$chain=j
      ps[[j]] = ptmp
    }

    cat('sending evaluations for ',pp,' in (', lb+(ub-lb)*pad/2,',',lb+(ub-lb)*(1-pad/2),')\n')
    rs = mcf$mylbapply(ps,mopt_obj_wrapper)

    rr1 <- evaluateParameters(ps,mcf)

    cat('got ', nrow(rr1), ' values\n')

    if (nrow(rr1)>0) {
      rr1$param = pp
      rr = rbind(rr,rr1)
    }
    cat('done with ',pp,'(',nn,'/',length(mcf$param_to_sample),')\n')
    nn = nn+1

  }

save(rr,file='slices.dat')

res = list()
res$p.start = p2
res$v.start = maxres
res$slices  = rr
return(res)

}

