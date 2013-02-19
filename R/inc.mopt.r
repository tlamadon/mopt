require(plyr)
require(MSBVAR)

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
#' @export
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
                                #This is a vector of names: c('delta', 'b')
 cf$run            = 0
 cf$shock_var      = 0.1
 cf$moments_to_use = c()
 cf$moments.data   = c()
 cf$moments.sd     = c()
 cf$np_shock       = 1
 cf$save_freq      = 25
 cf$initial_value = p	

 param.descript = data.frame()
 for (n in names(p)) {
   param.descript = rbind(param.descript,data.frame(param=n, lb=NA , ub=NA))
 }
 rownames(param.descript) <- names(p)
 cf$pdesc = param.descript

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

#' @export
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
        #  Vals is a list (later data.frame) of returned values and moments from each chain. If one chain returns NA, then ICONV deletes the whole line. However, the program crashes if all chains return NA value.   
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
  param.descript = data.frame()
  for (n in names(p)) {
    param.descript = rbind(param.descript,data.frame(param=n, lb=NA , ub=NA))
  }
  rownames(param.descript) <- names(p)

  param.descript[param.descript$param=='rho','lb'] = 0
  param.descript[param.descript$param=='rho','ub'] = 1

  return(param.descript)
}

#' @export
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
    
    if(Sys.info()[['sysname']]=='Windows') {
      cl <- makeCluster(detectCores())
      
      # setting up the slaves
      eval(parse(text = paste("clusterEvalQ(cl,setwd('",cf$wd,"'))",sep='',collapse=''))) 
      eval(parse(text = paste("clusterEvalQ(cl,source('",cf$source_on_nodes,"'))",sep='',collapse=''))) 
      
      # adding the normal lapply
      cf$mylapply = function(a,b) { return(parLapply(cl,a,b))}
      
      # adding the load balanced lapply
      cf$mylbapply = function(a,b) { return(clusterApplyLB(cl,a,b))}
      
      cf$N = length(cl)
    } else{
      cf$mylapply  = mclapply;
      cf$mylbapply = mclapply;
      cf$N=detectCores()      
    }
  } else {
    cf$mode = 'serial'
    cat('[mode=serial] NOT USING MPI !!!!! \n')    
    cf$mylapply  = lapply;
    cf$mylbapply = lapply;
    cf$N=3
  }
  return(cf)
}

#' this is the main function, it will run the
#' optimizer in parallel if called with MPI=TRUE
#' it will search for the minimum
#' @export
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
  # I find theta, breaks, nu, and kk are used algo.wl.r, acc is used to update sample var.
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

shockp <- function(name,value,shocksd,cf) {
  sh = rnorm(1,0,shocksd)
  # update value
  return(fitMirror( value * (1 + sh/100) ,
                              LB = cf$pdesc[name,'lb'],
                              UB = cf$pdesc[name,'ub']))
}

jumpParams.normalAndMirrored <- function(p,shocksd,VV,params.desc) {
  sh = rmultnorm(1,rep(0,nrow(VV)),VV) * shocksd

  # update value for each param
  for (pp in colnames(sh)) {
    p[[pp]] = fitMirror( p[[pp]] + sh[,pp] ,
                              LB = params.desc[pp,'lb'],
                              UB = params.desc[pp,'ub'])  
  }

  return(p)
}

list2df <- function(ll) {
 return(ldply(ll,function(l){ return(data.frame(rbind(unlist(l))))}))
}





