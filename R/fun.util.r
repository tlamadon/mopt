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

checkEvalStructure <- function(val,cf) {
  if (!("p" %in% names(val))) {
    cat('error, p is not in the return value of the objective function\n')
  } else {
    missing = setdiff(names(cf$initial_value),names(val$p))
    if (length(missing)>0) cat(paste(missing, collapse=', '), ' are missing from p in return value of objective function\n')
  }
  missing = setdiff(c('time','value','chain'),names(val))
  if (length(missing)>0) cat(paste(missing, collapse=', '), ' are missing from return value of objective function\n')
}

#' @export
evaluateParameters <- function(ps,cf,balance=FALSE) {
    
    #save evaluations to file
    #save(ps,file='lasteval.dat')

    cat('Sending parameter evaluations...\n')
	if (cf$mode=='mpi'){
		vals <- parLapply(cf$cl,ps,mopt_obj_wrapper,objfunc = cf$objfunc)
	} else if (balance) {
		vals = cf$mylbapply(ps,mopt_obj_wrapper,objfunc = cf$objfunc)
	} else {
	    vals = cf$mylapply(ps,mopt_obj_wrapper,objfunc = cf$objfunc)
	}
    cat('done\n')

    # process the return values 1 by 1
    rr = data.frame()
    for (val in vals) {

      # check the status
      rd = data.frame()
      if (!('status' %in% names(val))) next;
      if (val$status<0) {
        print(paste('error: ',val$error,'\n'))
        next
      }
      checkEvalStructure(val,cf)

      # get param value back, and the chain
      pt = val$p
      pt$chain = NULL
      names(pt) <- paste('p',names(pt),sep='.')
      rd  = data.frame(pt)

      rd$value  = val$value
      rd$status = val$status
      rd$time   = val$time
      rd$chain  = val$chain

      # collect the addititonal infos
      if (length(val$infos)>0) {
        rd.infos = data.frame(val$infos)
        colnames(rd.infos) <- paste('d',colnames(rd.infos),sep='.')
        rd = cbind(rd,rd.infos)
      }

      # get the moments
      if (!is.null(val$sm) & nrow(val$sm)>0) {
        rd.sm = val$sm$value
        names(rd.sm) = paste('m',val$sm$names,sep='.')
        rd.sm = data.frame(as.list(rd.sm))
        rd    = cbind(rd,rd.sm)
      }
      rr = rbind(rr,rd)
    }
    cat(sprintf('evaluations are finished [%d/%d]\n',nrow(rr),length(ps)))

    # return the evaluations
    rr$i     = cf$i
    rr$run   = cf$run
    return(rr)
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

#' @export
shockp <- function(name,value,shocksd,cf) {
  sh = rnorm(1,0,shocksd)
  # update value
  return(fitMirror( value * (1 + sh/100) ,
                              LB = cf$pdesc[name,'lb'],
                              UB = cf$pdesc[name,'ub']))
}

#' @export
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
