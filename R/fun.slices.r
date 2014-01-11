
#' Compute the objective function on a grid of params and show simulated moments
#' 
#' Will vary each parameter independently in a chosen range and report the value 
#' of the resulting simulated moments in relation to the moments in the data.
#' Can be used to construct a heuristic identification argument. Basically it can 
#' be seen which parameter affects which dimension of the model output, i.e. which 
#' simulated moment.
#' @param mcf object of class mopt
#' @param ns number of points in each dimension to evaluate
#' @param pad from bounds of parameter ranges. e.g. p in [0,1], avoid 0 and 1 with pad>0.
#' @param file \code{/path/to/your/file}
#' @return list with info and a data.frame \code{slices} summarizing all information of the exercise: parameter values, 
#' simulated moments, data moments. Input to \code{\link{plot.slices(slices)}}.
#' @export
#' @example examples/example-slices.r
compute.slices <- function(mcf,ns=30,pad=0.1,file="est.slices.RData") {

  # reading configuration
  # =====================
  pdesc = mcf$pdesc
  p     = mcf$initial_value

  # storing the time
  # =============================
  last_time = as.numeric(proc.time()[3])

  p2 = mcf$initial_value #predict(mcf,'p.all')

  cat('evaluate objective function at starting point\n')
  maxres =  MOPT_OBJ_FUNC(p2)
  cat(sprintf('%d nodes, %d parameters, %d points per grid \n',mcf$N,length(mcf$params_to_sample),ns))

  rr = data.frame() 
  pp = mcf$params_to_sample[1]
  nn=1

  for (pp in mcf$params_to_sample) {

    # we create a range, let's span the entire domain
    lb = mcf$pdesc[pp,'lb']
    ub = mcf$pdesc[pp,'ub']
    prange = seq( lb+(ub-lb)*pad/2,  lb+(ub-lb)*(1-pad/2) ,l=ns)

    # we put the values in a list of parameters
	# TODO watch out as some of the cluster functions take vectors
	# while others take lists.
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
    
	if (mcf$mode =='mpi'){
		rs = clusterApplyLB(cl=mcf$cl,x=ps,fun=mopt_obj_wrapper,objfunc=mcf$objfunc)
	} else {
		rs = mcf$mylbapply(ps,mopt_obj_wrapper,objfunc=mcf$objfunc)
	}

    rr1 = data.frame()
    for ( jj in 1:length(rs) ) {
      if (is.atomic(rs[[jj]])) next;
      if (is.null(rs[[jj]]))   next;
      if (rs[[jj]]$status==-1) next;

      sm = rs[[jj]]$sm
      sm$names = paste('sm',sm$names,sep='.')
      sm = daply(sm, .(names), function(d) {d$value[1]  })

      # get everything but the submoments
      rr1.tmp = data.frame( parvalue=prange[jj], 
                            t(unlist(ps[[jj]])) , 
                            t(unlist( rs[[jj]][setdiff(names(rs[[jj]]),c('sm','p','chain'))  ] ) ),
                            t(sm))
      rr1 = rbind(rr1, rr1.tmp);
    }

    cat('got ', nrow(rr1), ' values\n')

    if (nrow(rr1)>0) {
      rr1$param = pp
      rr = rbind(rr,rr1)
    }
    cat('done with ',pp,'(',nn,'/',length(mcf$params_to_sample),')\n')
    nn = nn+1

  }

  save(rr,mcf,file=file)
  res = list()
  res$p.start = p2
  res$v.start = maxres
  res$slices  = rr
  return(res)
}

#' plot.slices
#'
#' generates plots for each moments/parameter combinations
#' depends on output of function compute.slices(), which saves
#' a dataframe with plotting data and the used mopt config into a 
#' a file \code{est.slices.RData}
#' @param file path/to/est.slices.RData if not in getwd()
#' @param outpath path to directory where to save plots
#' @param type string indicating file type to save plot as. currently png and pdf only.
#' @export
plot.slices <- function(file=NULL,outpath='',type="png") {
  # we want to keep all submoments, value, param and param_value

	if (is.null(file)) {
		load('est.slices.RData')
	} else {
		load(file)
	}
  
  rr$conv=as.numeric(rr$status)>0
  
  rr.m          = melt(rr,id=c('param','param_value','conv'))
  rr.m          = subset(rr.m,rr.m$variable=='value' | str_detect(rr.m$variable,'sm'))
  rr.m$variable = gsub('sm.','',rr.m$variable)
  rr.m$from     = 'model'
  rr.m$value    = as.numeric(rr.m$value)

  rr.m=data.table(rr.m)

  # check if we have got the right number of moments in mcf as
  # in rr.m
  
  # this data.frame holds the initial values
  # of the parameters
  init.param.data = data.frame(value = unlist(mcf$initial_value) , param=names(mcf$initial_value))

  # we subset it to the ones we sampled over
  init.param.data = subset(init.param.data,param %in% unique(rr.m$param))

  for (pp in unique(rr.m$variable)) {
    if (pp == 'value') next;

    gp <- ggplot(subset(rr.m,variable==pp & from=='model')) + geom_point(aes(x=param_value,y=value,color=conv),size=1) + 
      geom_line(aes(x=param_value,y=value,group='param'),size=0.3) + 
      geom_hline(aes(yintercept=value),data=subset(mcf$data.moments,moment==pp),linetype=2) +
      geom_vline(aes(xintercept=value),data=init.param.data,linetype=2,color='red') +
      facet_wrap(~param,scales='free_x',ncol=3) +
      scale_y_continuous(paste('value of',pp))+ scale_x_continuous('varying parameter') + theme_bw()
    #print(gp)
    if (type=="png") ggsave(paste(outpath,'plot_ParamVsMoment_',pp,'.png',sep=''),width=10.6, height=5.93)
    if (type=="pdf") ggsave(paste(outpath,'plot_ParamVsMoment_',pp,'.pdf',sep=''),width=10.6, height=5.93)
  }

  pp ='value'
  gp <- ggplot(subset(rr.m,variable==pp & from=='model')) + geom_point(aes(x=param_value,y=value,color=conv),size=1) + 
      geom_line(aes(x=param_value,y=value,group='param'),size=0.3) + 
      geom_vline(aes(xintercept=value),data=params.data,linetype=2,color='red') +
      facet_wrap(~param,scales='free_x',ncol=3) +
      scale_y_continuous('objective function')+ scale_x_continuous('varying parameter') + theme_bw()
    #print(gp)
    if (type=="png") ggsave(paste(outpath,'plot_ParamVsMoment_',pp,'.png',sep=''),width=10.6, height=5.93)
    if (type=="pdf") ggsave(paste(outpath,'plot_ParamVsMoment_',pp,'.pdf',sep=''),width=10.6, height=5.93)
 }



#' Compute the objective function on a grid of params and show custom model output
#' 
#' Essentially the same as \code{\link{compute.slices}}, but does not report simulated 
#' moments but other model output. Useful for model output that is multidimensional. 
#' It's a simplified version of \code{\link{compute.slices}} in that it does not further
#' process the model output: it return a list with nests "parameter name", "value of parameter",
#' "model output".
#' For example instead of reporting the mean of a certain statistic, this function can
#' return a matrix or a higher dimensional array. Say you want to return the life-cycle
#' profile of a certain model variable x. This will be a vector of length N, where N is 
#' the number of periods in the model. The user has to design the MOPT_OBJ_FUN in such a way
#' that it returns the required output. There are 2 requirements for what \code{MOPT_OBJ_FUN} has to return.
#' First it has to be a list, second, the list needs components "status" (indicating whether a particular evaluation
#' is valid in some sense) and "output", which contains your custom model output.
#' @param mcf object of class mopt
#' @param ns number of points in each dimension to evaluate
#' @param pad from bounds of parameter ranges. e.g. p in [0,1], avoid 0 and 1 with pad>0.
#' @param file \code{/path/to/your/file.RData}
#' @return list by paramter name, parameter value index, containing the value of the parameter vector and a list \code{data} containing
#' your custom model output.
#' @export
#' @example examples/example-slices2.r
compute.slices2 <- function(mcf,ns=30,pad=0.1,file="est.slices.RData") {

  # reading configuration
  # =====================
  pdesc = mcf$pdesc
  p     = mcf$initial_value

  # storing the time
  # =============================
  last_time = as.numeric(proc.time()[3])

  p2 = mcf$initial_value #predict(mcf,'p.all')

  cat('evaluate objective function at starting point\n')
  maxres =  MOPT_OBJ_FUNC(p2)
  cat(sprintf('%d nodes, %d parameters, %d points per grid \n',mcf$N,length(mcf$params_to_sample),ns))

  rr = data.frame() 
  pp = mcf$params_to_sample[1]
  nn=1

  out <- vector("list",length(mcf$params_to_sample))
  names(out) <- mcf$params_to_sample

  for (pp in mcf$params_to_sample) {

    # we create a range, let's span the entire domain
    lb = mcf$pdesc[pp,'lb']
    ub = mcf$pdesc[pp,'ub']
    prange = seq( lb+(ub-lb)*pad/2,  lb+(ub-lb)*(1-pad/2) ,l=ns)

    # we put the values in a list of parameters
	# TODO watch out as some of the cluster functions take vectors
	# while others take lists.
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
    
	if (mcf$mode =='mpi'){
		rs = clusterApplyLB(cl=mcf$cl,x=ps,fun=mopt_obj_wrapper_custom,objfunc=mcf$objfunc)
	} else {
		rs = mcf$mylbapply(ps,mopt_obj_wrapper_custom,objfunc=mcf$objfunc)
	}

	# get model output for each parameter value
	out[[pp]] <- list()
    for ( jj in 1:length(rs) ) {
      if (is.null(rs[[jj]]))   next;
      if (rs[[jj]]$status==-1) next;
	  out[[pp]][[jj]] <- list()
	  out[[pp]][[jj]]$pars <- ps[[jj]]
	  out[[pp]][[jj]]$data <- rs[[jj]]$output
    }

    cat('got ', length(out[[pp]]), ' values\n')

    cat('done with ',pp,'(',nn,'/',length(mcf$params_to_sample),')\n')
    nn = nn+1

  }

  save(out,mcf,file=file)
  return(out)
}
