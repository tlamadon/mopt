#' Compute the objective function on a grid
#' 
#' Will vary each parameter independently from the 
#' starting value and store all values. By default
#' it will run for all parameters, otherwise it uses
#' the argument list
#' @export
#' @example examples/example-slices.r
compute.slices <- function(mcf,ns=30,pad=0.1) {

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

  save(rr,file='est.slices.dat')
  res = list()
  res$p.start = p2
  res$v.start = maxres
  res$slices  = rr
  return(res)
}

#' plot.slices
#'
#' generates plots for each moments/parameter combinations
#' @param path path where to put the plots
#' @param p list of ALL parameters
#' @param mcf mopt configuration of a model
#' @export
plot.slices <- function(p,mcf,path='') {
  # we want to keep all submoments, value, param and param_value
  if (!file.exists('est.slices.dat')){
	  cat('cannot find file est.slices.dat. you must call compute.slices first')
	  return(NULL)
  }
  load('est.slices.dat')
  nn = names(rr)
  rr$conv=rr$status>0
  
  rr.m = melt(rr,id=c('param','param_value','conv'))
  rr.m = subset(rr.m,rr.m$variable=='value' | str_detect(rr.m$variable,'sm'))
  rr.m$variable = gsub('sm.','',rr.m$variable)
  rr.m$from = 'model'

  rr.m=data.table(rr.m)
  params.data = subset(data.frame(value = unlist(p) , param=names(p)),param %in% unique(rr.m$param))
  for (pp in unique(rr.m$variable)) {
    if (pp == 'value') next;

    gp <- ggplot(subset(rr.m,variable==pp & from=='model')) + geom_point(aes(x=param_value,y=value,color=conv),size=1) + 
      geom_line(aes(x=param_value,y=value),size=0.3) + 
      geom_hline(aes(yintercept=value),data=subset(mcf$data.moments,moment==pp),linetype=2) +
      geom_vline(aes(xintercept=value),data=params.data,linetype=2,color='red') +
      facet_wrap(~param,scales='free_x',ncol=3) +
      scale_y_continuous(paste('value of',pp))+ scale_x_continuous('varying parameter') + theme_bw()
    #print(gp)
    ggsave(paste(path,'plot_ParamVsMoment_',pp,'.png',sep=''),width=10.6, height=5.93)
  }

  pp ='value'
  gp <- ggplot(subset(rr.m,variable==pp & from=='model')) + geom_point(aes(x=param_value,y=value,color=conv),size=1) + 
      geom_line(aes(x=param_value,y=value),size=0.3) + 
      geom_vline(aes(xintercept=value),data=params.data,linetype=2,color='red') +
      facet_wrap(~param,scales='free_x',ncol=3) +
      scale_y_continuous('objective function')+ scale_x_continuous('varying parameter') + theme_bw()
    #print(gp)
    ggsave(paste(path,'plot_ParamVsMoment_',pp,'.png',sep=''),width=10.6, height=5.93)
 }
