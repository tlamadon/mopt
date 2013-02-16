#' Compute the objective function on a grid
#' 
#' Will vary each parameter independently from the 
#' starting value and store all values. By default
#' it will run for all parameters, otherwise it uses
#' the argument list
#' @export
#' @family algos
algo.slices <- function(mcf,ns=30,pad=0.2) {

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