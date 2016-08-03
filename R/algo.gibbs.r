#' Metropolis Hasting MCMC chain updating
#' @export
#' @family algos
algo.gibbs <- function(rd,param_data,niter,cf,pdesc,priv) {

  # at each step, we compute many values for one given 
  # parameter, we then draw one of this values according
  # to the pseudo likelihood
  ps=list()
  params_to_sample  = cf$params_to_sample
  params_to_sample2 = paste('p',cf$params_to_sample,sep='.')
  N                 = cf$N
  
  rr = data.frame()
  # if param_data is empty then we accept all
  if (nrow(param_data)==0) {
    cf$tempering = 1
    cf$current_param = "none"
    priv = cf
    priv$N=1
    return(list(ps = computeInitialCandidates(N,cf) , evals=rd,priv=priv ))
  }

  # we draw the new parameter value according to the posterior probabilities
  # that we just computed
  new_param_i = sample.int(nrow(rd),1,prob=exp(-rd$value))
  best = rd[rd$chain==new_param_i,]
  best$chain=0
  rr = rbind(rr,best)

  # we display the evaluations
  if ("current_param"%in%names(priv)) {
    sel = rep(" ",N)
    sel[new_param_i]="*"
    flog.info("evaluations for %s",priv$current_param)
    v = sprintf("%+3.3f%s",as.numeric(rd[paste("p.",priv$current_param,sep="")][,1]),sel)
    flog.info("x= %s", paste(v,collapse="  "))
    v = sprintf("%+3.3f%s",as.numeric(rd["value"][,1]),sel)
    flog.info("f= %s",paste(v,collapse="  "))
    v = sprintf("%+3.3f%s",exp(-rd$value)/sum(exp(-rd$value)),sel)
    flog.info("p= %s",paste(v,collapse="  "))
  }
  
  # we randomly select a parameter
  param = sample(cf$params_to_sample, 1)
  grid = sort(cf$pdesc[param,'lb'] + (cf$pdesc[param,'ub']-cf$pdesc[param,'lb'])* ((seq(1/N,1,l=N) + rnorm(1)) %% 1))
  priv$current_param = param
  flog.info("nest param to evaluate: %s",param)

  # the next evaluations are just a uniform sequence from the bounds
  val_new = rd[rd$chain==new_param_i,]
  for (c.current in 1:N) {
    # store the realisation
    rr = rbind(rr,rd[rd$chain==c.current,])
    
    # compute new parameter
    p.new = cf$initial_value
    p.new$chain = c.current
    p.new[params_to_sample] = val_new[1,params_to_sample2]    
    p.new[[param]]          = grid[c.current]
    ps[[c.current]]         = p.new
  }
  
  v = sprintf("%+3.3f",grid)
  flog.info("support= %s",paste(v,collapse="  "))

   return(list(ps = ps , evals = rr, priv=priv ))
}