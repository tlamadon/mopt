#' Simple aglorithm which compute a slice, then takes the minimum
#' @export
#' @family algos
algo.slicemax <- function(rd,param_data,niter,cf,pdesc,priv) {

  # at each step, we compute many values for one given 
  # parameter, we then draw one of this values according
  # to the pseudo likelihood
  ps=list()
  params_to_sample  = cf$params_to_sample
  params_to_sample2 = paste('p',cf$params_to_sample,sep='.')
  N                 = cf$N
  rd$run=cf$i
  
  rr = data.frame()
  # if param_data is empty then we accept all
  if (nrow(param_data)==0) {
    cf$tempering = 1
    cf$current_param = "none"
    priv = cf
    priv$N=1
    return(list(ps = computeInitialCandidates(N,cf) , evals=rd,priv=priv ))
  }

  if ("current_param" %in% names(priv)) {
    curp = priv$current_param
    curp.list = priv$curp.list
  } else {
    curp="none"
    curp.list = rep(0,length(cf$params_to_sample))
    names(curp.list) = cf$params_to_sample
  }  

  if ("current_best_par" %in% names(priv)) {
    current_best_par = priv$current_best_par
    current_best_val = priv$current_best_val
  } else {
    current_best_val = Inf
    current_best_par = rd[1,]
  }  
  
  # select the best new value
  new_param_i = rd$chain[which.min(rd$value)]
  # compare this to the previous evaluation
  if (rd$value[rd$chain==new_param_i] < current_best_val) {
    flog.info("best new evaluation is better, we use it as new center (%f > %f)",current_best_val,rd$value[rd$chain==new_param_i])
    current_best_par = rd[rd$chain==new_param_i,]
    current_best_val = current_best_par$value
  } else {
    flog.info("best new evaluation is not better, we keep old one")
  }
  
  best = current_best_par
  best$chain=0
  rr = rbind(rd,best)

  # we display the evaluations
  if ("current_param"%in%names(priv)) {
    sel = rep(" ",N)
    sel[new_param_i]="*"
    flog.info("evaluations for %s",priv$current_param)
    v = rep(0,N)
    v[rd$chain] = as.numeric(rd[paste("p.",priv$current_param,sep="")][,1])
    flog.info("x= %s", paste(sprintf("%+3.3f%s",v,sel),collapse="  "))
    
    v[rd$chain] = as.numeric(rd["value"][,1])
    flog.info("f= %s", paste(sprintf("%+3.3f%s",v,sel),collapse="  "))
  }
  
  # we randomly select a parameter from the least evaluates once
  I = which(curp.list <= min(curp.list)+1)
  param = sample(names(curp.list)[I], 1)
  grid = sort(cf$pdesc[param,'lb'] + (cf$pdesc[param,'ub']-cf$pdesc[param,'lb'])* ((seq(1/N,1,l=N) + rnorm(1)) %% 1))

  priv$current_param = param
  priv$current_best_par = current_best_par
  priv$current_best_val = current_best_val
  curp.list[[param]] = curp.list[[param]] + 1
  priv$curp.list = curp.list
  
  flog.info("next param to evaluate: %s",param)

  # the next evaluations are just a uniform sequence from the bounds
  val_new = current_best_par
  for (c.current in 1:N) {
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