

#' compute the simulated score from the 
#' the chain. This is the matrix of the derivative of the 
#' moments with respect to the parameters around the best parameter
#' value.
#' @export
score.mopt_env <- function(me) { 


  mm =  melt(me,'mp')

  # use the best value
  ip0 = which.min(mm$objvalue)

  # compute the weight as the square distance from p0
  p0 <- as.list(mm[ip0,me$cf$params_to_sample,with=FALSE])
  dp0 = data.frame(i = 1:nrow(mm),p0)
  dp0$i <- NULL
  mm$weight = apply(  (as.matrix(mm[,me$cf$params_to_sample,with=FALSE]) -    as.matrix(dp0))^2 / (as.matrix(dp0))^2 , 1,mean)


  fmula = paste(me$cf$params_to_sample,collapse="+")
  M = matrix(0,length(me$cf$params_to_sample) , length(  me$cf$moments_to_use))

  rsq = 1:length(me$cf$moments_to_use)
  for (i in  1:length(me$cf$moments_to_use)) {
    ss     = summary(lm(formula(paste( me$cf$moments_to_use[[i]] , " ~ ",fmula)),mm,weights=1/(1+10000*mm$weight)))
    M[,i]  = coef(ss)[me$cf$params_to_sample,"Estimate"]
    rsq[i] = ss$r.squared 
  }

  rownames(M) <- me$cf$params_to_sample
  colnames(M) <- me$cf$moments_to_use

  return(M)
}


#' takes the best value, compute the score and
#' the standard error using the sandwich formula
#' @export
estimate.mopt_env <- function(me) {

  # get moments and compute diagonal matrix
  mms = as.list(me$cf$data.moments$sd)
  names(mms) = me$cf$data.moments$moment

  S  = score.mopt_env(me)

  VV = S %*% (diag( 1/as.numeric(mms[  me$cf$moments_to_use ] ))/1000^2) %*% t(S)

  mm = melt(me,'p')
  JJ = var(mm[temp==1,me$cf$params_to_sample,with=FALSE])

  SE = JJ %*% VV %*% JJ

  return(SE)

}
