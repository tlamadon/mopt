#' generates different plots from a chain
#' @export
plot.mopt_env <- function(me,what='na',select='untempered',taildata=0) {

  mopt = me$cf

  desc = list(
    'pdistr'='parameter posterior distribution',
    'ptime' ='parameters over time',
    'mdistr'='moment disitribution',
    'vtime' ='objective value over time',
    'mtime' ='moment over time',
    'mtable' = 'moment table',
    'pmv'   = 'parameter mean/variance per chain')
  graphs = names(desc)

  if (length(setdiff(what,graphs))>0) {
    cat("outputs available : \n")
    for (g in graphs) {
      cat(sprintf("%6.6s: %s \n",g,desc[[g]]))
    }
    return()
  }

  param_data = data.table(me$param_data)

  if (select=='untempered') {
    param_data = subset(param_data,chain %in% which(me$priv$chain.states$tempering==1))
  }

  mnames = grep('submoments.',names(param_data),value=TRUE)
  datamnames = grep('submoments.data',names(param_data),value=TRUE)
  mnames = setdiff(mnames, datamnames)

  if ( taildata>0 )  param_data = tail(param_data,taildata);

  # reporting the distributions of the parameters
  if ('pdistr' %in% what) {
    
    # we also want to append the best param value
    I = which(param_data$value == min(param_data$value,na.rm=TRUE))
    best_data = param_data[I,]

    gp <- ggplot(melt(param_data,measure.vars=paste('p',mopt$params_to_sample,sep='.'),id=c()),aes(x=value)) +
        geom_vline(aes(xintercept=value),data = melt(best_data,measure.vars=mopt$params_to_sample,id=c()),color='blue',linetype=2,size=0.3) +
        stat_density(fill='blue',alpha=0.5) + facet_wrap(~variable,scales='free')
    print(gp)
  }

  # reporting the evolution of the objective value
  if ('vtime' %in% what) {
    param_data$iter = 1:nrow(param_data)
    gp <- ggplot(param_data,aes(x=iter,y=value,color=chain)) + geom_line()
    print(gp)
  }

  # graph all moments
  if ('mdistr' %in% what) {
    quartz()
    gdata = melt(rename(param_data,c(value='objvalue')),measure.vars=mnames,id=c('objvalue'))
    gdata$variable = gsub('submoments.','',gdata$variable)
    save(gdata,file='tmp.dat')
    gp <- ggplot(gdata,aes(x=value)) + 
          geom_density(fill='blue',alpha=0.4,size=0) + 
          geom_vline(aes(xintercept=value),
                     data=data.frame(variable=mopt$data.moments$moment,
                                     value   =mopt$data.moments$value),
                     linetype=2,color='red')+
          facet_wrap(~variable,scales='free')
    print(gp)
  }

  # graph link between parameters and moments
  if ('pmpoints' %in% what) {
    quartz()
    gdata = melt(rename(param_data,c(value='objvalue')),measure.vars=mopt$params_all,id=c('objvalue'))
    gp <- ggplot(gdata,aes(x=value,y=objvalue)) + geom_point() + facet_wrap(~variable,scales='free')
    print(gp)
  }

  # link between submoments and parameters
  if ('pmreg' %in% what) {
    RHS = paste('log(',mopt$params_to_sample,')',sep= '',collapse=' + ')
    rr = data.frame()
    # we also want to weight the observations by how close they are to the optimal
    param_data$lmw = 1/(0.01+(param_data$value - min(param_data$value)))

    for ( ms in c(mnames,'log(value)') ) {
      fit = lm(paste(ms,'~',RHS),param_data,weights=lmw)
      r = model2frame(fit)
      r$dep = ms
      rr = rbind(rr,r)
    }
    rr = subset(rr,!variable %in% c('BIC','(Intercept)'))
    rr$dep = gsub('submoments.','',rr$dep)

    ggt <- ggtable(dep ~ variable) + ggt_cell_regression(rr,list(value='value',sd='sd',pval='pval')) +
           ggt_order('variable',mopt$params_to_sample) +
           ggt_order('dep',gsub('submoments.','',mnames))
    ggt$params$resize=0.7
    print(ggt);
  }

  # moment table 
  if ('pmreg' %in% what) {
    RHS = paste('log(',mopt$params_to_sample,')',sep= '',collapse=' + ')
    rr = data.frame()
    # we also want to weight the observations by how close they are to the optimal
    param_data$lmw = 1/(0.01+(param_data$value - min(param_data$value)))

    for ( ms in c(mnames,'log(value)') ) {
      fit = lm(paste(ms,'~',RHS),param_data,weights=lmw)
      r = model2frame(fit)
      r$dep = ms
      rr = rbind(rr,r)
    }
    rr = subset(rr,!variable %in% c('BIC','(Intercept)'))
    rr$dep = gsub('submoments.','',rr$dep)

    ggt <- ggtable(dep ~ variable) + ggt_cell_regression(rr,list(value='value',sd='sd',pval='pval')) +
           ggt_order('variable',mopt$params_to_sample) +
           ggt_order('dep',gsub('submoments.','',mnames))
    ggt$params$resize=0.6
    print(ggt);
  }

  if ('ptime' %in% what) {
    #param_data = data.frame(me$param_data)
    gdata = melt(rename(param_data,c(value='objvalue')),measure.vars=paste('p',me$cf$params_to_sample,sep='.'),id=c('i','chain'))
    gp <- ggplot(gdata,aes(x=i,color=chain,group=chain,y=value)) + 
      geom_point(size=0.5) + 
      geom_line(size=0.5) + 
      facet_wrap(~variable,scales='free') +
      theme_bw()
    print(gp)
  }

  if ('mtime' %in% what) {
    quartz()
    param_data$t = c(1:nrow(param_data))
    gdata = melt(rename(param_data,c(value='objvalue')),measure.vars=mnames,id=c('t','chain'))
    gdata$variable = gsub('submoments.','',gdata$variable)
    gp <- ggplot(gdata,aes(x=t,y=value)) + 
          geom_point(size=0.5) + 
          geom_hline(aes(yintercept=value,x=NULL,y=NULL),
                     data=data.frame(value = mopt$moments.data,
                     variable=names(mopt$moments.data)),
                     linetype=2,color='red')+
          facet_wrap(~variable,scales='free')
    print(gp)
  }

  if ('mtable' %in% what) {
    source('~/git/Utils/R/ggtable2.r')
    # get the optimal value of the submoments
    i = which.min(param_data$value)
    rr = t(param_data[i,])
    rr = data.frame(from='model',moment = rownames(rr),value = c(rr),sd=NA)
    rr = subset(rr,str_detect(rr$moment,'submoments.model'))
    rr$moment = str_replace( rr$moment,'submoments\\.model\\.','')
   
    dd = mopt$data.moments
    dd = subset(dd,moment %in% rr$moment)
    dd$from = 'data'
    rr = rbind(rr,dd)

    rownames(rr) <- NULL
    rr$pval=NA

    save(rr,file='tmp2.dat')
    ggt <- ggtable(moment ~ from) + ggt_cell_regression(rr,list(value='value',sd='sd',pval='pval'))
    print(ggt)
  }

  if ('pmv' %in% what) {
    mm = melt(me)
    gg <- ggplot(mm[,list(mean(value),sd(value)),list(chain,variable)],aes(x=chain,y=V1,ymin=V1-V2,ymax=V1+V2)) + 
            geom_pointrange() + facet_wrap(~variable,scales='free') + theme_bw()
    print(gg)
  }

}

#' print method for mopt_env
#' @export
print.mopt_env <- function(me) {
  cat(sprintf(" %3.d chains / %5.d evaluations / %3.d parameters estimated" , 
          me$param_data[,length(unique(chain))], me$param_data[,max(t)] , length(me$cf$params_to_sample))) 
}


mopt.newconf <- function(name) {

filename_make = 'Makefile'

if (file.exists(cf$file_chain)) {
  cat('Makefile already exists in current folder, saving to MakeMpi.inc\n')
  filename_make = 'MakeMpi.inc'
}

# creating and saving the make file
STR = 'runmpi:
	qsub qsub_start_mpi.sh 

clean:
	rm -rf *.out *.rout

tail:
	tail -f ./mpi_PROJECTNAME.out
'  
STR = gsub('PROJECTNAME',name,STR)
cat(STR,file=filename_make)

}

#' return some versions of the parameters
#' @export
#' @param what can be 'p.all' the best parameter set as a list, 'p.sd' for 
#' sampled parameters with standard deviations based on coldest chain, 'm' for list of 
#' simulated and data moments next to each other
predict.mopt_env <- function(me,what='p.all',base='') {


  # first type, is to return the parameters with the highest value
  cf = me$cf
  param_data = data.frame(me$param_data)
  I = which(param_data$value == min(param_data$value,na.rm=TRUE))[[1]]

  params_to_sample  = cf$params_to_sample
  params_to_sample2 = paste('p',cf$params_to_sample,sep='.')
  param_data = data.table(param_data)
  cat(sprintf("value=%f evals=%i\n",param_data[I[1]]$value, nrow(param_data)))

  if (what=='p.all') {
    pres = cf$initial_value
    pres[params_to_sample] = param_data[I,params_to_sample2,with=FALSE]
    return(pres)
  } 

  if (what=='p.sd') {
    VV = sqrt(diag(cov(param_data[chain==1,params_to_sample2,with=FALSE])))
    p = c(param_data[I,params_to_sample2,with=FALSE])
    return(data.frame(name = cf$params_to_sample, value = unlist(p), sd = c(VV)))
  }

  if (what=='m') {
    mnames      = grep('m\\.',names(param_data),value=TRUE)
    sim.moments = data.frame(model = as.numeric(param_data[I,mnames,with=FALSE]), moment = str_replace(mnames,'m\\.',''))
    sim.moments = merge(sim.moments,cf$data.moments,by='moment')
    return(sim.moments)
  }

  if (what=='coda') {
    mm = melt(me)
    mcs = list()
    for (i in 1:5) {
      mmm = cast(mm[chain==i],t~variable)
      mcs[[i]] = as.mcmc(mmm)
    }
    return(mcmc.list(mcs))
  }

  return(p)
}

#' Load an existing mopt config
#' 
#' loads a mopt config from file. This is very useful
#' for on-the-fly analysis of results that are generated
#' on a remote machine, or to process results on your local machine.
#' Fetching remote files requires a publickey authentication setup
#' for ssh to the remote server. your connetion must
#' work without having to type a password. It's easy
#' to setup the ssh-agent and add a passphrase to an existing
#' key. Preferably set up a config file for ssh. see the references.
#' @references \url{https://help.github.com/articles/working-with-ssh-key-passphrases}
#' \url{http://nerderati.com/2011/03/simplify-your-life-with-an-ssh-config-file/}
#' @param filename optional local filename
#' @param remote optional. full \code{scp} path: username@your.remote.com:~/path/to/remote/file.dat
#' @param reload NULL
#' @export
#' @examples
#' \dontrun{
#' me <- mopt.load(remote="hpc:git/wagebc/Rdata/evaluations.educ1.dat")
#' predict(me,'p.sd')
#' }
mopt.load <- function(filename='',remote='',reload=NULL) {

  if (str_length(remote)>0) {
  # generate a local tmp file
    filename = tempfile()
    # get the file using scp
    system(paste('scp',remote,filename))
  } 

	if (!is.null(reload)) {
    filename = tempfile()
    # get the file using scp
    system(paste('scp',reload$remote,filename))
  } 

  env = new.env()
  load(filename,envir=env)
  class(env) <- 'mopt_env'
  # add time and append accept rejects
  env$param_data = ddply(env$param_data,.(chain),function(d) {d$t = 1:nrow(d);d})
  env$param_data = data.table(env$param_data)
  env$remote = remote
  return(env)
}

#' melts data with only sampled parameters
#' if set cols='m' or cols='p' extract the wide format with 
#' only parameters or only moments
#' @export
melt.mopt_env <- function(me,cols=NULL) {
  param_data = data.frame(me$param_data)
  param_data = data.table(rename(param_data,c(value='objvalue')))

  # compute accept recject
  setkey(param_data,chain,t)
  param_data[,i.l1:=param_data[J(chain,t-1),i][,i]]
  param_data[,acc := i!=i.l1]

  # merge in temperature
  param_data = ddply(param_data,.(chain), function(d) {
    d$temp = me$priv$chain.states$tempering[[d$chain[1]]]
    return(d)
  })

  if (!is.null(cols)) {

    if (cols=='p') {
      gdata = param_data[,c('objvalue','temp','chain','acc','t',paste('p',me$cf$params_to_sample,sep='.'))]
      colnames(gdata)  =  str_replace( colnames(gdata) ,'p\\.','')
      return(data.table(gdata))
    } 

    if (cols=='m') {
      gdata = param_data[,c('objvalue','temp','chain','t','acc',paste('m',me$cf$moments_to_use,sep='.'))]
      colnames(gdata)  =  str_replace( colnames(gdata) ,'m\\.','')
      return(data.table(gdata))
    }

    if (cols=='mp') {
      gdata = param_data[,c('objvalue','temp','chain','t','acc',paste('m',me$cf$moments_to_use,sep='.'),paste('p',me$cf$params_to_sample,sep='.'))]
      colnames(gdata)  =  str_replace( colnames(gdata) ,'p\\.','')
      colnames(gdata)  =  str_replace( colnames(gdata) ,'m\\.','')
      return(data.table(gdata))
    }
  }

  gdata = melt(param_data,measure.vars=paste('p',me$cf$params_to_sample,sep='.'),id=c('objvalue','t','chain','temp','acc'))
  gdata$variable = str_replace(gdata$variable,'p.','')
  return(data.table(gdata))
}




