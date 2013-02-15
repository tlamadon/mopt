




# this example from Robert Casella's book "introducing monte carlo methods with R" (will call that "Casella").
# http://www.amazon.com/gp/product/1441915753?ie=UTF8&tag=chrprobboo-20&linkCode=as2&camp=1789&creative=390957&creativeASIN=1441915753
# not sure how legal it is to use his code. I'll check before this is made public.

# example 5.6. A wiggly function h(x,y) on (R x R)
# ------------------------------------------------

# this function is defined everywhere: it returns a value for any (x,y) you give to it.
wiggly <- function(plotit=FALSE){
	h <- function(x,y) (x * sin(20*y) + y* sin(20*x))^2 * cosh(sin(10*x)*x) + (x * cos(10*y) - y * sin(10*x))^2 * cosh(cos(20*y)*y)
	if (plotit) {
		x <- y <- seq(-3,3,le=435)
		z <- outer(x,y,h)
		par(bg="wheat",mar=c(1,1,1,1))
		persp(x,y,z,theta=155,phi=30,col="red2",ltheta=-120,shade=.75,border=NA,box=FALSE)
	}
	return(h)
}


# mixture of normal distributions
# -------------------------------

# this is a simple function that creates the log likelihood value of a mixture of d normal variables
# for d = 2 there is a plotting method

mixnorm <- function(x,means=c(0,2.5),sds=c(1,1),weight=c(0.25,0.75),plotit=FALSE){
	# x is the vector you choose: x[d] is your current guess of the mean of component d of the mixture
	# returns the value of the likelihood and a plot if 2 dimensions.
	d <- length(x)
	stopifnot(length(sds)==d)
	stopifnot(length(means)==d)
	stopifnot(length(weight)==d)
	stopifnot(sum(weight)==1)

	# make random data to evaluate the likelihood on
	nrand <- 300	# hard coded parameter: n data points into each dimension
	rdat <- lapply(1:d, function(j) rnorm(n=nrand,mean=means[j],sd=sds[j]))
	rdat <- c(unlist(rdat))

	# define likelihood function
	like <- function(x,means,sds,rdat,weight){
		d <- length(x)
		r <- matrix(0,length(rdat),d)
		for (i in 1:d) r[,i] <- weight[i] * dnorm(rdat - x[i],mean=means[i],sd=sds[i])
		r <- sum( log( rowSums( r ) ) )
		return(r)
	}

	if (plotit & d==2){
		# plot contour of likelihood evaluated at random data.
		z <- seq(,length=100)
		z <- seq(-2,5,length=100)
		y <- seq(-3,2,length=100)
		pdata <- expand.grid(z,y)
		pmat <- matrix(apply(pdata,1,like,means,sds,rdat,weight),nrow=100,ncol=100)
		contour(z,y,pmat,nlevels=50,xlab="x",ylab="y")
	}
	return(like(x,means,sds,rdat,weight))
}

# test: mixnorm(c(0,2.5),plotit=TRUE)
# test: mixnorm(c(0,2.5,3),means=c(0,2.5,4),sds=c(1,1,1),weight=c(0.25,0.5,0.25),plotit=TRUE)






