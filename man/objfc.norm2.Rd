\name{objfc.norm2}
\alias{objfc.norm2}
\title{create an objective function which computes
the distance between simulated means and true means}
\usage{
  objfc.norm2(mu = c(0, 0), sigma = diag(2), ns = 50)

  objfc.norm2(mu = c(0, 0), sigma = diag(2), ns = 50)
}
\arguments{
  \item{norm_mean}{2d vector with the true mean}

  \item{sample_size}{number of random draws to compute
  distance}

  \item{norm_mean}{2d vector with the true mean}

  \item{sample_size}{number of random draws to compute
  distance}
}
\description{
  create an objective function which computes the distance
  between simulated means and true means

  create an objective function which computes the distance
  between simulated means and true means
}
\examples{
myfunc <- objfc.norm2(c(0,0))
myfunc(list(m1=0,m2=0))
myfunc <- objfc.norm2(c(0,0))
myfunc(list(m1=0,m2=0))
}

