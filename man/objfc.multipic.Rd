\name{objfc.multipic}
\alias{objfc.multipic}
\title{create an objective function n dimensions and
k local maximums, returns maximum and function}
\usage{
  objfc.multipic(n, k, x.max)
}
\description{
  create an objective function n dimensions and k local
  maximums, returns maximum and function
}
\examples{
x.max = runif(10)
myfunc <- objfc.multipic(10,25,x.max=x.max)
myfunc(x.max)
}

