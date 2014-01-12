\name{mopt.load}
\alias{mopt.load}
\title{Load an existing mopt config}
\usage{
  mopt.load(filename = "", remote = "")
}
\arguments{
  \item{filename}{optional local filename}

  \item{remote}{optional. full \code{scp} path: username}
}
\description{
  loads a mopt config from file. This is very useful for
  on-the-fly analysis of results that are generated on a
  remote machine, or to process results on your local
  machine. Fetching remote files requires a publickey
  authentication setup for ssh to the remote server. your
  connetion must work without having to type a password.
  It's easy to setup the ssh-agent and add a passphrase to
  an existing key. Preferably set up a config file for ssh.
  see the references.
}
\examples{
\dontrun{
me <- mopt.load(remote="hpc:git/wagebc/Rdata/evaluations.educ1.dat")
predict(me,'p.sd')
}
}
\references{
  \url{https://help.github.com/articles/working-with-ssh-key-passphrases}
  \url{http://nerderati.com/2011/03/simplify-your-life-with-an-ssh-config-file/}
}

