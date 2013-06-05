mopt
====

Parallel Derivative Free Moment Optimization

This is a package intended for people who want to estimate a model based on some distance criteria. The library is built around the square minimization of the distance between a set of moments. The value of the objective function as well as the value of the moments at each evaluations are stored for later analysis. 

The parallelisation overhead is relatively costly, and so I would not recommend using this library if solving your model for a parameter set is under __1 second__. My experience has been that a solving time of around __30 seconds__ is a good place to start.

### Features:

 - works with MPI and OpenMP, but also in searial mode
 - offers serial optimization see [serial optim example](https://github.com/tlamadon/mopt/blob/master/examples/example-serial-optim.r) using the [minqa](http://cran.r-project.org/web/packages/minqa/index.html) package.
 - offers several MCMC implementations (see list of algorithms)
   - see [bgp example](https://github.com/tlamadon/mopt/blob/master/examples/example-bgp.r) for [Likelihood-Free Parallel Tempering](http://arxiv.org/abs/1108.3423)  by  Baragatti Grimaud and Pommeret 
 - offers several reporting tools, but also exports to classic mcmc types to use in R
 - offers a function to compute slices of the objective function in orthogonal directions, see [example](https://github.com/tlamadon/mopt/blob/master/examples/example-slices.r).

### References:

 - [Approximate Bayesian computation](http://en.wikipedia.org/wiki/Approximate_Bayesian_computation)
 - [Monte Carlo Markov Chain](http://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo)
 - "Determining the density of states for classical statistical models: A random walk algorithm to produce a flat histogram "
Wang, F. and Landau, D. Physical Review E  64  056101  (2001)
 - ["An MCMC approach to classical estimation"](http://papers.ssrn.com/sol3/papers.cfm?abstract_id=420371) Chernozhukov and Hong.
 
### Install

    install.packages('devtools');require(devtools);install_github('mopt',user='tlamadon')

### Contributors

 - [Ran Gu](https://github.com/lionup)
 - [Thibaut Lamadon](https://github.com/tlamadon)
 - [Florian Oswald](https://github.com/floswald)
