\name{EM.starting.point}
\alias{EM.starting.point}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Randomly assigns a starting point for the EM algorithm}
\description{
  This function should be invisible to most users, and is part of our
  the fitting routine using the EM algorithm.
  Our maximum likelihood procedure uses an iterative algorithm called
  Expectation-Maximization. This requires a starting point, chosen at
  random. EM.starting point randomly assigns this starting point. 
}
\usage{
EM.starting.point(d, trait = "binary")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{The dataframe that needs to be initialized}
  \item{trait}{Can be either ``binary'' or ''eQTL''}
}
\value{
  Returns the input data frame with reasonable random starting values.
}
\references{ A. P. Dempster, N. M. Laird, and D. B. Rubin, "Maximum likelihood from incomplete data via the EM algorithm," J. Royal Stat. Soc., vol. 39, pp. 1--38, 1977.}
\author{Vincent Plagnol \email{vincent.plagnol@cimr.cam.ac.uk} and Chris Barnes \email{christopher.barnes@imperial.ac.uk} }
