\name{test.posterior}
\alias{test.posterior}
\title{ Checks posterior probabilities are monotonic. }
\usage{
test.posterior(frame, ncomp, samples.by.disease = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{frame}{ Posterior data frame. }
  \item{ncomp}{ Number of components. }
  \item{samples.by.disease}{ List containing samples split on disease status. }
}
\description{
	The posterior probability of belonging to a particular component should fall to zero monotonically
	as the signal increases or decreases away from the component mean. This function checks for 
	posterior distributions that do not have this property.
}
\value{
	Returns TRUE is the posterior is not monotonic.
}
\author{Vincent Plagnol \email{vincent.plagnol@cimr.cam.ac.uk} and Chris Barnes \email{christopher.barnes@imperial.ac.uk} }
