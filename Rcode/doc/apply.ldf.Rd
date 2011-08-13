\name{apply.ldf}
\alias{apply.ldf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Applies a canonical correlation transformation to the data}
\description{
  Applies a canonical correlation transformation to the combination of the raw
  signal intensities with an initial set of posterior probabilities.
}
\usage{
apply.ldf(full.signal, posterior)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{full.signal}{A matrix with the raw signal intensity. One row per
  data point or sample in the data, and one column for the probability
  of each call. The matrix MUST have row names.}
  \item{posterior}{A matrix of posterior distribution for the
    calls. This matrix must have row names that match the signal
    intensity. The ordering does not have to be the same as the matrix
    of signals but each data point in ``full.signal'' must have a
    corresponding set of posterior probabilities.}
}
\details{
  Do not forget to add row names to both matrices.
}
\value{
  A one-dimensional vector with the transformed canonical corelation
  transformed values.
}
\author{Vincent Plagnol \email{vincent.plagnol@cimr.cam.ac.uk} and Chris Barnes \email{christopher.barnes@imperial.ac.uk} }
\keyword{cluster}