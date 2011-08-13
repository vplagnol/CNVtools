\name{apply.pca}
\alias{apply.pca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Applies to the data a principal component analysis}
\description{
  A simple wrapper around the R function prcomp.
}
\usage{
apply.pca(matrix.signal)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{matrix.signal}{A matrix containing the raw calls. The rows are
    the samples and the columns are the SNPs.}
}
\value{
  A one dimensional vector, one value per sample: this is the first
  principal component. 
}
\author{Vincent Plagnol \email{vincent.plagnol@cimr.cam.ac.uk} and Chris Barnes \email{christopher.barnes@imperial.ac.uk} }
\note{
  The output vector is normalized to have a standard deviation of 1.
}
\seealso{prcomp}
\keyword{cluster}
