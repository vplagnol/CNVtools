\name{getparams}
\alias{getparams}
\title{Return mixture parameters}
\description{
  This function should be invisible to most users. Given the full expanded data frame,
  getparams returns the number of components, the copy number, the mixture model parameters for each batch, 
  the likelihood of the model and the p(disease|c).
}
\usage{
getparams(d)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{d}{Full expanded data frame. }
}
\value{
  \item{ns}{Number of batches}
  \item{nc}{Copy number of this model}
  \item{nind}{Number of individuals}
  \item{lnL}{log likelihood}
  \item{alpha}{Matrix. The mixture proportions}
  \item{mean}{Matrix. The mixture means}
  \item{var}{Matrix. The mixture variances}
  \item{pdc}{Matrix. The p(disease|c) }
}
\author{ Vincent Plagnol \email{vincent.plagnol@cimr.cam.ac.uk} and Chris Barnes \email{christopher.barnes@imperial.ac.uk} }


