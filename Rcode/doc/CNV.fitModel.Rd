\name{CNV.fitModel}
\alias{CNV.fitModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Fits a mixture of Gaussian to a set of one dimensional points.}
\description{
  This is the workhorse function, essentially an R wrapper around a lot
  of C code. It fits GLM models to the data.
}
\usage{
CNV.fitModel(ncomp, 
	     nind, 
	     hyp = "H0", 
	     data, 
	     logit.offset,
             design.matrix.mean, 
	     design.matrix.variance,
             design.matrix.disease, 
	     pi.model = 0, 
	     mix.model = 10,
             control = list(tol = 1e-05, max.iter = 3000, min.freq= 4))
}
\arguments{
  \item{ncomp}{integer, number of components to fit to the data}
  \item{nind}{integer, total number of data points}
  \item{hyp}{Hypothesis, can be either H0 or H1}
  \item{data}{The data frame containing the data, in an expanded form
    (one point per individual and copy number)}
  \item{logit.offset}{An option most users will not use. It sets an offset
    when fitting the logit model for the disease status. This is used to
    obtain a profile likelihood when the disease parameter beta varies.}
  \item{design.matrix.mean}{The design matrix that relate mean cluster
    locations with batch.copy numbers.}
  \item{design.matrix.variance}{The design matrix for the cluster variances.}
  \item{design.matrix.disease}{The design matrix for the disease model.}
  \item{pi.model}{0,1,2 fit disease, hetero and quantitative models respectively.}
  \item{mix.model}{Specifies model for the components. }		
  \item{control}{A list of parameters that control the behavior of the fitting.}
}
\details{
  The user is very unlikely to actually use that function which is meant
  as an internal routine, a wrapper around the C code of the package.
  This function is called by the more user friendly function CNVtest.binary.
}
\value{
  \item{data}{The input expanded data frame, but with the posterior
    probabilities estimated.}
  \item{status}{A marker of convergence}
}
\author{Vincent Plagnol \email{vincent.plagnol@cimr.cam.ac.uk} and Chris Barnes \email{christopher.barnes@imperial.ac.uk} }
\seealso{CNVtest.binary}

\keyword{cluster}

