\name{ExpandData}
\alias{ExpandData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Expands a CNV input data frame for the maximum likelihood routines}
\description{
	This function should be invisible to most users. The methods within CNV.fitModel require that the CNV data
	is expanded N times where N is the number of copies. This allows the use of Generalized Linear Models (GLM)
	in constraining the Gaussian mixture component locations and spreads to be functions of the copy number.  	
}
\usage{
ExpandData(batch, trait, names, signal, ncomp, association.strata = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{batch}{List of vectors, one vector per batch in the data. Because each
    element in the list corresponds to a batch, each element should be a
    vector with a unique values repeated as many times as the number of
    data point in the batch.}
  \item{trait}{List of vectors, one vector per batch in the data. Each element of
    the list can be either a vector of quantitative traits or a vector
    of 0 and 1 in a case/control framework}
  \item{names}{List of vectors, one vector per batch in the data containing names
    for each data point, typically individual IDs.}
  \item{signal}{List, one vector per batch in the data.}
  \item{ncomp}{Integer, number of components one wants to fit to the data}
  \item{association.strata}{Optional, a factor vector containing the strata when using a stratified test of association.}
}
\value{
  An expanded data frame needed for CNVfit.binary.
}
\author{Vincent Plagnol \email{vincent.plagnol@cimr.cam.ac.uk} and Chris Barnes \email{christopher.barnes@imperial.ac.uk} }

\seealso{ CNVtest.binary }

