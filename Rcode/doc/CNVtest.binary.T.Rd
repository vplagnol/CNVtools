\name{CNVtest.binary.T}
\alias{CNVtest.binary.T}
\title{CNV association testing using T distributions}
\description{
  Test for CNV association with binary trait (typically case control)
  using a mixture of T distributions.
}
\usage{
CNVtest.binary.T(signal,batch, sample = NULL, disease.status = NULL,
ncomp, n.H0 = 5, n.H1 = 0, output = "compact",
model.mean = "~ strata(batch, cn)",
model.var = "~ strata(batch, cn)",
model.disease = "~ cn",
beta.estimated = NULL, 
start.mean = NULL, 
start.var = NULL,
control = list(tol = 1e-05, max.iter = 3000, min.freq = 4)) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{signal}{The vector of intensity values, meant to be a proxy for
    the number of copies.}
  \item{batch}{Factor, that describes how the data points should be
    separated in batches, corresponding to different tehnologies to
    measure the number of DNA copies, or maybe different cohorts in a
    case control framework.}
  \item{sample}{Optional (but recommended). A character vector
    containing a name for each data point, typically the name of the
    individuals.}
  \item{disease.status}{In the case control situation a vector of 0 and
    1 indicating which individuals are controls or cases.}
  \item{ncomp}{Number of components one wants to fit to the data.}
  \item{n.H0}{Number of times the EM should be used to maximize the
    likelihood under the null hypothesis of no association, each time
    with a different random starting point. The run that maximizes the 
    likelihood is stored.}
  \item{n.H1}{Number of times the EM should be used to maximize the
    likelihood under the alternate hypothesis of association present, each time
    with a different random starting point. The run that maximizes the
    likelihood is stored.}
\item{output}{The default value, ``compact'', returns a data frame with one line per sample.
    Any other setting witll return a much bigger data frame with one line per individual and copy number.
    This long format is the one used by the underlying fitting algorithm and is only useful if one attempts
    to use CNVtools in a non standard manner.}
  \item{model.mean}{Formula that relates the location of the means for
    the clusters with the number of copies and the different batches if
    there are multiple batches. Should be on the following:
    ``~strata(cn)'' or ``strata(batch, cn)''.}
  \item{model.var}{A formula describing the variance model, as above. The
    default is the free variance model ``~ strata(cn, batch)'' but could
    also be ``~ 1'', ``~ strata(cn)'' or ``~ strata(batch)''.}
  \item{model.disease}{A formula that relates the number of copies with
    the case/control status. The default is a linear trend model ``~
    cn''. Note that this formula will only matter under the alternate
    hypothesis and has no effect under the null.}
  \item{beta.estimated}{Optional. It is used if one wants to fit the
    model for a particular value of the log odds parameter beta
    (essentially if one is interested in the profile likelihood).
    In this case the disease model should be set to ' ~ 1' and the model
    to 'H1'. It will then provide the best model assuming the
    value of beta (the log odds ratio parameter) provided by the user.}
  \item{start.mean}{Optional. A set of starting values for the
    means. Must be numeric and the size must match ncomp. This argument
    can also be a matrix if one wants to specify multiple starting points. When
    passing a matrix as argument, the number of columns should equal the
    number of components, and the number of rows must be greater than
    max(n.H0, n.H1). When in a row some numbers are missing, CNVtools will
    pick the starting points randomly (the default).}
  \item{start.var}{Optional. A set of starting values for the
    variances. Must be numeric and the size must match ncomp. Can also
    be a matrix (see start.mean for details).}
  \item{control}{A list of parameters that control the behavior of the
    fitting. min.freq is the minimum number of data points in a copy
    number class before the algorithm sets the frequency of this class
    to zero. In the presence of a very rare genotype group it might be
    useful to lower this threshold. Note, however, that estimating the
    variance if there are very few individuals in a class may not be
    possible, so setting options such as constant variances
    (i.e. model.var = ' ~1') might be sensible.} 
}
\value{
  \item{model.H0}{The parameters for the best fit under H0.}
  \item{posterior.H0}{The output dataframe with the estimate posterior
    distribution under H0 as well as the most likely call.}
  \item{status.H0}{A character that describes the status of the fit under H0.
    The possible values are 'C' (converged), 'M' (maximum iterations reached), 'P' (posterior distribution problem).
    Fits that don't return 'C' should be excluded.  }
  \item{model.H1}{The parameters for the best fit under H1.}
  \item{posterior.H1}{The output dataframe with the estimate posterior
    distribution under H1}
  \item{status.H1}{A character that describes the status of the fit under H1.
    The possible values are 'C' (converged), 'M' (maximum iterations reached), 'P' (posterior distribution problem).
    Fits that don't return 'C' should be excluded.}
}
\references{
  Finite Mixture Models (Wiley Series in Probability and Statistics),
  G. Mc Lachlan and David Peel
}
\author{
  Vincent Plagnol and Chris Barnes
}
\seealso{
  CNVtest.binary
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{htest}
\keyword{regression}% __ONLY ONE__ keyword per line
\keyword{cluster}
