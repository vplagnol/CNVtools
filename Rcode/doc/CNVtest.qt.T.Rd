\name{CNVtest.qt.T}
\alias{CNVtest.qt.T}
\title{Fits a mixture of Gaussian to CNV data}
\description{
  This function fits a mixture of T distributions to Copy Number Variant
  data to explore potential correlations between the copy number and a
  quantitative trait.}
\usage{
CNVtest.qt.T(signal, batch, sample = NULL, qt = NULL, ncomp, n.H0=5, n.H1=0, 
	   model.mean = '~ strata(cn)',
           model.var = '~ strata(cn)',
	   model.qt = '~ cn',
           beta.estimated = NULL,
           start.mean = NULL,
           start.var = NULL,
	   control=list(tol=1e-5, max.iter = 3000, min.freq=4) )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{signal}{The vector of intensity values, meant to be a proxy for the
    number of copies.}
  \item{batch}{Factor, that describes how the data points should be
    separated in batches, corresponding to different tehnologies to
    measure the number of DNA copies, or maybe different cohorts in a
    case control framework.}
  \item{sample}{Optional (but recommended). A character vector
    containing a name for each data point, typically the name of the
    individuals.}
  \item{qt}{Quantitative trait values.}
  \item{ncomp}{Number of components one wants to fit to the data.}
  \item{n.H0}{Number of times the EM should be used to maximize the
    likelihood under the null hypothesis of no association, each time
    with a different random starting point. The run that maximizes the 
    likelihood is stored.}
  \item{n.H1}{Number of times the EM should be used to maximize the
    likelihood under the alternate hypothesis of association present, each time
    with a different random starting point. The run that maximizes the
    likelihood is stored.}
  \item{model.mean}{Formula that relates the location of the means for
    the clusters with the number of copies and the different batches if
    there are multiple batches. The default is ``~ strata(cn)''
    that assumes a free model for the cluster locations for each copy
    number. For this T distribution model there is only one
    alternative: '' ~ strata(cn, batch)'' assumes free variances for each
    combination of copy number and batch.}
  \item{model.var}{A formula as above, but to model the variances. The
    default is the free variance model for each copy number ``~
    strata(cn)''. There are three alternative variance models for this T 
    distribution model: ``~ strata(cn,batch)'', ``~ strata(batch)'' or
    even `` ~ 1'' (constant variances for all batches and components).} 
  \item{model.qt}{A formula that relates the number of copies with
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
    means. Must be numeric and the size must match ncomp.}
  \item{start.var}{Optional. A set of starting values for the
    variances. Must be numeric and the size must match ncomp.}
  \item{control}{A list of parameters that control the behavior of the
    fitting.}
}
\value{
  \item{model.H0}{The parameters for the best fit under H0.}
  \item{posterior.H0}{The output dataframe with the estimate posterior
    distribution under H0 as well as the most likely call.}
  \item{status.H0}{A character that describes the status of the fit under H0.
    The possible values are 'C' (converged), 'M' (maximum iterations
    reached), 'P' (posterior distribution problem). Fits that don't
    return 'C' should be excluded.}
  \item{model.H1}{The parameters for the best fit under H1.}
  \item{posterior.H1}{The output dataframe with the estimate posterior
    distribution under H1}
  \item{status.H1}{A character that describes the status of the fit under H1.
    The possible values are 'C' (converged), 'M' (maximum iterations
    reached), 'P' (posterior distribution problem). Fits that don't
    return 'C' should be excluded.} 
}

\author{Vincent Plagnol <vincent.plagnol@cimr.cam.ac.uk> and Chris Barnes <christopher.barnes@imperial.ac.uk>}
\seealso{apply.pca
	 apply.ldf}
\keyword{htest}
\keyword{cluster}
