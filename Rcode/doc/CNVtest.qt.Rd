\name{CNVtest.qt}
\alias{CNVtest.qt}
\title{Fits a mixture of Gaussian to CNV data}
\description{
  This function fits a mixture of Gaussians to Copy Number Variant
  data to explore potential correlations between the copy number and a
  quantitative trait.}

\usage{
CNVtest.qt(signal, batch, sample = NULL, qt = NULL, ncomp, n.H0=5, n.H1=0, 
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
    number.'' ~ strata(cn, batch)'' assumes free variances for each
    combination of copy number and batch. More traditional model
    specifications such as ' ~ cn' are also possible, but will converge
    more slowly and might have numerical stability issues.}
  \item{model.var}{A formula as above, but to model the variances. The
    default is the free variance model for each copy number ``~
    strata(cn)'' and the same model specifications as model.means can be
    used.}
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
  \item{control}{A list of parameters that control the behavior of the fitting.}
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
    Fits that don't return 'C' should be excluded.  }
}

\author{Vincent Plagnol <vincent.plagnol@cimr.cam.ac.uk> and Chris Barnes <christopher.barnes@imperial.ac.uk>}
\seealso{apply.pca
	 apply.ldf}
\examples{
	#Load data for CNV for two control cohorts 
	data(A112)
	raw.signal <- as.matrix(A112[, -c(1,2)])
	dimnames(raw.signal)[[1]] <- A112$subject

	#Extract CNV signal using principal components
	pca.signal <- apply.pca(raw.signal)

	#Extract batch, sample
	sample <- factor(A112$subject)
	batches <- rep("ALL",length(sample))

	#Create a fake quantitative trait
	trait <- rnorm(length(sample),mean=9.0,sd=1.0)

	#Fit the CNV with a three component model
	fit.pca <- CNVtest.qt(signal = pca.signal, sample = sample, batch = batches, 
		   	      qt = trait, ncomp = 3, n.H0=3, n.H1=3,
			      model.qt = "~ cn")
			  
	if(fit.pca[['status.H0']] == 'C' && fit.pca[['status.H1']] == 'C'){
	   #Calculate the likelihood ratio
	   LR <- -2*(fit.pca$model.H0$lnL - fit.pca$model.H1$lnL)			  
	
	   #Calculate the pvalue. Has 1 dof since we fit a trend model
	   pvalue <- 1 - pchisq(LR,1)
	}
}
\keyword{htest}
\keyword{cluster}
