\name{CNVtest.binary}
\alias{CNVtest.binary}
\title{Fits a mixture of Gaussian to CNV data}
\description{
  This function fits a mixture of Gaussians to Copy Number Variant data,
  both under the null hypothesis of no association and under the
  alternate hypothesis that the CNV frequencies differ between cases and
  controls. 
}
\usage{
CNVtest.binary(signal, batch, sample = NULL, disease.status = NULL, ncomp, 
	       n.H0 = 5, n.H1 = 0, 
               output = 'compact', 
	       model.mean = "~ strata(batch, cn)", 
	       model.var = "~ strata(batch, cn)",
               model.disease ="~ cn", 
               association.test.strata = NULL,
	       beta.estimated = NULL, 
	       start.mean = NULL, 
               start.var = NULL, 
	       control = list(tol = 1e-05, max.iter = 3000,min.freq = 4))
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
  \item{model.mean}{Formula that describes the linear model for the
    location of the mean signal intensity. The default is ``~ strata(cn,
    batch)'', which means that the mean intensity can take any value
    for any combination of the variables ``cn'' (for copy number) and
    ``batch''. More traditional model description such as ' ~
    as.factor(cn)' for example are also possible, but are likely to be
    slower to fit and less numerically stable than the ``strata''
    notation, which should be preferred.}
  \item{model.var}{A formula as above, but to model the
    variances. Whenever possible and to maximise speed and stability the
    model should be specified using the strata command, for example
    ``strata(batch, cn)'' (the default), meaning that variances are free
    to take any value for each combination of the variables ``batch''
    and ``copy  number''.
    Alternatives such as `` ~ cn'', i.e. variance proportional to the
    number of copies are allowed but slower to fit, and less stable
    numerically.}
  \item{model.disease}{A formula that links the number of copies with
    the case/control status. The default is a logit linear trend model
    ``~ cn''. Note that this formula will only matter under the
    alternate hypothesis and has no effect under the null (model
    descriptions using the ``strata'' command are not allowed for this
    model).}
  \item{association.test.strata}{Optional factor providing the strata when 
   using a stratified test of association (typically, but not always, 
   these are geographic regions of origins of the samples).}
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
  \item{control}{A list of parameters that control the behavior of the fitting.  
   min.freq is the minimum number of data points in a copy number class before 
   the algorithm sets the frequency of this class to zero.In the presence of a 
   very rare genotype group it might be useful to lower this threshold. 
   Note, however, that estimating the variance if there are very few individuals 
   in a class may not be possible, so setting options such as constant variances 
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

\author{Vincent Plagnol <vincent.plagnol@cimr.cam.ac.uk> and Chris Barnes <christopher.barnes@imperial.ac.uk>}
\seealso{apply.pca
	 apply.lda}
\examples{
	#Load data for CNV for two control cohorts 
	data(A112)
	raw.signal <- as.matrix(A112[, -c(1,2)])
	dimnames(raw.signal)[[1]] <- A112$subject

	#Extract CNV signal using principal components
	pca.signal <- apply.pca(raw.signal)

	#Extract batch, sample and trait information
 	batches <- factor(A112$cohort)
	sample <- factor(A112$subject)
	trait <- ifelse( A112$cohort == '58C', 0, 1)

	#Fit the CNV with a three component model
	fit.pca <- CNVtest.binary(signal = pca.signal, sample = sample, batch = batches, 
		   	          disease.status = trait, ncomp = 3, n.H0=3, n.H1=3,
				  model.disease = "~ cn")
			  
	if(fit.pca[['status.H0']] == 'C' && fit.pca[['status.H1']] == 'C'){
	   #Calculate the likelihood ratio
	   LR <- -2*(fit.pca$model.H0$lnL - fit.pca$model.H1$lnL)			  
	
	   #Calculate the pvalue. Has 1 dof since we fit a trend model
	   pvalue <- 1 - pchisq(LR,1)
	}
}
\keyword{htest}
\keyword{cluster}
\keyword{regression}
