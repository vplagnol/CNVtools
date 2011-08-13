\name{CNVtest.select.model}
\alias{CNVtest.select.model}
\title{Select number of components in a CNV }
\description{
  This function fits mixtures of Gaussians to Copy Number Variant data,
  and uses the Bayesian Information Criteria (BIC) to select the most
  likely number of components.
}
\usage{
CNVtest.select.model(signal, batch, sample = NULL, n.H0 = 3, 
		     method="BIC",
		     v.ncomp = 1:6,
                     v.model.component = rep('gaussian',6),
                     v.model.mean = rep("~ strata(cn)",6),
                     v.model.var = rep("~1", 6),
		     control=list(tol=1e-5, max.iter = 500, min.freq=4) )

}
\arguments{
  \item{signal}{The vector of intensity values, meant to be a proxy for the
    number of copies.}
  \item{batch}{ Factor, that describes how the data points should be
    separated in batches, corresponding to different tehnologies to
    measure the number of DNA copies, or maybe different cohorts in a
    case control framework. }
  \item{sample}{ Character vector containing a name for each data point,
    typically the name of the individuals.}
  \item{n.H0}{ Number of times the EM should be used to maximize the
    likelihood and calculate the BIC for each different model. }
  \item{v.ncomp}{ 
	Model specification. Numeric vector specifying number of components to attempt. See discussion.
  }
  \item{v.model.component}{Character vector defining the mixture model. Can either be 'T' or 'gaussian'.}
  \item{v.model.mean}{Model specification. Character vector specifying 
    different models for the component means. See discussion.}
  \item{v.model.var}{
  	Model specification. Character vector specifying different models for the component means. See discussion.
  }
  \item{method}{ Either BIC or AIC }
  \item{control}{A list of parameters that control the behavior of the fitting.}
}
\details{
  The function fits the different models, specified by the vectors v.ncomp, v.model.mean, v.model.var, to the data contained in signal.
  The lengths of v.ncomp, v.model.mean, v.model.var must be equal. The function iterates through the length of these vectors and fits the models
  n.H0 times, keeping the fit with the highest likelihood. The BIC and AIC is calculated for each model, the lowest BIC/AIC indicates the 'best' model.	
  
  In the default model specification, first the data is fit with 1 component, mean model = "~ 1" and variance model = "~ 1". Next the data is fit
  with 2 components, mean model = "~ as.factor(cn)" and variance model "~ 1" etc. 
  
}
\value{
  A data structure containing information from the fitting of the different models specified.

  \item{model}{
   A list (length = number of model fit) containing the models specified by the user.
  }
  \item{BIC}{
   A vector (length = number of model fit) containing the BIC values for each model.
  }
  \item{AIC}{
   A vector (length = number of model fit) containing the AIC values for each model.
  }
  \item{status}{
   A vector (length = number of model fit) containing the status of every fit of every model.
  }
  \item{np}{
   A vector (length = number of model fit) containing the number of parameters of each model.
  }
  \item{posteriors}{
   A list (length = number of model fit) containing the best posterior distribution number of each model.
  }
  \item{selected}{
  The number of the best model. This will be the model that has the lowest BIC or AIC depending on which method was specified.
  }
}
\references{ Schwarz, G. (1978) "Estimating the Dimension of a Model", Annals of Statistics, 6, 461-464. }
\author{Vincent Plagnol \email{vincent.plagnol@cimr.cam.ac.uk} and Chris Barnes \email{christopher.barnes@imperial.ac.uk} }
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
	
	results <- CNVtest.select.model(signal = pca.signal, batch = batches, 
	      	 		        sample = sample, n.H0 = 3)

        # Best model - with the default model setting this is also
	# the number of components
	best_model <- results$selected

	# Look at the fit
	cnv.plot( results[['posteriors']][[best_model]] )
}

\keyword{htest}
\keyword{cluster}
