\name{cnv.plot}
\alias{cnv.plot}
\title{Plots posterior probabilty distributions}
\description{
  Makes formatted density plots from the posterior data frame(s) returned by CNVtest.binary 
}
\usage{
cnv.plot(posterior, hist.or.dens='histogram', batch = NULL, freq = NULL, ...)
}
\arguments{
  \item{posterior}{The posterior distribution obtained from the CNVtools
    fitting  algorithm, for example using CNVtest.binary}
  \item{hist.or.dens}{Either 'histogram' or 'density' to plot the data as an
    histogram or using a kernel density estimator}
  \item{batch}{character vector (usually of length 1, but not always), designing the batches one wants to plot.}
  \item{freq}{This argument is only relevant when hist.or.dens='histogram'
    (the default). It matches the argument freq of the hist
    function. With freq = FALSE frequencies, and not raw counts, are
    shown in the histogram.}
  \item{...}{Usual arguments passed to the hist function, including main or breaks for example.}
}
\author{ Vincent Plagnol \email{vincent.plagnol@cimr.cam.ac.uk} and Chris Barnes \email{christopher.barnes@imperial.ac.uk}  }
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
	
	cnv.plot(fit.pca[['posterior.H0']], batch = '58C', breaks = 30)			  	
		  
}
