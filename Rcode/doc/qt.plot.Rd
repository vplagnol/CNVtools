\name{qt.plot}
\alias{qt.plot}
\title{Makes signal vs trait plots and posterior probabilty distributions}
\description{
  Makes signal vs trait and formatted density plots from the data frame returned by CNVtest.qt 
}
\usage{
qt.plot(DataFrame.list, main='', hist.or.dens='histogram')
}
\arguments{
  \item{DataFrame.list}{The output obtained from the CNVtools fitting algorithm CNVtest.qt}
  \item{main}{Potential title for the graph}
  \item{hist.or.dens}{Either 'histogram' or 'density' to plot the data as an
    histogram or using a kernel density estimator}  
}
\author{ Vincent Plagnol \email{vincent.plagnol@cimr.cam.ac.uk} and Chris Barnes \email{christopher.barnes@imperial.ac.uk}  }
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
			  
	qt.plot(fit.pca)	  
}
