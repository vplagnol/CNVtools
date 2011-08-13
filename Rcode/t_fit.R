
CNVtest.binary.T <- function(signal,
                             batch,
                             sample = NULL,
                             disease.status = NULL,
                             ncomp,
                             n.H0=5,
                             n.H1=0,
                             output = 'compact',
                             model.mean = '~ strata(batch, cn)',
                             model.var = '~ strata(batch, cn)',
                             model.disease = '~ cn',
                             beta.estimated = NULL,
                             start.mean = NULL,
                             start.var = NULL,
                             control=list(tol=1e-5, max.iter = 3000, min.freq=4) ){
  

  model.nu <- model.var ## for now on use the same model for variances and nus
  start.nu <- NULL
  
  nind <- length(signal)
  batch <- factor(batch)
  if (length(batch) != nind) {stop("Specification of batches does not have the same length as the input signal\n")}

  cat("Attempting to cluster the data with", nind, "individuals and", ncomp, "components\n")
  if ((sd(signal) > 3) || (sd(signal) < 0.33)) {
    cat("Warning: The signal you provided to the algorithm has a standard deviation significantly different from 1. 
	To maximize the stability of the underlying clustering algorithm we recommend you normalize this signal to make sure that the standard deviation is indeed 1.\n")
  }


  if (is.null(sample)) {sample <- paste(batch, c(1:nind), sep='_')} else {
    if (length(sample) != nind) {stop("The sample names you have specified do not match the number of data points\n")}
  }

  if (is.null(disease.status) ){
    if(n.H1 > 0){ stop("Must specify disease.status under H1")}	  
    disease.status <- rep(0, nind )       ##arbitrary
  }
  
  sample.split  <- split(x = sample, f = batch)
  signal.split  <- split(x = signal, f = batch)
  disease.split <- split(x = disease.status, f = batch)
  batch.split   <- split(x = batch, f=batch)

  trait.sample.list <- list( sample[which(disease.status == 0)], sample[which(disease.status == 1)] )
  
  ncohort <- length(levels(batch))
  
  data.for.Ccode <- ExpandData (batch = batch.split, trait = disease.split, names = sample.split, signal = signal.split, ncomp = ncomp)


  ########### Now checking that the starting values are what they should be
  ###### Turning them into matrices to deal with the possibility of multiple starting points

  for (start.values in c('mean', 'var')) {
    arg.start <- get(paste('start', start.values, sep = '.'))
    
    if ( !is.null(arg.start) ) {   ##now may we have provided starting values for the mean
      if (!(class(arg.start) %in% c('numeric', 'matrix'))) stop("If you provided starting values for ", start.values, ", these have to be either numeric of matrix")

      if (class(arg.start) == 'numeric') my.l <- length(arg.start)
      if (class(arg.start) == 'matrix') my.l <- dim(arg.start)[2]

      if (my.l != ncomp) stop("The size of starting.values for ", start.values, " (", my.l, ") does not fit the number of components (", ncomp, ")")


      if (class(arg.start) == 'numeric') assign(x = paste('start', start.values, sep = '.'), value = matrix(data = rep(arg.start, max(n.H0, n.H1)), nrow = max(n.H0, n.H1), byrow = TRUE))
    }
  }




  
############ Now we can start the estimation process
  model.spec <- get.model.spec('T',
                               model.mean,
                               model.var,
                               model.nu,
                               ncomp = ncomp,
                               nbatch = ncohort)[[1]] 

  status.H0 <- ""
  res <- list()
  offset <- rep(0, dim(data.for.Ccode)[1] )
  design.matrix.disease <- model.matrix(data =data.for.Ccode,object = as.formula(model.disease))

  if (n.H0 > 0) {      ################################### First fitting under H0
    best.model.H0 <- 0
    best.lnL.H0 <- -Inf
    best.status.H0 <- ""    


    for(i in 1:n.H0 ){
      cat("Iteration", i, "under H0\n")
      data.for.Ccode <- EM.starting.point(data.for.Ccode)
      
      for (start.values in c('mean', 'var', 'nu')) {  ### plug the starting values in the matrix
        arg.start <- get(paste('start', start.values, sep = '.'))
        if (!is.null(arg.start)) {
          line.arg.start <- arg.start[i,]
          if (sum(is.na(line.arg.start) == 0)) data.for.Ccode[, start.values] <- line.arg.start [ data.for.Ccode$cn ]  
        }
      }


      
      final.frame <- CNV.fitModel(ncomp,
                                  nind,
                                  hyp = "H0",
                                  data.for.Ccode,
                                  logit.offset = offset,
                                  design.matrix.mean = matrix(0),
                                  design.matrix.variance = matrix(0),
                                  design.matrix.disease,                                  
                                  mix.model = model.spec,
				  control=control)
                                  

      if(final.frame$status == "F") lnL <- -Inf
      else lnL <- getparams(final.frame$data)$lnL
	      
      if( lnL > best.lnL.H0 ){
        best.status.H0 <- final.frame$status
        best.lnL.H0 <- lnL
        best.model.H0 <- final.frame$data
      }
    }
    res$model.H0 <- getparams(best.model.H0)
    if (output == 'compact') {
      res$posterior.H0 <- compact.data.frame(best.model.H0)
      res$posterior.H0 <- res$posterior.H0[ match(sample, res$posterior.H0$subject),]
      res$status.H0 = best.status.H0
      if( best.status.H0 == "C" && test.posterior(frame = res$posterior.H0, ncomp = ncomp) == TRUE ) res$status.H0 = "P"
    } else res$posterior.H0 <- best.model.H0

    # Set the status
   
  }

  


  if (n.H1 > 0)  {  ################################## Now fitting under H1
    best.model.H1 <- 0
    best.lnL.H1 <- -Inf
    best.status.H1 <- ""
    status.H1 <- ""
    if (!is.null(beta.estimated))    offset <- beta.estimated * data.for.Ccode$cn
    
    for(i in 1:n.H1 ){
      cat("Iteration", i, "under H1\n")
      if  ((i == 1) & (n.H0 > 0)) data.for.Ccode <- best.model.H0  else {
        data.for.Ccode <- EM.starting.point(data.for.Ccode)
        
        for (start.values in c('mean', 'var', 'nu')) {  ### plug the starting values in the matrix
          arg.start <- get(paste('start', start.values, sep = '.'))
          if (!is.null(arg.start)) {
            line.arg.start <- arg.start[i,]
            if (sum(is.na(line.arg.start) == 0)) data.for.Ccode[, start.values] <- line.arg.start [ data.for.Ccode$cn ]  
          }
        }
      }

      final.frame <- CNV.fitModel(ncomp,
                                  nind,
                                  hyp = "H1",
                                  data = data.for.Ccode,
                                  logit.offset = offset,
                                  design.matrix.mean = matrix(0),
                                  design.matrix.variance = matrix(0),
                                  design.matrix.disease,
                                  mix.model = model.spec,
                                  control=control)
      
      if(final.frame$status == "F") lnL <- -Inf
      else lnL <- getparams(final.frame$data)$lnL
      
      if( lnL > best.lnL.H1 ){
        best.status.H1 <- final.frame$status
        best.lnL.H1 <- lnL
        best.model.H1 <- final.frame$data
      }
      
    }     
    res$model.H1 <- getparams(best.model.H1)
    res$posterior.H1 <- compact.data.frame(best.model.H1)
    res$posterior.H1 <- res$posterior.H1[ match(sample, res$posterior.H1$subject),]  ## just some reordering

    # Set the status
    res$status.H1 <- best.status.H1 
    if(best.status.H1 == "C" && test.posterior(frame = res$posterior.H1, ncomp = ncomp, samples.by.disease = trait.sample.list) == TRUE ) res$status.H1 <- "P"
  }
  
  return(res)
}





get.model.spec <- function(model.component,
                           model.mean,
                           model.var,
                           model.nu,
                           design.matrix.mean = NULL,
                           design.matrix.variance = NULL,
                           ncomp,
                           nbatch)
{

        
  nparam.alpha <- (ncomp-1)
  model.spec <- 0
  if( model.component != 'gaussian' && !(model.component %in% c('t', 'T')) ){
    cat("Did not specify valid model.component - going with gaussians\n")
    model.component <- 'gaussian'
  }
  
  model.mean <- gsub(model.mean, pattern = ' ', replacement = '')
  model.var <- gsub(model.var, pattern = ' ', replacement = '')
  full.model <- c('~strata(batch,cn)', '~strata(cn,batch)')

  
  if( model.component == 'gaussian' ){
    model.spec <- 10
############# set number of parameters = (n in alpha) + (n in mean) + (n in vars)


    if (!is.null(design.matrix.mean)) nparam.mean <- ncol(design.matrix.mean)
    if (!is.null(design.matrix.variance)) nparam.var <- ncol(design.matrix.variance)
    
    
    if (model.mean %in% full.model)  nparam.mean <- nbatch*ncomp
    if (model.mean == '~strata(cn)') nparam.mean <- ncomp
    
    if( model.var == '~ 1' )            nparameter.var <- 1
    if (model.var == '~strata(batch)' ) nparameter.var <- nbatch
    if (model.var == '~strata(cn)' )    nparameter.var <- ncomp
    if (model.var %in% full.model)      nparameter.var <- nbatch*ncomp

    nparameter <- nparam.alpha + nparam.mean + nparam.var
  }
  
  if( model.component %in% c('T', 't') ){
    nparameter <- nparam.alpha

    model.nu <- gsub(model.nu, pattern = ' ', replacement = '')
    
########## mean first
    if( model.mean == '~strata(cn)' ){
      model.spec <- 2300
      nparameter <- nparameter + ncomp
    }
    else if( model.mean %in% full.model) {
      model.spec <- 2400
      nparameter <- nparameter + nbatch*ncomp
    }
    else{
      cat("Mean formula not recognized, going with ~ strata(batch, cn)\n")
      model.spec <- 2400
      nparameter <- nparameter + nbatch*ncomp
    }  

########## variance
    if( model.var == '~ 1' ){ 
      model.spec <- model.spec + 10
      nparameter <- nparameter + 1
    }
    else if (model.var == '~strata(batch)' ){
      model.spec <- model.spec + 20
      nparameter <- nparameter + nbatch
    }
    else if (model.var == '~strata(cn)' ){
      model.spec <- model.spec + 30
      nparameter <- nparameter + ncomp
    }
    else if (model.var %in% full.model) {
      model.spec <- model.spec + 40
      nparameter <- nparameter + ncomp*nbatch
    }	
    else{
      cat("Variance formula not recognized, going with ~ strata(batch, cn)\n")
      model.spec <- model.spec + 40
      nparameter <- nparameter + ncomp*nbatch		
    } 

########## nu
    if (model.nu == '~ 1'){
      model.spec <- model.spec + 1
      nparameter <- nparameter + 1
    }	
    else if (model.nu == '~strata(batch)' ){
      model.spec <- model.spec + 2
      nparameter <- nparameter + nbatch	
    }
    else if (model.nu == '~strata(cn)' ){
      model.spec <- model.spec + 3
      nparameter <- nparameter + ncomp
    }
    else if (model.nu %in% full.model) {
      model.spec <- model.spec + 4
      nparameter <- nparameter + nbatch*ncomp
    }	
    else{
      cat("constrained variance fitting for t not implemented, going with ~strata(batch,cn)\n")
      model.spec <- model.spec + 4
      nparameter <- nparameter + nbatch*ncomp
    } 
  }
  return( list(model.spec,nparameter) )
}











### batch will be factor
## signal will be numeric
## qt is quantitative trait
CNVtest.qt.T <- function(signal,
                         batch,
                         sample = NULL,
                         qt = NULL, ncomp,
                         n.H0=5,
                         n.H1=0, 
                         model.mean = '~ strata(cn)',
                         model.var = '~ strata(cn)',
                         model.qt = '~ cn',
                         beta.estimated = NULL,
                         start.mean = NULL,
                         start.var = NULL,
                         control=list(tol=1e-5, max.iter = 3000, min.freq=4) ){

  model.nu <- model.var ##for now, let us do things in that way
  start.nu <- NULL
  
  nind <- length(signal)
  batch <- factor(batch)
  if (length(batch) != nind) {stop("Specification of batches does not have the same length as the input signal\n")}
  if (length(qt) != nind) {stop("Specification of qt does not have the same length as the input signal\n")}
    
  if (is.null(sample)) {sample <- paste(batch, c(1:nind), sep='_')} else {
    if (length(sample) != nind) {stop("The sample names you have specified do not match the number of data points\n")}
  }

  if (is.null(qt) ){
	if(n.H1 > 0){ stop("Must specify quantitative trait under H1")}	  
  	
	# Create arbitrary disease.status for fitting under null
  	qt <- rep(0, nind )       
  }

  sample.split  <- split(x = sample, f = batch)
  signal.split  <- split(x = signal, f = batch)
  qt.split <- split(x = qt, f = batch)
  batch.split   <- split(x = batch, f=batch)

  ncohort <- length(levels(batch))
  
  data.for.Ccode <- ExpandData (batch = batch.split, trait = qt.split, names = sample.split, signal = signal.split, ncomp = ncomp)


  if ( !missing(start.mean) && !is.null(start.mean) ) {   ##now may we have provided starting values for the mean
    if ((length(start.mean) != ncomp) || (class(start.mean) != 'numeric'))  {stop("You provided invalid starting values for the mean")}
  } else {start.mean <- NULL}

  
  if ( !is.null(start.var) && !missing(start.var) ) {   ##now may we have provided starting values for the variances
    if ((length(start.var) != ncomp) || (class(start.var) != 'numeric'))  {stop("You provided invalid starting values for the variances")}
  } else {start.var <- NULL}

  if ( !is.null(start.nu) && !missing(start.nu) ) {   ##now may we have provided starting values for the nu
    if ((length(start.nu) != ncomp) || (class(start.nu) != 'numeric'))  {stop("You provided invalid starting values for the nu")}
  } else {start.nu <- NULL}

  design.matrix.mean <- as.matrix(0)
  design.matrix.variance <- as.matrix(0)
  design.matrix.qt <- model.matrix(data = data.for.Ccode, object = as.formula(model.qt))

  model.spec <- get.model.spec('T',model.mean,model.var,model.nu,design.matrix.mean,design.matrix.variance,ncomp,ncohort)[[1]] 	 

  status.H0 <- ""
  res <- list()
  offset <- rep(0, dim(data.for.Ccode)[1] )

  
################################### First fitting under H0
  if (n.H0 > 0) {
    best.model.H0 <- 0
    best.lnL.H0 <- -Inf
    best.status.H0 <- ""    

    
    for(i in 1:n.H0 ){
      cat("Iteration", i, "under H0\n")
      data.for.Ccode <- EM.starting.point(data.for.Ccode)
      if (!is.null(start.mean)) data.for.Ccode$mean = start.mean[ data.for.Ccode$cn ]
      if (!is.null(start.var)) data.for.Ccode$var = start.var[ data.for.Ccode$cn ]            
      if (!is.null(start.nu)) data.for.Ccode$nu = start.nu[ data.for.Ccode$cn ]   

      final.frame <- CNV.fitModel(ncomp, nind, hyp = "H0", data.for.Ccode, logit.offset = offset,
                                  design.matrix.mean = matrix(0),
                                  design.matrix.variance = matrix(0),
                                  design.matrix.qt,
				  pi.model = 2,
				  mix.model = model.spec,
				  control=control)
      
      if(final.frame[['status']] == "F") lnL <- -Inf
      else lnL <- getparams(final.frame[['data']])$lnL
      	  
      if( lnL > best.lnL.H0 ){
        best.status.H0 <- final.frame[['status']]
        best.lnL.H0 <- lnL
        best.model.H0 <- final.frame[['data']]
      }
    }
    res$model.H0 = getparams(best.model.H0)
    res$posterior.H0 = compact.data.frame(best.model.H0)
    res$posterior.H0 <- res$posterior.H0[ match(sample, res$posterior.H0$subject),]	    

    # Set the status
    res[['status.H0']] = best.status.H0
    if( best.status.H0 == "C" && test.posterior(res[['posterior.H0']], ncomp) == TRUE ) res[['status.H0']] = "P"  	   
 }

  
################################## Now fitting under H1

  if (n.H1 > 0)  {
    best.model.H1 <- 0
    best.lnL.H1 <- -1000000000
    best.status.H1 <- ""
    status.H1 <- ""
    if (!is.null(beta.estimated))    offset <- beta.estimated * data.for.Ccode$cn

    for(i in 1:n.H1 ){
      cat("Iteration", i, "under H1\n")
      if  ((i==1) & (n.H0 > 0)) {data.for.Ccode <- best.model.H0}
      else{
        data.for.Ccode <- EM.starting.point(data.for.Ccode)
        if (!is.null(start.mean)) data.for.Ccode$mean = start.mean[ data.for.Ccode$cn ]
        if (!is.null(start.var)) data.for.Ccode$var = start.var[ data.for.Ccode$cn ]            
      }

      final.frame <- CNV.fitModel(ncomp,
                                  nind,
                                  hyp = "H1",
                                  data = data.for.Ccode,
                                  logit.offset = offset,
                                  design.matrix.mean,
                                  design.matrix.variance,
                                  design.matrix.qt,
				  pi.model = 2,
				  mix.model = model.spec,
                                  control=control)
            
      if(final.frame[['status']] == "F") lnL <- -1000000000
      else lnL <- getparams(final.frame[['data']])$lnL
	
      if( lnL > best.lnL.H1 ){
        best.status.H1 <- final.frame[['status']]
	best.lnL.H1 <- lnL
        best.model.H1 <- final.frame[['data']]
      }
      
    }     
    res[['model.H1']] = getparams(best.model.H1)
    res[['posterior.H1']] = compact.data.frame(best.model.H1)

    # Set the status
    res[['status.H1']] <- best.status.H1	
    # This need to be fixed. Current method of determining "P" will not work for QT
  }
	  
  return(res)
}

