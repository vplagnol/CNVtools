######## final set of functions for the CNVtools package



EM.starting.point <- function(d, trait='binary') {
  ncohort <- length( levels(d$batch) )
  ncomp <- range(d$cn)[2]
  d$alpha <- c(1/ncomp)
  
  d$pr <- c(-10)
  d$posterior <- c(0)

  if(trait == "binary"){
    Ncontrol <- length( d[d$cn == 1 & d$d == "0", ]$subject )
    Ncase <- length( d[d$cn == 1 & d$d == "1", ]$subject )
    
    if(Ncase == 0)  {d$pdc <- 0} else {d$pdc <- Ncase/(Ncontrol + Ncase)}
  }
  else  {d$pdc <- mean(d$trait); d$rs  <- sd(d$trait)}


  
                                        # Initialize mu and var for each cohort seperately
  lev <- levels(d$batch)
  for (my.batch in lev){
    ends <- range(d$signal)
    r <- diff(ends)

    pf <- 0.2
    if (ncomp > 1) {interval <- r/(ncomp-1)} else {interval <- r}
    mean <- seq(from=ends[1],by=interval,length.out=ncomp)
    mean <- mean + c( runif(1:ncomp, -1*pf*r, +1*pf*r ) )
    vars <- pmax( runif(1:ncomp,0.,diff(range(d$signal)/(6*ncomp) ) ), 0.1)
    nus <- runif(ncomp,5,15) 	

    mean <- sort(mean)
    d$mean <- ifelse ( d$batch == my.batch, mean[ d$cn ], d$mean )
    d$var <- ifelse ( d$batch == my.batch, vars[ d$cn ], d$var )
    d$nu <- ifelse ( d$batch == my.batch, nus[ d$cn ], d$nu )
    
  }

  return(d)
}





ExpandData <- function(batch, trait, names, signal,ncomp, association.strata = NULL)  {
  
                                        # ns is number of cohorts
                                        # m is list of data, each entry being a vector of data
                                        # disease is disease status (0 or 1) in same form as data
                                        # ncomp is number of components
  ncohort <- length(batch)
  
  cn <- vector(mode = "numeric", length = 0)
  trait.f <- vector(mode = "numeric", length = 0)
  batch.f <- vector(mode = "character", length = 0)
  subject <- vector(mode = "character", length = 0)
  signal.f <- vector(mode = "numeric", length = 0)
  a.strata.f <- vector(mode = "numeric", length = 0)
  
  for(nc in 1:ncomp){	  #for each possible number of copies I suppose
    for(ns in 1:ncohort){
      
      nind <- length(signal[[ns]])
      cn        <- append(cn,         rep(nc, nind) )  
      trait.f   <- append(trait.f,    trait[[ns]] )
      batch.f   <- append(batch.f,    as.character(batch[[ns]]) )	
      subject   <- append(subject,    as.character(names[[ns]]) )
      signal.f  <- append(signal.f,   signal[[ns]] )
      if (!is.null(association.strata)) a.strata.f  <- append(a.strata.f,   association.strata[[ns]] ) else a.strata.f <- append(a.strata.f,rep(1, nind))
      
    }	      
  }
  DataFrame <- data.frame(batch=batch.f,
                          cn=cn,
                          trait=trait.f,
                          subject=subject,
                          signal=signal.f,
                          alpha=0,
                          mean=0,
                          var=0,
                          nu=0,
                          pr=0,
                          posterior=0,
                          pdc=0,
                          strats.var = as.integer(0),
                          strats.mean = as.integer(0),
                          association.strata = a.strata.f)
  return (DataFrame)
}

getparams <- function(d)
{
  p <- list()
  p[["ns"]] <- length( levels(d$batch) )
  p[["nc"]] <- range(d$cn)[2]
  p[["nind"]] <- dim(d)[1]/p[["nc"]]

  maxLike <- tapply(d$pr, d$subject, max)  ##takes the max likelihood
  p[["lnL"]] <- sum ( maxLike + log( tapply(exp(d$pr - maxLike[ d$subject]), FUN=sum, d$subject)) )
  
  p[["alpha"]] <- matrix(0,nrow=p$nc,ncol=p$ns)
  p[["mean"]] <- matrix(0,nrow=p$nc,ncol=p$ns)
  p[["var"]] <- matrix(0,nrow=p$nc,ncol=p$ns)
  p[["nu"]] <- matrix(0,nrow=p$nc,ncol=p$ns)
  p[["pdc"]] <- matrix(0,nrow=p$nc,ncol=p$ns)
  
  lev <- levels(d$batch)
  for(j in 1:p$ns) {
    for(i in 1:p$nc) {
      p$mean[i,j] <- .Call("get_first_match", nrow(d), as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), as.numeric(j), d$mean)
      p$alpha[i,j] <- .Call("get_first_match", nrow(d), as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), as.numeric(j), d$alpha)
      p$var[i,j] <- .Call("get_first_match", nrow(d), as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), as.numeric(j), d$var)
      p$nu[i,j] <- .Call("get_first_match", nrow(d), as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), as.numeric(j), d$nu)
      p$pdc[i,j] <- .Call("get_first_match", nrow(d), as.numeric(d$cn), as.numeric(d$batch), as.numeric(i), as.numeric(j), d$pdc)
    }	
  }
  return (p)
}




## Compact the full data.frame into a more practical format, from the exapanded version to something tighter
compact.data.frame <- function (full.frame) 
{   

  full.frame <- full.frame[order(full.frame$cn),]
  full.frame.mod <- do.call(rbind.data.frame, split(x = full.frame$posterior, f = full.frame$subject))
  names(full.frame.mod) <- paste('P', c(1:(dim(full.frame.mod)[2])), sep='')       
  full.frame.mod$cn <-apply(full.frame.mod, FUN=which.max, MAR=1)    
  full.frame.mod$subject <- row.names(full.frame.mod)                              
  full.frame <- subset( full.frame[ ,c('subject', 'batch', 'signal', 'trait')], full.frame$cn == 1)
  full.frame <- merge(full.frame, full.frame.mod)        
  return (full.frame)
}


####################### plots the fit, it takes the list structure output from CNVtest.*
qt.plot <- function(DataFrame.list, main='', hist.or.dens = 'histogram') {
                                        #Assumes number of cohorts is one
  posterior.H0 <- DataFrame.list$posterior.H0    	
  LR <- -2*(DataFrame.list$model.H0$lnL - DataFrame.list$model.H1$lnL)
  
  cohort.list <- unique(as.character(posterior.H0$batch))
  ncohort <- length(cohort.list)
  ncomp <- max(posterior.H0$cn)
  
  if (ncohort != 1) {stop("QT plotting supported for only one batch\n")}

  par(mfrow=c(1,3))
  
  onecohort <- subset(posterior.H0, posterior.H0$batch == cohort.list[1])
  
                                        # Plot raw signal vs qt    
  fit <- lm(posterior.H0$trait ~ posterior.H0$signal)
  p <- summary(fit)$coefficients[2,4]
  plot( posterior.H0$signal, posterior.H0$trait, xlab='Signal', ylab='Trait', main=paste("Signal vs Trait : p = ",signif(p,3)) ) 
  abline(fit,col="blue")	
  
                                        # Plot MAP vs qt 
  fit <- lm(posterior.H0$trait ~ posterior.H0$cn)
  p <- summary(fit)$coefficients[2,4]
  plot( posterior.H0$cn, posterior.H0$trait, xlab='MAP (H0)', ylab='Trait', main=paste("MAP vs Trait : p = ",signif(p,3)) ) 
  abline(fit,col="blue")

	# Plot full model      
  loc.main <- paste("Full model : LR = ",signif(LR,3))	
  
  if (hist.or.dens == 'density') {
    dens <- density(onecohort$signal)
    plot( dens, col='red', xlim=range(onecohort$signal), xlab='Signal', ylab='', main=loc.main)
    my.max <- max(dens$y)
  }
      
  if (hist.or.dens == 'histogram') {
    my.hist <- hist( onecohort$signal, col='red', xlim=range(onecohort$signal), xlab='Signal', ylab='', main=loc.main, breaks=30)
    my.max <- max(my.hist$counts)
  }
  
  onecohort <- onecohort[ order(onecohort$signal), ]
  col <- 1          
  for (i in c(1:ncomp)) {
    lines(onecohort$signal, my.max*onecohort[, paste("P", i, sep='')], type="l", col=col); 
    col <- col + 2
  }
  if (main != '') title(main = main, outer = TRUE, line=-1)
}


####################### plots the data, it takes the posterior frame



cnv.plot <- function (posterior, hist.or.dens = "histogram", batch = NULL, freq = NULL, ...)  {
  
    if (!is.null(batch))  posterior <- posterior[posterior$batch %in% batch, ]
    posterior <- posterior[order(posterior$signal), ]
    if (hist.or.dens == "density") {
        dens <- density(posterior$signal)
        plot(dens, ...)
        my.max <- max(dens$y)
    }
    if (hist.or.dens == "histogram") {
      my.hist <- hist(posterior$signal, freq = freq, ...)
      my.max <- max(my.hist$counts)
      if (!is.null(freq)) {if (freq == FALSE) my.max <- max(my.hist$density)}
    }
    col <- 1
    ncomp <- max(posterior$cn)
    for (i in c(1:ncomp)) {
      lines(posterior$signal, my.max*posterior[, paste("P", i, sep = "")], type = "l", col = col)
      col <- col + 2
    }
}



### the workhorse function
CNV.fitModel <- function(ncomp,
                         nind,
                         hyp = "H0",
                         data,
                         logit.offset, 
                         design.matrix.mean,
                         design.matrix.variance,
                         design.matrix.disease,
                         pi.model = 0,
			 mix.model = 10,
			 control=list(tol=1e-5, max.iter = 3000, min.freq=4)) {

  pi.mod <- as.integer(pi.model)
  mix.mod <- as.integer(mix.model)

  ###### Here we make sure that the numbering of the strata variance makes sense
  initial <- sort(unique(data$strats.var))
  data$strats.var <- as.integer(rank(initial)[ match (x = data$strats.var, table = initial) ])

  initial <- sort(unique(data$strats.mean))
  data$strats.mean <- as.integer(rank(initial)[ match (x = data$strats.mean, table = initial) ])

  initial <- sort(unique(data$association.strata))
  data$association.strata <- as.integer(rank(initial)[ match (x = data$association.strata, table = initial) ])
  
  if(mix.model < 10) {stop("Specification of mix.model incorrect\n")}
  res <- .Call("C_fitmodel", 
	       as.integer(ncomp), 
               as.integer(nind), 
	       hyp, 
	       data, 
	       logit.offset,
               as.matrix(design.matrix.mean[, -1]), 
               as.matrix(design.matrix.variance[,-1]), 
               as.matrix(design.matrix.disease[,-1]), 
               control,
               mix.mod,	
               pi.mod)

  
  new.data <- res[[1]]
  data$posterior <- new.data[,1]
  data$mean <- new.data[,2]
  data$var <- new.data[,3]
  data$pr <- new.data[,4]
  data$alpha <- new.data[,5]
  data$pdc <- new.data[,6]
  data$nu <- new.data[,7]	

  return (list(data = data, status = res[[2]]))
}



### selects number of components using the BIC
CNVtest.select.model <- function(signal, batch, sample = NULL, n.H0 = 3, 
				 method="BIC",
				 v.ncomp = 1:6,
				 v.model.component = rep('gaussian',6),
				 v.model.mean = rep("~ strata(cn)",6),
				 v.model.var = rep("~1", 6),
				 control=list(tol=1e-5, max.iter = 500, min.freq=4) )
{
  start.values <- list()
  v.model.nu <- v.model.var

  y <- c( length(v.ncomp), length(v.model.component), length(v.model.mean), length(v.model.var), length(v.model.nu) ) 
  if( sum( ifelse( y == y[1], 1, 0 ) ) != 5 ){
    stop("ncomp, model.mean, model.var, model.nu, must have same length")
  }

  bics <- vector(mode = "numeric", length = 0);
  aics <- vector(mode = "numeric", length = 0);
  
  nind <- length(signal)
  batch <- factor(batch)
  if (length(batch) != nind) {stop("Specification of batches does not have the same length as the input signal\n")}
  
  if (is.null(sample)) {sample <- paste(batch, c(1:nind), sep='_')} else {
    if (length(sample) != nind) {stop("The sample names you have specified do not match the number of data points\n")}
  }
  
  sample.split  <- split(x = sample, f = batch)
  signal.split  <- split(x = signal, f = batch)
  disease.split <- split(x = rep(0,nind), f = batch )
  batch.split   <- split(x = batch, f=batch)
  nbatches <- length(levels(batch))

  ncohort <- length(levels(batch))
  nmodel <- length(v.ncomp)

  res <- list()
  res[['model']] <- vector("list",nmodel)
  res[['BIC']] <- vector("numeric",nmodel)		  
  res[['AIC']] <- vector("numeric",nmodel)
  res[['status']] <- vector("character",nmodel)
  res[['np']] <- vector("numeric",nmodel)
  res[['posteriors']] <- vector("list",nmodel)	
  res[['models']] <- vector("list",nmodel)
  res[['selected']] <- 0

  for(m in c(1:nmodel) ){
    message(paste("Fitting model ",m))

    status <- ""  
    	
    ncomp <- v.ncomp[m]	
    model.component <- v.model.component[m]
    model.mean <- v.model.mean[m]
    model.var <- v.model.var[m]
    model.nu <- v.model.nu[m]
     
    data.for.Ccode <- ExpandData (batch = batch.split, trait = disease.split, names = sample.split, signal = signal.split, ncomp = ncomp,association.strata = NULL)
    design.matrix.disease <- model.matrix(data = data.for.Ccode, object = as.formula("~ 1"))

    if( model.component == 'gaussian'){  ##need to fix this part to deal with the strata stuff
      special <- c("strata")      
      nparameters <- ncomp - 1  ##this is for alpha
      my.model <- 10
      for (design in c('var', 'mean')) {

        my.formula <- as.formula(get(paste('model', design, sep = '.')))
        Terms <- terms(my.formula, special, data=data.for.Ccode)
        strats <- attr(Terms, "specials")$strata

        if (!is.null(strats)) {


          mm <- list()
          mm[[1]] <- as.name("model.frame")
          mm[[2]] <- Terms
          names(mm)[2] <- "formula"
          mm[['data']] <- data.for.Ccode
          mm <- as.call(c(as.list(mm), list(na.action=as.symbol("na.omit"))))
          mm <- eval(mm)
          
          temps <- untangle.specials(Terms, "strata", 1)
          data.for.Ccode[, paste('strats', design, sep = '.')] <- as.integer(strata(mm[, temps$vars], shortlabel=TRUE))
          nstrata <- length(unique(data.for.Ccode[, paste('strats', design, sep = '.')]))
          nparameters <- nparameters + nstrata
          

          if (nstrata == 1) assign(x = paste('design.matrix', design, sep = '.'), value = model.matrix(data = data.for.Ccode, object = as.formula(' ~ 1') ))
          else assign(x = paste('design.matrix', design, sep = '.'), value = matrix(0))
        } else {
          assign(x = paste('design.matrix', design, sep = '.'), value = model.matrix(data = data.for.Ccode, object = my.formula ))
          nparameters <- nparameters + ncol(get(paste('design.matrix', design, sep = '.')))
        }
      }
    } else {  ##T case
      design.matrix.mean <- as.matrix(0)
      design.matrix.var <- as.matrix(0)
      model.spec <- get.model.spec(model.component,model.mean,model.var,model.nu,design.matrix.mean,design.matrix.var,ncomp,ncohort)
      my.model <-  model.spec[[1]]
      nparameters <- model.spec[[2]]
    }
    

    
    best.model.H0 <- 0
    best.lnL.H0 <- -Inf
    best.status <- "F"    
    
    offset <- rep( 0, dim(data.for.Ccode)[1] )
    suc <- 0    

    start.mean <- NULL
    start.var <- NULL
    start.nu <- NULL


    
    ########## This is where I can specify some starting values for the mean
    mean.label <- paste('M', ncomp, sep = '')      
    if (mean.label %in% names(start.values)) start.mean <- start.values [[ mean.label ]]

    for(i in 1:n.H0 ){
      message(paste("\t Iteration",i))
      data.for.Ccode <- EM.starting.point(data.for.Ccode)
      
      
      for (start.v in c('mean', 'var', 'nu')) {  ### plug the starting values in the matrix
        arg.start <- get(paste('start', start.v, sep = '.'))
        if (!is.null(arg.start)) {
          line.arg.start <- arg.start[i,]
          if (sum(is.na(line.arg.start) == 0)) data.for.Ccode[, start.v] <- line.arg.start [ data.for.Ccode$cn ]  
        }
      }
      
      final.frame <- CNV.fitModel(ncomp,
                                  nind,
                                  hyp = "H0",
                                  data = data.for.Ccode,
                                  logit.offset = offset,
                                  design.matrix.mean,
                                  design.matrix.var,
                                  design.matrix.disease,
                                  mix.model = my.model, 
				  control=control
                                  )
      
      if( final.frame[['status']] == "F" ) next;     
      suc <- suc + 1

      lnL <- getparams(final.frame[['data']])$lnL  
      if( lnL > best.lnL.H0 ){
        best.lnL.H0 <- lnL
        best.model.H0 <- final.frame[['data']]
	best.status <- final.frame[['status']]
      }
    }
 
    if( suc > 0 ){
    	bic <- -2*best.lnL.H0 + nparameters*log(nind)
    	aic <- -2*best.lnL.H0 + 2*nparameters;
    	bics <- append(bics, bic)
    	aics <- append(aics, aic)	    


    	res[['np']][[m]] <- nparameters
    	res[['model']][[m]] <- list( ncomp=ncomp, model.mean=model.mean, model.var=model.var )
    	res[['BIC']][[m]] <- bic		
    	res[['AIC']][[m]] <- aic
        res[['lnL']][[m]] <- best.lnL.H0
    	res[['posteriors']][[m]] <- compact.data.frame(best.model.H0)
	res[['models']][[m]] <- getparams(best.model.H0)
	
	# check for posterior problem
	if( best.status == "C" && ncomp != 1 && test.posterior(frame = compact.data.frame(best.model.H0), ncomp = ncomp) == TRUE ) res[['status']][[m]] = "P"
	else res[['status']][[m]] = best.status

	# update with best model
    	if(method == "BIC"){  res[['selected']] <- order(bics)[1] }
    	else{ res[['selected']] <- order(aics)[1] }	
     }
  }  

  return( res )	
}



test.posterior <- function(frame, ncomp, samples.by.disease = NULL) 
{
  bad.fit <- FALSE
  comp <- c(1:ncomp)
  col <- paste('P', comp, sep='')
  frame$average <- apply(MAR = 1, frame[, col, drop = FALSE], FUN = function (x) {sum(x * comp)})

  for (my.id in levels(frame$batch)) {
    small.frame <- subset(frame, frame$batch == my.id)	   

    if( is.null(samples.by.disease) ){
	small.frame <- small.frame[ order(small.frame$signal),]
	di <- diff(small.frame$average)    	
	should.not.be <- subset(di, di < 0)
	if (sum(should.not.be) < -0.01)  bad.fit <- TRUE
    }
    else{	
    	# split again into cases and controls	
    	for ( t in c(1:2) ){     
	  smaller.frame <- subset(small.frame, small.frame$subject %in% samples.by.disease[[t]] )    

    	  smaller.frame <- smaller.frame[ order(smaller.frame$signal),]
    	  di <- diff(smaller.frame$average)
    	  should.not.be <- subset(di, di < 0)
    	  if (sum(should.not.be) < -0.01)  bad.fit <- TRUE
    	}
     }
  }
  return (bad.fit)
}

### batch will be factor
## signal will be numeric
## disease status is 0/1


CNVtest.binary <- function(signal,
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
                           association.test.strata = NULL,
                           beta.estimated = NULL,
                           start.mean = NULL,
                           start.var = NULL,
			   control=list(tol=1e-5, max.iter = 3000, min.freq=4) ){


  nind <- length(signal)
  batch <- factor(batch)
  if (length(batch) != nind) {stop("Specification of batches does not have the same length as the input signal\n")}

  if (!is.null(association.test.strata)) {
    if (length(association.test.strata) != nind) stop("Specification of strata for association test does not have the same length as the input signal\n")
    association.test.strata <- as.numeric(factor(association.test.strata))
  }
  
  message(paste("Attempting to cluster the data with", nind, "individuals and", ncomp, "components"))
  if ((sd(signal) > 3) || (sd(signal) < 0.33)) {
    warning("The signal you provided to the algorithm has a standard deviation significantly different from 1. 
	To maximize the stability of the underlying clustering algorithm we recommend you normalize this signal to make sure that the standard deviation is indeed 1.")
  }

  
  if (is.null(sample)) {sample <- paste(batch, c(1:nind), sep='_')} else {
    if (length(sample) != nind) {stop("The sample names you have specified do not match the number of data points\n")}
    if (sum ( make.unique(as.character(sample)) != as.character(sample) ) > 0) {
      warning("Sample names not unique. CNVtools will modify them for uniqueness.")
      sample <- make.unique(as.character(sample))
    }
  }



  if (is.null(disease.status) ){
    if(n.H1 > 0){ stop("Must specify disease.status under H1")}	    		
    disease.status <- rep(0, nind )       # Create arbitrary disease.status for fitting under null
  }
  
  sample.split  <- split(x = sample, f = batch)
  signal.split  <- split(x = signal, f = batch)
  disease.split <- split(x = disease.status, f = batch)
  batch.split   <- split(x = batch, f=batch)
  if (!is.null(association.test.strata))  association.test.strata.split <- split(x =  association.test.strata , f=batch) else association.test.strata.split <- NULL


    
  trait.sample.list <- list( sample[which(disease.status == 0)], sample[which(disease.status == 1)] )
  
  ncohort <- length(levels(batch))
  
  data.for.Ccode <- ExpandData (batch = batch.split, trait = disease.split, names = sample.split, signal = signal.split, ncomp = ncomp, association.strata = association.test.strata.split)


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


############ Now let us specify the design matrices, i.e. the actual model for mean and variances
  design.matrix.disease <- model.matrix(data = data.for.Ccode, object = as.formula(model.disease))
  design.matrix.mean <- NULL
  design.matrix.var <- NULL
  special <- c("strata")

  for (design in c('var', 'mean')) {
    
    my.formula <- as.formula(get(paste('model', design, sep = '.')))
    Terms <- terms(my.formula, special, data=data.for.Ccode)
    strats <- attr(Terms, "specials")$strata
    if (!is.null(strats)) {
      
      m <- list()
      m[[1]] <- as.name("model.frame")
      m[[2]] <- Terms
      names(m)[2] <- "formula"
      m[['data']] <- data.for.Ccode
      m <- as.call(c(as.list(m), list(na.action=as.symbol("na.omit"))))
      m <- eval(m)
      
      temps <- untangle.specials(Terms, "strata", 1)
      data.for.Ccode[, paste('strats', design, sep = '.')] <- as.integer(strata(m[, temps$vars], shortlabel=TRUE))
      nstrata <- length(unique(data.for.Ccode[, paste('strats', design, sep = '.')]))
      if (nstrata == 1) assign(x = paste('design.matrix', design, sep = '.'), value = model.matrix(data = data.for.Ccode, object = as.formula(' ~ 1') ))
      else assign(x = paste('design.matrix', design, sep = '.'), value = matrix(0))
      
    } else {
      assign(x = paste('design.matrix', design, sep = '.'), value = model.matrix(data = data.for.Ccode, object = my.formula ))
    }
  }
  model.spec <- 10



############ Now we can start the estimation process
  status.H0 <- ""
  res <- list()
  offset <- rep(0, dim(data.for.Ccode)[1] )
  

  if (n.H0 > 0) {      ################################### First fitting under H0
    best.model.H0 <- 0
    best.lnL.H0 <- -Inf
    best.status.H0 <- ""    

    
    for(i in 1:n.H0 ){
      message(paste("Iteration", i, "under H0"))
      data.for.Ccode <- EM.starting.point(data.for.Ccode)

      for (start.values in c('mean', 'var')) {  ### plug the starting values in the matrix
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
                                  design.matrix.mean,
                                  design.matrix.var,
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
    res$model.H0$nu <- 0

    
    if (output == 'compact') {
      res$posterior.H0 <- compact.data.frame(best.model.H0)
      res$posterior.H0 <- res$posterior.H0[ match(sample, res$posterior.H0$subject),]
      res$status.H0 = best.status.H0
      if( best.status.H0 == "C" && test.posterior(frame = res$posterior.H0, ncomp = ncomp) == TRUE ) res$status.H0 = "P"
    } else res$posterior.H0 <- best.model.H0


  }

  


  if (n.H1 > 0)  {  ################################## Now fitting under H1
    best.model.H1 <- 0
    best.lnL.H1 <- -Inf
    best.status.H1 <- ""
    status.H1 <- ""
    if (!is.null(beta.estimated))    offset <- beta.estimated * data.for.Ccode$cn
    
    for(i in 1:n.H1 ){
      message(paste("Iteration", i, "under H1"))
      if  ((i == 1) & (n.H0 > 0)) data.for.Ccode <- best.model.H0  else {
        data.for.Ccode <- EM.starting.point(data.for.Ccode)
        
        for (start.values in c('mean', 'var')) {  ### plug the starting values in the matrix
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
                                  design.matrix.mean,
                                  design.matrix.var,
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

    if (output == 'compact') {
      res$posterior.H1 <- compact.data.frame(best.model.H1)
      res$posterior.H1 <- res$posterior.H1[ match(sample, res$posterior.H1$subject),]  ## just some reordering
    } else {
      res$posterior.H1 <- best.model.H1
    }

    # Set the status
    res$status.H1 <- best.status.H1 
    if(best.status.H1 == "C" && test.posterior(frame = res$posterior.H1, ncomp = ncomp, samples.by.disease = trait.sample.list) == TRUE ) res$status.H1 <- "P"
    res$model.H1$nu <- 0
  }
  
  return(res)
}






### batch will be factor
## signal will be numeric
## qt is quantitative trait
CNVtest.qt <- function(signal, batch, sample = NULL, qt = NULL, ncomp, n.H0=5, n.H1=0, 
		       model.mean = '~ strata(cn)',
                       model.var = '~ strata(cn)',
		       model.qt = '~ cn',
                       beta.estimated = NULL,
                       start.mean = NULL,
                       start.var = NULL,
		       control=list(tol=1e-5, max.iter = 3000, min.freq=4) ){
 
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

  
############ Now let us specify the design matrices, i.e. the actual model for mean and variances
  design.matrix.qt <- model.matrix(data = data.for.Ccode, object = as.formula(model.qt))
  design.matrix.var <- NULL
  design.matrix.mean <- NULL
  special <- c("strata")

  for (design in c('var', 'mean')) {
    
    my.formula <- as.formula(get(paste('model', design, sep = '.')))
    Terms <- terms(my.formula, special, data=data.for.Ccode)
    strats <- attr(Terms, "specials")$strata
    if (!is.null(strats)) {
      
      m <- list()
      m[[1]] <- as.name("model.frame")
      m[[2]] <- Terms
      names(m)[2] <- "formula"
      m[['data']] <- data.for.Ccode
      m <- as.call(c(as.list(m), list(na.action=as.symbol("na.omit"))))
      m <- eval(m)
      
      temps <- untangle.specials(Terms, "strata", 1)
      data.for.Ccode[, paste('strats', design, sep = '.')] <- as.integer(strata(m[, temps$vars], shortlabel=TRUE))
      nstrata <- length(unique(data.for.Ccode[, paste('strats', design, sep = '.')]))
      if (nstrata == 1) assign(x = paste('design.matrix', design, sep = '.'), value = model.matrix(data = data.for.Ccode, object = as.formula(' ~ 1') ))
      else assign(x = paste('design.matrix', design, sep = '.'), value = matrix(0))
      
    } else {
      assign(x = paste('design.matrix', design, sep = '.'), value = model.matrix(data = data.for.Ccode, object = my.formula ))
    }
  }
  model.spec <- 10

  
################################### First fitting under H0
  status.H0 <- ""
  res <- list()
  offset <- rep(0, dim(data.for.Ccode)[1] )

  

  if (n.H0 > 0) {
    best.model.H0 <- 0
    best.lnL.H0 <- -Inf
    best.status.H0 <- ""    

    
    for(i in 1:n.H0 ){
      message(paste("Iteration", i, "under H0"))
      data.for.Ccode <- EM.starting.point(data.for.Ccode)
      if (!is.null(start.mean)) data.for.Ccode$mean = start.mean[ data.for.Ccode$cn ]
      if (!is.null(start.var)) data.for.Ccode$var = start.var[ data.for.Ccode$cn ]            
      
      final.frame <- CNV.fitModel(ncomp, nind, hyp = "H0", data.for.Ccode, logit.offset = offset,
                                  design.matrix.mean,
                                  design.matrix.var,
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
    res$model.H0$nu <- 0
 }

  
################################## Now fitting under H1

  if (n.H1 > 0)  {
    best.model.H1 <- 0
    best.lnL.H1 <- -1000000000
    best.status.H1 <- ""
    status.H1 <- ""
    if (!is.null(beta.estimated))    offset <- beta.estimated * data.for.Ccode$cn

    for(i in 1:n.H1 ){
      message(paste("Iteration", i, "under H1"))
      if  ((i==1) & (n.H0 > 0)) {data.for.Ccode <- best.model.H0}
      else{
        data.for.Ccode <- EM.starting.point(data.for.Ccode)
        if (!is.null(start.mean)) data.for.Ccode$mean = start.mean[ data.for.Ccode$cn ]
        if (!is.null(start.var)) data.for.Ccode$var = start.var[ data.for.Ccode$cn ]            
      }

      final.frame <- CNV.fitModel(ncomp, nind, hyp = "H1", data = data.for.Ccode, logit.offset = offset,
                                  design.matrix.mean,
                                  design.matrix.var,
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
    #if(  best.status.H1 == "C" && test.posterior(res[['posterior.H1']], ncomp) == TRUE ) res[['status.H1']] <- "P"	   
    res$model.H1$nu 
  }
	  
  return(res)
}



apply.ldf <- function(full.signal, posterior) 
{
  posterior <- posterior [  dimnames(full.signal)[[1]],]
  if (sum(dimnames(full.signal[[1]]) != dimnames(posterior)[[1]]) > 0) stop("Names are not matching")
  
  can <- cancor(x = full.signal, y = posterior)
  ldf <- as.numeric( full.signal %*%  can$xcoef[,1])
  return (ldf/sd(ldf))
}

apply.pca <- function (matrix.signal) 
{
  # return the mean if there are only two probes
  if( dim(matrix.signal)[2] <= 2 ){
    m <- apply(matrix.signal, MAR=1, FUN=mean)    
    return(m)    
  }
  pca <- prcomp(matrix.signal, scale=TRUE)$x[,1]; 
  return (pca/sd(pca))
}









getQualityScore <- function(posterior){

  Qvec <- c()
  for (coh in levels (posterior$batch)) {
    posterior2 <- subset(posterior, posterior$batch == coh)
    sds <- tapply(posterior2$signal, FUN=sd, INDEX = posterior2$cn)
    means <- tapply(posterior2$signal, FUN=mean, INDEX = posterior2$cn)
    freq <- table(posterior2$cn)/length(posterior2$cn)

    ########removes some missing categories
    sds <- subset(sds, table(posterior2$cn) > 4)
    means <- subset(means, table(posterior2$cn) > 4)
    freq <- subset(freq, table(posterior2$cn) > 4)

    
    l <- length(means)
    
    if (l == 1) {return (NA)} else {   ############an exception for the one component situation
      dmeans <- abs(means[ 1: (l-1) ] - means[ 2:l ])    ##the differences of means for pairs of adjacent clusters
      av.sds <- (freq[1:(l-1)] * sds [ 1:(l-1) ]  + freq[ 2:l ] * sds [ 2:l ])/ ( freq[1:(l-1)] + freq[ 2:l ])   #the average std dev for pairs of adjacent clusters
      weights <- freq[1:(l-1)]*freq[ 2:l ]   ###weights for each pair of clusters
      Q <- sum(weights*dmeans/av.sds)/sum(weights)    ##the quality score
      Qvec <- append(Qvec, values = Q )
    }
  }
  return (min(Qvec))
}







###############
test <- function() {
  dyn.load("../src/CNVtools.so")
  load("../CNVtools/data/A112.RData")

  raw.signal <- as.matrix(A112[, -c(1,2)])
  dimnames(raw.signal)[[1]] <- A112$subject
  
  mean.signal <- apply(raw.signal, MAR=1, FUN=mean)
  pca.signal <- apply.pca(raw.signal)
  
  pdf("mean_pca_signal.pdf", width=10, height=5)
  par(mfrow=c(1,2))
  hist(mean.signal, breaks=50)
  hist(pca.signal, breaks=50)  
  dev.off()

  
  batches <- factor(A112$cohort)
  sample <- factor(A112$subject)
  trait <- ifelse( A112$cohort == '58C', 0, 1)

  # Determine number of components
  # ncomp <- CNVtest.select.model(signal = pca.signal, batch = batches, sample = sample, n.H0 = 5)
  ncomp <- 3


  start.mean.pca <- c(-1, 0.1, 1.8)
  start.mean.ldf <- c(-2, -0.5, 0.9)
  
  res.pca <- CNVtest.binary ( signal = pca.signal, sample = names(pca.signal), batch = batches, disease.status = trait, ncomp = ncomp, n.H0 = 1, n.H1 = 0, model.var = '~ cn', start.mean = start.mean.pca)
    
  
  pdf("look_at_fit_pca.pdf", width=10, height=5)
  par(mfrow=c(1,2))
  cnv.plot(res.pca$posterior.H0, main='Example of PCA fit')
  dev.off()
  

  ############################## Now use the ldf transformation
  posterior <- as.matrix((res.pca$posterior.H0)[, c('P1','P2', 'P3')])
  dimnames(posterior)[[1]] <- (res.pca$posterior.H0)$subject
  ldf.signal <- apply.ldf(raw.signal, posterior)
  print("LDF has been computed")


  
  ############## And now do a proper association test using the ldf transformed data
  res.ldf <- CNVtest.binary ( signal = ldf.signal, sample = names(pca.signal), batch = batches, disease.status = trait, ncomp = ncomp, start.mean = start.mean.ldf)
  
  pdf("look_at_fit_ldf.pdf", width=10, height=5)
  par(mfrow=c(1,2))
  cnv.plot(res.ldf$posterior.H0, main='Example of ldf fit')
  dev.off()
	
  message(paste("ldf status : ", res.ldf$status,"LR :",-2*(res.ldf$model.H0$lnL - res.ldf$model.H1$lnL),sep=' '))
  
  return(res.ldf)
}



test <- function() {
  dyn.load("../src/CNVtools.so")
  source("t_fit.R")
  
  set.seed(0)
  
  load("../CNVtools/data/A112.RData")
  raw.signal <- as.matrix(A112[, -c(1,2)])
  dimnames(raw.signal)[[1]] <- A112$subject

  A112$signal <- apply.pca(raw.signal)
  data <- A112
  data$batch <- data$cohort
  data$disease.status <- ifelse (data$cohort == "NBS", 1, 0)
  
  ncomp <- 3
  n.H0 <- 1
  n.H1 <- 1
  model.mean <- "~ strata(batch,cn)"
  model.var <- "~ strata(cn)"
  model.disease <- "~ cn"

  results <- CNVtest.select.model(signal=data$signal, batch = data$batch, sample = data$sample, n.H0 = 3, method="BIC")
  
  my.pca.fit <- CNVtest.binary.T(signal = data$signal,
                                 batch = data$batch,
                                 sample = data$sample,
                                 disease.status = data$disease.status,
                                 ncomp = ncomp, 
                                 n.H0 = n.H0,
                                 n.H1 = n.H1,
                                 model.mean = model.mean,
                                 model.var = model.var,
                                 model.disease = model.disease)

   
  pdf("test.pdf")
  cnv.plot(my.pca.fit$posterior.H0, batch = '58C', col = 'red', breaks = 40)
  cnv.plot(my.pca.fit$posterior.H1, batch = '58C', col = 'red', breaks = 40)
  dev.off()


  return(my.pca.fit)
}

#mtest <- test()
