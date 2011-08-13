

get.data <- function() {
  if (!file.exists("../CNVtools/data")) system ("mkdir ../CNVtools/data")
  for (name in c('A112')) {
    data <- read.table(paste("../examples/", name, ".dat", sep=''), header=TRUE)
    data <- subset(data, data$cohort %in% c('58C', 'NBS'))
    assign(name, value=data)
    save(list=name, file=paste("../CNVtools/data/", name, ".RData", sep=''))
  }

}

main <- function() {
  source("CNVtools.r")
  source("t_fit.R")
  list.fun <- c('EM.starting.point', 'ExpandData', 'CNV.fitModel', 'CNVtest.binary', 'CNVtest.qt', 'getparams', 'compact.data.frame', 'cnv.plot', 'qt.plot', 'apply.ldf', 'apply.pca', 'test.posterior','CNVtest.select.model','getQualityScore')
  list.fun.T <- c('CNVtest.binary.T', 'CNVtest.qt.T', 'get.model.spec')


  package.skeleton(name="CNVtools", list=c(list.fun, list.fun.T), path='../', force=TRUE)
  if (!file.exists("../CNVtools/src")) dir.create(path = "../CNVtools/src")
  if (!file.exists("../CNVtools/inst")) dir.create(path = "../CNVtools/inst")
  if (!file.exists("../CNVtools/inst/doc")) dir.create(path = "../CNVtools/inst/doc")
  if (!file.exists("../CNVtools/inst/doc/fig")) dir.create(path = "../CNVtools/inst/doc/fig")
  
  system("cp  ../src/*.cpp ../src/*.h  ../src/*.c ../CNVtools/src/")
  system("cp doc/*.Rd ../CNVtools/man/")
  file.copy(from = "vignette/vignette.Rnw", to ="../CNVtools/inst/doc/CNVtools-vignette.Rnw", overwrite = TRUE)
  file.copy(from = c("doc/DESCRIPTION", "doc/NAMESPACE"), to =  "../CNVtools/", overwrite = TRUE)
  file.copy(from  = "zzz.R", to="../CNVtools/R/", overwrite = TRUE)
}


main()
get.data()

warnings()
