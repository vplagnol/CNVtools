\name{get.model.spec}
\alias{get.model.spec}
\title{
Get model specifications (internal function)
}
\description{
Internal function to parse the formulas and extract codes for the model, as well as number of parameters
}
\usage{
get.model.spec(model.component, model.mean, model.var, model.nu, design.matrix.mean, design.matrix.variance, ncomp, nbatch)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model.component}{
%%     ~~Describe \code{model.component} here~~
}
  \item{model.mean}{
%%     ~~Describe \code{model.mean} here~~
}
  \item{model.var}{
%%     ~~Describe \code{model.var} here~~
}
  \item{model.nu}{
%%     ~~Describe \code{model.nu} here~~
}
  \item{design.matrix.mean}{
%%     ~~Describe \code{design.matrix.mean} here~~
}
  \item{design.matrix.variance}{
%%     ~~Describe \code{design.matrix.variance} here~~
}
  \item{ncomp}{
%%     ~~Describe \code{ncomp} here~~
}
  \item{nbatch}{
%%     ~~Describe \code{nbatch} here~~
}
}
\value{
A list with two components: model code and number of parameters of the model.
}
\author{
Vincent Plagnol and Chris Barnes
}
