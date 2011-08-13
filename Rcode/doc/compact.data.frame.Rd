\name{compact.data.frame}
\alias{compact.data.frame}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compacts the expanded data frame format needed by our fitting
  procedure into more compact and user friendly version}
\description{
 Small internal routine returning a more compact and user friendly version of the output of the fitting algorithm. 
}
\usage{
compact.data.frame(full.frame)
}
\arguments{
  \item{full.frame}{An expanded data frame (one point per data point and
    per component in the fit, ie. 1,000 individuals fitted on three
    components would have 3,000 rows. }
}
\details{
  This function should be invisible to most users and is part of the EM
  fitting procedure.
}
\value{
  A data frame in a compact version, with one row per data point and one
  column for each component: P1, P2, P3 in the three component case for
  the probabilities for the calls to be equal to 1,2 or 3.
}
\author{Vincent Plagnol \email{vincent.plagnol@cimr.cam.ac.uk} and Chris Barnes \email{christopher.barnes@imperial.ac.uk} }
