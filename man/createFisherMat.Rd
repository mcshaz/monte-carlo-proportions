\name{createFisherMat}
\alias{createFisherMat}
\docType{package}
\title{
Create n x n matrix of p values for proportional outcome counts
}
\description{
Create n x n matrix of p values for proportional outcome counts
}
\usage{
createFisherMat(allocationsPerArm)
}
\arguments{
  \item{allocationsPerArm}{an integer for half the trial number}
}
\value{
A matrix where each cell represents the fisher exact p value for c(j, allocationsPerArm-j, i, allocationsPerArm-i)
}
\examples{
\dontrun{
createFisherMat(8)
}
}
