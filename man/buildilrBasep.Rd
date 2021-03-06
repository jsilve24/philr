% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weighted_ILR.R
\name{buildilrBasep}
\alias{buildilrBasep}
\title{Weighted ILR Contrast Matrix}
\usage{
buildilrBasep(W, p)
}
\arguments{
\item{W}{sequantial binary partition matrix (e.g., signary matrix; output of
\link{phylo2sbp})}

\item{p}{weights (should not be closed)}
}
\value{
matrix
}
\description{
Weighted ILR Contrast Matrix
}
\examples{
p <- seq(.1,1,by=.2)
tr <- named_rtree(5)
sbp <- phylo2sbp(tr)
buildilrBasep(sbp, p)
}
\references{
J. J. Egozcue, V. Pawlowsky-Glahn (2016) \emph{Changing the Reference
  Measure in the Simplex and its Weighting Effects}. Austrian Journal of
  Statistics 45(4):25-44
}
\author{
Justin Silverman (adapted from compositions::gsi.buildilrBase)
}
