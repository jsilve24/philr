#' Shift data to origin given by p
#'
#' Shift must be applied before transformation
#'
#' @param x closed compositional data matrix (or vector)
#' @param p weights (should not be closed)
#' @return shifted data matrix \code{y} (no closure is applied) rows are
#'   samples, columns are parts
#' @author Justin Silverman & J. J. Egozcue
#' @references J. J. Egozcue, V. Pawlowsky-Glahn. \emph{Changing the Reference
#'   Measure in the Simplex and its Weighting Effects}. Austrian Journal of
#'   Statistics, 2016
#' @export
shiftp <- function(x, p){
  check.zeroes(p, "weights (p)")
  if(is.vector(x)) x <- matrix(x, nrow=1)
  y <- x / outer(rep(1,nrow(x)), p)
  y
}

# as given in equation 2
gp.mean <- function(y,p){
  sp <- sum(p) # as given in text on page 4
  exp(1/sp*rowSums(log(y)%*%diag(p)))
}

#equation 6
inner.prodp <- function(y1,y2,p){
  sum(p*clrp(y1,p)*clrp(y2,p))
}

normp <- function(y,p){
  sqrt(inner.prodp(y,y,p))
}

# pg 5
# for a pair of samples
distp.pair <- function(y1,y2,p){
  sqrt(sum(p*(clrp(y1,p)-clrp(y2,p))^2))
}

#pg 5
# This function is for a matrix where each row is a sample
distp <- function(y,p){
  n <- dim(y)[1]
  tmp <- matrix(rep(NA,n),ncol = n, nrow = n)
  for (i in 1:n){
    for (j in 1:i) {
      if (i!=j){
        tmp[i,j] <- distp.pair(y[i,],y[j,],p)
      }
    }
  }
  return(as.dist(tmp))
}

#' Weighted CLR Transform
#'
#' @param y shifted data matrix (e.g., output of \link{shiftp})
#' @param p weights (should not be closed)
#' @return matrix
#' @export
#'
#' @author Justin Silverman
#' @references J. J. Egozcue, V. Pawlowsky-Glahn. \emph{Changing the Reference
#'   Measure in the Simplex and its Weighting Effects}. Austrian Journal of
#'   Statistics, 2016
clrp <- function(y,p){
  check.zeroes(p, 'weights(p)')
  check.zeroes(y, 'dataset')

  if(is.vector(y)) y <- matrix(y, nrow=1)
  # TODO: Requires compositions::clo
  y <- compositions::clo(y)
  log(y/gp.mean(y,p))
}

#' Weighted ILR Transform
#'
#' Calculated using weighted CLR transform (\link{clrp})
#'
#' @param y shifted data matrix (e.g., output of \link{shiftp})
#' @param p weights (should not be closed)
#' @param V weighted contrast matrix (e.g., output of \link{buildilrBasep})
#'
#' @return matrix
#' @export
#' @author Justin Silverman
#' @references J. J. Egozcue, V. Pawlowsky-Glahn. \emph{Changing the Reference
#'   Measure in the Simplex and its Weighting Effects}. Austrian Journal of
#'   Statistics, 2016
ilrp <- function(y,p,V){
  check.zeroes(p, 'weights(p)')
  check.zeroes(y, 'dataset')

  # if(is.vector(y)) y <- matrix(y, nrow=1) # run in clrp anyways
  clrp(y,p)%*%diag(p)%*%V
}


#' Weighted ILR Contrast Matrix
#'
#' @param W sequantial binary partition matrix (e.g., signary matrix; output of
#'   \link{phylo2sbp})
#' @param p weights (should not be closed)
#' @return matrix
#' @export
#' @author Justin Silverman (adapted from \link[compositions]{gsi.buildilrBase})
#' @references J. J. Egozcue, V. Pawlowsky-Glahn. \emph{Changing the Reference
#'   Measure in the Simplex and its Weighting Effects}. Austrian Journal of
#'   Statistics, 2016
buildilrBasep <- function(W,p){
    check.zeroes(p, 'weights (p)')

    W = as.matrix(W)
    nc = ncol(W)
    D = nrow(W)
    isPos = (W > 0)
    isNeg = (W < 0)
    nPos = colSums(isPos*p)
    nNeg = colSums(isNeg*p)
    normConst = sqrt(nPos*nNeg/(nPos+nNeg))
    negPortion = -1/nNeg*t(isNeg)*normConst
    posPortion = 1/nPos*t(isPos)*normConst
    W = negPortion + posPortion
    return(t(W))
}
