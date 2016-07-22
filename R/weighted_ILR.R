#' miniclo
#'
#' small function to close (aka normalize by proportions, aka total sum scaling)
#' a dataset to a constant \code{k} (usually taken to be 1). After closure the row sums
#' of the dataset should sum to \code{k}.
#'
#' @param c dataset to be closed
#' @param k closure constant
#'
#' @return matrix (if c is a vector or matrix) or data.frame (if c is a data.frame)
#' @export
#'
#' @examples
#' c <- matrix(c(1,2,3,1,2,3,1,2,3), nrow = 3, byrow=TRUE)
#' miniclo(c)
#' miniclo(c, k=2)
miniclo <- function(c,k=1){
  check.zeroes(k, 'closure constant(k)')
  if(is.vector(c)) c <- matrix(c, nrow=1)
  (c/rowSums(c))*k
}


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
#' @examples
#' p <- seq(.1,1,by=.2)
#' c <- t(rmultinom(10,100,c(.1,.6,.2,.3,.2))) + 0.65   # add a small pseudocount
#' x <- miniclo(c)
#' shiftp(x, p)
shiftp <- function(x, p){
  check.zeroes(p, "weights (p)")
  if(is.vector(x)) x <- matrix(x, nrow=1)
  y <- x / outer(rep(1,nrow(x)), p)
  y
}



#' Weighted Geometric Means of Rows
#'
#' Calculates weighted geometric mean (see references). Note
#' if \code{p=rep(1, nrow(y))} (default) then this is just the
#' geometric mean of rows.
#'
#' @param y shifted data matrix (e.g., output of \link{shiftp})
#' @param p weights (should not be closed)
#'
#' @return vector (weighted geometric mean of rows)
#' @references J. J. Egozcue, V. Pawlowsky-Glahn. \emph{Changing the Reference
#'   Measure in the Simplex and its Weighting Effects}. Austrian Journal of
#'   Statistics, 2016
#'
#' @seealso gp.colMeans
#' @examples
#' p <- seq(.1,1,by=.2)
#' c <- t(rmultinom(10,100,c(.1,.6,.2,.3,.2))) + 0.65   # add a small pseudocount
#' x <- miniclo(c)
#' y <- shiftp(x, p)
#' gp.rowMeans(y, p)
gp.rowMeans <- function(y, p=rep(1, nrow(y))){
  sp <- sum(p) # as given in text on page 4
  exp(1/sp*rowSums(log(y)%*%diag(p)))
}


#' Geometric Means of Columns
#'
#' Calculates geometric mean of columns. Does not calculate
#' WEIGHTED geometric means (vs. \link{gp.rowMeans})
#'
#' @param x matrix or vector
#'
#' @return vector (geometric mean of columns)
#'
#' @seealso gp.rowMeans
#' @examples
#' g.colMeans(rbind(c(2,4,4), c(2,4,4)))
g.colMeans <- function(x){
  if(is.vector(x)) x <- matrix(x, nrow=1)
  exp(colMeans(log(x)))
}


#equation 7
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
  return(stats::as.dist(tmp))
}

#' Weighted CLR Transform
#'
#' @param y shifted data matrix (e.g., output of \link{shiftp})
#' @param p weights (should not be closed)
#' @return matrix
#' @details Note that this function will close the dataset \code{y} to 1.
#' @export
#'
#' @author Justin Silverman
#' @references J. J. Egozcue, V. Pawlowsky-Glahn. \emph{Changing the Reference
#'   Measure in the Simplex and its Weighting Effects}. Austrian Journal of
#'   Statistics, 2016
#' @examples
#' p <- seq(.1,1,by=.2)
#' c <- t(rmultinom(10,100,c(.1,.6,.2,.3,.2))) + 0.65   # add a small pseudocount
#' x <- miniclo(c)
#' y <- shiftp(x, p)
#' clrp(y, p)
clrp <- function(y,p){
  check.zeroes(p, 'weights(p)')
  check.zeroes(y, 'dataset')

  if(is.vector(y)) y <- matrix(y, nrow=1)
  y <- miniclo(y)
  log(y/gp.rowMeans(y,p))
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
#' @examples
#' # Weights
#' p <- seq(.1,1,by=.2)
#'
#' # Shifted Composition
#' c <- t(rmultinom(10,100,c(.1,.6,.2,.3,.2))) + 0.65   # add a small pseudocount
#' x <- miniclo(c)
#' y <- shiftp(x, p)
#'
#' # Contrast Matrix
#' tr <- named_rtree(5)
#' sbp <- phylo2sbp(tr)
#' V <- buildilrBasep(sbp, p)
#'
#' ilrp(y, p, V)
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
#' @examples
#' p <- seq(.1,1,by=.2)
#' tr <- named_rtree(5)
#' sbp <- phylo2sbp(tr)
#' buildilrBasep(sbp, p)
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
