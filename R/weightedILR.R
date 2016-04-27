# as given in equation 2
gp.mean <- function(y,p){
  sp <- sum(p) # as given in text on page 4
  exp(1/sp*rowSums(log(y)%*%diag(p)))
}

# as given in equation 2

# requires compositions clo
clrp <- function(y,p){
  y <- compositions::clo(y)
  log(y/gp.mean(y,p))
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

##### Weighted ILR Transform #####
# eq. 9
# x is the dataset
# p is the weightings
# V is the contrast matrix
ilrp <- function(x,p,V){
  clrp(x,p)%*%diag(p)%*%V
}


# eq. 7 and templated from gsi.buildilrBase from
# compositions package
buildilrBasep <- function(W,p){
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
