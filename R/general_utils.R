
# Generally Helpful for Validation ----------------------------------------

# check for zeros and throw error if present
# target is just a name for the warning message
check.zeroes <- function(x, target){
  if (any(x==0)){
    warning(paste(target, "should not contain zeroes"))
  }
}

# convert vector to row vector.
vec_to_mat <- function(x){
  if (is.vector(x)) {
    n <- names(x)
    x <- matrix(x, nrow = 1)
    colnames(x) <- n
  }
  x
}

# ACCESSOR FUNCTIONS FROM NAMES TO NUMBERS --------------------------------

#' Convert between node/tip labels and integer node numbers
#'
#' Useful if you want to convert between node labels (\code{c}), tip labels
#' (\code{t}) and the internal integer number that identifies that node
#' (\code{nn}). Particularly for use with plotting libraries.
#'
#' @param tr object of type \code{phylo}
#' @param x vector of numerics or characters
#' @return vector
#' @name name_nodenumber_conversion
#' @examples
#' tr <- named_rtree(5)
#' name.to.nn(tr, 'n1')
#' name.to.nn(tr,c('n1','n2','t1'))
#' nn.to.name(tr, 1:9)
NULL

#' @rdname name_nodenumber_conversion
#' @export
nn.to.name <- function(tr, x){
  if(!is.numeric(x)) stop('node numbers must be numeric')
  labels <- c(tr$tip.label, tr$node.label)
  labels[x]
}

#' @rdname name_nodenumber_conversion
#' @export
name.to.nn <- function(tr, x){
  if(!is.character(x)) stop('node/tip names (x) should be a character vector')
  labels <- c(tr$tip.label, tr$node.label)
  match(x, labels)
}


# Generally Helpful for Plotting ------------------------------------------

#' Converts wide format ILR transformed data to long format
#'
#' Converts wide format ILR transformed data (see \code{\link{philr}}) to long format
#' useful in various plotting functions where long format data is required.
#'
#' @param df PhILR transformed data in wide format (samples by balances) (see \code{\link{philr}})
#' @param labels vector (of length \code{nrow(df)}) with labels to group samples by
#' @return \code{df} in long format with columns
#' \itemize{
#' \item sample
#' \item labels
#' \item coord
#' \item value
#' }
#' @export
#' @importFrom tidyr gather_
#' @examples
#' tr <- named_rtree(5)
#' df <- t(rmultinom(10,100,c(.1,.6,.2,.3,.2))) + 0.65   # add a small pseudocount
#' colnames(df) <- tr$tip.label
#'
#' df.philr <- philr(df, tree=tr, part.weights='uniform',
#'       ilr.weights='uniform', return.all=FALSE)
#' convert_to_long(df.philr, rep(c('a','b'), 5))
convert_to_long <- function(df, labels){
  coord.names <- colnames(df)
  df.long <- as.data.frame(df)
  df.long$sample <- rownames(df.long)
  df.long$labels <- labels
  df.long <- gather_(df.long, "coord", "value",
                     setdiff(names(df.long), c("labels", "sample")))
  return(df.long)
}


# Helpful for generating examples -----------------------------------------

#' Generate random tree with named internal nodes
#'
#' Internal nodes are named by numbering and adding the prefix 'n'. This
#' function is laregly for use in examples throughout this package.
#'
#' @inheritParams ape::rtree
#' @importFrom ape rtree makeNodeLabel
#' @export
#'
#' @return An object of class "phylo"
#' @examples
#' named_rtree(5)
named_rtree <- function(n){
  tr <- rtree(n)
  makeNodeLabel(tr, 'number', 'n')
}
