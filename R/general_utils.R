
# Generally Helpful for Validation ----------------------------------------

# check for zeros and throw error if present
# target is just a name for the warning message
check.zeroes <- function(x, target){
  if (any(x==0)){
    warning(paste(target, "should not contain zeroes"))
  }
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

# Main accessor of node number through coordinate name
# depricated - use name.to.nn
c.to.nn <- function(tr, c){
  .Deprecated('name.to.nn')
  return(which(tr$node.label==c)+ape::Ntip(tr))
}

# Main accessor of node number through tip name
# depricated - use name.to.nn
t.to.nn <- function(tr, t){
  .Deprecated('name.to.nn')
  return(which(tr$tip.label==t))
}




# Generally Helpful for Plotting ------------------------------------------

#' Converts wide format ILR transformed data to long format
#'
#' Converts wide format ILR transformed data (see \code{\link{philr}}) to long format
#' useful in various plotting functions
#' (e.g., \code{\link{plot_density_breakdown_wtree}})
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
#' @importFrom tidyr gather
#' @examples
#' data(CSS)
#' df <- CSS$otu.table + 0.65   # add a small pseudocount
#' tree <- CSS$phy.tree
#' df.philr <- philr(df, tree, part.weights='anorm.x.gm.counts',
#'                   ilr.weights='blw.sqrt', return.all=FALSE, n_cores=1)
#' head(convert_to_long(df.philr, CSS$sample.data$BODY_SITE))
convert_to_long <- function(df, labels){
  coord.names <- colnames(df)
  df.long <- as.data.frame(df)
  df.long$sample <- rownames(df.long)
  df.long$labels <- labels
  df.long <- gather(df.long, coord, value, -labels, -sample)
  return(df.long)
}


# Helpful for generating examples -----------------------------------------

#' Generate random tree with named internal nodes
#'
#' Internal nodes are named by numbering and adding the prefix 'n'. This
#' function is laregly for use in examples throughout this package.
#'
#' @inheritParams ape::rtree
#'
#' @return An object of class "phylo"
named_rtree <- function(n){
  tr <- rtree(n)
  makeNodeLabel(tr, 'number', 'n')
}
