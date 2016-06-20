
# Generally Helpful for Validation ----------------------------------------

# check for zeros and throw error if present
# target is just a name for the warning message
check.zeroes <- function(x, target){
  if (any(x==0)){
    warning(paste(target, "should not contain zeroes"))
  }
}

# ACCESSOR FUNCTIONS FROM NAMES TO NUMBERS --------------------------------

# Main accessor of node number through coordinate name
c.to.nn <- function(tr, c){
  return(which(tr$node.label==c)+ape::Ntip(tr))
}

# Main accessor of node number through tip name
t.to.nn <- function(tr, t){
  return(which(tr$tip.label==t))
}

# Main accessor of node or tip name thorugh node number
nn.to.name <- function(tr, nn){
  n <- ape::Ntip(tr)
  if (nn <= n)return(tr$tip.label[nn])
  return(tr$node.label[nn-n])
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
#' @examples
#' library(phyloseq)
#' data(CSS)
#' df <- t(otu_table(CSS))
#' df <- df + 0.65   # add a small pseudocount
#' tree <- phy_tree(CSS)
#' df.philr <- philr(df, tree, part.weights='anorm.x.gm.counts',
#'                   ilr.weights='blw.sqrt', return.all=FALSE, n_cores=1)
#' head(convert_to_long(df.philr, get_variable(CSS, 'BODY_SITE'))
convert_to_long <- function(df, labels){
  #TODO: expand to labels can be an entire dataframe of metadata to add on
  coord.names <- colnames(df)
  df.long <- as.data.frame(df)
  df.long$sample <- rownames(df.long)
  df.long$labels <- labels
  df.long <- tidyr::gather(df.long, coord, value, -labels, -sample)
  return(df.long)
}
