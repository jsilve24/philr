#' Calculate Branch Length Weightings for ILR Coordinates
#'
#' Calculates the weightings for ILR coordinates based on branch lenghts of
#' a phylogenetic tree via a few different methods (see details).
#' @aliases calc.blw
#' @param tree a \code{phylo} class tree object that is binary (see \code{\link[ape]{multi2di}})
#' @param method options include: (default) \code{'sum.children'} and \code{'mean.descendants'}
#' see details for more information.
#' @return vector of weightings for ILR coordinates produced via specified method.
#' @details
#' ILR balances built from a binary partition of a phylogenetic tree
#' can be imbued with branch length information. This function is helpful in
#' calculating those weightings.\cr \cr
#' There are a number of methods for calculating these weightings, the default
#' \code{'sum.children'} calculates the weighting for a given balance as the sum
#' of its two direct children's branch length. An alternative that has been as yet less
#' studied is \code{'mean.descendants'} to calculate the weighting for a given balance
#' as the sum of its two direct children's branch lengths PLUS for each child the average
#' distance from it to its descendant tips. \cr \cr
#' \emph{Note:} That some trees contain tips with branch lengths of zero length. This can result
#' in that tip being unreasonably downweighted, as such this function automatically
#' adds a small pseudocount to those tips with zero length (equal to the smallest non-zero)
#' branch length on the tree.
#' @author Justin Silverman
#' @export
#' @seealso \code{\link{philr}}
#' @examples
#' data(CSS)
#' tree <- CSS$phy.tree
#' calculate.blw(tree, method='sum.children')[1:10]
#' calculate.blw(tree, method='mean.descendants')[1:10]
calculate.blw <- function(tree, method='sum.children'){
    nTips = ape::Ntip(tree)

    # Note that some of the terminal branches of the tree have zero length.
    # In these cases will replace those zero values with
    # the minimum branch length (greater than 0) found in the tree.
    # Note this is like adding a psudocount to branch lengths.
    min.nonzero <- min(tree$edge.length[tree$edge.length>0 & !is.na(tree$edge.length)])
    # Find the edges that end in tips and have zero length
    tip.edges.zero <- (tree$edge[,2] <= nTips) & (tree$edge.length == 0)
    # print warning/note statement
    n.replace <- sum(tip.edges.zero)
    if (n.replace >0 ){
        warning(paste('Note: a total of', sum(tip.edges.zero), 'tip edges with zero length',
        'have been replaced with a small pseudocount of the minimum',
        'non-zero edge length (',min.nonzero,').',sep=" "))
        tree$edge.length[tip.edges.zero] <- min.nonzero
    }

    if (method=='sum.children')return(blw.sum.children(tree))
    else if (method=='mean.descendants')return(blw.mean.descendants.sum.children(tree))
}

# Now calculate Branch Length Weighting as the sum of a nodes
#  two (direct) child's weights
blw.sum.children <- function(tree){
    nTips =  nTips = ape::Ntip(tree)
    X <- phangorn::Children(tree, (nTips+1):(nTips+tree$Nnode))
    fun <- function(x, el)sum(el[x])
    EL <- numeric(max(tree$edge))
    EL[tree$edge[,2]] <- tree$edge.length
    res <- sapply(X, fun, EL)
    if(!is.null(tree$node.label)) names(res) <- tree$node.label
    res
}

# Calculates the sum of the children's nodes average distance to descendant tips
# Which should be a slightly better calculation than MDTT when we are
# primarily interested in weightings on ILR coordinates
blw.mean.descendants.sum.children <- function(tree){
  nTips = ape::Ntip(tree)
  X <- phangorn::Children(tree, (nTips+1):(nTips+tree$Nnode))

  # Children's average branch length to tips (zero for tips)
  BMD <- mean_dist_to_tips(tree)
  names(BMD) <- NULL
  BMD <- c(numeric(nTips), BMD)

  # Each child's edge length
  EL <- numeric(max(tree$edge))
  EL[tree$edge[,2]] <- tree$edge.length

  fun <- function(x, el, bmd)sum(el[x] + bmd[x])
  res <- sapply(X, fun, EL, BMD)
  if(!is.null(tree$node.label)) names(res) <- tree$node.label
  return(res)
}

#' Mean distance from internal nodes to descendant tips
#'
#' Calculates the mean distance from each internal node to its descendant tips
#'
#' @inheritParams calculate.blw
#' @details This is a function used by \code{\link{calculate.blw}} when
#'   \code{method='mean.descendants'}, there this function is called twice, once
#'   for each direct child of a given internal node and the results are summed
#'   for each node.
#' @export
mean_dist_to_tips <- function(tree){
    nTips = ape::Ntip(tree)

    dist <- ape::dist.nodes(tree)
    node.numbers <- (nTips+1):(nTips+tree$Nnode)
    chld.tips <- phangorn::Descendants(tree, (nTips+1):(nTips+tree$Nnode), "tips")

    # Turn dist into a nodes (rows) x tips (columns) matrix
    dist <- dist[node.numbers, 1:nTips]
    rownames(dist) <- 1:tree$Nnode

    # Calculate the average distance of a node to its tips
    fun <- function(n, d, c)mean(d[n,c[[n]]])
    res <- sapply(1:tree$Nnode, fun, dist, chld.tips)
    if(!is.null(tree$node.label)) names(res) <- tree$node.label
    return(res)
}


