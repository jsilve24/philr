#' Create Sequential Binary Partition from Phylogenetic Tree
#'
#' This function converts a binary phylogenetic tree to sequential binary
#' parition to be used to then build an ILR basis for compositional
#' metagenomic data.
#'
#' @param tr a \code{phylo} tree object with n leaves
#' @return a n by n-1 matrix of the sequential binary partition sign matrix
#' @details
#' The choice of orientation for a balance (i.e., which of the two descendant
#' clades of an internal node is in the numerator or denominator of the
#' log-ratio) is given by the default of the function phangorn::Children
#' and that choise is used consistently throughout the philr package.
#' @author Justin Silverman
#' @export
#' @seealso \code{\link{philr}}
#' @examples
#' tr <- named_rtree(5)
#' phylo2sbp(tr)
phylo2sbp <- function (tr){
    nTips <- ape::Ntip(tr)
    ch <- phangorn::Children(tr, (nTips+1):(nTips+tr$Nnode))
    ti <- phangorn::Descendants(tr,(nTips+1):(nTips+tr$Nnode), 'tips')

    # Only look for children if its an internal node
    idx <- function(x,t){
        if (x <= nTips) return(x)
        return(t[[x-nTips]])
    }

    bpVec <- function(c, t){
        tmp <- rep(0, nTips)
        tmp[idx(c[1],t)] <- 1
        tmp[idx(c[2],t)] <- -1
        return(tmp)
    }

    sbp <- sapply(ch, bpVec, ti)

    colnames(sbp) <- tr$node.label
    rownames(sbp) <- tr$tip.label

    return(sbp)
}
