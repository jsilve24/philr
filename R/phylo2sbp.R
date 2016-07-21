#' Create Sequential Binary Partition from Phylogenetic Tree
#'
#' This function converts a binary phylogenetic tree to sequential binary
#' parition to be used to then build an ILR basis for compositional
#' metagenomic data.
#'
#' @param tr a \code{phylo} tree object with n leaves
#' @param n_cores (Optional) integer specifying the number of cores to use. See Details.
#' @return a n by n-1 matrix of the sequential binary partition sign matrix
#' @details
#' Parallelization is done through \code{parallel} package using type "FORK".
#' Note parallelization is rarely needed, even for trees of upwards of 40,000 leaves.
#' @author Justin Silverman
#' @export
#' @seealso \code{\link{philr}}
#' @examples
#' tr <- named_rtree(5)
#' phylo2sbp(tr)
phylo2sbp <- function (tr, n_cores=1){
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

    if (n_cores > 1){
        cl <- parallel::makeCluster(n_cores, type='FORK')
        sbp <- parallel::parSapply(cl, ch, bpVec, ti)
        parallel::stopCluster(cl)
    } else {
        sbp <- sapply(ch, bpVec, ti)
    }

    colnames(sbp) <- tr$node.label
    rownames(sbp) <- tr$tip.label

    return(sbp)
}
