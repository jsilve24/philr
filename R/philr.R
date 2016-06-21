#' Data transformation and driver of PhILR.
#'
#' This is the main function for building the phylogenetic ILR basis, calculating the
#' weightings (of the parts and the ILR coordinates) and then transforming the data.
#' @aliases  build.phylo.ilr
#' @param df matrix of data to be transformed (samples are rows,
#' compositional parts are columns) - zero must be dealt with either with pseudocount,
#' multiplicative replacement, or another method.
#' @param sbp (Optional) give a precomputed sbp matrix \code{link{phylo2sbp}}
#' if you are going to build multiple ILR bases (e.g., with different weightings).
#' @param part.weights weightings for parts, can be a named vector with names
#' corresponding to \code{colnames(df)} otherwise can be a string, options include:
#' \describe{
#' \item{\code{'uniform'}}{(default) uses the uniform reference measure}
#' \item{\code{'gm.counts'}}{geometric mean of parts of df}
#' \item{\code{'anorm'}}{aitchison norm of parts of df}
#' \item{\code{'anorm.x.gm.counts'}}{\code{'anorm'*'gm.counts'}, often gives good results}
#' }
#' @param ilr.weights weightings for the ILR coordiantes can be a named vector with names
#' corresponding to names of internal nodes of \code{tree} otherwise can be a string,
#' options include:
#' \describe{
#' \item{\code{'uniform'}}{(default) no weighting of the ILR basis}
#' \item{\code{'blw'}}{sum of children's branch lengths}
#' \item{\code{'blw.sqrt'}}{square root of \code{'blw'} option}
#' \item{\code{'mean.descendants'}}{sum of children's branch lengths PLUS the sum of
#' each child's mean distance to its descendent tips}
#' }
#' @param return.all return all computed parts (e.g., computed sign matrix(\code{sbp} ),
#' part weightings (code{p}), ilr weightings (code{ilr.weights}), contrast matrix (\code{V}))
#' as a list (default=\code{FALSE}) in addition to in addition to returning the transformed data (\code{df.ilrp}).
#' If \code{return.all==FALSE} then only returns the transformed data (not in list format)
#' If \code{FALSE} then just returns list containing \code{df.ilrp}.
#' @param n_cores (Optional) integer specifying the number of cores to use. See Details.
#' @inheritParams calculate.blw
#' @details
#' This is a utility function that pulls together a number of other functions
#' in \code{philr}. The steps that are executed are as follows:
#' \enumerate{
#' \item Create sbp (sign matrix) if not given
#' \item Create parts weightings if not given
#' \item Shift the dataset with respect to the new reference measure (e.g., part weightings)
#' \item Create the basis contrast matrix from the sign matrix and the reference measure
#' \item Transform the data based on the contrast matrix and the reference measure
#' \item Calculate the specified ILR weightings and multiply each balance by the corresponding weighting
#' }
#' Note for both the reference measure (part weightings) and the ILR weightings, specifying \code{‘uniform’} will
#' give the same results as not weighting at all.\cr \cr
#' Parallelization is done through \code{parallel} package using type "FORK".
#' Note parallelization is rarely needed, even for trees of upwards of 40,000 leaves.
#' @author Justin Silverman
#' @export
#' @import compositions
#' @seealso \code{\link{phylo2sbp}} \code{\link{calculate.blw}}
#' @examples
#' data(CSS)
#' df <- CSS$otu.table + 0.65   # add a small pseudocount
#' tree <- CSS$phy.tree
#' df.philr <- philr(df, tree, part.weights='anorm.x.gm.counts',
#'                   ilr.weights='blw.sqrt', return.all=FALSE, n_cores=1)
#' df.philr[1:5,1:5]
philr <- function(df, tree, sbp=NULL,
                            part.weights='uniform', ilr.weights='uniform',
                            return.all=FALSE, n_cores=1){
  # Check for Zero values in df
  if (any(df == 0)){
    stop('Zero values must be removed either through use of pseudocount, multiplicative replacement or other method.')
  }

  # Create the sequential binary partition sign matrix
  if (is.null(sbp)){
    print('Building Sequential Binary Partition from Tree...')
    sbp <-  phylo2sbp(tree, n_cores)
  } else {
    if ( (nrow(sbp)!=ncol(df)) | (ncol(sbp)!=ncol(df)-1) ){
      stop("given sbp does not match dimentions required dimentions (e.g., ncol(df) x ncol(df)-1")
    }
  }
  sbp <- sbp[colnames(df), ] #Very Important line to avoid logical errors

  # Now need to create the weights on the parts
  if (is.null(part.weights)){part.weights <- 'uniform'}
  if (part.weights=='uniform'){
    p <- rep(1, ncol(df))
    names(p) <- rownames(sbp)
  } else if (part.weights=='gm.counts'){
    p <- geometricmeanCol(df)
    p <- p[rownames(sbp)]
  } else if (part.weights=='anorm'){
    p <- apply(acomp(t(df)), 1, norm.acomp)
  } else if (part.weights=='anorm.x.gm.counts'){
    gm.counts <- geometricmeanCol(df)
    gm.counts <- gm.counts[rownames(sbp)]
    anorm <- apply(acomp(t(df)), 1, norm.acomp)
    anorm <- anorm[rownames(sbp)]
    p <- gm.counts*anorm
  }

  # Make sure everything is lined up (else throw an error)
  if (!all.equal(rownames(sbp), colnames(df), names(p))) {
      stop("Rownames of SBP, Colnames of df, and names of part weights not equal!")
  }

  # shift the dataset with respect to the new reference measure
  df <- as(acomp(df) - acomp(p), 'matrix')

  # Now create basis contrast matrix
  print('Building Contrast Matrix...')
  V <- buildilrBasep(sbp, p)

  # Finally transform the df
  print('Transforming the Data...')
  df.ilrp <- ilrp(df, p, V)

  # Now calculate ILR Weightings
  if (is.character(ilr.weights)){
    print('Calculating ILR Weights...')
    if (ilr.weights=='blw'){
      ilr.weights <- calculate.blw(tree, method='sum.children')
    } else if (ilr.weights=='blw.sqrt'){
      ilr.weights <- sqrt(calculate.blw(tree, method='sum.children'))
    } else if (ilr.weights=='uniform'){
      ilr.weights <- rep(1, ncol(df.ilrp))
      names(ilr.weights) <- colnames(df.ilrp)
    } else if (ilr.weights=='mean.descendants'){
        ilr.weights <- calculate.blw(tree, method='mean.descendants')
    }
  } else { # Validate input of ilr.weights
    if (!is.numeric(ilr.weights) | length(ilr.weights) != ncol(df.ilrp)){
      stop("ilr.weights must be recognized string or numeric vector of length = ncol(df)-1")
    }
  }

  #TODO: Speed up by not computing if 'uniform'
  ilr.weights <- ilr.weights[colnames(df.ilrp)]
  tmp.names <- colnames(df.ilrp)
  df.ilrp <- df.ilrp %*% diag(ilr.weights)
  colnames(df.ilrp) <- tmp.names

  # Build return list
  if (return.all==FALSE)return(df.ilrp)
  if (return.all==TRUE){
    l.return = list()
    l.return[['df.ilrp']] <- df.ilrp
    l.return[['sbp']] <- sbp
    l.return[['p']] <- p      # the part weights
    l.return[['V']] <- V      # the contrast matrix
    l.return[['ilr.weights']] <- ilr.weights
  }
  return(l.return)
}
