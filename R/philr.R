#' Data transformation and driver of PhILR.
#'
#' This is the main function for building the phylogenetic ILR basis, calculating the
#' weightings (of the parts and the ILR coordinates) and then transforming the data.
#' @aliases  build.phylo.ilr
#' @param df \strong{matrix} of data to be transformed (samples are rows,
#' compositional parts are columns) - zero must be dealt with either with pseudocount,
#' multiplicative replacement, or another method.
#' @param sbp (Optional) give a precomputed sbp matrix \code{\link{phylo2sbp}}
#' if you are going to build multiple ILR bases (e.g., with different weightings).
#' @param part.weights weightings for parts, can be a named vector with names
#' corresponding to \code{colnames(df)} otherwise can be a string, options include:
#' \describe{
#' \item{\code{'uniform'}}{(default) uses the uniform reference measure}
#' \item{\code{'gm.counts'}}{geometric mean of parts of df}
#' \item{\code{'anorm'}}{aitchison norm of parts of df (after closure)}
#' \item{\code{'anorm.x.gm.counts'}}{\code{'anorm'} times \code{'gm.counts'}}
#' \item{\code{'enorm'}}{euclidean norm of parts of df (after closure)}
#' \item{\code{'enorm.x.gm.counts'}}{\code{'enorm'} times \code{'gm.counts'}, often gives good results}
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
#' @param return.all return all computed parts (e.g., computed sign matrix(\code{sbp}),
#' part weightings (code{p}), ilr weightings (code{ilr.weights}), contrast matrix (\code{V}))
#' as a list (default=\code{FALSE}) in addition to in addition to returning the transformed data (\code{df.ilrp}).
#' If \code{return.all==FALSE} then only returns the transformed data (not in list format)
#' If \code{FALSE} then just returns list containing \code{df.ilrp}.
#' @param abund_values A single character value for selecting the
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{assay}} to be
#'   used. Only used when \code{df} is object from this class. Default: "counts".
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
#' Note for both the reference measure (part weightings) and the ILR weightings, specifying \code{'uniform'} will
#' give the same results as not weighting at all. \cr \cr
#' Note that some of the prespecified part.weights assume \code{df} is given as
#' counts and not as relative abundances. Except in this case counts or relative
#' abundances can be given.
#'
#' The tree argument is ignored if the \code{df} argument is
#' \code{\link[phyloseq:phyloseq-class]{phyloseq}} or
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}} object.
#' 
#' These objects can include a phylogenetic tree. If the phylogenetic tree
#' is missing from these objects, it should be integrated directly in these
#' data objects before running \code{philr}. Alternatively, you can always
#' provide the abundance matrix and tree separately in their standard formats.
#'
#' If you have a
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{SummarizedExperiment}}
#' object, this can be converted into
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{TreeSummarizedExperiment}}
#' to incorporate tree information.
#'
#'
#' @return matrix if \code{return.all=FALSE}, if \code{return.all=TRUE} then a list is returned (see above).
#' @author Justin Silverman; S3 methods by Leo Lahti
#' @export 
#' @seealso \code{\link{phylo2sbp}} \code{\link{calculate.blw}}
#' @examples
#' # Prepare example data
#' tr <- named_rtree(5)
#' df <- t(rmultinom(10,100,c(.1,.6,.2,.3,.2))) + 0.65 # add a small pseudocount
#' colnames(df) <- tr$tip.label
#' philr(df, tr, part.weights='enorm.x.gm.counts',
#'                ilr.weights='blw.sqrt', return.all=FALSE)
#'
#' # Running philr on a TreeSummarizedExperiment object
#'
#' ## Prepare example data
#' library(mia)
#' library(tidyr)
#' data(GlobalPatterns, package="mia")
#'
#' ## Select prevalent taxa 
#' tse <-  GlobalPatterns %>% subsetByPrevalentTaxa(
#'                                detection = 3,
#'                                prevalence = 20/100,
#'                                as_relative = FALSE)
#'
#' ## Pick taxa that have notable abundance variation across sammples
#' variable.taxa <- apply(assay(tse, "counts"), 1, function(x) sd(x)/mean(x) > 3.0)
#' tse <- tse[variable.taxa,]
#'
#' # Collapse the tree
#' tree <- ape::keep.tip(phy = rowTree(tse), tip = rowLinks(tse)$nodeNum)
#' rowTree(tse) <- tree
#'
#' ## Add a new assay with a pseudocount 
#' assays(tse)$counts.shifted <- assay(tse, "counts") + 1
#' 
#' ## Run philr for TreeSummarizedExperiment object
#' ## using the pseudocount data
#' res.tse <- philr(tse, part.weights='enorm.x.gm.counts',
#'                ilr.weights='blw.sqrt', return.all=FALSE,
#'                abund_values="counts.shifted")
#'
#'
#' # Running philr on a phyloseq object
#' \dontrun{
#'   pseq <- makePhyloseqFromTreeSummarizedExperiment(tse)
#'   res.pseq <- philr(pseq, part.weights='enorm.x.gm.counts',
#'                ilr.weights='blw.sqrt', return.all=FALSE)
#' }
#'
philr <- function(df, tree=NULL, sbp=NULL, part.weights='uniform', ilr.weights='uniform', return.all=FALSE, abund_values="counts") { UseMethod("philr") }


#' @export 
philr.default <- function(df, ...){
  warning(paste("philr does not know how to handle object of class ",
  class(df),
  "and can only be used on classes TBA"))
}

#' @export 
philr.numeric <- function(df, tree, ...){
    philr.data.frame(df, tree, ...)
}

#' @export 
philr.phyloseq <- function(df, ...){

    otu.table <- t(phyloseq::otu_table(df))
    tree <- phyloseq::phy_tree(df)

    philr.data.frame(otu.table, tree=tree, ...)

}



#' @export
philr.TreeSummarizedExperiment <- function(df, abund_values="counts", ...){
    tree <- TreeSummarizedExperiment::rowTree(df)
    otu.table <- t(SummarizedExperiment::assays(df)[[abund_values]])
    philr.data.frame(otu.table, tree=tree, ...)
}

#' @export 
philr.matrix <- function(df, tree, ...){
    philr.data.frame(df, tree, ...)
}

#' @export
philr.array <- function(df, tree, ...){
    philr.data.frame(df, tree, ...)
}

#' @export 
philr.data.frame <- function(df, tree, sbp=NULL,
                            part.weights='uniform', ilr.weights='uniform',
                            return.all=FALSE, abund_values=NULL, ...){
  
    # Convert df to mat with warning
    df.name <- deparse(substitute(df))
    if (is.data.frame(df)) {
        st <- paste("Converting ", df.name, " to matrix from data.frame",
            " consider adding as.matrix(", df.name, ") to function call",
            " to control conversion manually", sep="")
        warning(st)
        df <- as.matrix(df)
    }

    # Convert vector input for df to matrix
    df <- vec_to_mat(df)

    # Check for Zero values in df
    if (any(df == 0)){
        stop('Zero values must be removed either through use of pseudocount,
            multiplicative replacement or other method.')
    }

    # Create the sequential binary partition sign matrix
    if (is.null(sbp)){
        message('Building Sequential Binary Partition from Tree...')
        sbp <-  phylo2sbp(tree)
    } else {
        if ( (nrow(sbp)!=ncol(df)) | (ncol(sbp)!=ncol(df)-1) ){
            stop("given sbp does not match dimentions required dimentions 
                (e.g., ncol(df) x ncol(df)-1")
        }
    }
    if (is.null(colnames(df))) {
        stop("Input data must have column names that align with phylogenetic tree to 
            prevent logical errors.")
    }
    sbp <- sbp[colnames(df), ] #Very Important line to avoid logical errors

    # Now need to create the weights on the parts
    if (is.null(part.weights)){part.weights <- 'uniform'}
    if (part.weights=='uniform'){
        p <- rep(1, ncol(df))
        names(p) <- rownames(sbp)
    } else if (part.weights=='gm.counts'){
        p <- g.colMeans(df)
        p <- p[rownames(sbp)]
    } else if (part.weights=='anorm'){
        p <- apply(miniclo(t(df)), 1, function(x) normp(x, rep(1, nrow(df))))
    } else if (part.weights=='anorm.x.gm.counts'){
        gm.counts <- g.colMeans(df)
        gm.counts <- gm.counts[rownames(sbp)]
        anorm <- apply(miniclo(t(df)), 1, function(x) normp(x, rep(1, nrow(df))))
        anorm <- anorm[rownames(sbp)]
        p <- gm.counts*anorm
    } else if (part.weights=='enorm') {
        enorm <- apply(miniclo(t(df)), 1, function(x) sqrt(sum(x^2)))
        enorm <- enorm[rownames(sbp)]
        p <- enorm
    } else if (part.weights=='enorm.x.gm.counts'){
        gm.counts <- g.colMeans(df)
        gm.counts <- gm.counts[rownames(sbp)]
        enorm <- apply(miniclo(t(df)), 1, function(x) sqrt(sum(x^2)))
        enorm <- enorm[rownames(sbp)]
        p <- gm.counts*enorm
    }

    # Make sure everything is lined up (else throw an error)
    if (!all.equal(rownames(sbp), colnames(df), names(p))) {
        stop("Rownames of SBP, Colnames of df, and names of part weights not equal!")
    }

    # shift the dataset with respect to the new reference measure
    df <- shiftp(miniclo(df), p)

    # Now create basis contrast matrix
    message('Building Contrast Matrix...')
    V <- buildilrBasep(sbp, p)

    # Finally transform the df
    message('Transforming the Data...')
    df.ilrp <- ilrp(df, p, V)

    # Now calculate ILR Weightings
    if (is.character(ilr.weights)){
        message('Calculating ILR Weights...')
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
            stop("ilr.weights must be recognized string or numeric vector of 
            length = ncol(df)-1")
        }
    }

    # TODO: Speed up by not computing if 'uniform'
    if (!is.null(colnames(df.ilrp))){
        ilr.weights <- ilr.weights[colnames(df.ilrp)]
    }
    # ilr.weights <- ilr.weights[colnames(df.ilrp)]
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

#' Inverse of PhILR Transform
#'
#' @param df.ilrp transformed data to which the inverse transform will be applied
#' @param tree (optional) to be used to build sbp and contrast matrix (see details)
#' @param sbp (optional) the sbp (sequential binary partition) used to build a
#' contrast matrix (see details)
#' @param V (optional) the contrast matrix (see details)
#' @param part.weights weightings for parts, can be a named vector with names
#' corresponding to \code{colnames(df)}. Defaults to 'uniform' (part.weights = 1,...,1)
#' @param ilr.weights weightings for the ILR coordiantes can be a named vector with names
#' corresponding to names of internal nodes of \code{tree}.
#' Defaults to 'uniform' (ilr.weights = 1,...,1)
#'
#' @details This is a utility function for calculating the inverse of the
#' \code{\link{philr}} transform. Note that at least one of the following
#' parameters must be specified (\code{tree}, \code{sbp}, or \code{V}).
#' @return a matrix of compositions (rows are samples, columns are parts),
#' function removes the effects of ilr weights, part weights, and unshifts
#' the composition.
#'
#' @author Justin Silverman
#' @export
#' @seealso \code{\link{philr}}
#'
#' @examples
#'  tr <- named_rtree(5)
#'  df <- t(rmultinom(10,100,c(.1,.6,.2,.3,.2))) + 0.65   # add a small pseudocount
#'  colnames(df) <- tr$tip.label
#'  d <- philr(df, tr, part.weights='enorm.x.gm.counts',
#'                 ilr.weights='blw.sqrt', return.all=TRUE)
#'  d.inverted <- philrInv(d$df.ilrp, V=d$V, part.weights = d$p,
#'                         ilr.weights = d$ilr.weights)
#'  all.equal(miniclo(df), d.inverted)
philrInv <- function(df.ilrp, tree=NULL, sbp=NULL, V=NULL, part.weights=NULL, ilr.weights=NULL){
    if (all(c(is.null(tree), is.null(sbp), is.null(V)))){
        stop('At least one of the parameters (tree, sbp, V) must be non-null')
    }

    # Convert vector input for df to matrix
    df.ilrp <- vec_to_mat(df.ilrp)

    # Inverse ILR -> y (the shifted composition - but closed)
    if (is.null(part.weights)){
        part.weights <- rep(1, ncol(df.ilrp)+1)
    }

    if (!is.null(V)) NULL
    else if (!is.null(sbp)) V <- buildilrBasep(sbp, part.weights)
    else V <- buildilrBasep(phylo2sbp(tree), part.weights)
    V <- V[,colnames(df.ilrp)] # make sure everything is lined up

    # Undo ILR weights - note: doing nothing (e.g., not matching criteria)
    # is equivalent to undoing uniform ilr.weights (1,...,1)
    if (!is.null(ilr.weights)) {
        ilr.weights <- ilr.weights[colnames(df.ilrp)] # ensure things are lined up
        df.ilrp <- df.ilrp / outer(rep(1,nrow(df.ilrp)), ilr.weights)
    }

    # Make sure everything is lined up (else throw an error)
    if (!all.equal(colnames(V), colnames(df.ilrp))) {
        stop("Colnames of V, Colnames of df!")
    }

    if (!all(part.weights==1)){
        if (!all.equal(rownames(V), names(part.weights))){
            stop("rownames of V and names of parts.weights not equal")
        }
    }

    y <- ilrpInv(df.ilrp, V)
    x <- miniclo(shiftpInv(y, part.weights))
    x
}

