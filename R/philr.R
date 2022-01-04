#' Data transformation and driver of PhILR
#'
#' This is the main function for building the phylogenetic ILR basis,
#' calculating the weightings (of the parts and the ILR coordinates) and then
#' transforming the data.
#' 
#' @aliases build.phylo.ilr
#' @param x \strong{matrix}, \strong{TreeSummarizedExperiment},
#'   or \strong{phyloseq} object of data to be transformed
#'   (samples are rows, compositional parts are columns) -
#'   zero must be dealt with either with pseudocount, multiplicative
#'   replacement, or another method. 
#' @param sbp (Optional) give a precomputed sbp matrix
#'   \code{\link{phylo2sbp}} if you are going to build
#'   multiple ILR bases (e.g., with different weightings).
#' @param part.weights weightings for parts,
#'   can be a named vector with names corresponding to
#'   \code{colnames(x)} otherwise can be a string, options include:
#'   \describe{
#'   \item{\code{'uniform'}}{(default) uses the uniform reference measure}
#'   \item{\code{'gm.counts'}}{geometric mean of parts of x}
#'   \item{\code{'anorm'}}{aitchison norm of parts of x (after closure)}
#'   \item{\code{'anorm.x.gm.counts'}}{\code{'anorm'} times \code{'gm.counts'}}
#'   \item{\code{'enorm'}}{euclidean norm of parts of x (after closure)}
#'   \item{\code{'enorm.x.gm.counts'}}{\code{'enorm'} times \code{'gm.counts'},
#'   often gives good results} }
#' @param ilr.weights weightings for the ILR coordiantes can be a named vector
#'   with names corresponding to names of internal nodes of \code{tree}
#'   otherwise can be a string, options include:
#'   \describe{ \item{\code{'uniform'}}{(default) no weighting of the ILR basis}
#'   \item{\code{'blw'}}{sum of children's branch lengths}
#'   \item{\code{'blw.sqrt'}}{square root of \code{'blw'} option}
#'   \item{\code{'mean.descendants'}}{sum of children's branch lengths PLUS the
#'   sum of each child's mean distance to its descendent tips} }
#' @param return.all return all computed parts (e.g., computed sign matrix(\code{sbp}),
#'   part weightings (code{p}), ilr weightings (code{ilr.weights}),
#'   contrast matrix (\code{V})) as a list (default=\code{FALSE})
#'   in addition to in addition to returning the transformed data
#'   (\code{.ilrp}). If \code{return.all==FALSE} then only returns
#'   the transformed data (not in list format) If \code{FALSE} then
#'   just returns list containing \code{x.ilrp}.
#' @param abund_values A single character value for selecting the
#'   \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{assay}}
#'   to be used. Only used when \code{x} is object from this class. Default: "counts".
#' @param pseudocount optional pseudocount added to observation matrix (`x`) to avoid
#'   numerical issues from zero values. Default value is 0 which has no effect
#'   (allowing the user to handle zeros in their own preffered way before calling philr).
#'   Values < 0 given an error. 
#' @inheritParams calculate.blw

#' @details This is a utility function that pulls together a number of
#' other functions in \code{philr}. The steps that are executed are as follows:
#' \enumerate{ \item Create sbp (sign matrix) if not given \item Create parts
#' weightings if not given \item Shift the dataset with respect to the new
#' reference measure (e.g., part weightings) \item Create the basis contrast matrix
#' from the sign matrix and the reference measure \item Transform the data based on the
#' contrast matrix and the reference measure \item Calculate the specified ILR
#' weightings and multiply each balance by the corresponding weighting } Note
#' for both the reference measure (part weightings) and the ILR weightings,
#' specifying \code{'uniform'} will give the same
#' results as not weighting at all. \cr \cr Note that some of the prespecified
#' part.weights assume \code{x} is given as counts and not as relative abundances.
#' Except in this case counts or relative abundances can be given.
#'
#' The tree argument is ignored if the \code{x} argument is
#' \code{\link[phyloseq:phyloseq-class]{assay}} or
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{assay}}
#' These objects can include a phylogenetic tree. If the phylogenetic tree
#' is missing from these objects, it should be integrated directly in these
#' data objects before running \code{philr}. Alternatively, you can always
#' provide the abundance matrix and tree separately in their standard formats.
#'
#' If you have a 
#' \code{\link[SummarizedExperiment:SummarizedExperiment-class]{assay}},
#' this can be converted into
#' \code{\link[TreeSummarizedExperiment:TreeSummarizedExperiment-class]{assay}},
#' to incorporate tree information.
#'
#' @return matrix if \code{return.all=FALSE}, if \code{return.all=TRUE} then
#' a list is returned (see above).
#' @author Justin Silverman; S3 methods by Leo Lahti
#' @export 
#' @seealso \code{\link{phylo2sbp}} \code{\link{calculate.blw}}
#' @examples
#' # Prepare example data
#' tr <- named_rtree(5)
#' x <- t(rmultinom(10,100,c(.1,.6,.2,.3,.2))) + 0.65 # add a small pseudocount
#' colnames(x) <- tr$tip.label
#' philr(x, tr, part.weights='enorm.x.gm.counts',
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
#' ## Run philr for TreeSummarizedExperiment object
#' ## using the pseudocount data
#' res.tse <- philr(tse, part.weights='enorm.x.gm.counts',
#'                  ilr.weights='blw.sqrt', return.all=FALSE,
#'                  abund_values="counts", pseudocount=1)
#'
#' # Running philr on a phyloseq object
#' \donttest{
#'     pseq <- makePhyloseqFromTreeSummarizedExperiment(tse, abund_values="counts")
#'     res.pseq <- philr(pseq, part.weights='enorm.x.gm.counts',
#'                       ilr.weights='blw.sqrt', return.all=FALSE, pseudocount=1)
#' }
#'
philr <- function(x, tree=NULL, sbp=NULL, part.weights='uniform', ilr.weights='uniform', return.all=FALSE, abund_values="counts", pseudocount=0) { UseMethod("philr") }


#' @export 
philr.default <- function(x, ...){
    warning(paste("philr does not know how to handle object of class ",
    class(x),
    "and can only be used on classes TBA"))
}

#' @export 
philr.numeric <- function(x, tree, ...){
    philr.data.frame(x, tree, ...)
}

#' @export 
philr.phyloseq <- function(x, ...){

    otu.table <- t(phyloseq::otu_table(x))
    tree <- phyloseq::phy_tree(x)

    philr.data.frame(otu.table, tree=tree, ...)

}


#' @export
philr.TreeSummarizedExperiment <- function(x, tree=NULL, sbp=NULL,
                            part.weights='uniform', ilr.weights='uniform',
                            return.all=FALSE, abund_values="counts", ...){

    tree <- TreeSummarizedExperiment::rowTree(x)
    otu.table <- t(SummarizedExperiment::assays(x)[[abund_values]])
    philr.data.frame(otu.table, tree=tree, ...)
}

#' @export 
philr.matrix <- function(x, tree, ...){
    philr.data.frame(x, tree, ...)
}

#' @export
philr.array <- function(x, tree, ...){
    philr.data.frame(x, tree, ...)
}

#' @export 
philr.data.frame <- function(x, tree, sbp=NULL,
                            part.weights='uniform', ilr.weights='uniform',
                            return.all=FALSE, abund_values=NULL, pseudocount=0, ...){  
    if (pseudocount < 0) stop("pseudocount must be >= 0")
    if (!is.numeric(pseudocount)) stop("pseudocount must be numeric")
    if (length(pseudocount) != 1) stop("pseudocount must be a single number")

    # Convert x to mat with warning
    x.name <- deparse(substitute(x))
    if (is.data.frame(x)) {
        st <- paste("Converting ", x.name, " to matrix from data.frame",
            " consider adding as.matrix(", x.name, ") to function call",
            " to control conversion manually", sep="")
        warning(st)
        x <- as.matrix(x)
    }

    # Convert vector input for x to matrix
    x <- vec_to_mat(x)

    # add pseudocount
    x <- x + pseudocount
    
    # Check for Zero values in x
    if (any(x == 0)){
        stop('Zero values must be removed either through use of pseudocount,
            multiplicative replacement or other method. For convenience,
            you can use the "pseudocount" argument to philr.')
    }

    # Create the sequential binary partition sign matrix
    if (is.null(sbp)){
        message('Building Sequential Binary Partition from Tree...')
        sbp <-  phylo2sbp(tree)
    } else {
        if ( (nrow(sbp)!=ncol(x)) | (ncol(sbp)!=ncol(x)-1) ){
            stop("given sbp does not match dimentions required dimentions 
                (e.g., ncol(x) x ncol(x)-1")
        }
    }
    if (is.null(colnames(x))) {
        stop("Input data must have column names that align with phylogenetic 
            tree to prevent logical errors.")
    }
    sbp <- sbp[colnames(x), ] #Very Important line to avoid logical errors

    # Now need to create the weights on the parts
    if (is.null(part.weights)){part.weights <- 'uniform'}
    if (part.weights=='uniform'){
        p <- rep(1, ncol(x))
        names(p) <- rownames(sbp)
    } else if (part.weights=='gm.counts'){
        p <- g.colMeans(x)
        p <- p[rownames(sbp)]
    } else if (part.weights=='anorm'){
        p <- apply(miniclo(t(x)), 1, function(x) normp(x, rep(1, nrow(x))))
    } else if (part.weights=='anorm.x.gm.counts'){
        gm.counts <- g.colMeans(x)
        gm.counts <- gm.counts[rownames(sbp)]
        anorm <- apply(miniclo(t(x)), 1, function(x) normp(x, rep(1, nrow(x))))
        anorm <- anorm[rownames(sbp)]
        p <- gm.counts*anorm
    } else if (part.weights=='enorm') {
        enorm <- apply(miniclo(t(x)), 1, function(x) sqrt(sum(x^2)))
        enorm <- enorm[rownames(sbp)]
        p <- enorm
    } else if (part.weights=='enorm.x.gm.counts'){
        gm.counts <- g.colMeans(x)
        gm.counts <- gm.counts[rownames(sbp)]
        enorm <- apply(miniclo(t(x)), 1, function(x) sqrt(sum(x^2)))
        enorm <- enorm[rownames(sbp)]
        p <- gm.counts*enorm
    }

    # Make sure everything is lined up (else throw an error)
    if (!all.equal(rownames(sbp), colnames(x), names(p))) {
        stop("Rownames of SBP, Colnames of x, and names of part weights not equal!")
    }

    # shift the dataset with respect to the new reference measure
    x <- shiftp(miniclo(x), p)

    # Now create basis contrast matrix
    message('Building Contrast Matrix...')
    V <- buildilrBasep(sbp, p)

    # Finally transform the x
    message('Transforming the Data...')
    x.ilrp <- ilrp(x, p, V)

    # Now calculate ILR Weightings
    if (is.character(ilr.weights)){
        message('Calculating ILR Weights...')
        if (ilr.weights=='blw'){
            ilr.weights <- calculate.blw(tree, method='sum.children')
        } else if (ilr.weights=='blw.sqrt'){
            ilr.weights <- sqrt(calculate.blw(tree, method='sum.children'))
        } else if (ilr.weights=='uniform'){
            ilr.weights <- rep(1, ncol(x.ilrp))
            names(ilr.weights) <- colnames(x.ilrp)
        } else if (ilr.weights=='mean.descendants'){
            ilr.weights <- calculate.blw(tree, method='mean.descendants')
        }
    } else { # Validate input of ilr.weights
        if (!is.numeric(ilr.weights) | length(ilr.weights) != ncol(x.ilrp)){
            stop("ilr.weights must be recognized string or numeric vector of 
            length = ncol(x)-1")
        }
    }

    # TODO: Speed up by not computing if 'uniform'
    if (!is.null(colnames(x.ilrp))){
        ilr.weights <- ilr.weights[colnames(x.ilrp)]
    }
    # ilr.weights <- ilr.weights[colnames(x.ilrp)]
    tmp.names <- colnames(x.ilrp)
    x.ilrp <- x.ilrp %*% diag(ilr.weights)
    colnames(x.ilrp) <- tmp.names

    # Build return list
    if (return.all==FALSE)return(x.ilrp)
    if (return.all==TRUE){
        l.return = list()
        l.return[['x.ilrp']] <- x.ilrp
        l.return[['sbp']] <- sbp
        l.return[['p']] <- p      # the part weights
        l.return[['V']] <- V      # the contrast matrix
        l.return[['ilr.weights']] <- ilr.weights
    }
    return(l.return)
}

#' Inverse of PhILR Transform
#'
#' @param x.ilrp transformed data to which the inverse transform will be applied
#' @param tree (optional) to be used to build sbp and contrast matrix
#' (see details)
#' @param sbp (optional) the sbp (sequential binary partition) used to build a
#' contrast matrix (see details)
#' @param V (optional) the contrast matrix (see details)
#' @param part.weights weightings for parts, can be a named vector with names
#' corresponding to \code{colnames(x)}. Defaults to 'uniform'
#' (part.weights = 1,...,1)
#' @param ilr.weights weightings for the ILR coordiantes can be a named vector
#' with names corresponding to names of internal nodes of \code{tree}.
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
#'  x <- t(rmultinom(10,100,c(.1,.6,.2,.3,.2))) + 0.65 # add small pseudocount
#'  colnames(x) <- tr$tip.label
#'  d <- philr(x, tr, part.weights='enorm.x.gm.counts',
#'                 ilr.weights='blw.sqrt', return.all=TRUE)
#'  d.inverted <- philrInv(d$x.ilrp, V=d$V, part.weights = d$p,
#'                         ilr.weights = d$ilr.weights)
#'  all.equal(miniclo(x), d.inverted)
philrInv <- function(x.ilrp, tree=NULL, sbp=NULL, V=NULL, part.weights=NULL, ilr.weights=NULL){
    if (all(c(is.null(tree), is.null(sbp), is.null(V)))){
        stop('At least one of the parameters (tree, sbp, V) must be non-null')
    }

    # Convert vector input for x to matrix
    x.ilrp <- vec_to_mat(x.ilrp)

    # Inverse ILR -> y (the shifted composition - but closed)
    if (is.null(part.weights)){
        part.weights <- rep(1, ncol(x.ilrp)+1)
    }

    if (!is.null(V)) NULL
    else if (!is.null(sbp)) V <- buildilrBasep(sbp, part.weights)
    else V <- buildilrBasep(phylo2sbp(tree), part.weights)
    V <- V[,colnames(x.ilrp)] # make sure everything is lined up

    # Undo ILR weights - note: doing nothing (e.g., not matching criteria)
    # is equivalent to undoing uniform ilr.weights (1,...,1)
    if (!is.null(ilr.weights)) {
        ilr.weights <- ilr.weights[colnames(x.ilrp)] # ensure things are lined up
        x.ilrp <- x.ilrp / outer(rep(1,nrow(x.ilrp)), ilr.weights)
    }

    # Make sure everything is lined up (else throw an error)
    if (!all.equal(colnames(V), colnames(x.ilrp))) {
        stop("Colnames of V, Colnames of x!")
    }

    if (!all(part.weights==1)){
        if (!all.equal(rownames(V), names(part.weights))){
            stop("rownames of V and names of parts.weights not equal")
        }
    }

    y <- ilrpInv(x.ilrp, V)
    x <- miniclo(shiftpInv(y, part.weights))
    x
}

