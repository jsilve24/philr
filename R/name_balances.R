#' Name a balance (coordinate) based on taxonomy
#'
#' For a given ILR balance (coordinate) assigns a name to the balance based on a
#' provided taxonomy table. This is useful for interpretation of the balances.
#'
#' @param tr an object of class \code{'phylo'}
#' @param tax a matrix/data.frame of taxonomy, rownames should correspond to \code{tr$tip.labels}
#' columns should be taxonomic levels (named) with increasing taxonomic resolution from left to right
#' (e.g., Phylum to the left of Genus).
#' @param coord the name of a balance/internal node on the tree (given as a string)
#' @param method currently only \code{'voting'} implemented. See Details.
#' @param thresh threshold for assignment of taxonomy to a given part of a balance
#' (must be greater than 0.5 if \code{method='voting'}; see details).
#' @param return.votes whether voting results by taxonomic level should be shown for \code{coord}. Note: this is helpful when
#' \code{name.balance} does not return a clear winner, as may be the case when a given \code{coord} represents more than one
#' taxonomic lineage. votes are returned as a list indexed by \code{colnames(tax)} Options include:
#' \describe{
#' \item{\code{NULL}}{(default) only returns the combined consensus name of the balance}
#' \item{\code{'up'}}{adds tallied votes for the 'up' node to the output list}
#' \item{\code{'down'}}{adds tallied votes for the 'down'node to the output list}
#' \item{\code{'self'}}{adds tallied votes for \code{coord} to the output list}
#' }
#' @return If \code{return.votes=NULL} returns a string of the form (ex. 'Genus_Bacteroides/Phylum_Firmicutes'). Otherwise returns
#' a list with the above string as 'name', see Arguments for \code{show.votes} for other optional returned items.
#' @details A bit of terminology:
#' \describe{
#' \item{coord}{this is the same as the names of the balances which should be the same as the names of the internal nodes of \code{tr}}
#' \item{'up'}{this is the child node of \code{coord} that is represented in the numerator of the \code{coord} balance.}
#' \item{'down'}{this is the child node of \code{coord} that is represented in the denominator of the \code{coord} balance}
#' }
#' The method \code{'voting'} assigns the name of the each part of a balance (e.g., numerator and denominator / each
#' child of \code{coord}) as follows:
#' \enumerate{
#' \item First Subset \code{tax} to contain only descendent tips of the given child of \code{coord}
#' \item Second At the finest taxonomic (farthest right of \code{tax}) see if any one taxonomic label
#' is present at or above \code{thresh}. If yes output that taxonomic label (at that taxonomic level) as the label
#' for that child of \code{coord}. If no then move to coarser taxonomic level (leftward) and repeat.
#' }
#' @author Justin Silverman
#' @export
#' @seealso \code{\link{philr}}
#' @examples
#' library(phyloseq)
#' data(CSS)
#' tr <- phy_tree(CSS)
#' tax  <- tax_table(CSS)
#' name.balance(tr, tax, 'n1')
#' name.balance(tr, tax, 'n34', return.votes=c('up','down'))
name.balance <- function(tr, tax, coord, method="voting", thresh=0.95, return.votes=NULL){
  if (method=="voting"){
    # Get tips in 'up' and 'down' subtree
    l.tips <- get.ud.tips(tr,coord)

    tax <- as(tax, "matrix")
    # Subset tax table based on above
    tax.up <- tax[l.tips[['up']],]
    tax.down <- tax[l.tips[['down']],]

    # Get Voted Consensus for up and down taxa (character strings)
    up.voted <- vote.annotation(tax.up, voting.threshold=thresh)
    down.voted <- vote.annotation(tax.down, voting.threshold=thresh)

    # Combine into a string and output
    name <- paste(up.voted,"/",down.voted,sep="")

    if (is.null(return.votes)){
      return(name)
    } else {
      res <- list('name'=name)
    }
    if ('up' %in% return.votes){
      res[['up.votes']] <- tally.votes(tax, l.tips[['up']])
    }
    if ('down' %in% return.votes){
      res[['down.votes']] <- tally.votes(tax, l.tips[['down']])
    }
    if ('self' %in% return.votes){
      res[['self.votes']] <- tally.votes(tax, unlist(l.tips))
    }
    return(res)
  }
  # In the future can extend to other methods of annotation/naming (other than just voting)
}


# Returns a list of the 'up' and 'down' subtree's root nodes
# e.g., the child node of a given coordinate
# nn is node number
get.ud.nodes <- function(tr,coord, return.nn=FALSE){
  nn <- name.to.nn(tr, coord) # get node number
  l.nodes <- list()
  child <- phangorn::Children(tr, nn)
  if (return.nn==TRUE){
    l.nodes[['up']] <- child[1]
    l.nodes[['down']] <- child[2]
  } else{
    l.nodes[['up']] <- nn.to.name(tr, child[1])
    l.nodes[['down']] <- nn.to.name(tr, child[2])
  }
  return(l.nodes)
}

# Returns a list of the 'up' and 'down' subtree's values is a vector of tip ids (corresponds
# to up and down used for sbp creation)
# Each value is the ID of a tip
get.ud.tips <- function(tr,coord){
  l.tips <- list()
  child <- phangorn::Children(tr, name.to.nn(tr,coord))
  if (length(child) > 2) stop("Tree is not soley binary.") #TODO: Bit of validation - consider better location
  l.tips[['up']] <- sapply(unlist(phangorn::Descendants(tr,child[1],type='tips')), function(x) nn.to.name(tr, x))
  l.tips[['down']] <- sapply(unlist(phangorn::Descendants(tr,child[2],type='tips')), function(x) nn.to.name(tr, x))
  return(l.tips)
}

# Find most concerved
# returns character
# NA are not counted but held against the winner in voting
# Candidate must have >= voting.threshold to be considered the winner
vote.annotation <- function(tax, voting.threshold=0.95){
  if (voting.threshold <= 0.5)stop('voting.threshold must be > 0.5 for unique winner')
  if (is(tax, "character")){ # e.g., is there only 1 voter here
    tmp.names <- names(tax)
    tax <- matrix(tax,nrow=1)
    colnames(tax) <- tmp.names
  }
  nr <- nrow(tax) # the number of tips
  nc <- ncol(tax) # the number of taxonomic ranks in table
  name <- NULL

  for (i in seq(nc,1)){ # evaluate in decreasing order of taxonomic resolution
    votes <- tax[,i]
    if (all(is.na(votes))) {next}
    votes <- votes[!is.na(votes)] # drop NA votes but hold against the total number
    winner <- sort(table(votes), decreasing=TRUE)[1]
    if (!is.na(winner) & (winner/nr >= voting.threshold)){ # Arbitrarily set threshold to 95% of votes
      name <- names(winner)
      # Try and append taxonomic rank to name if columns of tax table are labled
      if (!is.null(colnames(tax))) {name <- paste(colnames(tax)[i],'_',name,sep="")}
      break
    } # Else try the next level up
  }

  if (is.null(name)){ # If no consensus vote found (above threshold)
    name <- "Unclear_Lineage_Identity"
  }
  return(name)
}

# Given a list of tax IDs (e.g., rownames in tax table)
# print out a easy to read list showing whats present at each taxonomic level.
# This is really to be used when You don't get a result you like with
# vote.annotation()
tally.votes <- function(tax,ids){
  l.votes <- list()
  nc <- ncol(tax) # the number of taxonomic ranks in the table
  for (i in seq(nc,1)){
      votes <- tax[ids,i]
      votes <- votes[!is.na(votes)]
      rank <- ifelse(!is.null(colnames(tax)), colnames(tax)[i],i)
      l.votes[[rank]] <- table(votes)
  }
  return(l.votes)
}
