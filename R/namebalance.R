require(phangorn)

####### ACCESSOR FUNCTIONS FROM NAMES TO NUMBERS #########
# Main accessor of node number through coordinate name
c.to.nn <- function(tr, c){
    return(which(tr$node.label==c)+Ntip(tr))
}

# Main accessor of node number through tip name
t.to.nn <- function(tr, t){
    return(which(tr$tip.label==t))
}

# Main accessor of node or tip name thorugh node number
nn.to.name <- function(tr, nn){
    n <- Ntip(tr)
    if (nn <= n)return(tr$tip.label[nn])
    return(tr$node.label[nn-n])
}


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
#' @param thresh threshold for assignment of taxonomy to a given part of a balance (see details).
#' @details
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
name.balance <- function(tr, tax, coord, method=c("voting"), thresh=0.95){
  names <- c()
  if ("voting" %in% method){
    # Get tips in 'up' and 'down' subtree
    l.tips <- get.ud.tips(tr,coord)

    # Subset tax table based on above
    tax.up <- tax[l.tips[['up']],]
    tax.down <- tax[l.tips[['down']],]

    # Get Voted Consensus for up and down taxa (character strings)
    up.voted <- vote.annotation(tax.up, voting.threshold=thresh)
    down.voted <- vote.annotation(tax.down, voting.threshold=thresh)

    # Combine into a string and output
    name.voted <- paste(up.voted,"/",down.voted,sep="")
    names <- c(names,name.voted)
  }

  # Can annotate with multiple methods and it will return a vector of names for that
  # coordinate
  return(names)
}


# Returns a list of the 'up' and 'down' subtree's root nodes
# e.g., the child node of a given coordinate
# nn is node number
get.ud.nodes <- function(tr,coord){
  nn <- c.to.nn(tr, coord) # get node number
  l.nodes <- list()
  child <- Children(tr, nn)
  l.nodes[['up']] <- nn.to.name(tr, child[1])
  l.nodes[['down']] <- nn.to.name(tr, child[2])
  return(l.nodes)
}

# Returns a list of the 'up' and 'down' subtree's values is a vector of tip ids (corresponds
# to up and down used for sbp creation)
# Each value is the ID of a tip
get.ud.tips <- function(tr,coord){
  l.tips <- list()
  child <- Children(tr, c.to.nn(tr,coord))
  if (length(child) > 2) stop("Tree is not soley binary.") #TODO: Bit of validation - consider better location
  l.tips[['up']] <- sapply(unlist(Descendants(tr,child[1],type='tips')), function(x) nn.to.name(tr, x))
  l.tips[['down']] <- sapply(unlist(Descendants(tr,child[2],type='tips')), function(x) nn.to.name(tr, x))
  return(l.tips)
}



# Find most concerved
# returns character
# NA are not counted but held against the winner in voting
# Candidate must have >= voting.threshold to be considered the winner
vote.annotation <- function(tax, voting.threshold=0.95){
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

# Write this function to take function(tr, tax, coord, child=c('up','down'))
# if child=NULL then output voting for coord rather than its child
# Do this so its more in line with the call to name.balance.

#' Tally 'votes' for taxonomic label of a coordinate/internal node of tree
tally.votes <- function(tax,ids){
    l.votes <- list()
    # Assume its a matrix
    nc <- ncol(tax) # the number of taxonomic ranks in the table
    for (i in seq(nc,1)){
        votes <- tax[ids,i]
        votes <- votes[!is.na(votes)]
        rank <- ifelse(!is.null(colnames(tax)), colnames(tax)[i],i)
        l.votes[[rank]] <- table(votes)
    }
    return(l.votes)
}

# May want to pair print.votes with something like
# Assign identity, where the user can view the votes themselves
# and choose how they would like to call it.
