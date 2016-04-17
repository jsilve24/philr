library(parallel)

calc.distances.parallel <- function(phyloseq,df.gm.blw, blw, no_cores=7){
  
  # Initiate cluster
  cl <- makeCluster(no_cores, type='FORK')
  
  # Set up list input for parallelization
  inpt <- 1:7
  
  # Define function for parallelization
  prl <- function(x){
    if (x==1){ # Naive Euclidean
      return(distance(phyloseq, method='euclidean'))
    }
    if (x==2){ # Bray Curtis
      return(distance(phyloseq, method='bray'))
    }
    if (x==3){ # Unweighted UniFrac
      return(distance(phyloseq, method='unifrac'))
    }
    if (x==4){ # Weighted UniFrac
      return(distance(phyloseq, method='wunifrac'))
    }
    if (x==5){ # Jaccard
      return(distance(phyloseq, method='jaccard'))
    }
    if (x==6){ # [W]eighted, [P]hylogenetically aware 
      return(dist(df.gm.blw, method='euclidean'))
    }
    if (x==7){ # [W]eighted, [N]ot Phylogenically aware
      df.gm <- df.gm.blw %*% diag(1/blw)
      colnames(df.gm) <- colnames(df.gm.blw)
      return(dist(df.gm, method='euclidean'))
    }
  }
  
  # Now acctually launch the parallelization
  l.dist <- parLapply(cl, inpt, prl)
  names <- c('naive.euclidean','bray.curtis','unweighted.unifrac','weighted.unifrac',
                'jaccard', 'WP.ilr', 'WN.ilr')
  names(l.dist) <- names[inpt]
  stopCluster(cl)
  return(l.dist)
}

# Same as above but given a list of matricies/data.frames 
# only calculates using euclidean distance
calc.euclidean.distances.parallel  <- function(l.df, no_cores=8){
  
  # Initiate cluster
  cl <- makeCluster(no_cores, type='FORK')
  
  inpt <- names(l.df)
  
  d <- function(x){
    df <- l.df[[x]]
    return(dist(df, method='euclidean'))
  }
  
  # Now acctually launch the parallelization
  l.dist <- parLapply(cl, inpt, d)
  names(l.dist) <- inpt
  
  stopCluster(cl)
  return(l.dist)
}


# Calculate PCoA for the list of distance matricies created using
# the calc.distance function
calc.pcoa.parallel <- function(phyloseq, l.dist, no_cores=8){
  # Initiate cluster
  cl <- makeCluster(no_cores, type='FORK')

  l.pcoa <- parLapply(cl, names(l.dist), function(x) ordinate(phyloseq, 'PCoA', distance=l.dist[[x]]))
  names(l.pcoa) <- names(l.dist)
  
  stopCluster(cl)
  return(l.pcoa)
}


# l.dist - list of distance matricies from calc.distances
# classes - vector of classes for each sample
eval.clustering.parallel <- function(l.dist, classes, no_cores=8){
  # Initiate cluster
  cl <- makeCluster(no_cores, type='FORK')
  
  dist.methods <- names(l.dist)
  
  evcl <- function(dist.method){
    k = length(unique(classes))
    
    l.evcl <- list() # temp list to return data
    l.evcl[['method']] <- dist.method
    
    pamClust <- pam(l.dist[[dist.method]],k=k, cluster.only=TRUE)
    l.evcl[['table']] <- table(classes, pamClust)
    clustats <- cluster.stats(l.dist[[dist.method]], 
                              as.numeric(classes), 
                              pamClust,
                              compareonly=TRUE)
    l.evcl[['vi']] <- clustats$vi
    l.evcl[['corrected.rand']] <- clustats$corrected.rand
    return(l.evcl)
  }
  
  l.clustats <- parLapply(cl, dist.methods, evcl)
  stopCluster(cl)
  return(l.clustats)
}