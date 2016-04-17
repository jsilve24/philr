# df - ilr transformed data in long data.frame format
# coord.name - string (name of coordinate to plot)
# tax - taxonomy table
plot.sites.density <- function(df, coord.name, tax){
    df.filter <- filter(df, coord==coord.name)
    
    p <- ggplot(df.filter, aes(x=value, fill=longsite)) +
       geom_density() +
       facet_grid(longsite~.,scales='free_y') +
       xlab('') + ylab('') +
       theme(strip.text.y=element_text(angle=0), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
       panel.grid.minor=element_blank()) +
       guides(fill=FALSE) +
       ggtitle(paste(coord.name, name.balance(tree, tax, coord.name, method='voting'),sep=': '))
    return(p)
}

# coord.name - string (name of coordinate to plot)
# phylo - a phylogenetic tree of type phylo4 object
# p.ggtree - for large trees, consider precomputing the tree skeleton
#    and passing that in to save time (especially if you will be calling this plot function
#    a number of times
plot.tree.balance <- function(coord.name, phylo, p.ggtree=NULL){
  if (is.null(p.ggtree)){
    p.ggtree <- ggtree(phylo)
  }
  l.children <- get.ud.nodes(phylo, coord.name)
  p <- p.ggtree +
    geom_hilight(l.children[['up']], fill='darkgreen', alpha=0.6) +
    geom_hilight(l.children[['down']], fill='steelblue', alpha=0.6) 
  #+ geom_tippoint(aes(color=Phylum))
  return(p)
}


# df - ilr transformed data in long data.frame format
# coord - string (names of coordinate to plot)
# tax - taxonomy table (e.g., output of tax_table() phyloseq)
# phylo - a phylogenetic tree of type phylo4 object
plot.sites.density.wtree <- function(df, coord.name, tax, phylo, ...){
    # First make the density plot
    p.density <- plot.sites.density(df, coord.name, tax)
    
    # Then create the tree plot
    p.tree <- plot.tree.balance(coord.name, phylo, ...)
    
    # Then combine the plots and display
    multiplot(p.tree, p.density,ncol=2, widths = c(.3,1))
}