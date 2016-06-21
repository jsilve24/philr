#' Plot function to visualize a balance conditioned on a discrete variable
#'
#' 'Breaksdown' the distance between groups along a given balance.
#'
#' @param df PhILR transformed data (see \code{\link{philr}}) in long format
#' (see \code{\link{convert_to_long}})
#' @param coord.name the name of a balance/internal node on the tree (given as a string)
#' @param name.balance logical (default: \code{FALSE}) of whether title should include output from
#' call to \code{\link{name.balance}}
#' @param tr (Optional) though needed if \code{name.balance=TRUE}
#' @inheritParams name.balance
#' @return plot created with ggplot2
#' @details
#' It is helpful to convert \code{df} to long format with the \code{\link{convert_to_long}} function
#' as this function requires a specific format for the data. Specifically it requires a data frame with four
#' columns \code{c('sample','labels', 'coord', 'value')}, where labels is the categorical variable to
#' group the data based on and value is the value of a given sample along \code{coord}.
#' @export
#' @import ggplot2
#' @examples
#' data(CSS)
#' df <- CSS$otu.table + 0.65   # add a small pseudocount
#' tree <- CSS$phy.tree
#' df.philr <- philr(df, tree, part.weights='anorm.x.gm.counts',
#'                   ilr.weights='blw.sqrt', return.all=FALSE, n_cores=1)
#' df.philr.long <- convert_to_long(df.philr, CSS$sample.data$BODY_SITE)
#' plot_density_breakdown(df.philr.long, 'n3', tree, CSS$tax.table)
plot_density_breakdown <- function(df, coord.name, tr=NULL, tax,
                                   name.balance=TRUE){
  df.filter <- subset(df, coord==coord.name)

  p <- ggplot(df.filter, aes(x=value, fill=labels)) +
    geom_density() +
    facet_grid(labels~., scales='free_y') +
    xlab('') + ylab('') +
    theme(strip.text.y=element_text(angle=0), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          panel.grid.minor=element_blank()) +
    guides(fill=FALSE)



  if (name.balance == TRUE){
    p <- p + ggtitle(paste(coord.name, name.balance(tr, tax, coord.name, method='voting'),sep=': '))
  } else{
    p <- p + ggtitle(coord.name)
  }

  return(p)
}

#' Plot a Balance
#'
#' Visualize a balance in place on a phylogenetic tree (numerator of balance is in green, denominator in blue)
#'
#' @inheritParams plot_density_breakdown
#' @inheritParams name.balance
#' @param color.tax taxonomic level (e.g., column name in \code{tax}) to color the tree by
#' @param plot.tax (Optional) Plot bar/heatmap of taxonomy labels at specified level
#' (e.g., column name in \code{tax})
#' @param color.tax (Optional) passed to gheatmap \code{color}
#' @param ... pass other arguments to ggtree (e.g., \code{layout='fan'})
#' @details
#' \code{tax} only needs to be specified if \code{color.tax != NULL}
#' @return plot created with ggtree
#' @export
#' @examples
#' data(CSS)
#' df <- CSS$otu.table + 0.65   # add a small pseudocount
#' tree <- CSS$phy.tree
#' df.philr <- philr(df, tree, part.weights='anorm.x.gm.counts',
#'                   ilr.weights='blw.sqrt', return.all=FALSE, n_cores=1)
#' df.philr.long <- convert_to_long(df.philr, CSS$sample.data$BODY_SITE)
#' plot_balance('n7', tree)
#' plot_balance('n7', tree, layout='fan')
#' plot_balance('n7', tree, tax=CSS$tax.table, color.tax='Phylum')
#'
#' # plot with phyla bar/strip
#' plot_balance('n7', tree, CSS$tax.table, plot.tax=c('Phylum'))
plot_balance <- function(coord.name, tr, tax=NULL,
                         plot.tax=NULL, color.tax=NULL, ...){
  ggtree.installed <- require('ggtree')
  if (ggtree.installed){
    # if (!is.null(color.tax)){
    #   tax.info <- split(rownames(tax), c(tax[,color.tax]))
    #   tr.tax <- ggtree::groupOTU(tr, tax.info)
    #   p.ggtree <- ggtree::ggtree(tr.tax, aes(color=group))
    # } else {
      p.ggtree <- ggtree::ggtree(tr, ...)
    # }
    l.children <- get.ud.nodes(tr, coord.name, return.nn=TRUE)
    p <- p.ggtree +
      ggtree::geom_hilight(l.children[['up']], fill='darkgreen', alpha=0.6) +
      ggtree::geom_hilight(l.children[['down']], fill='steelblue', alpha=0.6)

    if (!is.null(plot.tax)){
      t <- as.data.frame(tax[,plot.tax], stringsAsFactors=FALSE)
      t[is.na(t)] <- '' # reformat Missing values for gheatmap
      p <- ggtree::gheatmap(p, t, color=color.tax, colnames=FALSE, width=0.05)
    }

    return(p)
  }
}


#' Plot Balance on Tree next to Density Breakdown
#'
#' Visualizes Density breakdown \code{\link{plot_density_breakdown}} and balance in place on a phylogenetic tree
#' \code{\link{plot_balance}} (numerator of balance is in green, denominator in blue) next to each other.
#'
#' @inheritParams plot_balance
#' @inheritParams plot_density_breakdown
#' @return nothing returned, plot called with \code{multiplot}
#' @export
#' @examples
#' data(CSS)
#' df <- CSS$otu.table + 0.65   # add a small pseudocount
#' tree <- CSS$phy.tree
#' df.philr <- philr(df, tree, part.weights='anorm.x.gm.counts',
#'                   ilr.weights='blw.sqrt', return.all=FALSE, n_cores=1)
#' df.philr.long <- convert_to_long(df.philr, CSS$sample.data$BODY_SITE)
#' plot_density_breakdown_wtree(df.philr.long, 'n7', CSS$tax.table, tree)
#' plot_density_breakdown_wtree(df.philr.long, 'n7', CSS$tax.table, tree, plot.tax='Phylum')
plot_density_breakdown_wtree <- function(df, coord.name, tax, tr, name.balance=TRUE,
                                         plot.tax=NULL, color.tax=NULL, ...){
  ggtree.installed <- require('ggtree')
  if (ggtree.installed){
    # First make the density plot
    p.density <- plot_density_breakdown(df, coord.name, tr, tax, name.balance)

    # Then create the tree plot
    p.tree <- plot_balance(coord.name, tr, tax,
                           plot.tax=plot.tax, color.tax=color.tax, ...)

    # Then combine the plots and display
    ggtree::multiplot(p.tree, p.density, ncol=2, widths = c(.3,1))
  }
}
