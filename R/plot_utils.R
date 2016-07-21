#' annotate_balance
#'
#' annotate a balance oriented with respect to the PhILR transform.
#' That is, you can specify labels for the numerator (\code{up}) and
#' denominator (\code{down}).
#'
#' @param tr phylo object
#' @param coord named internal node/balance to annotate
#' @param p ggtree plot (tree layer), if \code{NULL} then a new plot will be
#' created.
#' @param labels label for the numerator and denominator of the balance respectively
#' @param offset offset for bar (if \code{bar=TRUE}) from tips
#' @param offset.text offset of text from bar (if \code{bar=TRUE}) or from tips
#' (if \code{bar=FALSE})
#' @param bar logical, should bar for each clade be plotted
#' @param barsize width of bar (if \code{bar=TRUE})
#' @param barfill fill of bar
#' @param geom geom used to draw label (e.g., \code{'text'} or \code{'label'})
#' @param ... additional parameters passed to \code{geom_rect} and specified \code{geom}
#'
#' @return ggplot object
#' @importFrom ggplot2 annotate
#' @importFrom ggtree ggtree get_clade_position
#'
#' @export
#' @author Justin Silverman
#'
#' @examples
#' tr <- named_rtree(10)
#'
#' annotate_balance(tr, 'n4', size=7)
#' annotate_balance(tr, 'n4', size=7, barsize=0.04, barfill='darkgreen', offset.text=0.05, color='red')
#' annotate_balance(tr, 'n4', bar=FALSE, size=7)
#' annotate_balance(tr, 'n4', bar=TRUE, size=7, labels=c('Num', 'Denom'), offset.text=.3)
#' annotate_balance(tr, 'n4', bar=TRUE, geom='label', size=8, offset.text=0.1)
annotate_balance <- function(tr, coord, p=NULL, labels=c('+','-'), offset=0,
                             offset.text=0.03, bar=TRUE, barsize=0.01,
                             barfill='darkgrey', geom='text', ...){

  # Check two labels were given
  if (length(labels) !=2) stop("two labels must be specified")
  names(labels) <- c('up', 'down')

  # get node numbers of children
  ch.coord <- unlist(get.ud.nodes(tr, coord))   # get ud.nodes (to orient) - tests if coord is a tip
  ch.nn <- name.to.nn(tr, ch.coord)   # Convert to node numbers
  names(ch.nn) <- names(ch.coord)

  if (is.null(p)){
    p <- ggtree::ggtree(tr)
  }

  # Create Dataframe that contains location and dimentions of bar
  xmax <- get_clade_position(p, node=name.to.nn(tr, coord))[,'xmax']
  df.up <- get_clade_position(p, node=ch.nn['up'])
  df.down <- get_clade_position(p, node=ch.nn['down'])
  df <- rbind(df.up, df.down)
  rownames(df) <- c('up','down')
  df[,'xmin'] <- xmax + offset
  df[, 'xmax'] <- df[,'xmin'] + barsize
  df[,'ymin'] <- df[,'ymin'] + 0.2
  df[,'ymax'] <- df[,'ymax'] - 0.2

  if (bar){
    p <- p +
      annotate(geom = geom,
               x=df['up',]$xmax+offset.text,
               y=(df['up',]$ymax + df['up',]$ymin)/2,
               label=labels['up'], ...) +
      annotate(geom = geom,
               x=df['down',]$xmax+offset.text,
               y=(df['down',]$ymax + df['down',]$ymin)/2,
               label=labels['down'], ...) +
      annotate('rect', xmin=df['up',]$xmin, xmax=df['up',]$xmax,
               ymin=df['up',]$ymin, ymax=df['up',]$ymax, fill=barfill) +
      annotate('rect', xmin=df['down',]$xmin, xmax=df['down',]$xmax,
               ymin=df['down',]$ymin, ymax=df['down',]$ymax, fill=barfill)
  } else {
    p <- p +
      annotate(geom = geom,
               x=xmax+offset.text,
               y=(df['up',]$ymax + df['up',]$ymin)/2,
               label=labels['up'], ...) +
      annotate(geom = geom,
               x=xmax+offset.text,
               y=(df['down',]$ymax + df['down',]$ymin)/2,
               label=labels['down'], ...)
  }
  p
}
