## ------------------------------------------------------------------------
library(philr)
library(phyloseq)

## ------------------------------------------------------------------------
data(CSS)
df <- t(otu_table(CSS)) + 1
tree <- phy_tree(CSS)

## ------------------------------------------------------------------------
df.philr <- philr(df, tree, part.weights='anorm.x.gm.counts', ilr.weights='blw.sqrt')

## ------------------------------------------------------------------------
df.philr.long <- convert_to_long(df.philr, get_variable(CSS, 'BODY_SITE'))
plot_balance('n7', tree)
plot_density_breakdown(df.philr.long, 'n7', tax_table(CSS))

