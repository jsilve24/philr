## ------------------------------------------------------------------------
library(philr)

## ------------------------------------------------------------------------
data(CSS)
df <- CSS$otu.table + 1
tree <- CSS$phy.tree

## ------------------------------------------------------------------------
df.philr <- philr(df, tree, part.weights='anorm.x.gm.counts', ilr.weights='blw.sqrt')

## ------------------------------------------------------------------------
# df.philr.long <- convert_to_long(df.philr, get_variable(CSS, 'BODY_SITE'))
# plot_balance('n7', tree)
# plot_density_breakdown(df.philr.long, 'n7', tax_table(CSS))

