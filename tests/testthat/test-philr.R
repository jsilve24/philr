context("philr and philrInv functions")

test_that("philrInv inverses the philr transform", {
  tr <- named_rtree(5)
  df <- t(rmultinom(10,100,c(.1,.6,.2,.3,.2))) + 0.65   # add a small pseudocount
  colnames(df) <- tr$tip.label

  d <- philr(df, tr, part.weights='enorm.x.gm.counts',
             ilr.weights='blw.sqrt', return.all=TRUE)

  expect_equal(philrInv(d$df.ilrp, V=d$V, part.weights = d$p, ilr.weights = d$ilr.weights),
               miniclo(df))
})
