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

test_that("philr and philrInv handles vectors and default uniform arguments", {
  tr <- named_rtree(5)
  df <- c(1,4,1,22,2)
  names(df) <- tr$tip.label

  df.ilr <- philr(df, tr, return.all=F)

  expect_equivalent(philrInv(df.ilr, tr),
               miniclo(df))
})


test_that("philr and philrInv conserve distances of circle", {
  phy1 <- named_rtree(3)
  phy2 <- named_rtree(3)

  t <- seq(0, 2*pi, by = 0.1)
  x <- cos(t)
  y <- sin(t)
  circ1 <- cbind(x, y)
  colnames(circ1) <- c("n1", "n2")

  circ2 <- philr(philrInv(circ1, tree = phy1), tree = phy2)

  expect_equal(max(dist(circ1)), max(dist(circ2)))
})
