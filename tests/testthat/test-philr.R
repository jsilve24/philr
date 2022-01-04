context("philr and philrInv functions")

test_that("philrInv inverses the philr transform", {
  tr <- named_rtree(5)
  x <- t(rmultinom(10,100,c(.1,.6,.2,.3,.2))) + 0.65   # add a small pseudocount
  colnames(x) <- tr$tip.label

  d <- philr(x, tree=tr, part.weights='enorm.x.gm.counts',
             ilr.weights='blw.sqrt', return.all=TRUE)

  expect_equal(philrInv(d$x.ilrp, V=d$V, part.weights = d$p, ilr.weights = d$ilr.weights),
               miniclo(x))
})

test_that("philr and philrInv handles vectors and default uniform arguments", {
  tr <- named_rtree(5)
  x <- c(1,4,1,22,2)
  names(x) <- tr$tip.label

  x.ilr <- philr(x, tree=tr, return.all=FALSE)

  expect_equivalent(philrInv(x.ilr, tr),
               miniclo(x))
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


test_that("philr handles data.frame input with warning", {
  tr <- named_rtree(5)
  x <- c(1,4,1,22,2)
  names(x) <- tr$tip.label
  xxxx <- as.data.frame(t(as.data.frame(x)))

  expect_warning(philr(xxxx, tree=tr, return.all=FALSE), "xxxx")
})


test_that("pseudocount works as expected", {

  tr <- named_rtree(5)
  pseudo <- 0.65
  x <- t(rmultinom(10,100,c(.1,.6,.2,.3,.2)))
  colnames(x) <- tr$tip.label  
  x.pseudo <- x + 0.65   # add a small pseudocount
  
  d1 <- philr(x.pseudo, tree=tr, part.weights='enorm.x.gm.counts',
             ilr.weights='blw.sqrt', return.all=TRUE)

  d2 <- philr(x, tree=tr, part.weights='enorm.x.gm.counts',
             ilr.weights='blw.sqrt', return.all=TRUE, pseudocount=pseudo)

  d3 <- philr(x.pseudo+5, tree=tr, part.weights='enorm.x.gm.counts',
             ilr.weights='blw.sqrt', return.all=TRUE)
	     
  expect_equal(d1, d2)
  expect_false(identical(d1,d3))

})
