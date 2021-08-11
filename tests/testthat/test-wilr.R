context("Weighted ILR Calculations")

test_that('check.zeroes throws proper waring',{
  expect_warning(check.zeroes(c(1,1,0), 'x'), 'x should not contain zeroes')
})

test_that('shiftp function handles both matricies and vectors', {
  x1 <- c(1,2,3)
  x2 <- matrix(c(1,2,3,1,2,3,1,2,3), nrow = 3, byrow=T)
  p <- c(1, .1, .5)

  expect_equal(shiftp(x1, p), matrix(c(1, 20, 6), nrow=1))
  expect_equal(shiftp(x2, p), matrix(rep(c(1, 20, 6), 3), nrow=3, byrow=T))
})

test_that('shiftpInv function reverses Shiftp', {
  x1 <- matrix(c(1,2,3), nrow=1)
  x2 <- matrix(c(1,2,3,1,2,3,1,2,3), nrow = 3, byrow=T)
  p <- c(1, .1, .5)
  y1 <- shiftp(x1, p)
  y2 <- shiftp(x2, p)

  expect_equal(shiftpInv(y1, p), x1)
  expect_equal(shiftpInv(y2, p), x2)
})

test_that('miniclo function handles both matricies and vectors', {
  x1 <- c(2,4,4)
  x2 <- matrix(c(2,4,4,2,4,4,2,4,4), nrow = 3, byrow=TRUE)

  expect_equal(miniclo(x1), matrix(x1/sum(x1), nrow = 1))
  expect_equal(miniclo(x2), x2/rowSums(x2))
  expect_equal(miniclo(data.frame(x2)), data.frame(x2/rowSums(x2)))
})

test_that('BuildilrBasep returns correct results', {
  sbp <- rbind(c(1,0), c(-1,1), c(-1,-1))
  #colnames(sbp) <- c('n1', 'n2')
  #rownames(sbp) <- c('otu1', 'otu2', 'otu3')
  p <- c(1, 1, 0.5)

  expect_equal(buildilrBasep(sbp, p),
                    rbind(c(0.7745967,  0.0000000),
                          c(-0.5163978, 0.5773503),
                          c(-0.5163978, -1.1547005)),
                    tolerance=1e-7)
})

test_that('g.rowMeans returns correct', {
  x <- matrix(c(2,4,4,2,4,4,2,4,4), nrow = 3, byrow=TRUE)
  p <- rep(1, 3)

  expect_equal(g.rowMeans(x, p), rep(32^(1/3), 3))

  # Test defaults have not been changed
  expect_equal(g.rowMeans(x), rep(32^(1/3), 3))
})

test_that('g.colMeans returns correct', {
  x <- matrix(c(2,4,4,2,4,4,2,4,4), nrow = 3, byrow=TRUE)

  expect_equal(g.colMeans(x), c(2,4,4))
})

test_that('normp returns zero for neutral element', {
  x <- c(2,2,2)
  p <- c(1,1,1)
  expect_equal(normp(x, p), 0)
})

test_that('ilrpInv reverses effect of ilrp', {
  p <- c(1,1,1,.1,.5)
  tr <- named_rtree(5)
  sbp <- phylo2sbp(tr)
  V <- buildilrBasep(sbp, p)
  x <- t(rmultinom(10,100,c(.1,.6,.2,.3,.2))) + 0.65   # add a small pseudocount
  x <- miniclo(x)
  y <- shiftp(miniclo(x), p)
  y.star <- ilrp(y, p, V)

  y.inversed <- ilrpInv(y.star, V)

  expect_equivalent(y.inversed, miniclo(y))
  expect_equivalent(miniclo(shiftpInv(y.inversed, p)), x)
})
