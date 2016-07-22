context("Weighted ILR Calculations")

test_that('check.zeroes throws proper waring',{
  expect_warning(check.zeroes(c(1,1,0), 'x'), 'x should not contain zeroes')
})

test_that('Shiftp function handles both matricies and vectors', {
  x1 <- c(1,2,3)
  x2 <- matrix(c(1,2,3,1,2,3,1,2,3), nrow = 3, byrow=T)
  p <- c(1, .1, .5)

  expect_equal(shiftp(x1, p), matrix(c(1, 20, 6), nrow=1))
  expect_equal(shiftp(x2, p), matrix(rep(c(1, 20, 6), 3), nrow=3, byrow=T))
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

test_that('gp.rowMeans returns correct', {
  x <- matrix(c(2,4,4,2,4,4,2,4,4), nrow = 3, byrow=TRUE)
  p <- rep(1, 3)

  expect_equal(gp.rowMeans(x, p), rep(32^(1/3), 3))

  # Test defaults have not been changed
  expect_equal(gp.rowMeans(x), rep(32^(1/3), 3))
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
