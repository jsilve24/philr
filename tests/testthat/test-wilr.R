context("Weighted ILR Calculations")

test_that('check.zeroes throws proper waring',{
  expect_warning(check.zeroes(c(1,1,0), 'x'), 'x should not contain zeroes')
})

test_that('Shift function handles both matricies and vectors', {
  x1 <- c(1,2,3)
  x2 <- matrix(c(1,2,3,1,2,3,1,2,3), nrow = 3, byrow=T)
  p <- c(1, .1, .5)

  expect_equal(shiftp(x1, p), matrix(c(1, 20, 6), nrow=1))
  expect_equal(shiftp(x2, p), matrix(rep(c(1, 20, 6), 3), nrow=3, byrow=T))
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
