context('Branch Length Calculations')

library(ape)
tr.fixed <- read.tree(text="((t2:0.3708901014,((t1:0.5037344745,
                           t3:0.5299885271)n4:0.7873590044,
                           t5:0.4477706181)n3:0.8255127743)n2:0.1799265754,
                           t4:0.9563850679)n1;")

test_that('mean_dist_to_tips gives correct results on fixed tree', {
  mdtt <- mean_dist_to_tips(tr.fixed)

  expect_equal(length(mdtt), 4)
  expect_equal(mdtt[3], 1.018737, tolerance=3e-7, check.attributes=F)
  expect_equivalent(mdtt[4], 0.5168615)
})

test_that('calculate.blw handles tip with edge length zero and gives correct results', {
  tr <- read.tree(text="((t1:0.3187429151,t2:0):0.3409723381,
                       t3:0.04305356415);")

  expect_warning(calculate.blw(tr), 'Note: a total of 1') # throws warning

  # replaces with pseudocount of min
  expect_equal(suppressWarnings(calculate.blw(tr, method='sum.children')),
               c(0.3840259, 0.3617965),
               tolerance=1e-7)
  expect_equal(suppressWarnings(calculate.blw(tr, method='mean.descendants')),
               c(0.5649241, 0.3617965),
               tolerance=1e-7)
})
