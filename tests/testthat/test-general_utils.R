context('tree_node_id_conversion functions')

tr.fixed <- ape::read.tree(text="(((t1:0.9480410812,(t4:0.5734868464,
                           t3:0.5969958168)n4:0.5359996252)n3:0.2657330192,
                           t2:0.6697232309)n2:0.8426657049,t5:0.1043712671)n1;")
tr.rand <- ape::rtree(10)
tr.rand <- ape::makeNodeLabel(tr.rand, method="number", prefix='n')

test_that('all functions handle vectors on random tree', {
  # on vectors of length 1
  expect_length(nn.to.name(tr.rand, 1), 1)
  expect_length(name.to.nn(tr.rand, 'n1'), 1)

  # on longer vectors
  expect_length(nn.to.name(tr.rand, c(1, 2, 11)), 3) # nodes and tips
  expect_length(name.to.nn(tr.rand, c('n1', 'n2', 'n3', 't1')), 4) # nodes and tips
})

test_that('all functions give correst results on fixed tree', {
  expect_equal(nn.to.name(tr.fixed, c(1,2,8)), c("t1", "t4", "n3"))
  expect_equal(name.to.nn(tr.fixed, c("t1", "t4", "n3")), c(1,2,8))
})

test_that('get.ud.*** throws error if given a tip', {
  expect_error(get.ud.tips(tr.fixed, 't1'), "t1 is a tip")
  expect_error(get.ud.nodes(tr.fixed, 't1'), "t1 is a tip")
})

test_that("vec_to_mat correctly handles named vectors", {
  a <- c("a"=5, "b"=2)
  b <- vec_to_mat(a)
  expect_equal(colnames(b), names(a))
  expect_equal(class(b), "matrix")
})
