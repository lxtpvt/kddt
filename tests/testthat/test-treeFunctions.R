test_that("Get all ancestors of a node by its id", {
  expect_equal(getRids(5), c(5,2,1))
})

test_that("Get all node's ids of a complete tree by its depth", {
  expect_equal(getAllNidByLevel(3), c(1:15))
})
