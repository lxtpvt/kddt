test_that("criterions", {
  expect_equal(criterion(c(1,2,3),c(1,3,1),type=1), 5/3)
})
