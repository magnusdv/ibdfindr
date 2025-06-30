
test_that("findIBD is consistent on brothersX", {
  res = findIBD(brothersX, verbose = FALSE)
  expect_snapshot(res)
})

test_that("findIBD is consistent on cousinsDemo", {
  res = findIBD(cousinsDemo, verbose = FALSE)
  expect_snapshot(res)
})

test_that("findIBD w/thompson is consistent on cousinsDemo", {
  res = findIBD(cousinsDemo, thompson = TRUE, verbose = FALSE)
  expect_snapshot(res)
})
