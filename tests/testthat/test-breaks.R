test_that("integer_breaks always returns integers and includes 1", {
  for (m in c(1, 5, 9, 10, 100, 250, 1000, 2596, 12345)) {
    b <- integer_breaks(c(1, m))
    expect_type(b, "integer")
    expect_equal(b, as.integer(b))
    expect_equal(b[1], 1L)
    expect_lte(max(b), m)
  }
})

test_that("integer_breaks output for the spec'd ranges", {
  expect_equal(integer_breaks(c(1, 9)), c(1L, 5L, 9L))
  expect_equal(integer_breaks(c(1, 100)), c(1L, 20L, 40L, 60L, 80L, 100L))
  expect_equal(integer_breaks(c(1, 250)), c(1L, 50L, 100L, 150L, 200L, 250L))
  expect_equal(integer_breaks(c(1, 2596)), c(1L, 500L, 1000L, 1500L, 2000L, 2500L))
})

test_that("integer_breaks handles tiny ranges", {
  expect_equal(integer_breaks(c(1, 1)), 1L)
  expect_equal(integer_breaks(c(1, 2)), c(1L, 2L))
  expect_equal(integer_breaks(c(1, 3)), c(1L, 2L, 3L))
})
