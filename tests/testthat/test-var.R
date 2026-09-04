# pa_var input validation, including warnings on type-irrelevant arguments.

test_that("pa_var validates its inputs", {
  v <- pa_var("treat", "binary", p = 0.5)
  expect_s3_class(v, "pa_var")
  expect_error(pa_var("d", "binary"), "strictly between")
  expect_error(pa_var("d", "binary", p = 1), "strictly between")
  expect_error(pa_var("d", "binary", p = c(0.4, 0.6)), "strictly between")
  expect_error(pa_var("z", "normal", sd = 0), "sd > 0")
  expect_error(pa_var("", "binary", p = 0.5))
  expect_error(pa_var("d", "binary", p = 0.5, icc = 1.5), "icc")
})

test_that("pa_var warns when arguments irrelevant to the type are supplied", {
  expect_warning(pa_var("z", "normal", p = 0.3), "ignored for the normal")
  expect_warning(pa_var("d", "binary", p = 0.5, sd = 2), "ignored for the binary")
  expect_warning(pa_var("d", "binary", p = 0.5, mean = 1), "ignored for the binary")
  expect_silent(pa_var("z", "normal", mean = 2, sd = 3))
  expect_silent(pa_var("d", "binary", p = 0.5))
  expect_silent(pa_var("d", "binary", p = 0.5, icc = 0.6))
})
