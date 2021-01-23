library(testthat)
source("kumaraswamy-def.R")
set.seed(2718)

context("test dkumaraswamy and pkumaraswamy")

# as the tests for pkumaraswamy and dkumaraswamy overlap to a large
# extent we put some tests in a wrapper 

### test_kumaraswamy_general
# is a wrapper around some tests that are valid for both, the density and the
# distribution function
##
test_kumaraswamy_general <- function(type) {
  if (type == "density") {
    kumaraswamy <- dkumaraswamy
    name <- "dkumaraswamy"
    beta <- dbeta
  }
  if (type == "distribution") {
    kumaraswamy <- pkumaraswamy 
    name <- "pkumaraswamy"
    beta <- pbeta
  }

  
  test_that(paste(name, "works for sensical inputs: output matches
    the output of the corresponding beta distribution"), {
    expect_equivalent(
      kumaraswamy(v, a = 1, b = 1:20), 
      beta(v, shape1 = 1, shape2 = 1:20))
    expect_equivalent( 
      kumaraswamy(v, a = 1:20, b = 1), 
      beta(v, shape1 = 1:20, shape2 = 1))
    # complex recycling with NAs and negative a 
    expect_equivalent(
      kumaraswamy(c(NA, -Inf, seq(0.1, 0.99, by = 0.01), NA), 
        a = c(NA, 2, 1, NA, -1), b = c(NA, 1, NA)), 
      beta(c(NA, -Inf, seq(0.1, 0.99, by = 0.01), NA), 
        shape1 = c(NA, 2, 1, NA, -1), shape2 = c(NA, 1, NA)))
  })
  
  test_that(paste(name, "gives correct warnings"), {
    # the warnings for the logical 
    expect_warning(kumaraswamy(TRUE, 1, 1), "x is logical")
    expect_warning(kumaraswamy(0.5, TRUE, 2), "a is logical")
    expect_warning(kumaraswamy(0.5, 1, TRUE), "b is logical")
    expect_warning(kumaraswamy(0.5, 1, 1, FALSE), "min is logical")
    expect_warning(kumaraswamy(0.5, 1, 1, 0, TRUE), "max is logical")
    # warnings for NA_real_ (partial matching works with expect_warning)
    expect_warning(kumaraswamy(NA_real_, 1, 1), "x contains NA")
    expect_warning(kumaraswamy(0.5, NA_real_, 1), "a contains NA")
    expect_warning(kumaraswamy(0.5, 1, NA_real_), "b contains NA")
    expect_warning(kumaraswamy(0.5, 1, 1, NA_real_), "min contains NA")
    expect_warning(kumaraswamy(0.5, 1, 1, 0, NA_real_), "max contains NA")
    # warnings for Inf
    expect_warning(kumaraswamy(Inf, 1, 1), "x contains Inf")
    expect_warning(kumaraswamy(Inf, -Inf, 1), "a contains Inf")
    expect_warning(kumaraswamy(0.5, 1, Inf), "b contains Inf")
    expect_warning(kumaraswamy(0.5, 1, 1, -Inf), "min contains Inf")
    expect_warning(kumaraswamy(0.5, 1, 1, 0, Inf), "max contains Inf")
    # warnings for input of length 0 
    expect_warning(kumaraswamy(numeric(0), 1, 1), "x has length 0")
    expect_warning(kumaraswamy(0.5, numeric(0), 1), "a has length 0")
    expect_warning(kumaraswamy(0.5, 1, numeric(0)), "b has length 0")
    expect_warning(kumaraswamy(0.5, 1, 1, numeric(0)), "min has length 0")
    expect_warning(kumaraswamy(0.5, 1, 1, 0, numeric(0)), "max has length 0")
    # warnings for a or b <= 0
    expect_warning(kumaraswamy(0.5, -1, 1), 
      "a contains at least one value <= 0")
    expect_warning(kumaraswamy(0.5, 1, -1), 
      "b contains at least one value <= 0")
    expect_warning(kumaraswamy(0.5, 1, 2, 3, 2))
  })
  
  test_that(paste(name, "gives correct errors"), {
    # type is checked
    expect_error(dumaraswamy("hallo", 1, 1))
    expect_error(dumaraswamy(0.5, "hallo", 1))
    expect_error(kumaraswamy(0.5, 1, 1, "hallo"))
    # NULL is not allowed
    expect_error(kumaraswamy(NULL, 1, 1))
    expect_error(kumaraswamy(0.5, 1, NULL))
    expect_error(kumaraswamy(0.5, 1, 1, NULL))
    # lists are not allowed
    expect_error(kumaraswamy(list(0.5), 1, 1))
  })
  

  test_that(paste(name, "gives correct output for edge-cases"), {
    expect_equivalent(
      kumaraswamy(numeric(0), NA, 1), 
      numeric(0))
    expect_equivalent(
      kumaraswamy(1, numeric(0), 1), 
      numeric(0))
    expect_true(is.na(kumaraswamy(0.5, NA, 1)))
    expect_true(is.na(kumaraswamy(NA, 1, 2)))
    expect_true(is.na(kumaraswamy(1, 1, 1, NA)))
    expect_equivalent(
      kumaraswamy(0.5, -1, 1), 
      dbeta(0.5, -1, 1)
    )
    expect_equivalent(
      kumaraswamy(0.5, 1, -1), 
      dbeta(0.5, 1, -1)
    )
  })
}




### test dkumaraswamy  
v <- runif(10)
test_kumaraswamy_general("density")

v <- runif(10)
test_that("min and max arguments work properly for dkumaraswamy", {
  expect_equal(
    dkumaraswamy(0.5, 0.5, 0.5, min = 0, max = 1),
    dkumaraswamy(1.5, 0.5, 0.5, min = 1, max = 2))
  expect_equal(# here we use the transformation theorem
    dkumaraswamy(v, 0.3, 0.7), 
    dkumaraswamy(v * 6 - 3, 0.3, 0.7, min = -3, max = 3) * 6
  )
})


test_that("dkumaraswamy gives correct warnings for log", {
  expect_warning(dkumaraswamy(0.5, 1, 1:2), 
    "at least one of a, b, min or max have length > 1. 
    All parameters, as well as x are recycled to the maximal length of 
      x, a, b, min and max")
  expect_warning(dkumaraswamy(0.5, 1, 1, c(0, 0)), 
    "at least one of a, b, min or max have length > 1. 
    All parameters, as well as x are recycled to the maximal length of 
      x, a, b, min and max")
  expect_warning(dkumaraswamy(0.5, 1, 1, 0, 1, c(TRUE, FALSE)), 
    "log has length > 1 and only first element is used")
  expect_warning(dkumaraswamy(0.5, 1, 1, 0, 1, 1), 
    "log is numeric and is coerced to logical")
})


test_that("dkumaraswamy is a density ", {
  expect_true(# integrates to 1 on the support 
    abs(integrate(dkumaraswamy, lower = 0, upper = 1, a = 1, 
      b = 2)[["value"]] - 1) < sqrt(.Machine[["double.eps"]]))
  # is 0 elsewhere 
  expect_true(all(dkumaraswamy(seq(-10, -0.01, by = 0.01), 1, 1) == 0))
  expect_true(all(dkumaraswamy(seq(1.01, 11, by = 0.01), 2, 3) == 0))
  # density is positive on the support 
  expect_true(all(dkumaraswamy(seq(0.01, 0.99, by = 0.01), 1, 2) > 0))
  expect_true(all(dkumaraswamy(c(-Inf, Inf), 1, 1) == c(0, 0)))
})

test_that("log argument works for dkumaraswamy", {
  expect_error(dkumaraswamy(0.5, 1, 1, log = NA))
  expect_error(dkumaraswamy(0.5, 1, 1, log = logical(0)))
  expect_error(dkumaraswamy(0.5, 1, 1, log = NULL))
  expect_error(dkumaraswamy(0.5, 1, 1, log = c(NA, TRUE)))
  expect_warning(dkumaraswamy(0.5, 1, 1, log = c(TRUE, FALSE)))
  expect_warning(dkumaraswamy(0.5, 1, 1, log = 0))
  expect_error(dkumaraswamy(0.5, 1, 1, log = "hallo"))
})

test_that("testing", {
  expect_equivalent(
    dkumaraswamy(0.5, 1, 1, min = c(0, 1), max = c(1, -1)), 
    c(1, NA_real_)
  )
})


# test that integrating the density concides with the values from pkumaraswamy
test_that("distribution function and density conincide", {
  expect_equal(
    integrate(dkumaraswamy, lower = 0, upper = 0.3, a = 1, b = 2)[["value"]], 
    pkumaraswamy(0.3, 1, 2), 
    tol = 0.01)
})

### testing pkumaraswamy 

test_kumaraswamy_general("distribution")

v <- runif(10)
test_that("min and max arguments work properly for pkumaraswamy", {
  expect_equal(# translations don't matter
    pkumaraswamy(0.5, 0.5, 0.5, min = 0, max = 1),
    pkumaraswamy(1.5, 0.5, 0.5, min = 1, max = 2))
  expect_equal(# can also stretch the support + translate
    pkumaraswamy(v, 0.3, 0.7), 
    pkumaraswamy(v * 6 - 3, 0.3, 0.7, min = -3, max = 3) 
  )
})

x1 <- seq(0.1, 0.9, by = 0.01)
test_that("distribution function is monotonous and 0/1 for x -> -Inf/Inf", {
  expect_numeric(pkumaraswamy(x1, a = 2, b = 3), sorted = TRUE)
  expect_true(all(pkumaraswamy(-100:0, a = 2, b = 3) == 0))
  expect_true(all(pkumaraswamy(1:100, a = 3, b = 10) == 1))
  expect_equivalent(
    pkumaraswamy(c(-Inf, Inf), 1, 2), 
    c(0, 1))
})

test_that("logical arguments are tested correctly for pkumaraswamy", {
  expect_error(pkumaraswamy(0.5, 1, 1, lower.tail = NA))
  expect_error(pkumaraswamy(0.5, 1, 1, lower.tail = logical(0)))
  expect_error(pkumaraswamy(0.5, 1, 1, lower.tail = NULL))
  expect_error(pkumaraswamy(0.5, 1, 1, lower.tail = c(NA, TRUE)))
  expect_warning(pkumaraswamy(0.5, 1, 1, lower.tail = c(TRUE, FALSE)))
  expect_warning(pkumaraswamy(0.5, 1, 1, lower.tail = 0))
  expect_error(pkumaraswamy(0.5, 1, 1, lower.tail = "hallo"))
  
  expect_error(pkumaraswamy(0.5, 1, 1, log.p = NA))
  expect_error(pkumaraswamy(0.5, 1, 1, log.p = logical(0)))
  expect_error(pkumaraswamy(0.5, 1, 1, log.p = NULL))
  expect_error(pkumaraswamy(0.5, 1, 1, log.p = c(NA, TRUE)))
  expect_warning(pkumaraswamy(0.5, 1, 1, log.p = c(TRUE, FALSE)))
  expect_warning(pkumaraswamy(0.5, 1, 1, log.p = 0))
  expect_error(pkumaraswamy(0.5, 1, 1, log.p = "hallo"))
})

test_that("lower.tail works for pkumaraswamy", {
  expect_equal(
    pkumaraswamy(x1, 0.4, 0.2, lower.tail = FALSE),
    1 - pkumaraswamy(x1, 0.4, 0.2))
})

test_that("log.p works for pkumaraswamy", {
  expect_equal(
    pkumaraswamy(x1, 0.4, 0.2, log.p = TRUE), 
    log(pkumaraswamy(x1, 0.4, 0.2)))
})


