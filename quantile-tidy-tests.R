library(testthat)
set.seed(31415926)

source("quantile-tidy-def.R")

# the testing might seem somewhat complicated at first but the idea is simple
# in ordered to test many argument combinations to also check for errors
# that result in very specific combinations of arguments we used
# expand.grid to create these argument lists 

# Then we test each all of the rows of the data-frames
# If an error throws, the corresponding row number is given so 
# one can immediately see which comination of arguments caused the error 

quantile_default <- stats:::quantile.default

### helper_equivalent
# checks equivalence of quantile_tidy and quantile.default
helper_equivalent <- function(x, probs = p, na.rm = FALSE, names = TRUE, 
  type = 7) {
  expect_equivalent( 
    quantile_tidy(x, probs, na.rm, names, type),
    quantile_default(x, probs, na.rm, names, type))
}

# this utility function makes it better readable to loop over the rows of
# dataframes that contain arguments for which we want to test 
helper_loop <- function(arguments, message = "",
  fun = helper_equivalent) {
  for (i in 1:NROW(arguments)) {
    args <- arguments[i,]
    # we need this apply to reduce the list of lists to a list 
    args <- lapply(args, FUN = function(x) x[[1]])
    test_that(paste(message, "row no.", i), {
      do.call(fun, args = args)
    })
  }
}

# now we define the arguments that will be used for the testing 
x <- rnorm(20)
x_complex <- complex(real = x, imaginary = rnorm(20))
p <- c(0, runif(8), 1)
x_date <- as.Date("2020-01-01, 1918-02-03, 1950-02-3")
x_ordered <- ordered(x)


# in the quantile.default function there are some bugs with respect to names 
# that are fixed in the quantile_tidy function, therefore only allow 
# names = TRUE 
args_numeric_1 <- expand.grid(
  "x" = list(x, c(x, NA), c(x, NA, Inf, -Inf), as.matrix(x), as.array(x), 
    c(x, rep(Inf, times = 10)), c(x, rep(1, times = 40))), 
  "probs" = list(p, c(p, NA), NULL, c(TRUE, FALSE), NA), 
  "na.rm" = list(TRUE), 
  "names" = list(TRUE), 
  "type" = as.list(1:9)
)

# if probs does not contain any NAs we can allow TRUE and FALSE for names
args_numeric_2 <- expand.grid(
  "x" = list(x, c(x, NA), c(x, NA, Inf, -Inf), as.matrix(x), as.array(x)), 
  "probs" = list(p, NULL, c(TRUE, FALSE)), 
  "na.rm" = list(TRUE), 
  "names" = list(TRUE, FALSE), 
  "type" = as.list(1:9)
)

args_complex <- expand.grid(
  "x" = list(x_complex, c(x_complex, Inf + 3i), c(x_complex, -Inf)), 
  "probs" = list(p, c(p, NA)), 
  "na.rm" = list(TRUE, FALSE), 
  "type" = as.list(1:9)
)

args_null <- expand.grid(
  "x" = list(NULL, 1:10), 
  "probs" = list(seq(0, 01, 0.25), NULL), 
  "names" = list(TRUE, FALSE), 
  "na.rm" = list(TRUE, FALSE), 
  "type" = as.list(1:9)
)

args_ord_dat_char <- expand.grid(
  "x" = list(x_date, c(x_date, NA, NA), x_ordered, 
    letters, c(letters, NA)), 
  "probs" = list(p, c(p, NA), NULL, c(TRUE, FALSE)), 
  "na.rm" = list(TRUE), 
  "names" = list(TRUE), 
  "type" = as.list(1, 3)
)

args_logicals <- expand.grid(
  "x" = x, 
  "probs" = seq(0, 1, 0.25), 
  "na.rm" = list(0, 1 + 1i), 
  "names" = list(TRUE, c(-1, 1))
)

# for these values it is actually relevant that we add the fuzz/eps when 
# defining the weights/positions 
args_rounding_issues <- expand.grid(
  "x" = list(c(1, 1e20)), 
  "probs" = 
    list(c(1 / 2 + 2 * .Machine$double.eps, 1 / 2 + 1 * .Machine$double.eps)), 
  "type" = as.list(1:9))


helper_loop(args_numeric_1)
helper_loop(args_numeric_2)
helper_loop(args_complex)
helper_loop(args_null)
helper_loop(args_ord_dat_char)
helper_loop(args_logicals)
helper_loop(args_rounding_issues)

test_that("correct outputs for x or probs == NULL", {
  helper_equivalent(x = NULL, type = 1)
  helper_equivalent(x = 1:10, probs = NULL, type = 3)
  helper_equivalent(x = NULL, type = 7)
})

test_that("breaks when x contains NAs but na.rm = FALSE", {
  expect_error(quantile_tidy(c(NA, 1:10), na.rm = FALSE))
})

test_that("probs limits works better than in the default function", {
  expect_error(
    quantile_tidy(x = 1:10, probs = 1 + 101 * .Machine$double.eps)
  )
  expect_error(
    quantile_tidy(x = 1:10, probs = -99 * .Machine$double.eps), NA
  )
  expect_error(
    quantile_default(x = 1:10, probs = -.Machine$double.eps)
  )
  expect_error(
    quantile_default(x = 1:10, probs = c(NA, -.Machine$double.eps)), NA
  )
  expect_equivalent(
    quantile_tidy(
      x = x, 
      probs = c(-50 * .Machine$double.eps, 1 + 50 * .Machine$double.eps)), 
    quantile(x = x, probs = c(0, 1))
  )
})



test_that("does not fail for probs c(NA, 0.5) and names == FALSE
  like quantile.default does", {
  expect_error(
    quantile_tidy(x, probs = c(NA, 0.5), names = FALSE), NA)
  expect_error(
    quantile_default(x, probs = c(NA, 0.5), names = FALSE)
  )
})
  
test_that("throws error if x is character or ordered factor and type is not
  in c(1L, 2L) as opposed to quantile.default", {
    expect_error(
      quantile_tidy(x = letters, type = 2)
    )
    expect_error(
      quantile_tidy(x = x_ordered, type = 4)
    )
    
    expect_error(
      quantile_tidy(x = letters, type = 1), NA
    )
    expect_error(
      quantile_tidy(x = x_ordered, type = 3), NA
    )
    
  })

# for type 4 a = 0,  = 1, therefore positions = probs * n 
# and lower_positions = floor(probs * n + 4 * .Machine$double.eps)
# and the weights that are below 4 * .Machine$double.eps are rounded to 0 
test_that("adding the fuzz when defining the weights works", {
  helper_equivalent(
    x = c(1, 1e20), 
    probs = c(1 / 2 + 2 * .Machine$double.eps, 1 / 2 + 1 * .Machine$double.eps), 
    type = 4)
  helper_equivalent(
    x = c(1, 1e20), 
    probs = c(1 / 2 + 2 * .Machine$double.eps, 1 / 2 + 1 * .Machine$double.eps), 
    type = 7)
})

test_that("lists are not accepted", {
  expect_error(quantile_tidy(x = list(1:3)))
  expect_error(quantile_tidy(1:10, probs = list(0,1)))
})

test_that("quantile_tidy handles types that are not in 1...9 correctly (and
  different than quantile.default)", {
    expect_equivalent(
      quantile_tidy(x, type = 6 + 0.9 * sqrt(.Machine$double.eps)), 
      quantile_tidy(x, type = 6)
    )
    expect_error(quantile_tidy(x, type = 7.9))
})

test_that("type works as desired", {
  expect_error(quantile_tidy(x = 1:10, type = c(1, 2)))
})

test_that("na.rm and names must be coercible to logical", {
  expect_error(quantile_tidy(1:10, na.rm = "a"))
  expect_error(quantile_tidy(1:10, names = "b"))
  expect_warning(quantile_tidy(1:10, names = 100))
  expect_warning(quantile_tidy(1:10, names = c(TRUE, FALSE)))
  expect_warning(quantile_tidy(1:10, na.rm = c(FALSE, TRUE)))
  expect_equivalent(
    quantile_tidy(1:10, probs = p, names = 1:3),
    quantile_tidy(1:10, probs = p,  names = TRUE)
  )
})

test_that("na.rm works", {
  expect_equivalent(
    quantile_tidy(c(rep(NA_real_, times = 10), x), na.rm = TRUE), 
    quantile_tidy(x, na.rm = FALSE)
  )
  expect_warning(quantile_tidy(NA_real_, na.rm = TRUE))
})
  
test_that("works if lower_positions == upper_positions ", {
  expect_equivalent(
    quantile_tidy(x, probs = 1:20 / 20, type = 4), 
    quantile(x, probs = 1:20 / 20, type = 4))
})

test_that("na.rm and names cannot be NA", {
  expect_error(quantile_tidy(1:10, na.rm = as.logical(NA)))
  expect_error(quantile_tidy(1:10, names = as.logical(NA)))
})

  
