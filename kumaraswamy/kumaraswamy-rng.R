### load packages
library(checkmate)

### general notes
# this file contains two rng functions to sample from the kumaraswamy
# distribution:
# https://en.wikipedia.org/wiki/Kumaraswamy_distribution
## parameters
# - a, b are shape parameters
# - min and max define the support 


### rkumaraswamy1
# generates pseudo-random numbers from the kumaraswamy distribution 
# by inverse transform sampling (utilizing runif(0,1))
## inputs
# - parameters for the kumaraswamy distribution: a, b, min, max
# - n : sample size
## values
# a numeric vector of length n containing the generated pseudo-random sample
##
rkumaraswamy1 <- function(n, a, b, min = 0, max = 1) {
  checked_inputs <- check_rkumaraswamy(n = n, a = a, b = b, min = min, 
    max = max)
  do.call(rkumaraswamy1_main, checked_inputs)
}


### rkumaraswamy2
# generates pseudo-random numbers from the kumaraswamy distribution
# it uses the relationship that X^(1/a) is kumaraswamy(a, b) distributed if X is 
# beta(1, b) distributed
## inputs
# - parameters for the kumaraswamy distribution: a, b, min, max 
# - n : sample size
## values
# a numeric vector of length n containing the generated pseudo-random sample
##
rkumaraswamy2 <- function(n, a, b, min = 0, max = 1) {
  checked_inputs <- check_rkumaraswamy(n = n, a = a, b = b, min = min, 
    max = max)
  do.call(rkumaraswamy2_main, checked_inputs)
}





### rkumaraswamy1_main
# does the main calculation for the rkumaraswamy1 function without input
# checking
## inputs and values are the same as for rkumaraswamy1
##
rkumaraswamy1_main <- function(n, a, b, min, max) {
  # output length should always correspond to n, even if one of the parameters
  # is of length 0 
  if (min(lengths(list(a, b, min, max))) == 0) return(rep_len(NA_real_, n))
  
  is_na <- is.na(a) | is.na(b) | is.na(min) | is.na(max)
  invalid_par <- a <= 0 | b <= 0 | max <= min | is_na
  
  if (length(invalid_par) == 1) { # if all parameters have length 1
    n_valid <- n * !invalid_par # either n or 0
  } else {# in this case invalid_par has length n
    n_valid <- n - sum(invalid_par)
  }
  
  # need to watch out that recycled parameters are correctly aligned with 
  # the sample from the uniform distribution
  
  p <- rep_len(NA_real_, length.out = n)
  p[!invalid_par] <- runif(n_valid)
  qkumaraswamy_main(p, a, b, min, max)
}

### rkumaraswamy2_main
# does the main calculation for the rkumaraswamy2 function without input
# checking
## inputs and values are the same as for rkumaraswamy2
rkumaraswamy2_main <- function(n, a, b, min, max) {
  # output length should always correspond to n, even if one of the parameters
  # is of length 0 
  if (min(lengths(list(a, b, min, max))) == 0) return(rep_len(NA_real_, n)) 
  if (n == 0) return(numeric(0))
  
  is_na <- is.na(a) | is.na(b) | is.na(min) | is.na(max)
  invalid_par <- a <= 0 | b <= 0 | max <= min | is_na
  # invalid_par is either of length 1 (if all parameters have length 1)
  # or of length n
  x <- rep_len(NA_real_, n)
  
  if (all(is.na(invalid_par))) return(rep(NA_real_, n))
  
  x[!invalid_par] <- (rbeta(n = n, shape1 = 1, shape2 = b))[!invalid_par]
  x[!invalid_par] <- (x ^ (1 / a) * (max - min) + min)[!invalid_par]
  x
}


### qkumaraswamy_main
# calculates the quantiles for the kumaraswamy distribution without input 
# checking. The function is written only to be used within rkumaraswamy1
# where the quantile function is required for the inverse transform sampling
## inputs
# parameters for the kumaraswamy distribution: a, b, min, max
# p : the probability values for which the quantiles are evaluated
## values
# vector containing the quantiles for p 
##
qkumaraswamy_main <- function(p, a, b, min, max) {
  if (length(p) == 0) return(numeric(0))
  quantiles <- (1 - (1 - p) ^ (1 / b)) ^ (1 / a)
  range <- max - min
  quantiles * range + min
}

### check_rkumaraswamy
# checks the inputs for the pseudo-random-number generators for the 
# kumaraswamy distribution
# in addition to that it recycles/shortenes parameters with length > 1 to 
# length n
## inputs
# see rkumaraswamy1 and rkumaraswamy2
## values
# a named list containing all the checked (and potentially modified) inputs 
##
check_rkumaraswamy <- function(n, a, b, min, max) {
  n <- check_n(n)
  checked_parameters <- check_parameters(n = n, a = a, b = b, min = min, 
    max = max, rng = TRUE)
  append(checked_parameters, list("n" = n))
}


### check_parameters
# this checks whether the parameters a, b, min and max are correctly specified
# (taking recycling into consideration)
# parameters with length > 1 are recycled to length n to avoid potential 
# recycling issues later 
## inputs:
# parameters for the kumaraswamy distribution: a, b, min, max
# n : length of x (for d/pkumaraswamy) or n (for rkumaraswamy1/2)
# rng : logical indicating whether we test the inputs for the rng or the 
# p/dkumaraswamy function (depends on how parameters are recycled, i.e. 
# either to the maximal length of x, a, b, min, max or n in the rng case 
## values
# a named list containing a, b, min, max in which all parameters with 
# length > are recycled to n_max 
##
check_parameters <- function(n, a, b, min, max, rng = FALSE) {
  # note that the ifelse version of the following
  # line causes an error if n is numeric(0) 
  # consider: ifelse(TRUE, numeric(0), 1) (return is NA)
  n_max <- if (rng) n else max(n, lengths(list(a, b, min, max)))
  
  # in the rng case it can happen that n is numeric(0) 
  if (length(n_max) == 0) n_max <- 0
  
  checked_a_b <- check_a_b(a = a, b = b) 
  checked_min_max <- check_min_max(min = min, max = max, n_max = n_max, 
    rng = rng)
  
  parameters <- append(checked_min_max, checked_a_b)
  
  # in case one parameter is of length 0 output only contains NA anyway
  # and therefore no further processing required 
  if (min(lengths(parameters)) == 0) return(parameters)
  
  is_par_vectorized <- any(lengths(parameters) > 1)
  
  if (is_par_vectorized && !rng) {
    warning("at least one of a, b, min or max have length > 1. 
    All parameters, as well as x are recycled to the maximal length of 
      x, a, b, min and max")
  }
  if (is_par_vectorized && rng) {
    warning("at least one of a, b, min or max have length > 1.
    All parameters are recycled/ shortened to length n")
  }
  
  ### THIS PART IS VERY IMPORTANT 
  # we recycle all parmaters with length > 1 to length n_max
  # this is an effective way to solve recycling issues later
  # for parameters with length 1 we let R internal recycling do the job
  
  parameters <- sapply(parameters, simplify = FALSE,
    FUN = function(x) if (length(x) > 1) rep_len(x, n_max) else x) 
  
  parameters
}


### check_min_max
## NOTE: earlier versions of this function used the least common denominator
# to optimize input checking for vectorized min and max. This made the 
# code more efficient but also more compliated and was therefore removed 

# Is used to check max and min. Note that a simple max > min check does not 
# suffice in case their only common divisor is 1 and one of them has length 
# > 1 (and both of them length > 1)
# To which length they are recycled depends 
# on whether we have an rng or density/distribution function. 
# rng-case: parameters are recycled to length n
# density/distribution-case: all parameters are recycled to the maximal length 
# of x, a, b, min, max 

## inputs
# n_max - the length to which a, b and min and max will be recycled
# max - max parameter from kumaraswamy
# min - min parameter from kumaraswamy
# rng - is the function used within the rng-functions
## outputs
# list containing min and max. If both have length > 1 they are recycled
# to length n_max (if only one has length > 1 there are no issues with a simple
# max > min test)
##
check_min_max <- function(min, max, n_max, rng = FALSE) {
  check_numeric_arg(min, "min")
  check_numeric_arg(max, "max")
  
  length_min <- length(min)
  length_max <- length(max)
  
  # we only need to be 
  # careful when testing max > min when both of them are of length > 1
  ## this can be made more efficient with least common multiple but KIS
  if (length_min > 1 & length_max > 1) {
    min <- rep_len(min, length.out = n_max)
    max <- rep_len(max, length.out = n_max)
  }
  
  invalid_max_min <- any(max <= min, na.rm = TRUE)
  
  if (invalid_max_min & !rng) {
    warning("max <= min for at least one value for the recycled vectors")
  }
  # slightly different warning in the rng case because there parameters can 
  # be shortened
  if (invalid_max_min & rng) {
    warning("max <= min for at least one value for the recycled/ shortened
      vectors")
  }
  
  list("min" = min, "max" = max)
}

### check_ab 
# checks whether a and b are correctly specified and gives informative
# warnings 
## inputs
# a, b : parameters for the kumaraswamy distribution
## values
# a list with a and b
##
check_a_b <- function(a, b) {
  check_numeric_arg(a, "a")
  check_numeric_arg(b, "b")
  
  if (any(a <= 0, na.rm = TRUE)) warning("a contains at least one value <= 0")
  if (any(b <= 0, na.rm = TRUE)) warning("b contains at least one value <= 0")
  
  list("a" = a, "b" = b)
}


### check_logical_arg
# this function checks whether an element can be reasonably be interpreted as 
# logical and gives tailored warning messages. It gives back the converted
# TRUE / FALSE of length 1 if possible, and stops otherwise 
# SIDE-NOTE:
# for the warning messages to work properly, the argument arg must be 
# a variable and not an expression (e.g. check_logical_arg(c(1,0)) gives
# weird warning messages but a <- c(1,0); check_logical_arg(a) works fine)
## inputs
# arg : the argument value
## value
# logical of length 1 
## 
check_logical_arg <- function(arg) {
  # an early exit which will cover almost all use-cases 
  if (test_logical(arg, len = 1, any.missing = FALSE)) {
    return(arg)
  }
  
  assert(
    check_logical(arg), 
    check_numeric(arg), 
    combine = "or"
  )
  assert_true(!is.na(arg[1])) # this ensures that arg has length > 1
  # AND that first entry is non NA. 
  # numeric(0) case is also covered: consider is.na(numeric(0)[1])
  
  arg_name <- substitute(arg) # for warning-message
  
  if (length(arg) > 1) {
    warning(paste(arg_name, "has length > 1 and only first element is used"))
    arg <- arg[1]
  }
  
  if (is.numeric(arg)) {
    arg <- as.logical(arg)
    warning(paste(arg_name, "is numeric and is coerced to logical"))
  }
  arg
}

### check_n
# checks the input n for rkumaraswamy1/2
# n must be numeric of length 0 or 1 and as.integer(n) == n must be TRUE
## inputs
# n : the value n as passed to rkumaraswamy1/2
## values
# n (potentially) converted to an interger
##
check_n <- function(n) {
  assert_numeric(n, any.missing = FALSE, max.len = 1, lower = 0, finite = TRUE)
  
  if (length(n) == 0) { # this is the case where n = numeric(0)
    warning("n has length 0")
    n <- 0 
    return(n)
  } 
  
  n_int <- as.integer(n)
  assert_true(n_int == n)
  n_int
}

### check_numeric_arg
# checks whether can be properly interpreted as a numeric arg
## inputs:
# arg : the argument to be checked
# arg_name : the name of the argument that is checked, if NULL
# substitute(arg) is used as the name (does not work when check_numeric_arg
# is nested inside other functions)
## values
# arg
##
check_numeric_arg <- function(arg, arg_name = NULL) {
  if (is.null(arg_name)) arg_name <- substitute(arg)
  assert(
    check_logical(arg), 
    check_numeric(arg), 
    combine = "or"
  )
  if (is.logical(arg)) warning(paste(arg_name, "is logical"))
  if (any(is.na(arg))) warning(paste(arg_name, "contains NAs and they are
    not removed"))
  if (length(arg) == 0) warning(paste(arg_name, "has length 0"))
  if (any(is.infinite(arg))) warning(paste(arg_name, "contains Inf"))
  arg
}