### load packages
library(checkmate)

### general notes
# this file contains the density and distribution function for the
# kumaraswamy distribution:
# https://en.wikipedia.org/wiki/Kumaraswamy_distribution
## parameters
# - a, b are shape parameters
# - min and max define the support 
## the implementation mostly follows the implementations of other distributions
# such as the normal or the beta distribution in e.g. the way vectorization
# of parameters (not only of x) works. 
# this allows to simulate situations in which parameters change conveniently
# it comes at the price of a higher complexity of input checking and some
# the main difference between this implementation and standard R implementations
# of distributions is that explicit input checking is conducted

### dkumaraswamy
# density function for the kumaraswamy distribution
# all arguments are vectorizable other than log (only first value is relevant)
# x, a, b, min and max are all recycled to their maximal length 
# density at the boundaries is not 0 to avoid possible problems during
# ML estimation 
## inputs:
# - parameters for the kumaraswamy distribution: a, b, min, max
# - x : vector where the density is to be evaluated
# - log : if TRUE the log density is returned
## value
# (log-) densities of the kumaraswamy distribution 
##
dkumaraswamy <- function(x, a, b, min = 0, max = 1, log = FALSE) {
  checked_inputs <- check_dkumaraswamy(x = x, a = a, b = b, min = min, 
    max = max, log = log)
  do.call(dkumaraswamy_main, checked_inputs)
}

### pkumaraswamy 
# distribution function for the kumaraswamy distribution
# All arguments but lower.tail and log.p are vectorizable (only first value
# is relevant)
# x, a, b, min, max all all recycled to their maximal length 
## inputs
# - parameters for the kumaraswamy distribution: a, b, min, max
# - x : a numeric vector where the distribution is to be evaluated
# - lower.tail : logical value indicating whether the return is F(x) or 1 - F(x)
# - log.p : logical value indicating whether the log of the probability
# should be returned
## values
# numeric vector containing (log of) F(x) (or 1 - F(x))
##
pkumaraswamy <- function(x, a, b, min = 0, max = 1, lower.tail = TRUE, 
  log.p = FALSE) {
  checked_inputs <- check_pkumaraswamy(x = x, a = a, b = b, min = min, 
    max = max, lower.tail = lower.tail, log.p = log.p)
  do.call(pkumaraswamy_main, checked_inputs)
}




### dkumaraswamy_main
# does the main calculations for dkumaraswamy without checking inputs
## inputs and values are as in dkumaraswamy
##
dkumaraswamy_main <- function(x, a, b, min, max, log) {
  # early exit in case there is a numeric(0)
  if (min(lengths(list(x, a, b, min, max))) == 0) return(numeric(0))

  # three different cases of assigning values: 
  # 1. invalid/NA parameters
  # 2. x has valid parameters and is outside the support
  # 3. x has valid parameters and is inside the support 
  
  ## NOTE
  # x, a, b, min and max are already either of length n_max or 
  # of length 1, so no recycling issues 
  
  ## initialize densities vector
  n_max <- max(lengths(list(x, a, b, min, max)))
  densities <- rep_len(NA_real_, n_max) 
  
  info <- get_assignment_info(x = x, a = a, b = b, min = min, max = max)
  x_scaled <- info[["x_scaled"]] # x_scaled = (x - min)/(max - min)
  is_ok <- info[["is_ok"]] # is_ok means no NA and valid parameters 
  outside_support <- info[["outside_support"]] # x is outside support
  # and has valid parameters
  range <- info[["range"]] # range = max - min
  
  densities[outside_support] <- 0
  
  densities[!outside_support & is_ok] <- 
    ((a * b * x_scaled ^ (a - 1) * (1 - x_scaled ^ a) ^ (b - 1)) * 
      1 / range)[!outside_support & is_ok]
  
  if (log) return(log(densities))
  densities
}

### pkumaraswamy_main
# does the main calculations for pkumaraswamy without checking inputs
## inputs and values are as in pkumaraswamy
##
pkumaraswamy_main <- function(x, a, b, min, max, lower.tail, log.p) {
  # early exit in case there is a numeric(0)
  if (min(lengths(list(x, a, b, min, max))) == 0) return(numeric(0))
  # now the other case in which all have length > 0 
  
  # three different cases of assigning values: 
  # 1. invalid/NA parameters
  # 2. x has valid parameters and is outside the support
  # 3. x has valid parameters and is inside the support 
  
  ## NOTE
  # x, a, b, min and max are already either of length n_max or 
  # of length 1, so no recycling issues 
  
  ## initialize probabilities vector
  n_max <- max(lengths(list(x, a, b, min, max)))
  probabilities <- rep_len(NA_real_, n_max) 
  
  info <- get_assignment_info(x = x, a = a, b = b, min = min, max = max)
  x_scaled <- info[["x_scaled"]] # (x - min)/(max - min)
  is_ok <- info[["is_ok"]] # valid parameter and not NA
  outside_support <- info[["outside_support"]] # outside support & valid 
  # parameter
  
  # assing values below min and above max depending on lower.tail
  c <- ifelse(lower.tail, 1, 0)
  probabilities[is_ok & (x_scaled < 0)] <- (1 - c)
  probabilities[is_ok & (x_scaled > 1)] <- c
  
  # assign values for valid parameters inside the support
  if (lower.tail) {
    probabilities[!outside_support & is_ok] <- 
      (1 - (1 - x_scaled^a)^b)[!outside_support & is_ok]
  } else {
    probabilities[!outside_support & is_ok] <- 
      ((1 - x_scaled^a)^b)[!outside_support & is_ok]
  }

  if (log.p) return(log(probabilities))
  
  probabilities
}

### check_dkumaraswamy
# checks the inputs for the density function of the kumaraswamy distribution
# parameters are (in case all of them are of length > 0) recycled to the same
# length (in case one of them is 0 this is not required because output 
# only contains NA)
# log is shortened to length 1 
## inputs:
# all the values that are passed to the dkumaraswamy function
## value: 
# a named list containing all the checked inputs 
##
check_dkumaraswamy <- function(x, a, b, min, max, log) {
  check_numeric_arg(x) # x is never modified 
  checked_parameters <- check_parameters(a = a, b = b, min = min, max = max, 
    n = length(x))
  log <- check_logical_arg(log)
  # x is included in the returned list for a clean do.call at the top-level
  append(checked_parameters, list("log" = log, "x" = x))
}

### check_pkumaraswamy
# checks the inputs for the kumaraswamy distribution function
# parameters are (in case all of them are of length > 0) recycled to the same
# length (in case one of them is 0 this is not required because output 
# only contains NA)
# log is shortened to length 1
## inputs:
# all the values that are passed to the pkumaraswamy function
## value: 
# a named list containing all the checked inputs 
##
check_pkumaraswamy <- function(x, a, b, min, max, lower.tail, log.p) {
  check_numeric_arg(x) # x is never modified 
  checked_parameters <- check_parameters(a = a, b = b, min = min, max = max, 
    n = length(x))
  lower.tail <- check_logical_arg(lower.tail)
  log.p <- check_logical_arg(log.p)
  # x is included in the returned list for a clean do.call at the top-level
  append(checked_parameters, list("lower.tail" = lower.tail, "log.p" = log.p, 
    "x" = x))
}

### get_assignment_info
# this function is a simple helper that is used when dealing with NAs and values 
# outside the support to assign the density/ distribution values 
#### ATTENTION
# at this point the recycling issues come into play. it is  important, that in 
# case 
##
get_assignment_info <- function(x, a, b, min, max) {
  n_max <- max(lengths(list(x, a, b, min, max)))
  range <- max - min
  
  is_na <- is.na(x) | is.na(a) | is.na(b) | is.na(min) | is.na(max) # is_na
  # already has output length
  invalid_par <- a <= 0 | b <= 0 | max <= min
  is_ok <- !(is_na | invalid_par)
  
  ## deal with valid parameters
  x_scaled <- (x - min) / range
  
  outside_support <- is_ok & !(x_scaled >= 0 & x_scaled <= 1)
  
  list(
    "is_ok" = is_ok,
    "outside_support" = outside_support, 
    "x_scaled" = x_scaled,
    "range" = range,
    "n_max" = n_max
  )
}


### check_parameters
# this checks whether the parameters a, b, min and max are correctly specified
# (taking recycling into consideration)
# parameters with length > 1 are recycled to length n to avoid potential 
# recycling issues later 
## inputs:
# parameters for the kumaraswamy distribution: a, b, min, max
# n : length of x (for d/pkumaraswamy) or n (for rkumaraswamy1/2)
# rng : logical indicating whether we test the inputs for the the rng or the 
# p/dkumaraswamy function (depends on how parameters are recycled, i.e. 
# either to the maximal length of x, a, b, min, max or n in the rng case 
## values
# a named list containing a, b, min, max in which all parameters with 
# length > are recycled to n_max 
##
check_parameters <- function(n, a, b, min, max, rng = FALSE) {
  # note that the ifelse version of the following
  # line causes an error if n is numeric(0) 
  # consider: ifelse(TRUE, numeric(0), 1) (return is of length 1)
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
# code - also at later stages - more complicated and was simplified. 

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
# list containing min and max. Parameters with length > 1 are recycled
# to length n_max
##
check_min_max <- function(min, max, n_max, rng = FALSE) {
  check_numeric_arg(min, "min")
  check_numeric_arg(max, "max")
  
  length_min <- length(min)
  length_max <- length(max)
  
  # we only need to be careful when testing max > min when both of them are of 
  # length > 1
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
# this function checks whether an element can be reasonably interpreted as 
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




