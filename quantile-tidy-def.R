
library(checkmate)


### quantile_tidy
# a tidy reimplementation for the quantile.default function
# see quantile-improved.pdf for a detailed description
## inputs
# @param x - the data from which to estimate the quantiles
# @param probs - the probability vectors for which the quantiles are estimated
# @param type - which algorithm is used for the estimation of quantiles 
# @param na.rm - should NAs be removed from x 
# @param names - should the output have the probabilities, to which the 
# quantile corresponds, as a name
## values
# a (named) vector containing the quantiles
##
quantile_tidy <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE,
  names = TRUE, type = 7L, ...) {
  checked_inputs <- check_inputs(x, probs, na.rm, names, type)
  
  for (z in c("x", "probs", "na.rm", "names", "type", "probs_original", 
    "probs_NA", "x_levels")) {
    assign(z, checked_inputs[[z]])
  }
  
  # calculate quantiles
  if (length(x) && length(probs)) {
    quantiles <- quantile_main(x = x, probs = probs, type = type)
  } else {
    quantiles <- rep(NA_real_, times = length(probs))
  }
  
  # modify the output and return it 
  modify_quantiles(quantiles = quantiles, x_levels = x_levels, probs = probs, 
    probs_NA = probs_NA, probs_original = probs_original, names = names)
}

#-----functions for the quantile estimation-------------------------------------

### quantile_main
# calculates the quantiles as defined in probs for the vector x according to 
# the type. Type 7 is special because of back-compatibility 
## inputs 
# @param x - the data from which to estimate the quantiles
# @param probs - the probability vectors for which the quantiles are estimated
# @param type - which algorithm is used for the estimation of quantiles 
## values
# a vector containing the computed quantiles 
##
quantile_main <- function(x, probs, type) {
  if (type == 7L) return(quantile_type_7(x, probs))
  quantile_general(x, probs, type)
}

### quantile_type_7
# this works a bit different than the other quantile algorithms 
# for the purpose of back-compatibility. It calculates the quantiles defined
# in probs for the x vector according to type 7 
## inputs
# @param x - the data from which to estimate the quantiles
# @param probs - the probability vectors for which the quantiles are estimated
## value
# a vector containing the computed quantiles 
##
quantile_type_7 <- function(x, probs) {
  length_x <- length(x)
  positions <- 1 + (length_x - 1) * probs
  lower_positions <- floor(positions)
  upper_positions <- ceiling(positions)
  relevant_positions <- unique(c(lower_positions, upper_positions))
  x <- sort(x, partial = relevant_positions)
  
  # initializ quantiles as the 'lower neighbour'
  quantiles <- x[lower_positions]
  
  # require_weighting are the ones for which the position lies between 
  # lower_positions and upper_positions
  require_weighting <- x[upper_positions] != quantiles
  
  # weights are in (0,1) 
  weights <- (positions[require_weighting] - lower_positions[require_weighting])
  # more efficient than
  # weights <- (positions - lower_positions)[require_weighting]
  
  # in those cases Inf cannot cause a problem because the weights are 
  # in (0,1)
  quantiles[require_weighting] <- (1 - weights) * quantiles[require_weighting] + 
    weights * x[upper_positions[require_weighting]]
  
  quantiles
}

### quantile_general
# calculates quantiles for all types 1 ... 9 
## inputs 
# @param x - the data from which to estimate the quantiles
# @param probs - the probability vectors for which the quantiles are estimated
# @param type - which algorithm is used for the estimation of quantiles 
## values
# the estimated quantiles
##
quantile_general <- function(x, probs, type) {
  length_x <- length(x)
  # get positions and lower positions
  positions <- get_positions(probs, n = length_x, type)
  # how are the lower and upper neighbours weighted 
  weights <- get_weights(probs, positions, type)
  # which positions are relevant for the partial sort
  relevant_positions <- get_relevant_positions(weights, positions, length_x)
  # calcualte quantiles
  get_quantile(x, weights, positions, relevant_positions, length_x)
}

### get_weights
# defines the weights by calling the get_weights_discrete functions for type 
# 1,2,3 and calling the get_weights_conntinuous function for types 4...9
## inputs
# @param probs - the probability vectors for which the quantiles are estimated
# @param positions - list with the exact position and the lower position
# (lower position corresponds to k in the help page and exact position 
# to a + probs * (n + 1 - a - b))
# @param type - which algorithm is used for the estimation of quantiles 
## value
# a vector of length probs with values in [0,1]
##
get_weights <- function(probs, positions, type) {
  if (type <= 3) return(get_weights_discrete(probs, positions, type))
  get_weights_continuous(probs, positions, type)
}

### get_weights_discrete
# defines weights for type 1,2,3
get_weights_discrete <- function(probs, positions, type) {
  positions_exact <- positions[["exact"]]
  positions_lower <- positions[["lower"]]
  switch(type, 
    positions_exact > positions_lower, 
    ((positions_exact > positions_lower) + 1) * 1 / 2, 
    (positions_exact != positions_lower) | ((positions_exact %% 2L) == 1L))
}

### get_weights_continuous 
# defines weights fot type 4, ..., 9
get_weights_continuous <- function(probs, positions, type) {
  weights <- positions[["exact"]] - positions[["lower"]]
  # set those weights that are almost 0 to 0 
  eps <- 4 * .Machine$[["double.eps"]]
  almost_zero <- abs(weights) < eps
  weights[almost_zero] <- 0
  
  weights
}


### get relevant positions
# determines which positions are relevant for the partial sort
get_relevant_positions <- function(weights, positions, n) {
  lower_positions <- positions[["lower"]]
  upper_positions <- lower_positions + 1
  relevant_lower <- weights < 1 
  relevant_upper <- weights > 0
  
  relevant_positions <- c(lower_positions[relevant_lower], 
    upper_positions[relevant_upper])
  
  # some positions can be outside the range as described 
  # above 
  in_range <- relevant_positions > 0 & relevant_positions <= n
  
  # we have to take care of the boundary cases
  # we actually only need 1 in case 0 or -1 is in relevant_positions and we
  # only need n in case n + 1/ n + 2 is in relevant_positions

  unique(c(1, relevant_positions[in_range], n))
}


### get_a_b
# get the values a b for discrete and continuous quantiles
get_a_b <- function(type) {
  switch(type,
    {a <- 0; b <- 1}, 
    {a <- 0; b <- 1}, 
    {a <- -0.5; b <- 1.5},
    {a <- 0; b <- 1}, 
    {a <- b <- 0.5}, 
    {a <- b <- 0}, 
    {a <- b <- 1}, 
    {a <- b <- 1 / 3}, 
    {a <- b <- 3 / 8})
  list("a" = a, "b" = b)
}

## get_positions
# get the exact and lower positions that are relevant for determining
# the relevant positions for the partial sort and the weights 
get_positions <- function(probs, n, type) {
  intermediate <- get_a_b(type)
  a <- intermediate[["a"]] 
  b <- intermediate[["b"]]
  
  # actually we could also add the fuzz in case of a discrete quantile
  # this is somewhat arbitrary in the original quantile function
  eps <- 0 
  if (type > 3) eps <- 4 * .Machine[["double.eps"]]
  
  exact_positions <- a + probs * (n + 1 - a - b)
  lower_positions <- floor(exact_positions + eps)
  list(
    "exact" = exact_positions, 
    "lower" = lower_positions
  )
}



### get_quantile
# caculates the quantiles, once all relevant information was extracted
get_quantile <- function(x, weights, positions, relevant_positions, n) {
  lower_positions <- positions[["lower"]]
  x <- sort(x, partial = relevant_positions)
  # we extend x to the left and right, because for some algorithms the 
  # positions can become -1, 0, n + 1, n + 2
  x <- c(x[1L], x[1L], x, x[n], x[n])
  
  # Here we have to follow the code from the stats.quantile function 
  # exactly as otherwise we encounter different results in very specific 
  # situations, e.g. if x is complex and contains Infinity 
  
  # the + 2L is relevant because of the potential -1, 0, n + 1, n + 2 in
  # the positions, it ensures that the quantiles for positions <= 0 
  # are assigned the smallest value in x and quantiles for positions >= n + 1
  # are assigned the largest value in x
  
  quantiles <- x[lower_positions + 2L]
  quantiles[weights == 1] <- x[lower_positions + 3L][weights == 1]
  
  require_weighting <- 
    (0 < weights) & (weights < 1) & 
    x[lower_positions + 2L] != x[lower_positions + 3L]
  
  # note that the if is important to include the if condition as otherwise
  # we get problems for non-numeric x, (even empty character vector cannot
  # be multiplied)
  if (any(require_weighting)) {
    quantiles[require_weighting] <- ((1 - weights) * x[lower_positions + 2L] + 
        weights  * x[lower_positions + 3L])[require_weighting]
  }
  
  quantiles
}


#------functions for input checking---------------------------------------------

### check_inputs 
# checks whether the inputs of the quantile function are of admissible type and
# transforms them to standard form
## inputs:
# - all the arguments from the quantile function except for ...
## value
# a list with 
# x - x without NAs
# x_levels - levels of x 
# probs - probs without NAs 
# probs_original - original probs 
# probs_NA - logical vector indicating which of the original probs is NA
# type - type as an integer (if type was integerish before)
# na.rm - na.rm coerced to logical and only first element
# names - names coerced to logicak and only first element
##
check_inputs <- function(x, probs, na.rm, names, type) {
  type <- check_type(type)
  na.rm <- check_logical_arg(na.rm)
  names <- check_logical_arg(names)
  checked_probs <- check_probs(probs)
  checked_x <- check_x(x, type , na.rm )
  
  list(
    "x" = checked_x[["x"]], 
    "x_levels" = checked_x[["x_levels"]],
    "probs" = checked_probs[["probs"]], 
    "probs_original" = checked_probs[["probs_original"]], 
    "probs_NA" = checked_probs[["probs_NA"]],
    "type" = type, 
    "names" = names
  )
}

### check_logical_arg
# this function checks whether an element can be reasonably interpreted as 
# logical and gives tailored warning messages and gives back the converted
# TRUE / FALSE of length 1 if possible, and stops otherwise 
# for the warning messages to work properly, the argument arg must be 
# a variable and not an expression (e.g. check_logical_arg(c(1,0)) gives
# weird warning messages but a <- c(1,0); check_logical_arg(a) works fine)
# only to be used in this context !!!
## inputs
# @param arg - the argument value
## value
# logical of length 1 
##
check_logical_arg <- function(arg) {
  arg_name <- substitute(arg) # for warning/ stop-message
  # an early exit which will cover almost all uses of the function
  if (test_logical(arg, len = 1, any.missing = FALSE)) {
    return(arg)
  }
  
  # otherwise we require arg to be atomic with length > 0
  # and to be convertible to type logical 
  assert(
    check_atomic(arg), 
    check_true(length(arg) > 0), 
    combine = "and"
  )
  
  if (length(arg) > 1) {
    warning(paste(arg_name, "has length > 1 and only first element is used"))
    arg <- arg[1]
  }
  
  assert_true(!is.na(arg))
  
  # now we deal with the case in which the first element of args is non-logical
  # but also non-NA
  if (!is.logical(arg)) {
    arg <- as.logical(arg)
  }
  
  if (!is.na(arg)) { # if the as.logical was not successfull arg is now NA 
    warning(paste(arg_name, "was coerced to logical"))
  } else {
    stop(paste(arg_name, "could not be meaningfully interpreted as logical"))
  }
  
  arg
}


### check_type
# this function checky whether type is an integerish value in 1:9 and converts
# it if possible and gives the converted value back, and breaks otherwise
## inputs
# type - the value which is to be checked
## value
# an integer in 1:9
##
check_type <- function(type) {
  assert_integerish(type, len = 1, any.missing = FALSE, lower = 1, upper = 9)
  if (!is.integer(type)) {
    type_original <- type
    type <- as.integer(round(type)) 
    
    if (type_original != type) {
      warning(paste0("type is almost an integer (with a tolerance of ", 
        sqrt(.Machine[["double.eps"]]), ") and was rounded to the next integer"))
    }
  }
  type
}

### check_probs
# checks whether all values in probs are within (-eps, 1 + eps), converts
# values in (-eps, 0) to 0, values in (1, 1 + eps) to 1, 
# breaks if values lie outside (-eps, 1 + eps)
## arguments:
# @param probs the vector that supposedly contains probabilities 
## values
# the checked probs with values in (-eps, 0) roundedto 0 and values in 
# (1, 1 + eps) rounded to 1
##
check_probs <- function(probs) {
  assert(
    check_numeric(probs), 
    check_logical(probs),
    check_null(probs), 
    combine = "or"
  )
  
  # early exit if probs is NULL
  if (is.null(probs)) {
    warning("probs is NULL")
    return(
      list(
        "probs" = NULL, 
        "probs_original" = NULL, 
        "probs_NA" = logical(0))
    )
  }
  
  eps <- 100 * .Machine[["double.eps"]]
  probs_original <- probs
  probs_NA <- is.na(probs)
  probs <- probs[!probs_NA]
  probs_in_range <- probs >= -eps & probs <= 1 + eps
  
  if (!all(probs_in_range)) {
    stop("probs contains values outside of [0,1]")
  }
  
  # correct values
  probs <- pmax(0, pmin(probs, 1))
  
  if (!length(probs)) {
    warning("probs only contains NAs")
    # here we could also introduce an early exit, this would make the code
    # more complicated without much gain as this case rarely happens 
    # in addition to that one might also want to get the other warnings
    # with respect to x
  }
  
  list(
    "probs" = probs, 
    "probs_original" = probs_original, 
    "probs_NA" = probs_NA)
}

## check_x
# checks whether x is specified correctly, which depends on the type and na.rm
check_x <- function(x, type, na.rm) {
  ## early exit if x is NULL
  if (is.null(x)) {
    warning("x is NULL")
    return(
      list(
        "x" = NULL, 
        "x_levels" = NULL
      )
    )
  }
  
  assert(
    check_numeric(x), 
    check_complex(x),
    check_logical(x), 
    check_character(x), 
    check_date(x),
    check_factor(x, ordered = TRUE),
    combine = "or"
  )
  
  if ((is.factor(x) || is.character(x) || is(x,"Date")) && 
      !(type %in% c(1L, 3L))) {
    stop("ordered factors, Dates and characters are only allowed for type 1 or 
      3")
  }
  
  x_levels <- levels(x) # if x is not a factor this is NULL
  
  x_NA <- is.na(x)
  if (any(x_NA) && !na.rm) {
    stop("x contains NA but na.rm was specified to FALSE")
  }
  if (any(x_NA)) {
    x <- x[!x_NA]
    warning("NAs in x were removed")
  }
  if (!length(x)) {
    warning("x only contains NA values")
  }
  
  return(list(
    "x" = x, 
    "x_levels" = x_levels
  ))
}


### format_quantile_names
# here we have to be carful with the probs includes NAs and 
# names == TRUE bug in the original function
format_quantile_names <- function(probs) {
  digits <- max(2L, getOption("digits"))
  # display probs as percents, therefore multiply by 100
  probs <- probs * 100
  if (!length(probs)) return(character(0))
  
  # function depends on length of probs for efficiency 
  if (length(probs) < 100) {
    paste0(formatC(probs, format = "fg", width = 1, digits = digits), "%")
  } else {
    paste0(format(probs, trim = TRUE, digits = digits), "%")
  }
}

### modify_quantiles
# deal with factors and NAs
modify_quantiles <- function(quantiles, x_levels, probs, probs_NA, probs_original, 
  names) {
  if (!is.null(x_levels)) { # in case x is a ordered factor
    quantiles <- factor(quantiles, levels = seq_along(x_levels),
      labels = x_levels, ordered = TRUE)
  }
  
  # in case probs contains NA values we have to extend the returned 
  # quantile-vector with the NAs
  if (any(probs_NA)) {
    quantiles_including_NA <- numeric(length(probs_original))
    quantiles_including_NA[!probs_NA] <- quantiles
    quantiles_including_NA[probs_NA] <- NA_real_
    quantiles <- quantiles_including_NA
  }
  
  # we have to check for the length of probs as well, otherwise 
  # we get different results than stats:::quantile in case probs is NULL
  # and names is TRUE 
  if (names & length(probs)) {
    names(quantiles)[!probs_NA] <- format_quantile_names(probs)
    names(quantiles)[probs_NA] <- ""
  }
  quantiles
}




