library(checkmate)
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
# this function checks whether an element can be meaninfully interpreted as 
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
## arguments_
# type - the value which is to be checked
## value: an integer 
check_type <- function(type) {
  assert_integerish(type, len = 1, any.missing = FALSE, lower = 1, upper = 9)
  if (!is.integer(type)) {
    type_original <- type
    type <- as.integer(round(type)) 
    
    if (type_original != type) {
      warning(paste0("type is almost an integer (with a tolerance of ", 
        sqrt(.Machine$double.eps), ") and was rounded to the next integer"))
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
  
  eps <- 100 * .Machine$double.eps
  probs_original <- probs
  probs_NA <- is.na(probs)
  probs <- probs[!probs_NA]
  probs_in_range <- probs >= -eps & probs <= 1 + eps
  
  if (!all(probs_in_range)) {
    stop("probs contains values outside of [0,1]")
  }
  
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
    stop("ordered factors, Dates and characters are only allowed for type 1 or 3")
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


# here we have to be carful with the probs includes NAs and 
# names == TRUE bug 
# note that ... can not overwrite arguments like digits
# in addition to that we can also speciy getOption within the function 
# body, as it is a global Option that does not depend on the environment
# (defaults are evaluated in the environment from which the function 
# is called, as opposed to assignments in the function body)
format_quantile_names <- function(probs) {
  digits <- max(2L, getOption("digits"))
  probs <- probs * 100
  if (!length(probs)) return(character(0))
  
  if (length(probs) < 100) {
    paste0(formatC(probs, format = "fg", width = 1, digits = digits), "%")
  } else {
    paste0(format(probs, trim = TRUE, digits = digits), "%")
  }
}

modify_quantiles <- function(quantiles, x_levels, probs, probs_NA, probs_original, 
  names) {
  if (!is.null(x_levels)) { # in case x is a ordered factor
    quantiles <- factor(quantiles, levels = seq_along(x_levels),
      labels = x_levels, ordered = TRUE)
  }
  
  # in case probs contains NA values we have to extend the returned 
  # quantile-vector 
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
