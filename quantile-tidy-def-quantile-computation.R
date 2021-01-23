### quantile_main
# calculates the quantiles as defined in probs for the vector x according to 
# the type. Type 7 is special because of back-compatibility 
## inputs 
# @param x - the data from which to estimate the quantiles
# @param probs - the probability vectors for which the quantiles are estimated
# @param type - which algorithm is used for the estimation of quantiles 
## values
# a vector containing the computed quantiles 

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

quantile_type_7 <- function(x, probs) {
  length_x <- length(x)
  positions <- 1 + (length_x - 1) * probs
  lower_positions <- floor(positions)
  upper_positions <- ceiling(positions)
  relevant_positions <- unique(c(lower_positions, upper_positions))
  x <- sort(x, partial = relevant_positions)
  
  quantiles <- x[lower_positions]
  
  # require_weighting <- positions > lower_positions & x[upper_positions] != quantiles
  require_weighting <- x[upper_positions] != quantiles
  
  weights <- (positions[require_weighting] - lower_positions[require_weighting])
  # probably slightly more efficient than
  # weights <- (positions - lower_positions)[require_weighting]
  
  # in those cases Inf does not cause a problem because we cannot multiply
  # Inf * 0 
  quantiles[require_weighting] <- (1 - weights) * quantiles[require_weighting] + 
    weights * x[upper_positions[require_weighting]]
  
  quantiles
}

### quantile_general
# calculates quantiles for all type 1 ... 9 (type 7 works slightly different
# than in quantile_type_7)
## inputs 
# @param x - the data from which to estimate the quantiles
# @param probs - the probability vectors for which the quantiles are estimated
# @param type - which algorithm is used for the estimation of quantiles 

quantile_general <- function(x, probs, type) {
  length_x <- length(x)
  positions <- get_positions(probs, n = length_x, type)
  weights <- get_weights(probs, positions, type)
  relevant_positions <- get_relevant_positions(weights, positions, length_x)
  quantiles <- get_quantile(x, weights, positions, relevant_positions, length_x)
  quantiles
}

### get_weights
# defines the weights by calling the get_weights_discrete functions for type 1,2,3
# and calling the get_weights_conntinuous function for types 4...9
## inputs
# @param probs - he probability vectors for which the quantiles are estimated
# @param positions - the value k as defined in the pseudocode 
# @param type - which algorithm is used for the estimation of quantiles 
## value
# a vector of length probs with values in [0,1]
get_weights <- function(probs, positions, type) {
  if (type <= 3) return(get_weights_discrete(probs, positions, type))
  get_weights_continuous(probs, positions, type )
}

### get_weights_discrete
# defines weights for type 1,2,3
get_weights_discrete <- function(probs, positions, type) {
  positions_exact <- positions[["exact"]]
  positions_lower <- positions[["lower"]]
  weights <- switch(type, 
    positions_exact > positions_lower, 
    ((positions_exact > positions_lower) + 1) * 1/2, 
    (positions_exact != positions_lower) | ((positions_exact %% 2L) == 1L))
  weights
}

### get_weights_continuous 
# defines weights fot type 4, ..., 9
get_weights_continuous <- function(probs, positions, type) {
  weights <- positions[["exact"]] - positions[["lower"]]
  # set those weights that are almost 0 to 0 
  eps <- 4 * .Machine$double.eps
  almost_zero <- abs(weights) < eps
  weights[almost_zero] <- 0
  
  weights
}



get_relevant_positions <- function(weights, positions, n) {
  lower_positions <- positions[["lower"]]
  upper_positions <- lower_positions + 1
  relevant_lower <- weights < 1 
  relevant_upper <- weights > 0
  
  relevant_positions <- c(lower_positions[relevant_lower], 
    upper_positions[relevant_upper])
  
  in_range <- relevant_positions > 0 & relevant_positions <= n
  
  # we have to take care of the boundary cases
  # we actually only need 1 in case 0 is in relevant_positions and we
  # only need n in case n + 1 is in relevant_positions
  # but we can also simply add 1 and 0 in all cases (what's more efficient
  # probably depends on the case)
  
  relevant_positions <- unique(c(1, relevant_positions[in_range], n))
  
  relevant_positions
}

# as we can unify the calculation of positions of both - the discrete and the 
# continuous 

get_a_b <- function(type) {
  switch(type,
    {a <- 0; b <- 1}, 
    {a <- 0; b <- 1}, 
    {a <- -0.5; b <- 1.5},
    {a <- 0; b <- 1}, 
    {a <- b <- 0.5}, 
    {a <- b <- 0}, 
    {a <- b <- 1}, 
    {a <- b <- 1/3}, 
    {a <- b <- 3/8})
  list("a" = a, "b" = b)
}

get_positions <- function(probs, n, type) {
  intermediate <- get_a_b(type)
  a <- intermediate[["a"]] 
  b <- intermediate[["b"]]
  
  # actually we could also add the fuzz in case of a discrete quantile
  # this is somewhat dubious in the original quantile function
  eps <- 0 
  if (type > 3) eps <- 4 * .Machine$double.eps
  
  exact_positions <- a + probs * (n + 1 - a - b)
  lower_positions <- floor(exact_positions + eps)
  list(
    "exact" = exact_positions, 
    "lower" = lower_positions
  )
}




# this function is the same for both discrete and continuous
# what it does: it determines which positions are to be passed to the 
# partial argument of the sort function

# get_relevant_positions
# for both - discrete and continuous quantiles - the 


get_quantile <- function(x, weights, positions, relevant_positions, n) {
  lower_positions <- positions[["lower"]]
  x <- sort(x, partial = relevant_positions)
  # we extend x to the left and right, because for some algorithms the 
  # positions can become -1, 0, n + 1, n + 2
  x <- c(x[1L], x[1L], x, x[n], x[n])
  
  # Here we have to follow the code from the stats.quantile function 
  # exactly as otherwise we encounter different results in very specific situations
  # , e.g. if x is complex and contains Infinity 
  # quantiles[weights == 0] <- x[lower_positions + 2L][weights == 0]
  quantiles <- x[lower_positions + 2L]
  quantiles[weights == 1] <- x[lower_positions + 3L][weights == 1]
  

  require_weighting <- 
    (0 < weights) & (weights < 1) & x[lower_positions + 2L] != x[lower_positions + 3L]
  # note that the if is important to include the if condition as otherwise
  # we get problems for non-numeric x, (even empty character vector cannot
  # be multiplied)
  if (any(require_weighting)) {
    quantiles[require_weighting] <- ((1 - weights) * x[lower_positions + 2L] + 
        weights  * x[lower_positions + 3L])[require_weighting]
  }
  
  quantiles
}
