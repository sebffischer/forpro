#  File src/library/stats/R/quantile.R
#  Part of the R package, https://www.R-project.org
#
#  Copyright (C) 1995-2017 The R Core Team
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  https://www.R-project.org/Licenses/


quantile <- function(x, ...) UseMethod("quantile")

quantile.POSIXt <- function(x, ...)
  .POSIXct(quantile(unclass(as.POSIXct(x)), ...), attr(x, "tzone"))

# @param x is the sample on the basis of which the quantiles will be estimated
# @param probs is a numeric vector specifying which quantiles are to be 
# estimated
# @param na.rm specifies whether NAs should be removed
# @param names specifies whether the output should be a named vector (the names
# being the corresponding quantiles) or not
# @param type specifies which of the 9 possible algorithms will be used
# for the quantile estimation
# @param ... possible additional arguments that are not used however 
quantile.default <-
  function(x, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE,
    type = 7, ...)
  {
    #----------------this block is only for factors----------------------------#
    
    # here we check that if x is a factor it is also a ordered factor 
    # and that the specified type is either 1 or 3. non-orderd factors are not
    # allowed because there is no sensible sorting that can be done 
    # and therefore no distribution function in the usual sense is defined, 
    # this implies that it also does not make sense to estimate quantiles. 
    # we require the type to be 1 or 3 for ordered factors because these are 
    # the only types for which no multiplication is ever neccessary. 
    # (for the other types it depends on the exact values of x and p so  
    # allowing them for ordered factors would mean that the function would 
    # work inconsistently, meaning sometimes it would give an error and 
    # sometimes not

    if(is.factor(x)) {
      if(is.ordered(x)) {
        if(!any(type == c(1L, 3L)))
          stop("'type' must be 1 or 3 for ordered factors")
      } else
        stop("factors are not allowed")
      lx <- levels(x) # for ordered factors this will be needed later as the 
      # levels will be also assigned to the output, i.e. the quantile vector
    } else lx <- NULL
    
    #--------------------------------------------------------------------------#
    
    #-----------------------deal with NAs in x---------------------------------#
    
    if (na.rm) # remove NAs if na.rm == TRUE
      x <- x[!is.na(x)]
    else if (anyNA(x)) # function breaks if data contains any NAs but user 
      # specified na.rm = FALSE (which is the default)
      stop("missing values and NaN's not allowed if 'na.rm' is FALSE")
    
    #---------deal with NAs in probs and values that are not probabilities-----#
    
    # eps: tolerance margin for which values in probs
    # are still accepted (i,e. [0 - eps, 1 + eps]) which is important when 
    # the input probs is effected by rounding errors
    eps <- 100*.Machine$double.eps # assign tolerance value 
    # in case there is at least one value in probs that is non-NA and outside 
    # the interval defined above the function breaks
    if (any((p.ok <- !is.na(probs)) & (probs < -eps | probs > 1+eps)))
      stop("'probs' outside [0,1]") # break if the values of probs are not
    # probabilities (within a reasonable margin of error)
    # p.ok is a logical vector indicating which of the values in probs are 
    # non-NA
    
    n <- length(x) # define the length of the vector (where NAs are already
    # removed)
    if(na.p <- any(!p.ok)) { # in case there is at least one NA in probs
      o.pr <- probs # save the original probs vector as specified by the user
      # this will be needed later to modify the output 
      probs <- probs[p.ok] # remove the NA-elements in probs
      probs <- pmax(0, pmin(1, probs)) # set values in [-eps,0] to 0
      # and values in [1, 1 + eps] to 1
      # probs is now a vector without NAs and all values lie within [0,1]
    }
    # note that in case probs values are in (-eps, 0) or (1, 1 + eps) 
    # and probs does not contain NAs those values are not corrected
    
    #--------------------------------------------------------------------------#
    
    #-------------------the actual algorithm for the quantiles-----------------#
    
    np <- length(probs) # number of elements in probs without the NAs
    if (n > 0 && np > 0) { # in case there is at least one non-NA element in 
      # x and probs respectively 
      if(type == 7) {# in the default case do
        # this works a bit differently than the other continuous quantiles
        # in order to ensure back-compatibility
        
        # note: I also use variable names as specified in the quantile
        # help page where the general form of the quantile function is defined    
        # this means that y is the same variable as h 
        
        # we need the variable index to define y and to find j and j + 1 
        index <- 1 + (n - 1) * probs 
        # lo corresponds to j in the help page and hi[i] will either correspond 
        # to j + 1 (in case is.integer(index[i]) == FALSE) and otherwise
        # it holds lo[i] == hi[i] (then we also don't need x[j + 1]) as
        # y == 0
        lo <- floor(index) # j
        hi <- ceiling(index) # j + 1 (special case: y == 0: hi == j)
        
        # we only do a partial sorting for efficiency 
        # Note that if i is in partial, the i-th order
        # statistic of the partially sorted x will coincide with the i-th 
        # order statistic of the fully sorted x and in addition to that
        # all the values x[1], ..., x[i-1] of the partially sorted vector are
        # smaller than x[i] and x[i + 1], ..., x[n] are larger 
        # note that the partial argument is different for type 7 compared
        # to the other types, the reason for this is that calculation of quantiles
        # for the other types makes use of x[j] != x[j + 1] so we really need
        # hi <- lo + 1, instead of hi <- ceiling(index) as defined here 
        # in addition to that we do not need to add 1 and n to the partial sort
        # argument because index is always in [1, n]
        x <- sort(x, partial = unique(c(lo, hi))) 
        # now we initialize all quantiles as x[j] (j from help-page)
        qs <- x[lo]
        
        # in i we store the indices that still require updating, i.e. 
        # Q(p) is not  x[j] but (1 - y) * x[j] + y * x[j] with y != 0
        
        # note that the way i is defined here is unncessarily complicated as:
        # x[hi][i] != qs[i] => index[i] > index[lo]
        # so we would be good with i <- which(x[hi] != qs)
        i <- which(index > lo & x[hi] != qs) # '!=' for '>' working w/ complex
        # h is the y as defined in the help-page and let's us weight
        # x[j] and x[j + 1]
        h <- (index - lo)[i] # > 0	by construction
        ##	    qs[i] <- qs[i] + .minus(x[hi[i]], x[lo[i]]) * (index[i] - lo[i])
        ##	    qs[i] <- ifelse(h == 0, qs[i], (1 - h) * qs[i] + h * x[hi[i]])
        # now we update the values for which y (i.e. h) is in (0,1]
        qs[i] <- (1 - h) * qs[i] + h * x[hi[i]]
      } else { # in case type is not 7 
        if (type <= 3) {
          ## Types 1, 2 and 3 are discontinuous sample qs.
          # nppm will be used to get the indices j and the weights 
          nppm <- if (type == 3) n * probs - .5 # n * probs + m; m = -0.5
          else n * probs          # m = 0
          # j is the j as defined in the help page
          j <- floor(nppm)
          # h again corresponds to the y in the description and is the vector
          # that contains how x[j] and x[j + 1] are to be weighted 
          # to estimate the quantiles 
          h <- switch(type, # update weights according to type
            # a great visualization for these weights can be found in the 
            # Sample quantiles in statistical packages paper that is also 
            # referenced in the help page 
            (nppm > j),		# type 1
            ((nppm > j) + 1)/2, # type 2
            (nppm != j) | ((j %% 2L) == 1L)) # type 3
        } else {
          ## Types 4 through 9 are continuous sample qs.
          # all the types define the nppm (that let's us define j and h, which 
          # corresponds to y in the help page)
          # via the same formula: a + probs * (n + 1 - a - b)
          # but with different values for a and b, so here we assign the values
          # a and b according to the type 
          switch(type - 3, # type - 3 means that 4 becomes 1 etc. 
            {a <- 0; b <- 1},    # type 4
            a <- b <- 0.5,   # type 5
            a <- b <- 0,     # type 6
            a <- b <- 1,     # type 7 (unused here)
            a <- b <- 1 / 3, # type 8
            a <- b <- 3 / 8) # type 9
          ## need to watch for rounding errors here#
          fuzz <- 4 * .Machine$double.eps # fuzz is used to account for 
          # rounding errors 
          nppm <- a + probs * (n + 1 - a - b) # n*probs + m
          # as mentioned above nppm let's us define j and y/h
          
          # if the calculation of nppm is affected by rounding errors this could
          # mean that floor(nppm) as calculated by machine arithmacy is 
          # one integer larger than when calculated exactly. 
          j <- floor(nppm + fuzz) # m = a + probs*(1 - a - b)
          # again j corresponds to the j in the help page
          
          # define the weights 
          h <- nppm - j
          # set those weights that are almost 0 to 0 
          # if adding the fuzz makes a difference for the floor function 
          # in the lines above 
          # then the weight h is also set to zero. This improves the condition 
          # of the problem (see condition number for multiplication)
          if(any(sml <- abs(h) < fuzz)) h[sml] <- 0
        }
        # j and j + 1 are defined as in the help page. 
        # note that we need to include the 1, because if we e.g. set type = 3
        # and probs = 0 j = c(-1) and j + 1 = c(0) but the first element 
        # has to be sorted
        # analogously we need to include n e.g. for the case type == 6 and probs 
        # = 1 as then j == n + 1 and j + 1 == n + 2, however the estimated 
        # quantile is the n-th order statistic and would otherwise not be sorted 
        # correctly
        x <- sort(x, partial =
            unique(c(1, j[j>0L & j<=n], (j+1)[j>0L & j<n], n))
        )
        # this concatenation is done for the reason as just described 
        x <- c(x[1L], x[1L], x, x[n], x[n])
        ## h can be zero or one (types 1 to 3), and infinities matter
        ####        qs <- (1 - h) * x[j + 2] + h * x[j + 3]
        ## also h*x might be invalid ... e.g. Dates and ordered factors
        qs <- x[j+2L] # initialize quantiles as the "lower neighbour"
        # note that it is not sufficient to only assign the values for which 
        # h == 0 because we later combare x[j + 2L] != x[j + 3L] for the indices
        # 0 < h < 1. 
        # Doing in that non-obvious manner improves how infinities are dealt 
        # with if we don't do it like that we will e.g. encounter problems if x 
        # is complex and contains Inf as e.g. (Inf + 2i) * 0.5 is Inf+NaNi 
        qs[h == 1] <- x[j+3L][h == 1] # update values for which h == 1
        # now we update values for which h is in (0,1) and the lower and upper 
        # neighbour do not coincide (then no weighting is required anyway + 
        # some special cases for Infinities)
        # therefore we calculate a logical vector other that gives the indices
        # that still require an update 
        other <- (0 < h) & (h < 1) & (x[j+2L] != x[j+3L]) # '!=' for '<' in 
        # complex case
        # (there is no complete ordering on the complex numbers)

        # note that it is important to have the if(any(other)) because if we 
        # don't, it becomes a problem for ordered factors or characters as for 
        # those no + and * is defined (even for emtpy characters/factors)
        if(any(other)) qs[other] <- ((1-h)*x[j+2L] + h*x[j+3L])[other]
      }
    } else { # this is the case in which length(probs) == 0 or length(x) == 0
      qs <- rep(NA_real_, np) # define quantiles as a vector of real NAs 
      # where the length is the length of probs (with NAs already removed)
    }
    
    #--------------------------------------------------------------------------#
    
    #---------------------------now we deal with factors again-----------------#
    # this has to be done because the 'factor structure' of x was not passed
    # to qs (once x <- c(x[1], x[1], x, x[n], x[n])) x is no longer a factor
    # because it get's unclassed 
    # a factor 
    if(is.character(lx)) # in case x is not a factor lx is NULL and therefore
      # not a character 
      # a ordered factor is created so that the output corresponds to the
      # class of the input 
      # as qs is in 1, ..., length(lx) the levels attribute has to be set to 
      # seq_along(lx) and the labels are assigned according to the 
      # originallevels of x 
      qs <- factor(qs, levels = seq_along(lx), labels = lx, ordered = TRUE)
    
    #--------------------------------------------------------------------------#
    
    
    #------------------now we assign the names in case names == TRUE-----------#
    if(names && np > 0L) { # we format the probs into character vectors and
      # numbers are display using %, i.e. 1 corresponds to "100%" see
      # documentation of format_perc for details
      names(qs) <- format_perc(probs)
    }
    
    #--------------------------------------------------------------------------#
    
    
    #---------------deal with NAs in probs and return quantiles----------------#
    if(na.p) { # do this more elegantly (?!)
      # we assign the quantile values for those probs that are non NA
      o.pr[p.ok] <- qs
      # we assign the name "" to the NAs in probs
      names(o.pr) <- rep("", length(o.pr)) # suppress <NA> names
      # we assign the (formatted) names of qs to the corresponding values
      names(o.pr)[p.ok] <- names(qs)
      # return the named vector
      o.pr
    } else qs 
    #--------------------------------------------------------------------------#
  }

##' Formatting() percentages the same way as quantile(*, names=TRUE).
##' Should be exported
##' (and format.pval() moved to stats; both documented on same page)

# this function formats a numeric vector ínto a character vector that displays 
# the numbers as percentages 
# @param x - a numeric vector
# @param digits - how many decimal places
# @param probability - if TRUE 1 corresponds to 100%, otherwise to 1%
# @param use.fC - logical that indicates whether formatC or format is used
# by the way the function is written it seems that formatC is better for 
# small vectors but format scales better with the length of x 
# @param ... can be passed to format

format_perc <- function(x, digits = max(2L, getOption("digits")),
  probability = TRUE, use.fC = length(x) < 100, ...)
{
  if(length(x)) { # if length >= 1 this is is TRUE
    if(probability) x <- 100 * x # in case x is a probability vector
    # 1 will be displayed as "100%"
    
    # now x (the probability vector) is trasformed into a character vector and
    # x is display as percentages
    
    paste0(if(use.fC) ## formatC is slow for long x
      # format == "fg" means that digits indicates how many signficiant 
      # digits are to be printed
       formatC(x, format = "fg", width = 1, digits=digits)
      # formats an x for pretty printing
      # format(c(1,10), trim = FALSE) evaluates to c(" 1", "10")
      # format(c(1,10), trim = TRUE) evaluates to c("1", "10")
      else format(x, trim = TRUE, digits=digits, ...), "%")
    # if length(x) == 0 we return an empty character vector
  } else character(0)
}

# simply a wrapper that returns the difference between the 0.25 and 0.75 
# quartiles (IQR probably means inter-quartile range)
IQR <- function (x, na.rm = FALSE, type = 7)
  diff(quantile(as.numeric(x), c(0.25, 0.75), na.rm=na.rm, names = FALSE,
    type = type))