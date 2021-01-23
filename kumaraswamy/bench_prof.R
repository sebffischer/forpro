library(bench)
library(profvis)

source("kumaraswamy-rng.R")

set.seed(1804)

n = 1e5

#----input checking------------------------------------------------------------#

a_scalar <- 0.5
b_scalar = 0.3
min_scalar = -1
max_scalar = 2

bench_input_checking_scalar <- bench::mark(
  check_rkumaraswamy(n, a_scalar, b_scalar, min_scalar, max_scalar)
)

  
a_vector <- rep_len(rexp(1e2), n)
b_vector <- rep_len(rexp(1e2), n)
min_vector <- rep_len(rnorm(1e2), n)
max_vector <- rep_len(rnorm(1e2), n)

bench_input_checking_vector <- bench::mark(
  check_rkumaraswamy(n, a_vector, b_vector, min_vector, max_vector)
)



#----actual calculations-------------------------------------------------------#
  
bench_rkum_scalar <- bench::mark(
  rkumaraswamy1_main(n, a_scalar, b_scalar, min_scalar, max_scalar),
  rkumaraswamy2_main(n, a_scalar, b_scalar, min_scalar, max_scalar), 
  check = FALSE
)

#names(bench_rkum_scalar)[1] <- "function"
#bench_rkum_scalar[,1] <- c("rkumaraswamy1_main", "rkumaraswamy2_main")
#bench_rkum_scalar[, 2:9] <- round(bench_rkum_scalar[, 2:9], 2)



bench_rkum_vector <- bench::mark(
  rkumaraswamy1_main(n, a_vector, b_vector, min_vector, max_vector),
  rkumaraswamy2_main(n, a_vector, b_vector, min_vector, max_vector), 
  check = FALSE
)

names(bench_rkum_vector)[1] <- "function"
bench_rkum_vector[,1] <- c("rkumaraswamy1_main", "rkumaraswamy2_main")

profvis::profvis(
  {
  rkumaraswamy1_main(n, a_scalar, b_scalar, min_scalar, max_scalar)
  rkumaraswamy2_main(n, a_scalar, b_scalar, min_scalar, max_scalar)}
)

### save results 

# save(list = c("bench_input_checking_scalar", 
              "bench_input_checking_vector", 
              "bench_rkum_scalar", 
              "bench_rkum_vector"), 
  file = "benchmark_results.RData")