library(microbenchmark)
library(survival)
library(glmnet)

Rcpp::sourceCpp("src/subvec_mult.cpp")

n <- 1e6
u <- rnorm(n)
v <- rnorm(n)

microbenchmark::microbenchmark(
 v1 = subvec_mult_v1(u, v, n),
 v2 = subvec_mult_v2(u, v, n),
 v3 = subvec_mult_v3(u, v, n),
 v4 = subvec_mult_v4(u, v, n),
 v5 = subvec_mult_v4(u, v, n)
)



