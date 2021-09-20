
library(microbenchmark)

Rcpp::sourceCpp("src/pass_by_ref.cpp")

x <- matrix(rnorm(1e5), nrow = 1e3, ncol = 1e2)

microbenchmark(
 gb = global_bind(x),
 av = advanced_constructor(x)
)
