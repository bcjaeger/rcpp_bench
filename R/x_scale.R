
library(microbenchmark)

Rcpp::sourceCpp("src/x_scale.cpp")

x <- matrix(rnorm(1e5), nrow = 1e3, ncol = 1e2)
w <- rnorm(1e3)

bm = microbenchmark(
 fun_1 = x_scale_wtd_1(x[, ], w),
 fun_2 = x_scale_wtd_2(x[, ], w),
 fun_R = scale(x), # just for reference
 times = 100
)

print(bm)
