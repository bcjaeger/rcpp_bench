
library(microbenchmark)

Rcpp::sourceCpp("src/x_scale.cpp")

x <- matrix(rnorm(1e7), nrow = 1e5, ncol = 1e2)
w <- abs(rnorm(1e5))

ncols = 80
nrows = 50000

.cols <- sort(sample(ncol(x), ncols))
.rows <- sort(sample(nrow(x), nrows))

..cols <- as.integer(.cols-1)
..rows <- as.integer(.rows-1)

smat = x_scale_wtd(x, w, ..cols, ..rows)
iter = x_scale_wtd_iter(x, w, ..cols, ..rows)


bm = microbenchmark(
 smat = x_scale_wtd(x, w, ..cols, ..rows),
 iter = x_scale_wtd_iter(x, w, ..cols, ..rows),
 times = 100,
 unit = 'relative'
)

print(bm)
