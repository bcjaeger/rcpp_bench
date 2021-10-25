
library(microbenchmark)
Rcpp::sourceCpp("src/submat_mult.cpp")

ncol = 100
nrow = 10000

ncol_sub = 4
nrow_sub = 1000

x <- matrix(rnorm(ncol*nrow), ncol=ncol, nrow=nrow)

head(x)

r_rows <- sample(nrow, nrow_sub)
c_rows <- r_rows - 1

r_cols <- sample(ncol, ncol_sub)
c_cols <- r_cols - 1

beta <- rnorm(n = ncol_sub)

bm <- microbenchmark(
 r = x[r_rows, r_cols] %*% beta,
 iter_2loop = iter_2loop(x, c_rows, c_cols, beta),
 submat = submat(x, c_rows, c_cols, beta),
 submat2 = submat2(x, c_rows, c_cols, beta),
 iter_bigout = iter_bigout(x, c_rows, c_cols, beta),
 iter_bigbeta = iter_bigbeta(x, c_rows, c_cols, beta),
 unit = 'relative',
 times = 500
)

print(bm)


