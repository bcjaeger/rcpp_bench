

source("src/join_versus_prealloc.cpp")

library(microbenchmark)

n_row = 1000
n_col = 10

microbenchmark(
 test0 = arma_resize(n_row = n_row, n_col = n_col),
 test1 = arma_join(n_row = n_row, n_col = n_col),
 test2 = arma_prealloc(n_row = n_row * 10, n_col = n_col, max_row = n_row),
 test3 = arma_join_list(n_row = n_row, n_col = n_col),
 unit = 'relative'
)

# Maximum number of nodes in a binary tree of height h is 2^h – 1.
# if min_obs is 10 and we have 160 obs -> 16 leaves -> height is 4




