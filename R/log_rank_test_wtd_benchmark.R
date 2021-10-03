

library(microbenchmark)
library(survival)
library(glmnet)

Rcpp::sourceCpp("src/log_rank_test_wtd.cpp")

.pbc <- pbc[complete.cases(pbc), ]
.pbc <- .pbc[order(.pbc$time), ]
.pbc$status[.pbc$status > 0] <- .pbc$status[.pbc$status > 0] - 1

x <- as.matrix(.pbc[, -c(1,2,3,6), drop = FALSE])
y <- Surv(.pbc$time, .pbc$status)

wts <- rep(1, nrow(x))

group <- 1+rbinom(nrow(x), 1, 2/4)
g <- group - 1

log_rank_test_wtd(y, wts, g)
survdiff(y ~ group)$chisq

wts <- sample(1:5, size = nrow(.pbc), replace = TRUE)

log_rank_test_wtd(y, wts, g)

.pbc_wtd <- .pbc[ rep(seq(nrow(.pbc)), times = wts), ]
y_wtd <- Surv(.pbc_wtd$time, .pbc_wtd$status)

group_wtd <- rep(group, times = wts)

survdiff(y_wtd ~ group_wtd)$chisq
