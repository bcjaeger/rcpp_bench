
library(microbenchmark)
library(survival)
library(glmnet)

# Rcpp::sourceCpp("src/log_rank_test.cpp")


.pbc <- pbc[complete.cases(pbc), ]
.pbc <- .pbc[order(.pbc$time), ]
.pbc$status[.pbc$status > 0] <- .pbc$status[.pbc$status > 0] - 1

x <- as.matrix(.pbc[, -c(1,2,3,6), drop = FALSE])
y <- Surv(.pbc$time, .pbc$status)

# x <- orsf2::flchain_x
# x <- x[, - which(colnames(x)=='sexM')]
# y <- orsf2::flchain_y
# y[1:3, 1] <- y[1:3, 1]+1e-4

# random_samp <- sort(sample(nrow(x), nrow(x) * 0.95))
# x <- x[random_samp, ]
# y <- y[random_samp, ]

time <- y[, 1]
status <- y[, 2]

# set.seed(10)

group <- 1+rbinom(nrow(x), 1, 3/4)

rho <- 0

tt <- function(n = nrow(y), rho){

 ngroup <- 2

 obs <- rep(0, ngroup)
 exp <- rep(0, ngroup)
 var <- rep(0, ngroup*ngroup)
 risk <- rep(0, ngroup)
 kaplan <- rep(0, nrow(y))

 survdiff2(n, ngroup, rho, time, status,
           group, obs, exp, var, risk, kaplan)

 var <- matrix(var, ncol=2)

 df <- exp > 0

 temp2 <- ((obs - exp)[df])[-1]
 vv <- (var[df,df])[-1,-1, drop=FALSE]

 chi <- sum(solve(vv, temp2) * temp2)
 chi

}

yy <- as.matrix(y)
gg <- group-1

bm <- microbenchmark(
 tt = tt(n = nrow(y), rho = rho),
 # srv = survdiff(Surv(time, status) ~ group, rho = rho)$chisq,
 v0 = log_rank_test_v0(yy, gg),
 v1 = log_rank_test_v1(yy, group, rho = rho),
 v2 = log_rank_test_v2(yy, group, rho = rho),
 v3 = log_rank_test_v3(yy, gg),
 times = 10000,
 unit = 'relative'
)

print(bm)

