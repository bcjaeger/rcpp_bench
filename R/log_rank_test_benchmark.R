
library(microbenchmark)
library(survival)
library(glmnet)

Rcpp::sourceCpp("src/log_rank_test.cpp")

.pbc <- pbc[complete.cases(pbc), ]
.pbc <- .pbc[order(.pbc$time), ]
.pbc$status[.pbc$status > 0] <- .pbc$status[.pbc$status > 0] - 1

x <- as.matrix(.pbc[, -c(1,2,3,6), drop = FALSE])
y <- Surv(.pbc$time, .pbc$status)

# x <- orsf2::flchain_x
# x <- x[, - which(colnames(x)=='sexM')]
# y <- orsf2::flchain_y
# y[1:3, 1] <- y[1:3, 1]+1e-4

# vet <- survival::veteran
#
# x <- as.matrix(vet[,c(1,5,6,7,8)])
# y <- as.matrix(vet[,c('time', 'status')])
# y <- Surv(time = vet$time, event = vet$status)
# x <- x[order(y[,'time']), ]
# y <- y[order(y[,'time']), ]

# random_samp <- sort(sample(nrow(x), nrow(x) * 0.95))
# x <- x[random_samp, ]
# y <- y[random_samp, ]

time <- y[, 1]
status <- y[, 2]

# set.seed(10)

group <- 1+rbinom(nrow(x), 1, 2/4)

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
        times = 1000
)

print(bm)


# mdl <- coxph(y~x)
#
# n <- nrow(x)
#
# lc_vals <- predict(mdl)
#
# lc_sort <- order(lc_vals)
# n_events <- n_obs <- 0
#
# for(i in seq(1, n)){
#  n_events = n_events + y[lc_sort[i], 'status']
#  n_obs = n_obs + 1
#  cat("n_events:", n_events, " ")
#  cat("n_obs:", n_obs, "\n")
#  if(n_events >= 5 && n_obs >= 10) break
# }
#
# ii <- i
#
# lc_index_lwr <- lc_sort[i]
#
# n_events <- n_obs <- 0
# gg <- rep(0, nrow(x))
#
# for(i in seq(n, 1)){
#  n_events = n_events + y[lc_sort[i], 'status']
#  n_obs = n_obs + 1
#  gg[lc_sort[i]] = 1
#  cat("n_events:", n_events, " ")
#  cat("n_obs:", n_obs, "\n")
#  if(n_events >= 5 && n_obs >= 10) break
# }
#
# total_steps <- i - ii
#
# v4 = log_rank_test_v4(yy, gg)
#
# index  = v4[,1]
# g_risk = v4[,2]
# n_risk = v4[,3]
# deaths = v4[,4]
#
# temp2 = g_risk / n_risk
# temp1 = deaths * temp2 * (n_risk-deaths) / (n_risk-1)
#
# expected = deaths * temp2
#
# V = temp1 * (1 - temp2)
#
# o = sum(status * gg)
# e = sum(expected)
# v = sum(V, na.rm = T)
#
# (o-e)^2 / v
# survdiff(Surv(time, status) ~ gg, rho = rho)$chisq
#
# i = i-1
#
# n_cps = 10
# i_step <- round(total_steps / (n_cps))
#
# library(tibble)
# ggdat <- tibble()
#
# while(i > ii){
#
#  gg[lc_sort[i]] <- 1
#
#  for(j in seq(nrow(v4), 1)) {
#   if (index[j] + 1 == lc_sort[i]) {
#    break
#   } else if (index[j] + 1 > lc_sort[i]) {
#    j = j+1
#    break
#   }
#  }
#
#  k <- seq(j, nrow(v4))
#
#  # g_risk[k] = g_risk[k] + 1
#
#  o = o + status[lc_sort[i]]
#  e = e + sum(deaths[k]/n_risk[k])
#
#  temp2[k] = temp2[k] + 1 / n_risk[k]
#  temp1[k] = temp1[k] + deaths[k] * (n_risk[k]-deaths[k]) / (n_risk[k] * (n_risk[k]-1))
#
#  V[k] = temp1[k] * (1 - temp2[k])
#
#  if((i-lc_index_lwr) %% i_step == 0){
#
#   v = sum(V, na.rm = T)
#
#   bcj = (o-e)^2 / v
#   tt = survdiff(Surv(time, status) ~ gg, rho = rho)$chisq
#
#   cat("i: ", i, "; ",
#       "j: ", j,  "; ",
#       "diff: ", round(bcj - tt, 3),
#       '\n',
#       sep = '')
#
#
#   ggdat <- rbind(ggdat, tibble(i = i, bcj = bcj, tt = tt))
#
#  }
#
#
#  i <- i - 1
#
# }
#
# library(ggplot2)
#
# ggplot(ggdat) +
#  aes(x = i, y = bcj) +
#  geom_line()
#







