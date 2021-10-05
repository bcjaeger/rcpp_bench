
library(microbenchmark)
library(survival)
library(glmnet)

Rcpp::sourceCpp("src/log_rank_test_multi.cpp")


.pbc <- pbc[complete.cases(pbc), ]
.pbc <- .pbc[order(.pbc$time), ]
.pbc$status[.pbc$status > 0] <- .pbc$status[.pbc$status > 0] - 1

x <- as.matrix(.pbc[, -c(1,2,3,6), drop = FALSE])
y <- Surv(.pbc$time, .pbc$status)

x <- orsf2::flchain_x
x <- x[, - which(colnames(x)=='sexM')]
y <- orsf2::flchain_y
y[1:3, 1] <- y[1:3, 1]+1e-4

mdl <- coxph(y~x)

summary(mdl)

XB <- predict(mdl)

out = lrt_multi_v1(y, XB, 10, 2)

library(tidyverse)
as_tibble(out) |>
 set_names('x', 'y') |>
 ggplot(aes(x=x, y=y)) +
 geom_point()
