

library(microbenchmark)
library(survival)
library(glmnet)

Rcpp::sourceCpp("src/cox_fit_submat.cpp")


.pbc <- pbc[complete.cases(pbc), ]
.pbc <- .pbc[order(.pbc$time), ]

.pbc$status[.pbc$status > 0] <- .pbc$status[.pbc$status > 0] - 1

x <- as.matrix(.pbc[, -c(1,2,3,6), drop = FALSE])
y <- Surv(.pbc$time, .pbc$status)

x <- orsf2::flchain_x
x <- x[, - which(colnames(x)=='sexM')]
y <- orsf2::flchain_y
y[1:3, 1] <- y[1:3, 1]+1e-4



# x <- x[,2, drop=F]

set.seed(1)

weights <- sample(1:5, size = nrow(x), replace = TRUE)

iter_max = 2

coxph_args <- list(x = x,
                   y = y,
                   strata = NULL,
                   init = rep(0, ncol(x)),
                   control = coxph.control(iter.max = iter_max),
                   offset = rep(0, nrow(x)),
                   weights = weights,
                   method = 'breslow',
                   rownames = NULL,
                   resid = FALSE,
                   nocenter = 0)

ncols = 5
nrows = 2000

.cols <- sort(sample(ncol(x), ncols))
.rows <- sort(sample(nrow(x), nrows))

..cols <- as.integer(.cols-1)
..rows <- as.integer(.rows-1)

orsf <- newtraph_cph_submat(x[],
                            y,
                            weights,
                            ..cols,
                            ..rows,
                            method = 0,
                            eps = 1e-9,
                            iter_max = iter_max,
                            rescale = TRUE)

orsf2 <- newtraph_cph_iit(x[],
                          y,
                          weights,
                          ..cols,
                          ..rows,
                          method = 0,
                          eps = 1e-9,
                          iter_max = iter_max,
                          rescale = TRUE)



surv <- coxph.fit(x = x[.rows, .cols],
                  y = y[.rows, ],
                  strata = NULL,
                  init = rep(0, ncols),
                  control = coxph.control(iter.max = iter_max),
                  offset = rep(0, nrows),
                  weights = weights[.rows],
                  method = 'breslow',
                  rownames = NULL,
                  resid = FALSE,
                  nocenter = 0)

# verify that orsf and surv get the same coefficients
print(max(abs(orsf[, 1] - surv$coefficients)))
# verify that orsf and surv get the same standard error
print(max(abs(orsf[, 2] - sqrt(diag(surv$var)))))


bmark <- microbenchmark(

        orsf = newtraph_cph_submat(x[,],
                                   y,
                                   weights,
                                   ..cols,
                                   ..rows,
                                   method = 0,
                                   eps = 1e-9,
                                   iter_max = iter_max,
                                   rescale = TRUE),

        orsf2 = newtraph_cph_iit(x[,],
                                 y,
                                 weights,
                                 ..cols,
                                 ..rows,
                                 method = 0,
                                 eps = 1e-9,
                                 iter_max = iter_max,
                                 rescale = TRUE),

        surv = coxph.fit(x = x[.rows, .cols],
                         y = y[.rows, ],
                         strata = NULL,
                         init = rep(0, ncols),
                         control = coxph.control(iter.max = iter_max),
                         offset = rep(0, nrows),
                         weights = weights[.rows],
                         method = 'breslow',
                         rownames = NULL,
                         resid = FALSE,
                         nocenter = 0),

        unit = 'relative',

        times = 1000

)

print(bmark)

