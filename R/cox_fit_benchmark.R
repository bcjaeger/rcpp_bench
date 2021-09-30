

library(microbenchmark)
library(survival)
library(glmnet)

Rcpp::sourceCpp("src/cox_fit.cpp")


.pbc <- pbc[complete.cases(pbc), ]
.pbc <- .pbc[order(.pbc$time), ]

.pbc$status[.pbc$status > 0] <- .pbc$status[.pbc$status > 0] - 1

x <- as.matrix(.pbc[, -c(1,2,3,6), drop = FALSE])
y <- Surv(.pbc$time, .pbc$status)

# x <- orsf2::flchain_x
# x <- x[, - which(colnames(x)=='sexM')]
# y <- orsf2::flchain_y
# y[1:3, 1] <- y[1:3, 1]+1e-4



# x <- x[,2, drop=F]

set.seed(1)

weights <- sample(1:5, size = nrow(x), replace = TRUE)

iter_max = 10

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

glmnet_args <- list(x = x,
                    y = y,
                    family = 'cox',
                    weights = weights)


orsf <- newtraph_cph(x[, , drop = FALSE],
                     y,
                     weights,
                     method = 0,
                     eps = 1e-5,
                     iter_max = iter_max,
                     rescale = TRUE)

surv <- do.call(coxph.fit, coxph_args)

# verify that orsf and surv get the same coefficients
print(max(abs(orsf[, 1] - surv$coefficients)))
# verify that orsf and surv get the same standard error
print(max(abs(orsf[, 2] - sqrt(diag(surv$var)))))


bmark <- microbenchmark(

    orsf = newtraph_cph(x = x[, , drop = FALSE],
                        y = y[, ],
                        weights,
                        method = 0,
                        eps = 1e-5,
                        iter_max = iter_max,
                        rescale = TRUE),

    surv = suppressWarnings(do.call(coxph.fit, coxph_args)),

    glmnet = do.call(glmnet, glmnet_args),

    glmnet_cv = do.call(cv.glmnet, glmnet_args),

    unit = 'relative',

    times = 100

)

print(bmark)

