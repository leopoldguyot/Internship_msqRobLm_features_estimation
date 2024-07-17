####---- Load packages and data ----####

library("scp")
library("msqrob2")
library("BiocParallel")
library("tidyverse")

####---- Utility functions ----####

.createMockBalanced <- function(k, runs) {
    out <- split(runs, runs)
    out <- lapply(out, function(x) {
        index <- 1:length(x)
        x[sample(index)] <- rep(1:k, length.out = length(index))
        x
    })
    out <- S4Vectors::unname(do.call(c, out))
    as.factor(out)
}

####---- Load data ----####

sce <- readRDS("data/leduc_mock.rds")
set.seed(1234)
sce$Mock <- .createMockBalanced(2, sce$Set)

sub_assay <- sce["EATTEFSVDAR_2",] # non fitted with msqrob

sub_assay <- msqrob(
    sub_assay,
    formula = ~ -1 + Channel + Set + Mock,
    ridge = FALSE, robust = FALSE,
    overwrite = TRUE)
debug(msqrobLm)
msqrobLm(assay(sub_assay), ~ -1 + Channel + Set + Mock, colData(sub_assay))


custom_msqlm <- function (y, formula, data, robust = TRUE, maxitRob = 5)
{
    myDesign <- model.matrix(formula, data)
    models <- apply(y, 1, function(y, design) {
        obs <- is.finite(y)
        type <- "fitError"
        model <- list(coefficients = NA, vcovUnscaled = NA,
                      sigma = NA, df.residual = NA, w = NA)
        if (sum(obs) > 0) {
            X <- design[obs, , drop = FALSE]
            X <- X[, colMeans(X == 0) != 1, drop = FALSE]
            y <- y[obs]
            colnames_orig <- colnames(design)
            if (robust) {
                mod <- try(MASS::rlm(X, y, method = "M", maxit = maxitRob),
                           silent = TRUE)
                if (!is(mod, "try-error")) {
                    type <- "rlm"
                }
            }
            else {
                mod <- try(lm.fit(X, y))
                if ((!is(mod, "try-error")) & mod$rank == ncol(X)) {
                    type <- "lm"
                }
            }
            if (type == "rlm") {
                w <- mod$w
                sigma <- sqrt(sum(mod$w * mod$resid^2)/(sum(mod$w) -
                                                            mod$rank))
                df.residual <- sum(mod$w) - mod$rank
                if (df.residual < 2L)
                    type <- "fitError"
            }
            if (type == "lm") {
                w <- NULL
                sigma <- sqrt(sum(mod$residuals^2/mod$df.residual))
                df.residual <- mod$df.residual
                if (df.residual < 2L)
                    type <- "fitError"
            }
            if (type != "fitError") {
                coef <- rep(NA, length(colnames_orig))
                names(coef) <- colnames_orig
                coef[names(mod$coef)] <- mod$coef
                vcovUnscaled <- matrix(NA, nrow = length(colnames_orig),
                                       ncol = length(colnames_orig))
                rownames(vcovUnscaled) <- colnames(vcovUnscaled) <- colnames_orig
                vcovUnscaled[names(mod$coef), names(mod$coef)] <- msqrob2:::.vcovUnscaled(mod)
                model <- list(coefficients = mod$coef, vcovUnscaled = .vcovUnscaled(mod),
                              sigma = sigma, df.residual = df.residual,
                              w = w)
            }
        }
        .StatModel(type = type, params = model, varPosterior = as.numeric(NA),
                   dfPosterior = as.numeric(NA))
    }, design = myDesign)
    hlp <- limma::squeezeVar(var = vapply(models, getVar, numeric(1)),
                             df = vapply(models, getDF, numeric(1)))
    for (i in seq_len(length(models))) {
        mydf <- hlp$df.prior + getDF(models[[i]])
        models[[i]]@varPosterior <- as.numeric(hlp$var.post[i])
        models[[i]]@dfPosterior <- as.numeric(mydf)
    }
    return(models)
}
