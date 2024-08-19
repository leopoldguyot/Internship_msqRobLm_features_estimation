#### Packages

library(scp)
library(msqrob2)
library(limma)
library(microbenchmark)
library(dplyr)

#### Functions

msqrobLm_fix <- function (y, formula, data, robust = TRUE, maxitRob = 5)
{
    myDesign <- model.matrix(formula, data)
    models <- apply(y, 1, function(y, design) {
        obs <- is.finite(y)
        type <- "fitError"
        model <- list(coefficients = NA, vcovUnscaled = NA,
                      sigma = NA, df.residual = NA, w = NA)
        if (sum(obs) > 0) {
            X <- design[obs, , drop = FALSE]
            qrX <- qr(X)
            X <- X[, qrX$pivot[seq_len(qrX$rank)], drop = FALSE] # change
            # => old version
            # X <- X[, colMeans(X == 0) != 1, drop = FALSE]
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
                model <- list(coefficients = mod$coef, vcovUnscaled = msqrob2:::.vcovUnscaled(mod),
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

msqrobLm_wrapper <- function(sce){
    msqrob2::msqrobLm(assay(sce),
             ~ 1 + Channel + Set + Mock,
             colData(sce),
             robust = FALSE)
}

msqrobLm_fix_wrapper <- function(sce){
    msqrobLm_fix(assay(sce),
             ~ 1 + Channel + Set + Mock,
             colData(sce),
             robust = FALSE)
}

limma_wrapper <- function(sce){
    y <- assay(sce)
    design <- model.matrix(~ 1 + Channel + Set + Mock, colData(sce))
    limma::lmFit(y, design)
}

scplainer_wrapper <- function(sce){
    scp::scpModelWorkflow(sce,
                          formula = ~ 1 + Channel + Set + Mock)
}


### Data loading

sce <- readRDS("data/leduc_mock.rds")
set.seed(1234)
sce$Mock <- .createMockBalanced(2, sce$Set)

sce1000 <- sce[sample(1:nrow(sce), 1000),]
colData(sce1000) <- droplevels(colData(sce1000))
sce2000 <- sce[sample(1:nrow(sce), 3000),]
colData(sce2000) <- droplevels(colData(sce2000))
sce5000 <- sce[sample(1:nrow(sce), 5000),]
colData(sce5000) <- droplevels(colData(sce5000))
sce10000 <- sce[sample(1:nrow(sce), 10000),]
colData(sce10000) <- droplevels(colData(sce10000))
sce18000 <- sce[sample(1:nrow(sce), 18000),]
colData(sce18000) <- droplevels(colData(sce18000))

### Benchmarking
time <- Sys.time()

bench_res <- microbenchmark(msqrobLm_wrapper(sce2000),
                            msqrobLm_fix_wrapper(sce2000),
                            limma_wrapper(sce2000),
                            scplainer_wrapper(sce2000),
                            msqrobLm_wrapper(sce5000),
                            msqrobLm_fix_wrapper(sce5000),
                            limma_wrapper(sce5000),
                            scplainer_wrapper(sce5000),
                            msqrobLm_wrapper(sce10000),
                            msqrobLm_fix_wrapper(sce10000),
                            limma_wrapper(sce10000),
                            scplainer_wrapper(sce10000),
                            msqrobLm_wrapper(sce18000),
                            msqrobLm_fix_wrapper(sce18000),
                            limma_wrapper(sce18000),
                            scplainer_wrapper(sce18000),
                            times = 3)
bench_res <- bench_res %>%
    as.data.frame() %>%
    mutate(Method = sub("_wrapper.*", "", expr),
           featureNumber = as.numeric(sub(".*sce([0-9]+).*", "\\1", expr)))

write.csv(bench_res, file = "data_output/benchmark_results.csv")


# for only the estimable peptides
colData(sce) <- droplevels(colData(sce))
test <- msqrobLm_wrapper(sce)
modelled <- unlist(lapply(test, function(x) x@type == "lm"))

subset_modelled <- sce[modelled,]
colData(subset_modelled) <- droplevels(colData(subset_modelled))

bench_res_modelled <- microbenchmark(msqrobLm_wrapper(subset_modelled),
                                     msqrobLm_fix_wrapper(subset_modelled),
                                     limma_wrapper(subset_modelled),
                                     scplainer_wrapper(subset_modelled),
                                     times = 5)
bench_res_modelled <- bench_res_modelled %>%
    as.data.frame() %>%
    mutate(Method = sub("_wrapper.*", "", expr),
           featureNumber = as.numeric(sub(".*sce([0-9]+).*", "\\1", expr)))

write.csv(bench_res_modelled, file = "data_output/benchmark_results_modelled.csv")

print(time - Sys.time())
