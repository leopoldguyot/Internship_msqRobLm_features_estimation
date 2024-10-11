library("scp")
library("msqrob2")
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

profvis::profvis({
    msqrob(
        sce[1:4000,], formula = ~ -1 + Set + Mock,
        ridge = FALSE, robust = FALSE
    )
})
