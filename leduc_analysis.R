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

####---- Run msqrob LM analysis ----####

BiocParallel::register(BiocParallel::MulticoreParam(4)) ## 4 workers
sce <- msqrob(
    sce[1:10,], formula = ~ -1 + Set + Mock, ## only fixed effects
    ridge = FALSE, robust = FALSE
)

###----- Run msqrob with random effect ----####

sce <- msqrob(
    sce,
    formula = ~ -1 + Channel + (1|Set) + Mock,
    modelColumnName = "random_model",
    ridge = FALSE, robust = FALSE)

####---- Run scplainer analysis ----####

sce <- scpModelWorkflow(
    sce, formula = ~ 1 + Channel + Set + Mock # -1
)

####---- Compare number of fits ----####

## Number of successful fits by msqrob
msq_feat <- sum(sapply(rowData(sce)$msqrobModels, function(x) x@type != "fitError"))
## Number of succesful fits by random msqrob
ran_feat <- sum(sapply(rowData(sce)$random_model, function(x) x@type != "fitError"))


## Number of successful fits by scplainer
scpl_feat <- sum(scpModelFilterNPRatio(sce, filtered = FALSE) >= 1)

tot_feat <- nrow(sce)

msq_feat
scpl_feat
ran_feat
## These two numbers should be about the same...


# investigate

not_mod <- lapply(rowData(sce)$msqrobModels, function(x){
    any(is.na(getCoef(x)))
})
length(not_mod)-sum(unlist(not_mod)) # 2889 feat modelled

rowData(sce)$not_modelled <- unlist(not_mod)
rowData(sce)$scplainer_modelled <- scpModelFilterNPRatio(sce, filtered = FALSE) >= 1

modelled <- sce[!rowData(sce)$not_modelled]

non_modelled <- sce[rowData(sce)$not_modelled]

non_modelled_only_msqrob <- sce[rowData(sce)$not_modelled & rowData(sce)$scplainer_modelled]

non_modelled_scplainer <- sce[!rowData(sce)$scplainer_modelled]

y <- assay(non_modelled_scplainer[1,])

rownames(y)
scp:::scpModelFitList(sce, filtered = F)[rownames(y)][[1]]@coefficients
scpModelFilterNPRatio(sce, filtered = FALSE)[rownames(y)]
y <- as.data.frame(y)
y_long <- pivot_longer(y, cols = everything(), names_to = "sample", values_to = "Intensity")
col_data <- colData(sce)
col_data$sample <- rownames(col_data)
col_data <- as_tibble(col_data)
y_long <- left_join(y_long, col_data, by = "sample")

lm(Intensity ~ 1 + Channel + Set + Mock, data = y_long) # differences between lm and scplainer



y2 <- assay(non_modelled_only_msqrob[1,])

scp:::scpModelFitList(sce)[rownames(y2)][[1]]@coefficients # NP ratio > 1
rowData(sce)[rownames(y2),]$msqrobModels # fitError

y2 <- as.data.frame(y2)
y2_long <- pivot_longer(y2, cols = everything(), names_to = "sample", values_to = "Intensity")

y2_long <- left_join(y2_long, col_data, by = "sample")

lm(Intensity ~ 1 + Channel + Set + Mock, data = y2_long) # lm reussi
# differences between lm and scplainer

sub_assay <- sce[rownames(y2),]


sce <- msqrob(
    sub_assay,
    formula = ~ -1 + Channel + Set + Mock,
    ridge = FALSE, robust = FALSE,
    overwrite = TRUE)
