# ----1. load packages and data----
library(dplyr)
library(janitor)
library(r.jive)

#setup
paired     <- readRDS("data/clean/paired_hmo_metagenome.rds")
microbiome <- readRDS("data/clean/microbiome_clean.rds")  

hmo_cols <- setdiff(names(readRDS("data/clean/hmo_percent_unique.rds")),
                    c("time", "name", "enrollment_id", "secretor",
                      "diversity", "evenness", "days", "days_from_61", "sum"))

# ---- 2. link samples----
paired <- paired |> mutate (sample_clean = make_clean_names(sample_name))
sum(paired$sample_clean %in%  names(microbiome))

# ----3. build the matrices ----
tax_cols <- c("domain","phylum","class","order","family","genus","species","ogu")

micro_mat <- as.matrix(microbiome[, paired$sample_clean])
rownames(micro_mat) <- microbiome$ogu
colnames(micro_mat) <- paired$enrollment_id

hmo_mat <- t(as.matrix(paired[, hmo_cols]))     # t() flips to HMOs x samples
colnames(hmo_mat) <- paired$enrollment_id

identical(colnames(micro_mat), colnames(hmo_mat))   # must be TRUE
dim(hmo_mat); dim(micro_mat)

# ---- 4. filter and transform ----
keep <- rowMeans(micro_mat > 0) >= 0.10
micro_mat <- micro_mat[keep, ]
clr <- function(m) {
  pseudo <- min(m[m > 0]) / 2   # small number added so log() never sees a zero
  l <- log(m + pseudo)
  sweep(l, 2, colMeans(l))      # subtract each sample's average log value
}
hmo_clr   <- clr(hmo_mat)
micro_clr <- clr(micro_mat)

#---- 5. run JIVE ----
#set.seed(123)
#jive_fit <- jive(list(HMO = hmo_clr, Microbiome = micro_clr), method = "perm")

# ---- 6. results ----
class(jive_fit)
jive_fit$rankJ
jive_fit$rankA
summary(jive_fit)
showVarExplained(jive_fit)

saveRDS(jive_fit, "data/clean/jive_fit.rds")
