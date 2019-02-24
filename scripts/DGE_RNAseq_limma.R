#rm(list = ls())
# ------------------------------------------------------------------------------
# Loading source files and libraries
source(here::here("scripts","DGE_functions.R"))
#callr::rscript(here::here("scripts","DGE_RNAseq_limma.R"))

projectlibraries <-
  c(
    "PerformanceAnalytics",
    "edgeR",
    "limma",
    "stringr",
    "RColorBrewer",
    "tidyverse",
    "data.table",
    "rlang",
    "here",
    "callr"
  )

# call the function to check status, else, install or update or load packages
ProjectLibraries(projectlibraries)




# To remove packages
# lapply(projectlibraries, remove.packages, character.only = TRUE)
# lapply(projectlibraries, remove.packages)
# ------------------------------------------------------------------------------
# setwd
# setwd("~/R_Cruz_Nov18")

# ------------------------------------------------------------------------------
# Read the RDS file generated from the RNAseq pipleline
# pipeline output

# rawcounts <- read.table("cruz_Nov18_counts.txt",header=T,row.names=1)
# pipeline.RDS <-
#   readRDS(file = "Input_RNAseq_Shiny.rds")
# pipeline.RDS$rawcounts <- rawcounts

# pipeline.RDS <-
#   readRDS(file = here("raw_data/output_RNAseq_pipeline.rds"))
pipeline.RDS <-
  read_rds(path = here("raw_data/output_RNAseq_pipeline.rds"))

# pipeline.RDS <-
#   read_rds(path = here("raw_data/test_pipeline.rds"))

# laneinfotags <- c(
#   "_L001",
#   "_L002"
# )
# checkLanetags <- paste(laneinfotags, collapse = "|")
#
# gsubtagsRowCols <- function(x) {
#   gsub(checkLanetags, "", colnames(x))
#   #gsub(checkLanetags, "", rownames(x))
# }
#
# looplist <- lapply(pipeline.RDS, gsubtagsRowCols)
#
# # remove the added strings
# for (i in seq_along(looplist)) {
#   colnames(pipeline.RDS[[i]]) <- unlist(looplist[i])
# }
#
# head(pipeline.RDS)

pipeline.RDS <- lapply(pipeline.RDS, function(df) as_tibble(df))
pipeline.RDS$featurecounts_summary <- pipeline.RDS$featurecounts_summary
rawcounts <- pipeline.RDS$featurecounts

# ------------------------------------------------------------------------------
# variable inputs from project to project
SampleCategory <- read_delim(here("raw_data/sample.txt"),
                             delim = "\t",trim_ws = T)
# # trimspaces
# SampleCategory <-
#   SampleCategory %>%
#   rownames_to_column("sample_full_id") %>%
#   map_df(., ~
#   trimws(.))

GeneNames <-
  read_delim(here("raw_data/gene_names.txt"),delim = "\t",trim_ws = T)
# GeneNames <-
#   GeneNames %>%
#   rownames_to_column("ensemble_id") %>%
#   map_df(., ~ trimws(.))

# this has all the sampletypes, may be need more appropriate naming
samples.stage <- factor(SampleCategory$Base)


# ------------------------------------------------------------------------------
# DEG analysis codes
# The rawcounts should have rownames and columns as numeric
# if not, then do these steps
# Moreover, we need a matrix for calculations of DEG
rawcounts_num <- rawcounts[, -1]
rawcounts_num <- data.matrix(rawcounts_num)
rownames(rawcounts_num) <- rawcounts$sample_id

isexpr <-
  rowSums(cpm(rawcounts_num) > 2) >= dim(rawcounts_num)[2] / 2
x <- rawcounts_num[isexpr, ]

# Rawcounts columnames
rawcounts.colnames <- colnames(rawcounts_num)


# create design matrix
design <- model.matrix(~ 0 + samples.stage)
# colnames(design) <- c("CZ10FF", "CZ1E", "CZ5G", "CZwt")
colnames(design) <- levels(samples.stage)
y <- voom(x, design, plot = TRUE)

# merge genenames
voom <- y$E

# convert this into dataframe, as it is easier for merging later

voom <- data.frame(ensembl_id = rownames(voom), voom, row.names = NULL)

# Tibble is a better way to work with df. Gives more info to understand
voom <- voom %>% as_tibble() %>% mutate_if(is.factor, as.character)

# merge based on Ensemble id
voom_names <- merge(GeneNames, voom, by = "ensembl_id")
# rownames(voom_names) <- voom_names$Row.names
# voom_names$Row.names <- NULL

# Get combinations with paired comparisions
combinations <-
  combn(levels(samples.stage), 2, FUN = paste, collapse = "-")

# matrix is based on samples.stage, so no need to enter comparisions manually
contrast.matrix <-
  makeContrasts(
    contrasts = combinations,
    levels = design
  )
# contrast.matrix <-
#   makeContrasts(CZ10FF - CZwt,
#                 CZ1E - CZwt,
#                 CZ5G - CZwt,
#                 CZ10FF - CZ1E,
#                 CZ10FF - CZ5G,
#                 CZ1E - CZ5G,
#                 levels = design)

# Calculate the linear fit
fit <- lmFit(y, design)
fit2 <- eBayes(contrasts.fit(fit, contrast.matrix))

# This variable controls the way comparision column name should be
comparisions <-
  as.vector(combn(levels(samples.stage), 2, FUN = paste, collapse = "_vs_"))

datalist <- list()
sigdatalist <- list()
for (i in 1:length(comparisions)) {
  # The colnames we would like to append with sample comparision name
  colsTorename <- c(col1 = "logFC", col2 = "adj.P.Val")
  # for each coef, append the comparision and LogFc/p.val
  names(colsTorename) <- c(
    paste0(comparisions[i], "_logFC"),
    paste0(comparisions[i], "_adj.P.Val")
  )

  # Extract toptable from the linear model fit

  # -- without filtering
  df <- topTable(fit2, coef = index(comparisions)[i], n = Inf) %>%
    rownames_to_column("ensembl_id") %>%
    rename(!!colsTorename)
  # mapt them to voom_df that has rawcounts
  voom_df <- voom_names
  datalist[[i]] <-
    left_join(df, voom_df, by = "ensembl_id")

  # -- with filtering for significance
  sigdf <-
    topTable(fit2, coef = index(comparisions)[i], n = Inf) %>%
    rownames_to_column("ensembl_id") %>%
    filter(adj.P.Val <= 0.05) %>%
    filter(logFC >= 1 | logFC <= -1) %>%
    rename(!!colsTorename)

  # do we need to merge with voom_df ?
  sigdatalist[[i]] <- sigdf
}

names(datalist) <- comparisions
names(sigdatalist) <- comparisions
datalist <- lapply(datalist, function(df) as_tibble(df))
sigdatalist <- lapply(sigdatalist, function(df) as_tibble(df))

# Colnames to select in the final dataframe
LogFc.colnames <- paste0(comparisions, sep = "_logFC")
adj.P.Val.colnames <- paste0(comparisions, sep = "_adj.P.Val")
# hypens need to be converted to .
rawcounts.colnames <- gsub("-", ".", rawcounts.colnames)
colnames.toselect <- c(
  "ensembl_id",
  "gene_name",
  "AveExpr",
  rawcounts.colnames,
  LogFc.colnames,
  adj.P.Val.colnames
)

# merge the list of dataframes by ensemble ID column
complete_df <- reduce(datalist, full_join,
  by = "ensembl_id",
  suffix = c("_lstx", "_lsty")
)

colnames(complete_df) <- make.unique(gsub(
  "_lstx|_lsty", "",
  colnames(complete_df)
),
sep = "_"
)
final_df <- complete_df %>% dplyr::select(colnames.toselect)

# Push deg analysis data to the RDS file
pipeline.RDS$deg_analysis <- final_df
# Push DEG analysis into the pipeline RDS object
# save the file RDS file
# saveRDS(
#   pipeline.RDS,
#   file = "pipeline_with_DEG.RDS.rds",
#   ascii = FALSE,
#   version = NULL,
#   compress = TRUE,
#   refhook = NULL
# )
