source(here::here("scripts","DGE_RNAseq_limma.R"))
#callr::rscript(here::here("scripts","DGE_RNAseq_limma.R"))


# rename columns by adding prefixes
# These prefixes customize which columns should be selected for plots
# Total_reads = Total_pairs * 2
if (!is_empty(pipeline.RDS$insert_size_metrics)) {
  pipeline.RDS$hisat2_metrics$Total_pairs <- pipeline.RDS$hisat2_metrics$Total_pairs * 2
}


colnames(pipeline.RDS$hisat2_metrics) <-
  c(
    "xaxis_sample_id",
    "yaxis_UniquelyMapped",
    "yaxis_MultiMapped",
    "yaxis_UnMapped",
    "norm_Total_reads",
    "Overall_alignment_rate"
  )

colnames(pipeline.RDS$Infer_metrics) <-
  c(
    "xaxis_sample_id",
    "yaxis_Percent_reads_failed_to_determine",
    "yaxis_Percentage_of_reads_forward_mapped",
    "yaxis_Percentage_of_reads_reverse_mapped"
  )

# coloumn names changed with a prefix, so it is easy to grep for plot
# for single stranded we do not have insert_size_metrics table
if (!is_empty(pipeline.RDS$insert_size_metrics)) {
  colnames(pipeline.RDS$insert_size_metrics) <-
    c(
      "xaxis_sample_id",
      "MEDIAN_INSERT_SIZE",
      "MEDIAN_ABSOLUTE_DEVIATION",
      "yaxis_MEAN_INSERT_SIZE",
      "STANDARD_DEVIATION",
      "READ_PAIRS",
      "PAIR_ORIENTATION"
    )

}

# no coloumns need to be changed
colnames(pipeline.RDS$junctions) <-
  c(
    "xaxis_sample_id",
    "yaxis_Known_Splicing_Junctions",
    "yaxis_Novel_Splicing_Junctions",
    "yaxis_Partial_Novel_Splicing_Junctions",
    "norm_Total_Splicing_Junctions"
  )

# Total reads in Ribosomal contamination may vary from actual total reads
# Hence, sum the columns rowwise
pipeline.RDS$ribo <- pipeline.RDS$ribo %>%
  mutate(norm_Total_Reads = select(., Ribosomal, Not_ribosomal, Unmapped) %>%
    rowSums(na.rm = TRUE))
colnames(pipeline.RDS$ribo) <-
  c(
    "xaxis_sample_id",
    "yaxis_Ribosomal",
    "yaxis_Non-ribosomal",
    "yaxis_Unmapped",
    "Percent_contamination",
    "norm_Total_Reads"
  )

# Total_reads is part of read-distribution
colnames(pipeline.RDS$read_distribution) <-
  c(
    "xaxis_sample_id",
    "Total_Reads",
    "norm_Total_Tags",
    "Total_Assigned_Tags",
    "CDS_Exons",
    "5'UTR_Exons",
    "3'UTR_Exons",
    "Introns",
    "TSS_up_10kb",
    "TES_down_10kb",
    "Total_Unassigned_Tags",
    "yaxis_Exon_percent",
    "yaxis_Intron_percent",
    "yaxis_Intergenic_percent",
    "yaxis_Exon_read_perKb",
    "yaxis_Intron_read_perKb",
    "yaxis_up/downstream_reads_perKb"
  )

# Create the ID Column
# pipeline.RDS$featurecounts_summary <-
#  pipeline.RDS$featurecounts_summary %>%
#  rownames_to_column("sample_id")

# Add the total_reads column

pipeline.RDS$featurecounts_summary$norm_Total_Reads <- pipeline.RDS$hisat2_metrics$norm_Total_reads

colnames(pipeline.RDS$featurecounts_summary) <-
  c(
    "xaxis_sample_id",
    "yaxis_Assigned",
    "yaxis_Unassigned_Ambiguity",
    "yaxis_Unassigned_NoFeatures",
    "Unassigned_Secondary",
    "Unassigned_Unmapped",
    "norm_Total_Reads"
  )

# Change sample_id in featurecounts with ensemble ID
colnames(pipeline.RDS$featurecounts)[1] <- c("xaxis_ensembl_id")

# Volcano Plot / MA plot
colnames(pipeline.RDS$deg_analysis) <-
  gsub(
    "ensembl_id",
    "xaxis_ensembl_id",
    colnames(pipeline.RDS$deg_analysis)
  )
colnames(pipeline.RDS$deg_analysis) <-
  gsub(
    "gene_name",
    "xaxis_gene_name",
    colnames(pipeline.RDS$deg_analysis)
  )

prefixLogFC <-
  colnames(pipeline.RDS$deg_analysis)[grepl("_logFC", colnames(pipeline.RDS$deg_analysis))]
logFcindex <- match(prefixLogFC, names(pipeline.RDS$deg_analysis))
colnames(pipeline.RDS$deg_analysis)[logFcindex] <-
  paste("yaxis_", colnames(pipeline.RDS$deg_analysis[, logFcindex]), sep = "")

prefixP.Val <-
  colnames(pipeline.RDS$deg_analysis)[grepl("_adj.P.Val", colnames(pipeline.RDS$deg_analysis))]
P.Valindex <- match(prefixP.Val, names(pipeline.RDS$deg_analysis))
colnames(pipeline.RDS$deg_analysis)[P.Valindex] <-
  paste("yaxis_", colnames(pipeline.RDS$deg_analysis[, P.Valindex]), sep = "")



# Heatmap
matcindexs <-
  colnames(pipeline.RDS$deg_analysis) %in% rawcounts.colnames
colnames(pipeline.RDS$deg_analysis)[matcindexs] <-
  paste("cpm_", rawcounts.colnames, sep = "")

# create the input_parameters dataframe

adddedstrings <- c("xaxis_", "yaxis_", "norm_", "cpm_")
checkvars <- paste(adddedstrings, collapse = "|")


# Recreate the cleaner input for the Rshiny input

rshiny_input <- list()
rshiny_input$hisat_alignment <- pipeline.RDS$hisat2_metrics
rshiny_input$duplication_metrics <- pipeline.RDS$duplication_metrics
rshiny_input$featurecounts_summary <-
  pipeline.RDS$featurecounts_summary
rshiny_input$featurecounts <- pipeline.RDS$featurecounts
rshiny_input$ribosomal_contamination <- pipeline.RDS$ribo
rshiny_input$infer_metrics <- pipeline.RDS$Infer_metrics
#if (!is_empty(pipeline.RDS$insert_size_metrics)) {
  rshiny_input$insert_size_metrics <- pipeline.RDS$insert_size_metrics
#}

rshiny_input$read_distribution_percent <- pipeline.RDS$read_distribution %>%
  dplyr::select(contains("xaxis_"), contains("_percent"), contains("norm_"))
rshiny_input$read_distribution_perkb <- pipeline.RDS$read_distribution %>%
  dplyr::select(contains("xaxis_"), contains("_perKb"))
rshiny_input$splicing_junctions <- pipeline.RDS$junctions
rshiny_input$volcanto_plot <- pipeline.RDS$deg_analysis %>%
  dplyr::select(AveExpr, contains("xaxis_"), contains("yaxis_"))
rshiny_input$ma_plot <- pipeline.RDS$deg_analysis %>%
  dplyr::select(AveExpr, contains("xaxis_"), contains("yaxis_"))
rshiny_input$heatmap <- pipeline.RDS$deg_analysis %>%
  dplyr::select(AveExpr,contains("xaxis_"), contains("yaxis_"), contains("cpm_"))
rshiny_input$mds_plot <- pipeline.RDS$deg_analysis %>%
  dplyr::select(xaxis_ensembl_id, contains("cpm_"))

# MA plot  - per comparision
# --  x-axis : Normalized mean
# -- y-axis : log2 foldchange
# -- pvalue to adjust

# Volcano Plot - per comparision
# -- x-ais: Log2 fold change
# -- y axis: p-value

# Heatmap Plot
# -- x-ais: average cpm counts
# -- y axis: gene names
# pvaue to be controlled

# MDS Plot
# -- x-ais: sample_ids
# -- y axis: cpm

# Clean formatting the data
rshiny_input <- lapply(rshiny_input, as.tibble)
rshiny_input <- lapply(rshiny_input, function(ldf)
  mutate_if(ldf, is.factor, function(x) as.numeric(as.character(x))))

 # rshiny_input <-
 #   lapply(rshiny_input, function(ldf)
 #     mutate_if(ldf, is.numeric, funs(round(., 2))))

# sprintf(c(10.5334,1000), fmt = '%#.3f')
# round(c(10.5334,1000), 3)

# Working -----
# Create a parameter file for Rshiny
# Sidebar Menu names

CreateInputParameters <- function(df, match) {
  matchnames <- colnames(df)[grepl(match, colnames(df))]
  if (match %in% "adj.P.Val") {
    if_else(is_empty(matchnames[grepl("adj.P.Val", matchnames)]) == T, F, T)
  } else {
    matchnames[grepl(match, matchnames)] %>% gsub(match, "", .)
  }
}


lsdtdf <- rshiny_input
# The function does not accept multiple strings, hence added a new col
# for cpm_y_axis. So that, cpm counts are grepped.


rshiny_input$input_parameters <- tibble(
  name = names(lsdtdf),
  x_axis = sapply(lsdtdf, CreateInputParameters, "xaxis_"),
  y_axis = sapply(lsdtdf, CreateInputParameters, "yaxis_"),
  cpm_y_axis = sapply(lsdtdf, CreateInputParameters, "cpm_"),
  normalization_axis = sapply(lsdtdf, CreateInputParameters, "norm_"),
  pval_adjust = sapply(lsdtdf, CreateInputParameters, "adj.P.Val"),
  data_index = sapply(lsdtdf, function(x) if_else(nrow(x) > 0, 1 , 0))

)
# Collapsing values in inputparemters dataframe
# adding the cpm values into the y_axis
# push the cpm counts into the y_axis

# For heatmap we need pvalue, and cpm and as well as
rshiny_input$input_parameters$y_axis$heatmap <-
  c(
    rshiny_input$input_parameters$y_axis$heatmap,
    rshiny_input$input_parameters$cpm_y_axis$heatmap
  )
rshiny_input$input_parameters$y_axis$mds_plot <-
  rshiny_input$input_parameters$cpm_y_axis$mds_plot


#  droppping the cpm_y_XIS
rshiny_input$input_parameters <-
  rshiny_input$input_parameters %>%
  dplyr::select(-c(cpm_y_axis))


# ------------------------------------------------------------------------------
# # Clean up the added prefix strings across all dataframes
adddedstrings <- c("xaxis_", "yaxis_", "norm_", "cpm_")
checkvars <- paste(adddedstrings, collapse = "|")

GsubOverlist <- function(x) {
  gsub(checkvars, "", colnames(x))
}
zz <- lapply(rshiny_input, GsubOverlist)

# remove the added strings
for (i in seq_along(zz)) {
  colnames(rshiny_input[[i]]) <- unlist(zz[i])
}

# check colnames before submitting
lapply(rshiny_input, colnames)



# Create a final output for Rshiny save final output for Rshiny
# save the file RDS file
# saveRDS(
#   rshiny_input,
#   file = here("output_data/input_data.rds"),
#   ascii = FALSE,
#   version = NULL,
#   compress = TRUE,
#   refhook = NULL
# )

write_rds(path = here("output_data/input_data.rds"),
          rshiny_input,compress = "xz")

