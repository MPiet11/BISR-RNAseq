library(tidyverse)
library(data.table)

# read multiple text files and get their path
list_of_files <- list.files(
  path = paste(getwd(), "/rseqc_reports/", sep = ""), recursive = TRUE,
  pattern = "\\.txt", full.names = T
)

# # Read all the files and create a FileName column to store filenames
#list_of_DT <- lapply(list_of_files, fread)
list_of_DT <- lapply(list_of_files, function(x)
  read_delim(x, delim = "\t",trim_ws = T))


find.list <- list(".txt", "rseqc_")
find.string <- paste(unlist(find.list), collapse = "|")
names(list_of_DT) <- gsub(find.string, replacement = "",
                          x = list.files(paste(getwd(), "/rseqc_reports/", sep = ""), pattern = "\\.txt", full.names = FALSE))

# Some df has ID column as X1, change it to sample_id ----
input_data <- list_of_DT
input_data <- lapply(input_data, function(x) setNames(x, sub("X1", "sample_id", names(x))))

# Create hisat2_metrics ----
input_data$hisat2_metrics <-
  merge(input_data$hisat2_uniq_multi_unmap,
    input_data$hisat2_overall_alignment,
    by = "sample_id"
  ) %>% as_tibble()

# Create new df after calculating percent ribosomal contamination ----
input_data$ribo$Percent_contamination <-
  round(
    (input_data$ribo$Ribosomal / (input_data$ribo$Ribosomal + input_data$ribo$Not_ribosomal)) * 100,
    3
  )


# defining the dataframe
df <- input_data$featurecounts_summary

# featurecounts summary, all rows with zero, exclude the id columns ----
input_data$featurecounts_summary <- input_data$featurecounts_summary[rowSums(abs(input_data$featurecounts_summary[, -1])) != 0, ]

# custom function to transpose while preserving names
transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column() %>%
    tibble::as_tibble()
  return(t_df)
}

# using the function
input_data$featurecounts_summary <- transpose_df(input_data$featurecounts_summary)
colnames(input_data$featurecounts_summary) <- input_data$featurecounts_summary[1,]
input_data$featurecounts_summary <- input_data$featurecounts_summary[-1,]

# Matrix holds only one class unlike dataframe
# Hence, even numeric are coerced into character
# Hence, need to convert chr cols & !sample_id
# into numeric

input_data$featurecounts_summary <- input_data$featurecounts_summary %>%
  mutate_if(function(col) is_character(col) &
              !all(col == .$sample_id), as.numeric)

# Read distribution calculate respective percentages ----
dt <- input_data$read_distribution_P12

PercentExons <- function(x) {
  round(((rowSums(x[, c("CDS_Exons", "5'UTR_Exons", "3'UTR_Exons")] / rowSums(
    x
  ))) * 100), 3)
}
PercentIntrons <- function(x) {
  round(((rowSums(x[, c("Introns")] / rowSums(
    x
  ))) * 100), 3)
}
PercentIntergenic <- function(x) {
  round(((rowSums(x[, c("TSS_up_10kb", "TES_down_10kb", "Total_Unassigned_Tags")] /
    rowSums(
      x
    ))) * 100), 3)
}
vector_KB <- input_data$read_distribution_P3
exon_length <- rowSums(vector_KB[, c("CDS_Exons", "5'UTR_Exons", "3'UTR_Exons")])
intron_length <- rowSums(vector_KB[, c("Introns")])
tss_length <- rowSums(vector_KB[, c("TSS_up_10kb", "TES_down_10kb")])

dt %>% mutate(
  Exon_percent = PercentExons(dt[, 5:11]),
  Intron_percent = PercentIntrons(dt[, 5:11]),
  Intergenic_percent = PercentIntergenic(dt[, 5:11]),
  Exon_read_perKb = round(rowSums(dt[, 5:7]) / (exon_length / 1000), 3),
  Intron_read_perKb = round(rowSums(dt[, 8]) / (intron_length / 1000), 3),
  `up/downstream_reads_perKb` = round(rowSums(dt[, 9:10]) / (tss_length / 1000), 3)
) -> input_data$read_distribution



## Save an RDS list object
write_rds(path = paste(getwd(),'/raw_data/output_RNAseq_pipeline.rds', sep = ""),
          input_data,compress = "xz")
# saveRDS(input_data, file = paste(getwd(),'/raw_data/output_RNAseq_pipeline.rds', sep = ""), ascii = FALSE, version = NULL, compress = TRUE)
# Save an RDS list object
#input_data <- readRDS(paste(getwd(), "/../data/clean/Input_RNAseq_Shiny.rds"))


