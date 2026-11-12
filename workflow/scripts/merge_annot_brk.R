args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)

brk1_file <- args[1]
brk2_file <- args[2]
info_file <- args[3]
merged_file <- args[4]

brk1_annot <- read_tsv(brk1_file,
             col_names = TRUE,
             show_col_types = FALSE)
brk2_annot <- read_tsv(brk2_file,
             col_names = TRUE,
             show_col_types = FALSE)
info_field <- read_tsv(info_file,
             col_names = TRUE,
             show_col_types = FALSE)

if (nrow(brk1_annot) == 0 || nrow(brk2_annot) == 0 || nrow(info_field) == 0) {

  colnames(brk1_annot) <- paste(colnames(brk1_annot), "brk1",  sep="_")
  colnames(brk2_annot) <- paste(colnames(brk2_annot), "brk2",  sep="_")

  df_merged <- cbind(brk1_annot, brk2_annot, info_field)

  write_tsv(df_merged, merged_file)

} else {

  brk1_annot <- read_tsv(brk1_file,
               col_names = TRUE,
               show_col_types = FALSE)
  brk2_annot <- read_tsv(brk2_file,
               col_names = TRUE,
               show_col_types = FALSE)
  info_field <- read_tsv(info_file,
               col_names = TRUE,
               show_col_types = FALSE)

  df_merged <- brk1_annot %>%
    full_join(
      brk2_annot, by = "name", suffix= c("_brk1", "_brk2")
      ) %>%
    full_join(
      info_field, by = "name"
      ) %>%
    arrange(
      seqnames_brk1, start_brk1, end_brk1
      ) %>%
    select(
      -width_brk1,-strand_brk1,-width_brk2,-strand_brk2
      )

  write_tsv(df_merged, merged_file)

}