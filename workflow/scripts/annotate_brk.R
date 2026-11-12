args <- commandArgs(trailingOnly = TRUE)

library(rtracklayer)
library(GenomicFeatures)
library(GenomicRanges)
library(tidyverse)
library(openxlsx)

gencode_annot_file <- args[1]
fi <- args[2]
fo <- args[3]

df <- read_tsv(gencode_annot_file,
               col_names = TRUE,
               show_col_types = FALSE)

gencode <- makeGRangesFromDataFrame(
  df,
  seqnames.field = "seqnames",
  start.field    = "start",
  end.field      = "end",
  strand.field   = "strand",
  keep.extra.columns = TRUE
)

annotate_intervals <- function(query, annotation) {

  hits <- findOverlaps(query, annotation)
  
  annot_overlaps <- cbind(
    as.data.frame(query[queryHits(hits)]),
    as.data.frame(mcols(annotation[subjectHits(hits)]))
  )

  nohit <- setdiff(seq_along(query), queryHits(hits))
  annot_nohit <- as.data.frame(query[nohit])

  if (nrow(annot_nohit) != 0) {

    for (col in colnames(mcols(gencode))) {
      annot_nohit[[col]] <- NA   # or "." if you prefer
    }

  } else {

    message("annot_nohit has zero intervals.")
  } 

  annot_all <- rbind(annot_overlaps, annot_nohit)

  collapse_cols <- c(
    "type", 
    "gene_id", 
    "gene_type", 
    "gene_name",
    "tag",
    "transcript_id",
    "transcript_type",
    "transcript_name",
    "exon_number",
    "exon_id"
    )

  # reorder by interval order
  annot_all_collapsed <- annot_all %>%
    mutate(ID = paste(paste(seqnames,start, sep=":"), end, sep="-")) %>%
    group_by(ID) %>%
    summarise(
      # keep all other columns as first value
      across(-all_of(collapse_cols), first),
      
      # collapse only selected columns
      across(all_of(collapse_cols), ~ {
        paste(., collapse = ",")
      }),
      
      .groups = "drop"
  ) 

  return(annot_all_collapsed)
}

brk <- import(fi)
brk_annot <- annotate_intervals(brk, gencode)

write_tsv(brk_annot, fo)

