args <- commandArgs(trailingOnly = TRUE)

library(rtracklayer)
#library(GenomicFeatures)
#library(GenomicRanges)
library(tidyverse)

# ---------- PARSE GTF FILE INTO BED ----------
gtf_file <- args[1] 
bed_file <- args[2]

gtf <- import(gtf_file)
df <- as.data.frame(gtf) %>%
  filter(type != "transcript") %>%
  select(seqnames, start, end, width, strand, type, gene_id, gene_type, gene_name, tag, transcript_id, transcript_type, transcript_name, exon_number, exon_id)

write.table(df, bed_file, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)



