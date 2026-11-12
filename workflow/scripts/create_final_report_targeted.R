library(xlsx)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

fusion_500kbp <- args[1] 
fusion_1Mbp <- args[2] 
snpIndel <- args[3] 
flt3 <- args[4]
report <- args[5]

fusion_500kbp <- fread(fusion_500kbp)
fusion_500kbp$sample <- gsub(".fastq", "", fusion_500kbp$sample) 

fusion_1Mbp <- fread(fusion_1Mbp)
fusion_1Mbp$sample <- gsub(".fastq", "", fusion_1Mbp$sample) 

snpIndel <- fread(snpIndel)

flt3 <- fread(flt3)

write.xlsx(fusion_500kbp, report, sheetName = "FUSION_500kbp", row.names=FALSE, append = FALSE)
write.xlsx(fusion_1Mbp, report, sheetName = "FUSION_1Mbp", row.names=FALSE, append = TRUE)
write.xlsx(snpIndel, report, sheetName = "ROI", row.names=FALSE, append = TRUE)
write.xlsx(flt3, report, sheetName = "FLT3", row.names=FALSE, append = TRUE)