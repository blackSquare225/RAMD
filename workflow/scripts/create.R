library(xlsx)
library(fread)

args <- commandArgs(trailingOnly=TRUE)

snpIndel <- args[1] 
flt3 <- args[2]
report <- args[3]

# Load data 
snpIndel <- fread(snpIndel)
flt3 <- fread(flt3)

# write excel report file
write.xlsx(snpIndel,file=report, sheetName = "ROI", row.names=FALSE, append = FALSE)
write.xlsx(flt3,file=report, sheetName = "FLT3", row.names=FALSE, append = TRUE)
