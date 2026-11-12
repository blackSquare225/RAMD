args = commandArgs(trailingOnly=TRUE)

library(openxlsx)
library(data.table)

vcf_inv <- args[1] 
vcf_bnd <- args[2] 
chr3_rearr <- args[3]
report <- args[4]

vcf_inv <- fread(vcf_inv)
vcf_bnd <- fread(vcf_bnd)
chr3_rearr <- fread(chr3_rearr)

wb <- createWorkbook()

addWorksheet(wb, "INV")
writeData(wb, "INV", vcf_inv)

addWorksheet(wb, "BND")
writeData(wb, "BND", vcf_bnd)

addWorksheet(wb, "chr3_rearr")
writeData(wb, "chr3_rearr", chr3_rearr)

saveWorkbook(wb, report, overwrite = TRUE)


