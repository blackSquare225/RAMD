args = commandArgs(trailingOnly=TRUE)

library(openxlsx)
library(data.table)

QDNAseq_50kbp <- args[1] 
QDNAseq_100kbp <- args[2] 
ichorCNA <- args[3]
report <- args[4]

QDNAseq_50kbp <- fread(QDNAseq_50kbp)
QDNAseq_100kbp <- fread(QDNAseq_100kbp)
ichorCNA <- fread(ichorCNA)

wb <- createWorkbook()

addWorksheet(wb, "QDNAseq_50kbp")
writeData(wb, "QDNAseq_50kbp", QDNAseq_50kbp)

addWorksheet(wb, "QDNAseq_100kbp")
writeData(wb, "QDNAseq_100kbp", QDNAseq_100kbp)

addWorksheet(wb, "ichorCNA_50kbp")
writeData(wb, "ichorCNA_50kbp", ichorCNA)

saveWorkbook(wb, report, overwrite = TRUE)


