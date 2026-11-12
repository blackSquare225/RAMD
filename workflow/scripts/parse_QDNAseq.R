args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(tidyverse)

fi <- args[1]
sampleName <- args[2] 
fo <- args[3]

df <- fread(fi) %>% 
	mutate(
		SAMPLE_NAME = sampleName, 
		CHROMOSOME = paste0("chr", CHROMOSOME),
		CN_STATUS = ifelse(LOG2_RATIO_MEAN < 0, "DEL", "DUP")  
		) %>% 
	filter(
		DATAPOINTS >= 10 & (LOG2_RATIO_MEAN < -0.1 | LOG2_RATIO_MEAN > 0.1)
		) %>%
	select(
		CHROMOSOME,START,STOP,CN_STATUS,DATAPOINTS, LOG2_RATIO_MEAN, SAMPLE_NAME
		)

write_tsv(df, fo, col_names = FALSE)