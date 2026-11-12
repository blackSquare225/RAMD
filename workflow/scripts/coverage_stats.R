args <- commandArgs(trailingOnly = TRUE)

library(tidyverse)
library(data.table)
library(R.utils)
library(gridExtra)

fi <- args[1]
plot <- args[2]
fo_tsv <- args[3]
fo_csv <- args[4]

cov_bed <- fread(fi) 

colnames(cov_bed) <- c("CHROM","START","END","GENE","COVERAGE")

exons <- c("RUNX1","STAG2","BCOR","EZH2","TP53","ZRSR2","SF3B1","U2AF1")

bed_for_plot  <- cov_bed %>%
	filter(GENE %in% exons) %>%
	group_by(GENE) %>% 
	mutate(regions = seq(1,length(GENE))) %>%
	ungroup()

sub <- unique(bed_for_plot$regions)
colors <- c(SF3B1 = "#FFDB6D", EZH2 = "#C4961A", TP53 = "#F4EDCA", RUNX1 = "#D16103", U2AF1 = "#C3D7A4", ZRSR2 = "#52854C", BCOR = "#4E84C4", STAG2 = "#293352")

options(warn = -1)
 p <- lapply(unique(bed_for_plot$GENE), function(gene) {
  ggplot(bed_for_plot[bed_for_plot$GENE == gene, ], aes(x = regions, y = COVERAGE, fill = regions)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, fill = colors[gene]) +
    ggtitle(paste(gene)) +
    theme_minimal() +
    theme(legend.position = "none") +
    scale_x_discrete(limits = sub)
})

# Arrange the plots in a grid
pdf(plot)
grid.arrange(grobs = p, ncol = 2)
dev.off()

cov_summary <- cov_bed %>%
	group_by(GENE) %>%
	mutate(MEAN_COV = round(mean(COVERAGE)),
		SD_COV = round(sd(COVERAGE))
		) %>%
	ungroup()

write_tsv(cov_summary, fo_tsv)
write_csv(cov_summary, fo_csv)



