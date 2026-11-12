library(data.table)
library(tidyverse)

panel_genes <- readLines("gene_panel.lst")

aml_ohsu_2018_mut <- fread("./aml_ohsu_2018/data_mutations.txt") %>%
		filter(Hugo_Symbol %in% panel_genes) %>%
		separate(HGVSc, into=c("ENST", "HGVSc"), sep=":") %>%
		mutate(DATASET = "aml_ohsu_2018") %>%
		select(Hugo_Symbol,HGVSc,HGVSp,DATASET) %>%
		distinct() %>%
		mutate(across(everything(), ~ na_if(.x, "")))


aml_ohsu_2022_mut <- fread("./aml_ohsu_2022/data_mutations.txt") %>%
		filter(Hugo_Symbol %in% panel_genes) %>%
		separate(HGVSc, into=c("ENST", "HGVSc"), sep=":") %>%
		mutate(DATASET = "aml_ohsu_2022") %>%
		select(Hugo_Symbol,HGVSc,HGVSp,DATASET) %>%
		distinct() %>%
		mutate(across(everything(), ~ na_if(.x, "")))

aml_target_gdc_mut <- fread("./aml_target_gdc/data_mutations.txt") %>%
	filter(Hugo_Symbol %in% panel_genes) %>%
	separate(HGVSc, into=c("ENST", "HGVSc"), sep=":") %>%
	mutate(DATASET = "aml_target_gdc") %>%
	select(Hugo_Symbol,HGVSc,HGVSp,DATASET) %>%
	distinct() %>%
	mutate(across(everything(), ~ na_if(.x, "")))

aml_tcga_gdc_mut <- fread("./aml_tcga_gdc/data_mutations.txt") %>%
	filter(Hugo_Symbol %in% panel_genes) %>%
	separate(HGVSc, into=c("ENST", "HGVSc"), sep=":") %>%
	mutate(DATASET = "aml_tcga_gdc") %>%
	select(Hugo_Symbol,HGVSc,HGVSp,DATASET) %>%
	distinct() %>%
	mutate(across(everything(), ~ na_if(.x, "")))

laml_tcga_pan_can_atlas_2018_mut <- fread("./laml_tcga_pan_can_atlas_2018/data_mutations.txt") %>%
	filter(Hugo_Symbol %in% panel_genes) %>%
	separate(HGVSc, into=c("ENST", "HGVSc"), sep=":") %>%
	mutate(DATASET = "laml_tcga_pan_can_atlas_2018") %>%
	select(Hugo_Symbol,HGVSc,HGVSp,DATASET) %>%
	distinct() %>%
	mutate(across(everything(), ~ na_if(.x, "")))

laml_tcga_pub_mut <- fread("./laml_tcga_pub/data_mutations.txt") %>%
	filter(Hugo_Symbol %in% panel_genes) %>%
	separate(HGVSc, into=c("ENST", "HGVSc"), sep=":") %>%
	mutate(DATASET = "laml_tcga_pub") %>%
	select(Hugo_Symbol,HGVSc,HGVSp,DATASET) %>%
	distinct() %>%
	mutate(across(everything(), ~ na_if(.x, "")))

laml_tcga_mut <- fread("./laml_tcga/data_mutations.txt") %>%
	filter(Hugo_Symbol %in% panel_genes) %>%
	separate(HGVSc, into=c("ENST", "HGVSc"), sep=":") %>%
	mutate(DATASET = "laml_tcga") %>%
	select(Hugo_Symbol,HGVSc,HGVSp,DATASET) %>%
	distinct() %>%
	mutate(across(everything(), ~ na_if(.x, "")))


merged_aml <- rbind(
	aml_ohsu_2018_mut, 
	aml_ohsu_2022_mut, 
	aml_target_gdc_mut,
	aml_tcga_gdc_mut,
	laml_tcga_pan_can_atlas_2018_mut,
	laml_tcga_pub_mut,
	laml_tcga_mut
	) %>%
	distinct(Hugo_Symbol,HGVSc,HGVSp, .keep_all=TRUE)

colnames(merged_aml) <- c("GENE_SYMBOL","MUTATION_CDS","MUTATION_AA","DATABASE")

write_tsv(merged_aml, "PUBLIC_AML_DATA.tsv")