library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)

fi <- args[1]
cosmic_db <- args[2]
fp_db <- args[3]
benign_db <- args[4]
vus_db <- args[5]
diagnostic_db <- args[6]
public_db <- args[7]
fo <- args[8]

# Load input file
df <- read_tsv(fi, show_col_types = FALSE)

# Remove population variants
popAF_cols <- c("EUR_AF","gnomADe_AF","gnomADe_NFE_AF","gnomADg_AF","gnomADg_NFE_AF")
is_popVars <- apply(df[, popAF_cols], 1, function(row) any(!is.na(row) & row >= 0.01))
df$label <- ifelse(is_popVars, "pop_vars", "")
df <- df %>% 
  filter(label != "pop_vars") %>%
  select(-label) 

# Filter variants according to some technical requirements PVAL cut-off to 0.001
df <- df %>% 
  mutate(
    FREQ = as.numeric(gsub("%", "", FREQ))
    ) %>%
  filter(FREQ >= 3) %>% # allele frequency cut-off to 3%
  filter(PVAL < 0.001) %>% # PVAL cut-off to 0.001
  filter(IMPACT %in% c("HIGH", "MODERATE")) %>% # ENSEMBL-VEP annotation
  separate(HGVSc, into = c("ENST", "MUTATION_CDS"), sep = ":", extra="warn")  

# Load databases
COSMIC_db <- read_tsv(cosmic_db, show_col_types = FALSE)
FP_db <- read_tsv(fp_db, show_col_types = FALSE)
BENIGN_db <- read_tsv(benign_db, show_col_types = FALSE)
VUS_db <- read_tsv(vus_db, show_col_types = FALSE)
DIAGNOSTIC_db <- read_tsv(diagnostic_db, show_col_types = FALSE)
PUBLIC_db <- read_tsv(public_db, show_col_types = FALSE)

# Separate the HGVSp field
df <- df %>% separate(HGVSp, into = c("ENSP", "MUTATION_AA"), sep = ":", extra = "warn")

# Helper function to annotate variants based on AA/CDS priority
annotate <- function(df, db, label) {
  # AA-level match
  df_aa <- df %>%
    filter(!is.na(MUTATION_AA)) %>%
    left_join(
      db %>% select(GENE_SYMBOL, MUTATION_AA) %>% mutate(match = 1),
      by = c("SYMBOL" = "GENE_SYMBOL", "MUTATION_AA" = "MUTATION_AA")
    )

  # CDS-level match
  df_cds <- df %>%
    filter(is.na(MUTATION_AA) & !is.na(MUTATION_CDS)) %>%
    left_join(
      db %>% select(GENE_SYMBOL, MUTATION_CDS) %>% mutate(match = 1),
      by = c("SYMBOL" = "GENE_SYMBOL", "MUTATION_CDS" = "MUTATION_CDS"),
    )

  # Combine and fill unmatched with 0
  bind_rows(df_aa, df_cds) %>%
    mutate(!!label := if_else(is.na(match), 0, 1)) %>%
    select(-match)
}

# Apply annotation functions
df_annot <- df %>%
  annotate(COSMIC_db, "COSMIC") %>%
  annotate(FP_db, "FP") %>%
  annotate(BENIGN_db, "BENIGN") %>%
  annotate(VUS_db, "VUS") %>%
  annotate(DIAGNOSTIC_db, "DIAGNOSTIC") %>%
  annotate(PUBLIC_db, "PUBLIC_AML_DATA") %>%
  distinct()

# Write final output
write_tsv(df_annot, fo)
