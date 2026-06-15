#This script will be used to create the colData dataframe required for deSEQ2
# the coefficients specify the design (variables being compared in the differential expression analysis)

## libraries ##
library(tidyverse)
library(stringr)
library(glue)

############# GLOBAL VARIABLES #############
OUT_DIR = "geo_gse118126"
META_TSV = file.path(OUT_DIR, "meta/SRP156347.metadata.tsv")
COUNT_FILE = file.path(OUT_DIR, "v2.1/counts/read_counts.txt")
NSKIP = 1

############################################

# extracting count data column/header info.
count_df = read_delim(COUNT_FILE, delim = "\t", skip = NSKIP, show_col_types = F)

count_mtrx = count_df  %>% 
            select(grep(".bam", names(count_df)))
colnames(count_mtrx) = str_extract(basename(colnames(count_mtrx)), pattern = "^[A-Za-z0-9]+")



## Creating colData aka col_map
meta_df = read_tsv(META_TSV, show_col_types = F)
experiment_samples = meta_df$experiment_title
genotype = regmatches(experiment_samples, regexpr("(?<=GSM[0-9]{7}: )[A-Za-z0-9]+", experiment_samples, perl = TRUE))
tissue = regmatches(experiment_samples, regexpr("(?<=\\([A-Z]{2}\\)) [A-Za-z]+", experiment_samples, perl = TRUE))
treatment = regmatches(experiment_samples, regexpr("(?<=treated with) [A-Za-z0-9]+", experiment_samples, perl = TRUE))
dpi = regmatches(experiment_samples, regexpr("[0-9]+(?= dpi)", experiment_samples, perl = TRUE))
replicate = str_extract(experiment_samples, pattern = "(?<=replicate)[1-9]+")


col_map = tibble("accession" = trimws(meta_df$run_accession), 
                 "sample_id" = trimws(meta_df$sample_accession),
                 "genotype" = trimws(genotype),
                 "tissue" = trimws(tissue),
                 "treatment" = trimws(treatment),
                 "dpi" = trimws(dpi),
                 "replicate" = trimws(replicate)
                )
col_map = col_map  %>% 
    mutate("coef1" = genotype,
           "coef2" = tissue,
           "coef3" = treatment,
           "coef4" = dpi
           )

## Ensure first column of colData equals the matrix colname
if (any(colnames(count_mtrx) != meta_df$run_accession)){
    cat("Matrix Colnames:", colnames(count_mtrx), "\n")
    cat("ColData Colanmes:", meta_df$study_accession, "\n\n")
    stop("Matrix column names do not align with col_map column1 names")
}

## write col_map to file.
  write.table(
    col_map,    
    file.path(OUT_DIR, "meta/colData.tsv"),
    col.names = TRUE,
    row.names = FALSE
  )

