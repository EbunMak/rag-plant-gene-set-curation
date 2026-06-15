#This script will be used to create the colData dataframe required for deSEQ2
## coefficients are used in the deSEQ2 parameter 'design' to specify the variables being compared for differential expression analysis

## libraries ##
library(tidyverse)
library(stringr)
library(glue)

############# GLOBAL VARIABLES #############
OUT_DIR = "geo_gse232243_24s"
META_TSV = file.path(OUT_DIR, "meta/SRP437013.metadata.tsv")
COUNT_FILE = file.path(OUT_DIR, "counts/read_counts.txt")
NSKIP = 1

############################################

# extracting count data column/header info.
count_df = read_delim(COUNT_FILE, delim = "\t", skip = NSKIP, show_col_types = F)

count_mtrx = count_df  %>% 
            select(grep(".bam", names(count_df)))
colnames(count_mtrx) = str_extract(basename(colnames(count_mtrx)), pattern = "^[A-Za-z0-9]+")

## Creating colData aka col_map
meta_df = read_tsv(META_TSV, show_col_types = F)
experiment_samples = meta_df$experiment_desc


replicate = regmatches(experiment_samples, regexpr("(?<=dpi_)[1-9]+", experiment_samples, perl = TRUE))
dpi = regmatches(experiment_samples, regexpr("(?<=_)[1-9](?=dpi)", experiment_samples, perl = TRUE))
strain = regmatches(experiment_samples, regexpr("(?<=Ta_).*(?=_[1-9]dpi)", experiment_samples, perl = TRUE))


col_map = tibble("accession" = trimws(meta_df$run_accession), 
                 "sample_id" = trimws(meta_df$sample_accession),
                 "strain" = trimws(strain),
                 "dpi" = trimws(dpi),
                 "replicate" = trimws(replicate)
                )
col_map = col_map  %>% 
    mutate("coef1" = strain,
           "coef2" = dpi
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

