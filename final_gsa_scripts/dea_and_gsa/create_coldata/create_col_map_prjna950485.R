#This script will be used to create the colData dataframe required for deSEQ2
## Creating col_map bor biosample download: prjna950485

## libraries ##
library(tidyverse)
library(stringr)
library(glue)

############# GLOBAL VARIABLES #############
OUT_DIR = "prjna950485_out"
META_TSV = file.path(OUT_DIR, "meta/PRJNA950485.metadata.tsv")
COUNT_FILE = file.path(OUT_DIR, "counts/read_counts.txt")
NSKIP = 1
############################################

# extracting count data column/header info.
count_df = read_delim(COUNT_FILE, delim = "\t", skip = NSKIP, show_col_types = F)

count_mtrx = count_df  %>% 
            select(grep(".bam", names(count_df)))
colnames(count_mtrx) = str_extract(basename(colnames(count_mtrx)), pattern = "^[A-Za-z0-9]+")


## This is a function that will remove characters that should not end up in the levels for deseq design
## Can Add to this as new ones come up.
remove_unsafe_characters = function(vec){

  new_vec = str_replace_all(vec, "%", "p_")
  return(new_vec)
}


## Creating colData aka col_map
meta_df = read_tsv(META_TSV, show_col_types = F)

exp_samples = meta_df$library_name
strain = str_extract(exp_samples, pattern = ".*_root")
strain = remove_unsafe_characters(strain)
treatment = str_extract(exp_samples, pattern = "(?<=root_).*(?=_)")
treatment = remove_unsafe_characters(treatment)
replicate = str_extract(exp_samples, pattern = "[1-9]$")


col_map = tibble("accession" = meta_df$run_accession, 
                 "sample_id" = meta_df$sample_accession,
                 "strain" = strain,
                 "treatment" = treatment,
                 "replicate" = replicate
                )
col_map = col_map  %>% 
    mutate("coef1" = strain,
           "coef2" = treatment
           )

## Ensure first column of colData equals the matrix colname
if (any(colnames(count_mtrx) != col_map$accession)){
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
