#!/usr/bin/env Rscript

####### LOAD LIBRARIES #########
#BiocManager::install("apeglm")
print("Loading Required Packages...")
suppressPackageStartupMessages({
library (DESeq2) #load the package
library(tximport)
library(readr)
library(tidyverse)
library(stringr)
library(argparse)
library(glue)
})
################################

############# DIRECTORIES AND FILES --> GLOBAL VARIABLES #######
parser = ArgumentParser(description = "Differential Expression Analysis for RNAseq Data ")

# optional arguments with defaults
parser$add_argument("--organism",
                    default = "hsapiens", 
                    type = "character",
                    help = "organism options: ['hsapiens', 'taestivum']"
)
parser$add_argument("--count_file",
                    default = "",
                    type = "character",
                    help = "File countaining raw RNA expression counts."  #counts should not be normalized.
)
parser$add_argument("--out_dir",
                    default = "",   #default output is where the count_file is located
                    type = "character",
                    help = "Output Directory for Files"
)
parser$add_argument("--nskip",
                    default = 1,     #skip the first line which contains a note, then you can read in matrix as a delimiter file.
                    type = "integer",
                    help = "Number of rows in the matrix count file to skip before getting to the column heading row. Only necessary if count file is provided"
)
parser$add_argument("--alpha",
                    default = 0.05,
                    type = "numeric",
                    help = "Critical Value (Default=0.05)"
)
parser$add_argument("--gene_out_type",
                    default = "refseq_gene",   #refseq_gene or refseq_tx
                    type = "character",
                    help = "Gene output type to be written to file. Options ['symbol', 'entrez', 'refseq_gene', 'refseq_tx']"
)
parser$add_argument("--raw_file_dir",
                    default = "",      #../../geo_data/wheat_test1/RAW/
                    type = "character",
                    help = "Directory containin raw salmon .sf files if provided"
)
parser$add_argument("--tx2gene_file",
                    default = "../../geo_data/wheat_tx2gene_convert.txt",   #
                    type = "character",
                    help = "file containing table mapping transcripts to genes and any other gene name formats"
)
parser$add_argument("--col_map",
                    default = "",   #
                    type = "character",
                    help = "file containing colData for DEseq. Must contain headers 'sampleid', 'coef[1-9]+', and 'replicate'. Can contain other informatory columns as well."
)
parser$add_argument("--ref_seq_cvt_file",
                    default = "",   #
                    type = "character",
                    help = "file containing information on converting between different refseq versions (should be tab delimited between version names)"
)
parser$add_argument("--cvt_direction",
                    default = "",   #
                    type = "character",
                    help = "Direction of conversion (can be '' for no conversion, or 'v2.1_v1.1' or 'v1.1_v2.1')"
)

args = parser$parse_args()

organism = args$organism
gene_out_type = args$gene_out_type
nskip = args$nskip
out_dir = args$out_dir
raw_file_dir = args$raw_file_dir
tx2gene_file = args$tx2gene_file
col_map = args$col_map
count_file = args$count_file
alpha = args$alpha
cvt_file = args$ref_seq_cvt_file
cvt_direction = args$cvt_direction

########################################

### Handle creation of output directory 
create_output_dir = function(out_dir){

  if (out_dir != ""){
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  } else {
    dir.create("out", recursive = TRUE, showWarnings = FALSE)
  }
  return(out_dir)
}


### This function reads in count matrix created from .bam files
read_in_count_file = function(count_file, nskip){

  count_df = read_delim(count_file, delim = "\t", skip = nskip, show_col_types = F)
  
  count_mtrx = count_df  %>% select(grep(".bam", names(count_df)))
  colnames(count_mtrx) = str_extract(basename(colnames(count_mtrx)), pattern = "^[A-Za-z0-9]+") #everything up do \.
  count_mtrx = as.matrix(count_mtrx)
  rownames(count_mtrx) = count_df$Geneid
  
  return(count_mtrx)
}


#### These functions are for if the count data comes from a Salmon alignment
# Point to all your .sf files and dload as a list (with nicer names)
read_salmon_files = function(raw_file_dir){
  
  salmon_files <- list.files(raw_file_dir, pattern = "\\.sf$", full.names = TRUE)
  names(salmon_files) <- gsub("_quant.sf", "", basename(salmon_files))
  return(salmon_files)
}


# get a transcript-to-gene mapping table
load_tx2gene_table = function(tx2gene_file, salmon_files=FALSE){
  
  #if a tx2gene file is given use it
  if (tx2gene_file != ""){
    tx2gene = read_delim(tx2gene_file, delim="\t", show_col_types = F)
	print("tx2gene file loaded")
  }
  #otherwise just map all transcripts to the base RefSeq gene   
  else if (salmon_files){
    warning("No transcript to gene file provided. Creating tx2gene table for RefSeq transcript types only.")
    tx2gene = read_delim(salmon_files[[1]], delim="\t", show_col_types = F) %>% 
      dplyr::select(Name) %>% 
      mutate("gene_id" = sub("\\.\\d+$", "", Name)) %>%    # e.g. TraesCS1A03G0004200.1 -> TraesCS1A03G0004200
      rename("tx_id" = "Name")
  } else {
    stop("No tx2gene table created")
  }
  return(tx2gene)
}


#load transcript abundance data and combine from all raw files.
import_transcript_abundance_data = function(salmon_files, tx2gene, txOut = FALSE){

  if (txOut){
    txi <- tximport(salmon_files, type = "salmon", txOut = txOut)  #keep analysis at transcript level
  } else {
    txi <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene)  #convert transcripts to respective genes
  }
  print("Transcript Data Successfully Imported")
  return(txi)
}


read_in_col_map = function(col_map){

  cat("Loading ColData")
  col_map_df = read_delim(col_map, show_col_types = F)
  col_map_df = col_map_df  %>% mutate(across(starts_with("coef"), ~factor(.x, levels = sort(unique(.x))))) #convert coefficients to factors
  
  cat("\nFactor Levels:\n")
  print(col_map_df %>% select(starts_with("coef")) %>% lapply(levels))
  cat("\n")

  return(col_map_df)
}

# use this to determine the interaction order for the design matrix based on the number of replicates in the sample
# it assumes replicates are labelled in the column 'replicate' and coeficients are in 'coef[1-9]+' column
# interaction order will only go as high as ^3 and is determined by the number of replicates for each biological state.
get_interaction_order = function(col_map_df){
  
  cat("Getting interaction order for design...\n")
  coefs_df = dplyr::select(col_map_df, grep("coef", ignore.case = TRUE, names(col_map_df))) 
  n_coefs = ncol(coefs_df)

  n_levels <- apply(coefs_df, 
                    MARGIN = 2, 
                    FUN = function(x) length(unique(x)))
  min_reps = min(table(col_map_df$replicate))
  total_samples = prod(n_levels)*min_reps
  
  # Number of parameters for a given interaction order
  n_params = function(order) {
    1 + sum(sapply(1:order, function(k) choose(n_coefs, k)))
  }
  
  cat(sprintf("n coefficients: %d\n", n_coefs))
  cat(sprintf("Total samples: %d\n\n", total_samples))
  
  # Only check up to order 2 maximum (^2)
  max_order <- min(n_coefs, 2)
  feasible_order <- 1
  
  for (order in 1:max_order) {
    params <- n_params(order)
    min_needed <- params * min_reps
    feasible <- total_samples >= min_needed
    
    if (feasible) feasible_order <- order
  }
  
  cat(sprintf("\nRecommended interaction order: %d\n", feasible_order))
  return(feasible_order)
}

## This is to drop interactions where not all subcategories are present (leads to 0s in matrix)
## Used to determine which variables will have interactions in the 'design'
test_complete_interactions = function(col_map_df, test_pairs){

  count_n_reps = col_map_df %>%
    count(across(all_of(test_pairs))) %>%
    complete(!!!syms(test_pairs), fill = list(n = 0))

  if (0 %in% count_n_reps$n){
    return(FALSE)  #This combination of pairs has some missing combinations (cannot use as interaction)
  } else {
    return(TRUE)   #This combination of pairs has no missing combos, can use as interaction
  }
  
}


## Design created from parsing col_map file for all columns beginning with coef[1-9]+.
## design_type can be 'main_effects', 'pairwise', 'all'  (where all will do all interaction terms if there are > 2 coefficients)
create_design = function(col_map_df){
  
  int_order = get_interaction_order(col_map_df)  #gets the interaction order for the design based on number of coefficients and replicates.
  coefs = names(col_map_df)[grep("coef[1-9]+", ignore.case = T, names(col_map_df))]

  ## First order intearction
  if (int_order == 1){
    design = paste("~(", paste(coefs, collapse = " + "), ")")

  ## Second Order Interaction
  } else if (int_order == 2){
    coefs_for_design = c(coefs)   #start with 1st order terms
    all_combos = combn(coefs, 2, simplify = FALSE)  #gives all pairwise combinations of coef_cols

    for (i in 1:length(all_combos)){                 #need to check that there are no missing combinations
      no_zeros = test_complete_interactions(col_map_df, all_combos[[i]])
      if (no_zeros){
        coefs_for_design = c(coefs_for_design, paste(all_combos[[i]], collapse = "*"))
        }
    }
    design = paste("~(", paste(coefs_for_design, collapse = "+"), ")")

  ## Third Order Interaction
  } else if (int_order == 3){
    design = paste("~(", paste(coefs, collapse = "+"), ")^3")
  } else {
    stop(glue("Interaction order not determined.  Interaction Order:  {int_order}"))
  }
  return(as.formula(design))
}

### This is where the differential expression analysis will be performed, list of all result objects will be returned
## Each result element in the list is one pairwise comparison for each coefficient in the GLM
# shrinkage is applied to not overestimate fold changes for lowly expressed genes
perform_de_analysis = function(col_map, design, txi = NULL, count_mtrx = NULL, alpha = 0.05){

  cat("ColData Headers:\t", colnames(col_map), "\n")
  cat("Design:\t", as.character(design), "\n\n")
  if (!is.null(txi)){
    dds = DESeqDataSetFromTximport(txi = txi, colData = col_map, design = design)
  } else if (!is.null(count_mtrx)){
    dds = DESeqDataSetFromMatrix(countData = count_mtrx,  colData = col_map, design = design)
  }

  # Collapse if there are multiple technical replicates for the same specimen
  if (nrow(col_map) > length(unique(col_map$sample_id))){
    warning("There are multiple technical replicates per specimen, summing replicates for each specimen.")
    dds = collapseReplicates(dds, groupby = col_map$sample_id)  #collapse 
  }

  #differential expression analysis
  dds = DESeq(dds)
  res_names = resultsNames(dds)  #this is a list of the coefficients for the linear model (specified in design)

  all_res = list()
  for (n in res_names[2:length(res_names)]){
    print(glue("Getting Results for {n}"))
    #res = lfcShrink(dds, coef = n, type = "apeglm", returnList = T)  #apply shrinkage for lowly expressed genes
    res = results(dds, name = n, alpha = alpha)  #2-tailed p-values
    all_res[[n]] = res
  }
  return(all_res)
}



## Convert from refseq 1.1 Traes to 2.1 traes
convert_refseq = function(cvt_file, vec_to_cvt, cvt_direction, out_dir){

  cat("Old vector:", vec_to_cvt[1:5], "\n")
  cvt_df = read_delim(cvt_file, show_col_types = F)
  
  if (cvt_direction == "v2.1_v1.1"){
    lookup = setNames(object = cvt_df$v1.1, cvt_df$v2.1)  #object, names
  } else if (cvt_direction == "v1.1_v2.1"){
    lookup = setNames(object = cvt_df$v2.1, cvt_df$v1.1)  #headers of conversion file should be v1.1 and v2.1, etc.
  } else {
    stop(glue("refseq conversion direction must be 'v2.1_v1.1' or 'v1.1_v2.1', not {cvt_direction}"))
  }
  
  new_vec = character(length(vec_to_cvt))
  missing_ids = c()
  cat("Converting Refseq version IDs:  ", cvt_direction, "...\n")
  for (i in 1:length(vec_to_cvt)){

    if (!vec_to_cvt[i] %in% names(lookup)){
      missing_ids = c(missing_ids, vec_to_cvt[i])
      new_vec[i] = NA
    } else {
      new_vec[i] = lookup[vec_to_cvt[i]]
    }
  
  if (i %% 1000 == 0){
     cat(i, ":\tOld ID:",  vec_to_cvt[i], "\tNew ID:", lookup[vec_to_cvt[i]], "\n")  ## This is a progress updater
  }
  }

  cat("New vector:", new_vec[1:5], "\n")
  missing_filename = file.path(out_dir, paste0("missing_refseqs_", str_extract(cvt_direction, ".*(?=_)"), ".txt"))
  writeLines(missing_ids, missing_filename)
  cat(length(missing_ids),
   "missing IDs in refseq version conversion file.\nWritten to ", 
   missing_filename, "\n")

  return(new_vec)
}



## extract relevant data from results for ORA and GSEA
gather_gsea_ora_data = function(res, alpha=0.05){

  gsea_ranked = as_tibble(res, rownames = "gene") %>% 
    filter(!is.na(`stat`)) %>%
    arrange(desc(`stat`)) %>%    #order by ranked test statistic (keep all genes)
    dplyr::select("gene", "stat")

  ora_sig <- subset(res, padj < alpha) %>%  #Only keep significant genes for ORA
    as_tibble(rownames = "gene") %>% 
    filter(!is.na(padj)) %>%
    dplyr::select(gene)
  
  return(list("gsea_ranked"=gsea_ranked, "ora_sig"=ora_sig))
}


write_gsa_files = function(out_dir, filename, results, gsea_ranked, ora_sig, gene_out_type){

  ### GSEA FILE (2 columns needed)
  gsea_file = file.path(out_dir, paste0("GSEA_", filename, "_", gene_out_type, "_genes_of_interest.rnk"))
  write_delim(
    gsea_ranked,
    file = gsea_file,
    delim = "\t",
    col_names = F
  )

  ### ORA ##  
  ##### Input for ORA ##  --> #Just a list of genes of interest####
  ora_goi_file = file.path(out_dir, paste0("ORA_", filename, "_", gene_out_type, "_genes_of_interest.txt"))
  write.table(
    ora_sig,    #A list of all the significantly differentiated genes.
    ora_goi_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )

  ## Also need to create a reference gene file for ORA ##
  ora_reference_file = file.path(out_dir, paste0("ORA_", gene_out_type, "_reference_genes.txt"))
  write.table(
    rownames(results),  #all genes
    ora_reference_file,
    row.names = FALSE,
    col.names = FALSE,
    quote = FALSE
  )
  
  ## Maybe I should also write a file with the parameters used.??
  
  print(glue("Created File Output: {gsea_file}"))
  print(glue("Created File Output: {ora_goi_file}"))
  print(glue("Created File Output: {ora_reference_file}"))
}


#write all extra data to a file
write_meta_data = function(results, col_map, filename, out_dir, gene_out_type){
  
  meta_out = file.path(out_dir, "extra_data")
  dir.create(meta_out, recursive=TRUE, showWarnings=FALSE)

  # Write raw results and meta data to file
  results_with_rownames <- cbind("Gene" = rownames(results), results)
  full_results = file.path(meta_out, paste0("full_results_", filename, "_", gene_out_type, ".txt"))

  write.table(
    results_with_rownames,
    full_results,
    col.names = TRUE,
    row.names = FALSE
  )
  result_params = file.path(meta_out, paste0("result_params_", filename, ".txt"))
  write.table(
    mcols(results, use.names = TRUE) ,
    result_params,
    col.names = NA
  )
  print(glue("Created File Output: {full_results}"))
  print(glue("Created File Output: {result_params}"))
}

## Collect appropriate info and write to files.
parse_results = function(all_res, col_map, out_dir, gene_out_type){

  for (r in names(all_res)){
    output = gather_gsea_ora_data(all_res[[r]])  #function from sourced file
    gsea_ranked = output$gsea_ranked
    ora_sig = output$ora_sig
    rm(output)
    
    write_gsa_files(out_dir, r, all_res[[r]], gsea_ranked, ora_sig, gene_out_type)  #write files for gsa in webgestalt
    write_meta_data(all_res[[r]], col_map, r, out_dir, gene_out_type)
  }
}


## The head cheif function --> calling all functions
main_func = function(){
  
  out_dir = create_output_dir(out_dir)
  
  ### NEED TO DETERMINE THE FILE INPUT TYPE TO KNOW WHICH FUNCTIONS TO CALL ##
  ### If Input type is a raw .sf files, need txi import specific functions
  if (raw_file_dir != ""){   #If directory for RAW .sf files are provided.
    input_type = "tx" 
    salmon_files = read_salmon_files(raw_file_dir)
    tx2gene_full = load_tx2gene_table(tx2gene_file, salmon_files = salmon_files)  #Will use tx2gene file first before creating its own
    
    #Reduce table to the format required by the function Transcript --> Gene
    tx2gene_red = tx2gene_full %>% 
      dplyr::select("Transcript stable ID", "Gene stable ID")
    
    if (gene_out_type == "refseq_tx"){
      txOut = TRUE      #do DE analysis at transcript level
    } else {
      txOut = FALSE
    }
    txi = import_transcript_abundance_data(salmon_files, tx2gene_red, txOut)
    
    ## Convert 2.1 to 1.1 if needed.
    if (cvt_direction != ""){
      rownames(txi) = cvt_refseq(cvt_file, rownames(txi$counts), cvt_direction, out_dir)
      before_NA_removal = nrow(txi$counts)
      txi$counts = txi$counts[!is.na(rownames(txi$counts)), ]
      cat(before_NA_removal-nrow(txi$counts), " Rows containing GeneIDs as NA were dropped.\n")
      cat(nrow(txi$counts), "Rows Remaining in Count Matrix\n\n")
    }
    ### If input type is matrix of counted read mappings
    } else if (count_file != ""){
      input_type = "count_mtrx"
      count_mtrx = read_in_count_file(count_file, nskip)
    
    if (!is.integer(count_mtrx)){
      warning("Matrix countains fractional amounts: Rounding...")
      count_mtrx = round(count_mtrx, 0)  #round fractional counts (from fractional multimapping reads.)
    }

    if (cvt_direction != ""){
      rownames(count_mtrx) = convert_refseq(cvt_file, rownames(count_mtrx), cvt_direction, out_dir)
    }
    #Remove any rows where GeneIds are NA
    before_NA_removal = nrow(count_mtrx)
    count_mtrx = count_mtrx[!is.na(rownames(count_mtrx)), ]
    cat(before_NA_removal-nrow(count_mtrx), " Rows containing GeneIDs as NA were dropped.\n")
    cat(nrow(count_mtrx), "Rows Remaining in Count Matrix\n\n")
  }
  
  #Read in colData info
  col_map_df = read_in_col_map(col_map)
  
  ## Determine design of the linear model
  design = create_design(col_map_df)
  
  #Run DE analysis
  if (input_type == "tx"){
    all_res = perform_de_analysis(txi = txi, col_map = col_map_df, design = design, alpha = alpha)
  } else if (input_type == "count_mtrx"){
    all_res = perform_de_analysis(count_mtrx = count_mtrx, col_map = col_map_df, design = design, alpha = alpha)
  }

  # Extract valuable information from results
  parse_results(all_res, col_map_df, out_dir, gene_out_type)
}

if (!interactive()) {
  main_func()
  warnings()
  print("FINISHED")
}

