#!/usr/bin/env Rscript


#O######## LIBRARIES #########
print("Loading Packages...")
suppressPackageStartupMessages({
library(WebGestaltR)
library(glue)
library(stringr)
library(argparse)
})

options(warn = 1)

#### ARGUMENTS ##########
parser = ArgumentParser(description = "Run Gene Set Analysis for 'ORA' or 'GSEA' ")

# optional arguments with defaults
parser$add_argument("--method", 
                    default = "ORA", 
                    choices = c("ORA", "GSEA"),
                    help = "Enrichment method [default: ORA]"
                    )
parser$add_argument("--organism",
                    default = "hsapiens",
                    help = "Organism [default: hsapiens], options['hsapiens', 'taestivum']"
                    )
parser$add_argument("--out_dir",
                    default = "webgestalt_out",
                    help = "Path for output directory"
                    )
parser$add_argument("--goi_file",
                    default = "",
                    help = "File containing genes of interest.
                    If --method='ORA', file should be one column of genes and .txt
                    If --method='GSEA', file should be 2 columns of genes(ranked) \\t foldchange and .rnk"
                    )
parser$add_argument("--db_file",
                    default = "",
                    help = "A GMT file containing gene set database info."
                    )
parser$add_argument("--enrich_db_type",
                    default = "ensembl_gene_id",
                    help = "The ID format that the proteins are in for the database GMT file. "
                    )
parser$add_argument("--goi_type",
                    default = "ensembl_gene_id",
                    help = "The ID format that the proteins for the GOI file are in. "
                    )
parser$add_argument("--reference_type",
                    default = "ensembl_gene_id",
                    help = "The ID format that the proteins for the ORA reference file are in. "
                    )
parser$add_argument("--reference_file",
                    default = "",
                    help = "A .txt file containing a list of genes in 1 column.
                    Only required for the 'ORA' method."
                    )
parser$add_argument("--project_name",
                    default = "",
                    help = "Project name to name directory and files"
)
parser$add_argument("--minNum",
                    default = 3,
                    type = "integer",
                    help = "Minimum number of annotated genes"
)
parser$add_argument("--maxNum",
                    default = 1000,
                    type = "integer",
                    help = "Maximum number of annotated genes"
)
parser$add_argument("--threads",
                    default = 1,
                    type = "integer",
                    help = "Number of Threads for Running WebGestalt"
)
parser$add_argument("--perNum",
                    default = 1000,
                    type = "integer",
                    help = "Permutation Number for GSEA (default: 1000)"
)


### Access arguments ###
args = parser$parse_args()

enrich_method = args$method
organism = args$organism
out_dir = args$out_dir
goi_file = args$goi_file
enrich_db_file = args$db_file
enrich_db_type = args$enrich_db_type
goi_type = args$goi_type
reference_gene_file = args$reference_file
reference_type = args$reference_type
project_name = args$project_name
minNum = as.integer(args$minNum)
maxNum = args$maxNum
threads = as.integer(args$threads)
perNum = as.integer(args$perNum)
#################################
#print(glue("Working Directory: ", getwd()))

## Check enrich method is valid ##
if (!enrich_method %in% c("ORA", "GSEA")){
  stop(glue("enrich_method must be either 'ORA', 'GSEA'"))
}

if (!organism %in% c("hsapiens", "taestivum")){
  stop(glue("Organism must be 'hsapiens' or 'taestivum', not {organism}"))
}

### Check the file paths are valid ##
if (!file.exists(enrich_db_file)){
  stop(glue("File {enrich_db_file} not found."))
}
if (!file.exists(goi_file)){
    stop(glue("File {goi_file} not found."))
}
  
#check ORA specifics
if (enrich_method == "ORA"){
  if (!file.exists(reference_gene_file)){
    stop(glue("File {reference_gene_file} not found."))
  }
  #check goi file has one column and is .txt
} else if (enrich_method == "GSEA"){
  ## Check GSEA rank file is a .rnk file with 2 columns.
}


if (organism %in% c("taestivum")){
  cat("Webgestalt does not support organism: ", organism, ".    Converting to 'others'\n")
  organism = "others"   #taestivum not supported in WebGestalt
}

goi_lab = regmatches(basename(goi_file), regexpr("(?<=[GSEA|ORA]_).+(?=_genes_of_interest)", basename(goi_file), perl = TRUE))
print(glue("Performing {enrich_method} analysis for {goi_lab}"))

################ RUN WEBGESTALT #####################
if (enrich_method == "ORA"){
  
  #Create ORA specific directory
  output = file.path(out_dir, "ORA")
  if (!dir.exists(output)) {
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    print(glue("Creating Directory: {output}"))
  }

  web_results = WebGestaltR(enrichMethod = enrich_method,
              organism = organism,
              enrichDatabaseFile = enrich_db_file,   #this is our curated database GMT file  #currently a symbols file
              enrichDatabaseType = enrich_db_type,    #entrezgene, #kegg  #genesymbol  #entrezgene_protein-coding
              interestGeneFile = goi_file,
              interestGeneType = goi_type,  #using symbols.gmt
              referenceGeneFile =  reference_gene_file,    #"genome",   # n
              referenceGeneType = reference_type, #TYPE is not needed if 'organisms' is others.
              minNum = minNum,
              maxNum = maxNum,
               sigMethod = "fdr",    # fdr or 'top'  #default=fdr
               fdrMethod = "BH",   #default=BH
               fdrThr = 0.05, #alpha  #default=0.05
               outputDirectory = output,
               projectName = paste0(goi_lab, "_", project_name),
               listName = goi_lab,
               nThreads = threads,
              ## ORA specific parameters below ##
              useWeightedSetCover = TRUE,  #default=TRUE,
              useAffinityPropagation = FALSE,   #default=FALSE
              usekMedoid = TRUE,   #default=TRUE
              kMedoid_k = 25       #default=25
              )
  print("DONE")
  
} else if (enrich_method == "GSEA"){
  
  #Create GSA specific directory
  output = file.path(out_dir, "GSEA")
  if (!dir.exists(output)) {
    dir.create(output, recursive = TRUE)
    print(glue("Creating Directory: {output}"))
  }

  #run webgestalt analysis
  WebGestaltR(enrichMethod = enrich_method,
              organism = organism,
              enrichDatabaseFile = enrich_db_file,   #this is our curated database GMT file  #currently a symbols file
              enrichDatabaseType = enrich_db_type,  #entrezgene, #kegg  #genesymbol  #entrezgene_protein-coding
              interestGeneFile = goi_file,        #must be 2 columns for GSEA
              interestGeneType = goi_type,  #using symbols.gmt
              collapseMethod = "median", #could be mean, min, max, median, etc. 
              minNum = minNum,
              maxNum = maxNum,
              sigMethod = "fdr",   # fdr or 'top'
              fdrMethod = "BH",
              fdrThr = 0.05, #alpha
             # useWeightedSetCover = T,
              outputDirectory = output,
              projectName = paste0(goi_lab, "_", project_name),     #project_name
              listName = goi_lab,
              nThreads = threads,
              #GSEA specific
              perNum = perNum,     #number of permutations; default=1000
              gseaP = 1,          #exponential scaling factor of phenotype score, default=1   #LEARN MORE
              saveRawGseaResult = TRUE,
              gseaPlotFormat = "png"   #or 'svg'
  )
  print("DONE")
  cat("\n\n")
  
} else if (enrich_method == "NTA"){
  stop("'NTA' is not supported")
}



