## parse GSA results and make a summary csv of the gene set count information (significant gene set results for each project id)

###### Libraries ########
library(tidyverse)
library(glue)
#########################


############## SOME GLOBAL VARIABLES (DIRECTORY PATHS + FILENAMES) TO ACCESS DATA ###########
top_dir_list = c("prjna950485_out", "geo_gse118126", "geo_gse232243_24s")
webgestalt_dir = "webgestalt_out"
gmt_types_list = c("original_gmt", "reconstructed_gmt", "direct_prompt_gmt")
gsa_method = c("GSEA", "ORA")
coef1_only_pat = "Project_coef1_(?!.*coef[1-9]+.*)"   #comparing coef1 only
coef2_only_pat = "Project_coef2_(?!.*coef[1-9]+.*)"   #comparing coef2 only
coef3_only_pat = "Project_coef3_(?!.*coef[1-9]+.*)"   #comparing coef3 only
coef4_only_pat = "Project_coef3_(?!.*coef[1-9]+.*)"   #comparing coef4 only

coef1_coef2_pat = "Project_coef1.*coef2.*"            #comparing coef1 and coef2
coef1_coef3_pat = "Project_coef1.*coef3.*"            #comparing coef1 and coef2
coef1_coef4_pat = "Project_coef1.*coef4.*"            #comparing coef1 and coef4
coef2_coef3_pat = "Project_coef2.*coef3.*"
coef2_coef4_pat = "Project_coef2.*coef4.*"
coef3_coef4_pat = "Project_coef3.*coef4.*"

project_patterns = c("coef1" = coef1_only_pat, 
                     "coef2" = coef2_only_pat, 
                     "coef3" = coef3_only_pat,
                     "coef4" = coef4_only_pat,
                     "coef1_coef2" = coef1_coef2_pat,
                     "coef1_coef3" = coef1_coef3_pat,
                     "coef1_coef4" = coef1_coef4_pat,
                     "coef2_coef3" = coef2_coef3_pat,
                     "coef2_coef4" = coef2_coef4_pat,
                     "coef3_coef4" = coef3_coef4_pat
                     )

affinity_propagation_name = "enriched_geneset_ap_clusters_.*.txt"
weighted_set_cover_name = "enriched_geneset_wsc_topsets_.*.txt"
results_tab_name = "enrichment_results_.*.txt"
report_name = "Report_.*.html"
###########################################################################################


## Going to loop through these directories to collect some initial results

#First setting up a tibble that I can add count data too.
gene_set_data_list = list()


#1. Loop through GEO Sets (geo):
proj_counter = 0
for (geo in top_dir_list){

    #2. Loop Through GMT Types (gmt)
    for (gmt in gmt_types_list){

        #3. Loop Through GSA Method (method)
        for (method in gsa_method){
            project_path = file.path(geo, webgestalt_dir, gmt, method)

            #4. Loop Through Projects (pat)
            for (pat in project_patterns){
                project_list = list.files(project_path)
                pat_name = names(project_patterns[project_patterns == pat])
                proj_dir = project_list[grepl(pattern = pat, x = project_list, perl = TRUE)]

                if (length(proj_dir) == 0){
                    warning(glue("Project with coefficient: {pat_name}, not present for directory: {project_path}"))
                }
                #5. Maybe there are multiple coef1 files, etc. (p)
                for (p in proj_dir){

                    results_path = file.path(project_path, p)

                    #6. Now read in the info from the files in this path
                    results_list = list.files(results_path)

                    #.1 affinity propagation data
                    ap_file = results_list[grepl(pattern = affinity_propagation_name,
                                                 x = results_list, 
                                                perl = TRUE
                                                )]
                    if (length(ap_file) == 0){
                        print("No Significant AP Results")
                        ap_count = NA
                    } else {
                        ap_data = readLines(file.path(results_path, ap_file))
                        ap_count = length(ap_data)
                    }
                    
                    #.2 weighted set cover data
                    wsc_file = results_list[grepl(pattern = weighted_set_cover_name,
                                                 x = results_list, 
                                                perl = TRUE
                                                )]
                    if (length(wsc_file) == 0){
                        print("No Significant WSC Results")
                        wsc_count = NA
                        wsc_coverage = NA
                    } else {
                        wsc_data = readLines(file.path(results_path, wsc_file))
                        wsc_count = length(wsc_data)-1  #header is the coverage score
                        wsc_coverage = str_extract(wsc_data[1], pattern="\\d+\\.\\d+")
                    }
                    #.3 results table data
                    results_file = results_list[grepl(pattern = results_tab_name,
                                                 x = results_list, 
                                                perl = TRUE
                                                )]
                    if (length(results_file) == 0){
                        print("No Significant FDR Results")
                        fdr_count = NA
                    } else {
                        results_data = read_csv(file.path(results_path, results_file), show_col_types = F)
                        fdr_count = nrow(results_data)
                    }

                    #.4 html table data
                    html_report = results_list[grepl(pattern = report_name,
                                 x = results_list, 
                                 perl = TRUE
                                )]
                    if (length(html_report) == 0){
                        print("No html report generated")
                        html_link = NA
                    } else {
                        html_link = file.path("http://localhost:52330", results_path, html_report)
                    }

                    ## Maybe open up html file and take a look?
                    proj_counter = proj_counter+1   #Number of Projects We are Loop

                    7. #### Append new row/data to tibble. #######
                    new_row = tibble("geo" = geo,
                         "gmt" = gmt,
                         "method" = method,
                         "coef" = pat_name,
                         "project" = p,
                         "ap_count" = ap_count,
                         "wsc_count" = wsc_count,
                         "wsc_coverage" = wsc_coverage,
                         "fdr_count" = fdr_count,
                         "html_link" = html_link
                        )

                    gene_set_data_list[[proj_counter]] = new_row
                }
            }
        }
    }
}
proj_counter


gene_set_df = bind_rows(gene_set_data_list)  %>% arrange(geo, method, gmt, coef) 
#gene_set_df  %>% view()

write_csv(x = gene_set_df,
          file = "gsa_results/gsa_count_long.csv"
          )

#view(gene_set_df)


wide_gene_set_df = gene_set_df  %>% 
   # filter(method == "GSEA")  %>% 
    pivot_wider(names_from = gmt,
                 values_from = c(ap_count, wsc_count, fdr_count, wsc_coverage))  %>% 
    arrange(geo, method, coef) 

#wide_gene_set_df  %>% view()

write_csv(x = wide_gene_set_df,
          file = "gsa_results/gsa_count_wide.csv"
          )

