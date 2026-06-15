## Parse Enrichment rsults FDR, AP, and WSC and make a summary json of the gene sets and project names.
## This creates a JSON FILE with all the results

library(tidyverse)
library(glue)
library(jsonlite)

############## GLOBAL VARIABLES ##############
top_dir_list = c("prjna950485_out", "geo_gse118126", "geo_gse232243_24s")
webgestalt_dir = "webgestalt_out"
gmt_types_list = c("original_gmt", "reconstructed_gmt", "direct_prompt_gmt")
gsa_method = c("GSEA", "ORA")

coef1_only_pat = "Project_coef1_(?!.*coef[1-9]+.*)"
coef2_only_pat = "Project_coef2_(?!.*coef[1-9]+.*)"
coef3_only_pat = "Project_coef3_(?!.*coef[1-9]+.*)"
coef4_only_pat = "Project_coef3_(?!.*coef[1-9]+.*)"
coef1_coef2_pat = "Project_coef1.*coef2.*"
coef1_coef3_pat = "Project_coef1.*coef3.*"
coef1_coef4_pat = "Project_coef1.*coef4.*"
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
                     "coef3_coef4" = coef3_coef4_pat)

results_tab_name = "enrichment_results_.*.txt"
affinity_propagation_name = "enriched_geneset_ap_clusters_.*.txt"
weighted_set_cover_name = "enriched_geneset_wsc_topsets_.*.txt"
##############################################




gene_set_names = list()
proj_counter = 0

for (geo in top_dir_list){
    for (gmt in gmt_types_list){
        for (method in gsa_method){
            project_path = file.path(geo, webgestalt_dir, gmt, method)

            for (pat in project_patterns){
                project_list = list.files(project_path)
                pat_name = names(project_patterns[project_patterns == pat])
                proj_dir = project_list[grepl(pattern = pat, x = project_list, perl = TRUE)]

                if (length(proj_dir) == 0){
                    warning(glue("Project with coefficient: {pat_name}, not present for directory: {project_path}"))
                    next
                }

                for (p in proj_dir){
                    proj_counter = proj_counter + 1
                    results_path = file.path(project_path, p)
                    results_list = list.files(results_path)

                    key = glue("{geo}__{gmt}__{method}__{pat_name}__{p}")

                    entry = list(
                        geo     = geo,
                        gmt     = gmt,
                        method  = method,
                        coef    = pat_name,
                        project = p
                    )

                    # FDR gene sets
                    results_file = results_list[grepl(pattern = results_tab_name,
                                                      x = results_list, perl = TRUE)]
                    if (length(results_file) > 0){
                        results_data = read_tsv(file.path(results_path, results_file),
                                                show_col_types = FALSE)
                        entry$fdr_genesets = results_data$geneSet
                    } else {
                        entry$fdr_genesets = list()
                    }

                    # AP gene sets
                    ap_file = results_list[grepl(pattern = affinity_propagation_name,
                                                 x = results_list, perl = TRUE)]
                    if (length(ap_file) > 0){
                        ap_data = read_tsv(file.path(results_path, ap_file),
                                           show_col_types = FALSE, col_names = FALSE)
                        entry$ap_genesets = ap_data$X1
                    } else {
                        entry$ap_genesets = list()
                    }

                    # WSC gene sets
                    wsc_file = results_list[grepl(pattern = weighted_set_cover_name,
                                                  x = results_list, perl = TRUE)]
                    if (length(wsc_file) > 0){
                        wsc_data = readLines(file.path(results_path, wsc_file))
                        entry$wsc_genesets = wsc_data[-1]
                    } else {
                        entry$wsc_genesets = list()
                    }

                    gene_set_names[[proj_counter]] = entry
                }
            }
        }
    }
}

write_json(gene_set_names, path = "gsa_results/geneset_names.json", pretty = TRUE, auto_unbox = TRUE)
