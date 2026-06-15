##### Parse gsa summary JSON results ###
#get the counts for the shared gene sets, added gene sets, and removed gene sets


######### LIBRARIES #########   
library(tidyverse)
library(jsonlite)
library(glue)

############################

df = read_csv("gsa_results/gsa_count_long.csv", show_col_types = FALSE)

######################################
top_dir_list = c("prjna950485_out", "geo_gse118126", "geo_gse232243_24s")
gmt_types_list = c("original_gmt", "reconstructed_gmt", "direct_prompt_gmt")
gsa_method = c("GSEA", "ORA")
coef_list = unique(df$coef)
project_list = str_extract(df$project  %>% unique(), pattern = "Project_.*refseq_gene")
#########################################

data = read_json("gsa_results/geneset_names.json")

## Summarize Results And return as a row in a tibble
compare_results = function(result_list){

    ## COMPARE FOR WSC
    og_wcs = result_list$original$wcs
    rc_wcs = result_list$reconstructed$wcs
    # print(og_wcs)
    # cat("\n")
    # print(rc_wcs)
    # cat("\n")
    shared_wcs = length(intersect(og_wcs, rc_wcs))
    removed_wcs = length(setdiff(og_wcs, rc_wcs))
    added_wcs = length(setdiff(rc_wcs, og_wcs))
    union_wsc = length(union(og_wcs, rc_wcs))

    # COMPSRE FOR AP
    og_ap = result_list$original$ap
    rc_ap = result_list$reconstructed$ap
    union_ap = length(union(og_ap, rc_ap))
    shared_ap = length(intersect(og_ap, rc_ap))
    removed_ap = length(setdiff(og_ap, rc_ap))
    added_ap = length(setdiff(rc_ap, og_ap))
    
    #COMPARE FOR FDR
    og_fdr = result_list$original$fdr
    rc_fdr  = result_list$reconstructed$fdr
    shared_fdr = length(intersect(og_fdr, rc_fdr))
    removed_fdr = length(setdiff(og_fdr, rc_fdr))
    added_fdr = length(setdiff(rc_fdr, og_fdr))
    union_fdr = length(union(og_fdr, rc_fdr))

    #Join up Into one row
    new_row = tibble("shared_wsc" = shared_wcs,
                     "shared_ap" = shared_ap,
                     "shared_fdr" = shared_fdr,
                     "added_wsc" = added_wcs,
                     "added_ap" = added_ap,
                     "added_fdr" = added_fdr,
                     "removed_wsc" = removed_wcs,
                     "removed_ap" = removed_ap,
                     "removed_fdr" = removed_fdr,
                     "frac_wsc" = shared_wcs/union_wsc,
                     "frac_ap" = shared_ap/union_ap,
                     "frac_fdr" =  shared_fdr/union_fdr,
                     "n_og_wsc" = length(og_wcs),
                     "n_og_ap" = length(og_ap),
                     "n_og_fdr" = length(og_fdr),
                     "n_rc_wsc" = length(rc_wcs),
                     "n_rc_ap" = length(rc_ap),
                     "n_rc_fdr" = length(rc_fdr),
                     "union_ap" = union_ap,
                     "union_wsc" = union_wsc,
                     "union_fdr" = union_fdr
                     )

    return(new_row)
}


### Extract Related INFO FROM JSONS, COMPARING GMTs
filter_json = function(json_data, project, method){

    meta_data = list()
    compare_results = list()
    for (i in 1:length(json_data)){

        proj_shortened = str_extract(json_data[[i]]$project, pattern = "Project_.*refseq_gene")
        if ((proj_shortened == project) &&
          (json_data[[i]]$method == method)
         ){
            ## For meta data
            meta_data = list("geo" = json_data[[i]]$geo,
                             "method" = method,
                             "coef" = json_data[[i]]$coef,
                             "project" = project
                             )

            ## For actual results
            fdr_gs = json_data[[i]]$fdr_genesets  %>% unlist()
            ap_gs  = json_data[[i]]$ap_genesets  %>% unlist()  %>% trimws()
            wsc_gs = json_data[[i]]$wsc_genesets  %>% unlist()

            if (json_data[[i]]$gmt == "original_gmt"){
                compare_results[["original"]] = list("fdr" = fdr_gs, "ap" = ap_gs, "wcs" = wsc_gs)
            } else if (json_data[[i]]$gmt == "reconstructed_gmt"){
                compare_results[["reconstructed"]] = list("fdr" = fdr_gs, "ap" = ap_gs, "wcs" = wsc_gs)
            } else if (json_data[[i]]$gmt == "direct_prompt_gmt") {
                compare_results[["direct_prompt"]] = list("fdr" = fdr_gs, "ap" = ap_gs, "wcs" = wsc_gs)
            }
         }
    }
    return(list("meta_data" = meta_data, "results" = compare_results))
}

fix_vector_length = function(og_vec, rc_vec, dp_vec){

    og_vec = if (is.null(og_vec)) NA else og_vec
    rc_vec = if(is.null(rc_vec)) NA else rc_vec
    dp_vec = if(is.null(dp_vec)) NA else dp_vec


    # Find the max length
    max_len <- max(length(og_vec), length(rc_vec), length(dp_vec))

    # Pad each vector with NA to match max length
    length(og_vec) <- max_len
    length(rc_vec) <- max_len
    length(dp_vec) <- max_len

    df <- tibble(og = og_vec, rc = rc_vec, dp = dp_vec)
    return(df)
}


write_full_comparison_files = function(result_list, meta_list){

    base_filename = glue("{meta_list$geo}_{meta_list$method}_{paste(unlist(meta_list$coef), collapse='_')}_{meta_list$project}")
    ## WCS file
    og_wcs = result_list$original$wcs  %>% unlist()
    rc_wcs = result_list$reconstructed$wcs  %>% unlist()
    dp_wcs = result_list$direct_prompt$wcs  %>% unlist()

    wcs_df = fix_vector_length(og_vec = og_wcs, rc_vec = rc_wcs, dp_vec = dp_wcs)
    dir.create("gsa_results/wcs", showWarnings = FALSE) 
    filename_wcs = file.path("gsa_results", "wcs", glue(base_filename, "_wcs.tsv"))
    write_tsv(x = wcs_df, file = filename_wcs)

    #AP File
    og_ap = result_list$original$ap
    rc_ap = result_list$reconstructed$ap
    dp_ap = result_list$direct_prompt$ap
    ap_df = fix_vector_length(og_vec = og_ap, rc_vec = rc_ap, dp_vec = dp_ap)
    dir.create("gsa_results/ap", showWarnings = FALSE) 
    filename_ap = file.path("gsa_results", "ap", glue(base_filename, "_ap.tsv"))
    write_tsv(x = ap_df, file = filename_ap)

    og_fdr = result_list$original$fdr
    rc_fdr = result_list$reconstructed$fdr
    dp_fdr = result_list$direct_prompt$fdr
    fdr_df = fix_vector_length(og_vec = og_fdr, rc_vec = rc_fdr, dp_vec = dp_fdr)
    dir.create("gsa_results/fdr", showWarnings = FALSE) 
    filename_fdr = file.path("gsa_results", "fdr", glue(base_filename, "_fdr.tsv"))
    write_tsv(x = fdr_df, file = filename_fdr)

    cat("File written:  ", filename_wcs, "\n")
    cat("File written:  ", filename_ap, "\n")
    cat("File written:  ", filename_fdr, "\n")
}



all_rows = list()
row_num = 1
for (p in project_list){
    for (m in gsa_method){
        res = filter_json(json_data = data, project = p, method = m)
        meta = res$meta_data
        results = res$results

        ## NOW DETERMINE HOW THE GENESETS HAVE CHANGED Numerically
        new_row = compare_results(results)
        new_row$method = m
        new_row$geo = meta$geo
        new_row$project = p
        new_row$coef = paste(unlist(meta$coef), collapse = "_")
        all_rows[[row_num]] = new_row
        row_num = row_num + 1

        #WRITE FILES FOR COMPARISON
        write_full_comparison_files(results, meta)
    }
}

summary_df = bind_rows(all_rows)  %>% 
    select(geo, method, coef, project, n_og_wsc, n_og_ap, n_og_fdr, n_rc_wsc, n_rc_ap, n_rc_fdr,
            union_wsc, union_ap, union_fdr, shared_wsc, shared_ap, shared_fdr,
            added_wsc, added_ap, added_fdr, removed_wsc, removed_ap, removed_fdr,
            frac_wsc, frac_ap, frac_fdr)  %>% 
    arrange(geo, method, coef, project)
view(summary_df)

write_csv(summary_df, file = "gsa_results/gs_compare_summary.csv")

sample(1:12, size = 2, replace = F)
