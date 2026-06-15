## Plotting gsa results

library(tidyverse)
library(knitr)
library(kableExtra)


df = read_csv("gsa_count_long.csv", show_col_types = F)

## FIX GEO Accessions ##
df = df %>% 
  mutate(geo = case_when(geo == "geo_gse118126"  ~ "PRJNA484520",      #GEO
                         geo == "geo_gse232243_24s" ~ "PRJNA971477",   #24s
                         geo == "prjna950485_out" ~ "PRJNA950485"))   #bio


## function to pull out relevant results from count dataframe
filter_dataframe = function(method_list, coef_list, gmt_list, slice_list = NA, splice = TRUE){
  
  if (splice){
  df_filtered = df  %>% 
    filter(method %in% method_list)  %>% 
    filter(coef %in% coef_list)  %>% 
    filter(gmt %in% gmt_list)  %>% 
    group_by(geo, coef, gmt, method)  %>% 
    slice(slice_list) 
  } else {
    df_filtered = df  %>% 
      filter(method %in% method_list)  %>% 
      filter(coef %in% coef_list)  %>% 
      filter(gmt %in% gmt_list)  %>% 
      group_by(geo, coef, gmt, method)
  }
  
  return(df_filtered)
  
}


df_filtered = filter_dataframe(method_list = c("GSEA"),
                               coef_list = df$coef %>% unique(),
                               gmt_list = c("original_gmt", "reconstructed_gmt", "direct_prompt_gmt"),
                               splice = FALSE)


## PIVOT LONGER ##  --> adjust some of the naming for better plotting
df_filtered_longer = df_filtered  %>% 
  pivot_longer(cols = c(ap_count, wsc_count, wsc_coverage, fdr_count),
               names_to = "count_metric",
               values_to = "value"
  ) 
df_filtered_longer = df_filtered_longer %>% 
  mutate(value = case_when(is.na(value) ~ 0,
                           TRUE ~ value),
         redundancy_reduction = case_when(count_metric == "ap_count" ~ "AP",
                                          count_metric == "wsc_count" ~ "WSC",
                                          count_metric == "fdr_count" ~ "FDR",
                                          count_metric == "wsc_coverage" ~ "WSC Coverage"),
         gmt = case_when(gmt == "direct_prompt_gmt" ~ "Direct Prompt",
                         gmt == "original_gmt" ~ "Original",
                         gmt == "reconstructed_gmt" ~ "Reconstructed")
         )

df_filtered_longer$gmt = factor(df_filtered_longer$gmt, 
                                levels = c("Original", "Reconstructed", "Direct Prompt"))
df_filtered_longer$redundancy_reduction = factor(df_filtered_longer$redundancy_reduction,
                                                 levels = c("FDR", "AP", "WSC"))

#### Making Plots

### facetted by GEO
all_p = df_filtered_longer %>% 
  select(geo, gmt, method, coef, redundancy_reduction, value) %>% 
  filter(redundancy_reduction %in% c("AP", "WSC", "FDR"))  %>%
  ggplot(aes(x = redundancy_reduction, y = value, fill = gmt)) +
  geom_boxplot(outlier.shape = NA) +    #so outliers aren't plotted twice
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75, seed = 10),
    size = 1.5, alpha = 0.6
  ) +
  facet_wrap(~ geo) +
  labs(x = "Redundancy Reduction", y = "Gene Set Count", fill = "Annotation",
       #title = "Counting Significantly Enriched Gene Sets Reported for each Gene Set Enrichment Analysis"
       ) +
  scale_fill_brewer(palette = "Blues") +
  theme(
        axis.title.x = element_text(margin = margin(t = 10))
        )

print(all_p)

ggsave(filename = "plots/all_count_boxplot.pdf", 
       plot = all_p, 
       scale = 3,
       units = "cm",
       dpi = 600,
       width = 8,
       height = 6
       )


#### Getting a table of the Average Counts (Averaged Between Coefficients in the Design)
all_average = df_filtered_longer %>% 
  select(geo, gmt, method, coef, redundancy_reduction, value) %>% 
  group_by(geo, gmt, method, redundancy_reduction) %>% 
  reframe(average = mean(value),
          median = median(value))

all_average %>% na.omit() %>% view()
