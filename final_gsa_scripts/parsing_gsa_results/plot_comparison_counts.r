## ALL of these comparisons are between the original and reconstructed dataframes.
### make plots from gs_compare_summary.csv

#libraries
library(tidyverse)


# fix bio accession names
summary_df = read_csv("gs_compare_summary.csv", show_col_types = F)

summary_df = summary_df %>% 
  mutate(geo = case_when(geo == "geo_gse118126"  ~ "PRJNA484520",      #GEO
                         geo == "geo_gse232243_24s" ~ "PRJNA971477",   #24s
                         geo == "prjna950485_out" ~ "PRJNA950485")) %>%   #bio
  select(-c(frac_ap, frac_wsc, frac_fdr))

######### Plots of shared, added, and removed counts normalized to the original annotation size  counts/n_original
################## ADDING FRACTIONAL AMOUNTS FOR --> comparing n in original
summary_df = summary_df %>% 
  mutate(shared_n_frac_wsc = shared_wsc/n_og_wsc,
         shared_n_frac_ap = shared_ap/n_og_ap,
         shared_n_frac_fdr = shared_fdr/n_og_fdr,
         added_n_frac_wsc = added_wsc/n_og_wsc,
         added_n_frac_ap = added_ap/n_og_ap,
         added_n_frac_fdr = added_fdr/n_og_fdr,
         removed_n_frac_wsc = removed_wsc/n_og_wsc,
         removed_n_frac_ap = removed_ap/n_og_ap,
         removed_n_frac_fdr = removed_fdr/n_og_fdr,
  )

## fixing the labels for better plotting
n_frac_df = summary_df %>% 
  select(c(1:4, grep("_n_frac", names(summary_df))))

n_frac_df_long = n_frac_df %>% 
  pivot_longer(cols = 5:13,
               names_to = "metric",
               values_to = "value")

n_frac_df_long = n_frac_df_long %>% 
  mutate(redundancy_reduction = str_extract(metric, pattern = "(?<=_)[a-z]+$"),
         redundancy_reduction = toupper(redundancy_reduction),
         comparison = str_extract(metric, pattern = "[a-z]+(?=_)")
         )

n_frac_df_long$redundancy_reduction = factor(n_frac_df_long$redundancy_reduction, levels = c("FDR", "AP", "WSC"))


## Make plot
n_frac_counts = n_frac_df_long %>% 
  filter(method == "GSEA") %>% 
  ggplot(aes(x = redundancy_reduction, y = value, fill = comparison)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75, seed = 10),
    size = 1.5, alpha = 0.5
  ) +
  facet_wrap(~geo) +
  scale_fill_viridis_d() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Redundancy Reduction", 
       y = "Gene Set Count / Original Annotation Gene Set Count", fill = "Comparison") +
  #ggtitle("Comparing Enriched Gene Sets Between Original GO Annotation and LLM Reconstructed Annotation for GSEA") +
  theme(axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        )
print(n_frac_counts)

ggsave(plot = n_frac_counts,
       filename = "plots/n_frac_gs_comparison_counts.pdf", 
       dpi = 600,
       width = 8, 
       height = 6, 
       units = "cm",
       scale = 3.5
       )


### A summary table of average counts
n_frac_summary_table = n_frac_df_long %>% 
  filter(method == "GSEA") %>% 
  group_by(geo, method, comparison, redundancy_reduction) %>% 
  reframe(n = n(),
            median = median(value),
            mean = mean(value))

n_frac_summary_table %>% view()


############################################## ADDITIONAL PLOTS NOT USED IN PUBLICATION ###########################
#n_frac_summary_table %>% 
#  filter(comparison == "shared") %>% 
#  pull(median) %>% 
#  min()
#
#n_frac_summary_table %>% 
#  filter(comparison == "shared") %>% 
#  pull(median) %>% 
#  max()
#
#
#n_frac_summary_table %>% 
#  filter(comparison == "shared") %>% 
#  pull(median) %>% 
#  median()

#view(summary_df)




############################# PLOTS BEFORE FRACTIONAL AMOUNTS #################
#long_df = summary_df %>% 
#    pivot_longer(cols = names(summary_df)[grep("wsc|ap|fdr", names(summary_df))],
#                 names_to = "metric",
#                 values_to = "values"
#    ) %>% 
#  mutate(comparison = str_extract(metric, pattern = "^[^_]+"),
#         redundancy_reduction = str_extract(metric, pattern = "(?<=_).*$"))
#
#long_df$comparison = factor(long_df$comparison, levels = c("added", "removed", "shared"))
##view(long_df)
#
#
#### CREATE PLOT FOR GSA
#gsea_just_counts = long_df %>% 
#  filter(method == "GSEA") %>% 
#  filter(comparison %in% c("added", "removed", "shared")) %>% 
#  ggplot(aes(x = redundancy_reduction, y = values, fill = comparison)) +
#  geom_boxplot() +
#  #geom_jitter() +
#  facet_wrap(~geo) +
#  scale_fill_viridis_d() +
#  labs(x = "Redundancy Reduction", y = "Gene Set Count", fill = "Comparison",
#       title = "Comparing Enriched Gene Sets Between Original GO Database and LLM Reconstructed Database for GSEA") +
#  theme(axis.title.x = element_text(margin = margin(t = 10)),
#        axis.title.y = element_text(margin = margin(r = 10))
#  )
#print(gsea_just_counts)
#
#ggsave(plot = gsea_just_counts,
#       filename = "plots/gsea_gs_comparison_just_counts.pdf", 
#       #dpi = 600,
#       width = 8, 
#       height = 6, 
#       units = "cm",
#       scale = 3.5
#)
#
#
### CREATE PLOT FOR ORA
#ora_just_counts = long_df %>% 
#  filter(method == "ORA") %>% 
#  filter(comparison %in% c("added", "removed", "shared")) %>% 
#  ggplot(aes(x = redundancy_reduction, y = values, fill = comparison)) +
#  geom_boxplot() +
#  #geom_jitter() +
#  scale_fill_viridis_d() +
#  facet_wrap(~geo) +
#  labs(x = "Redundancy Reduction", y = "Gene Set Count", fill = "Comparison",
#       title = "Comparison of Differences in Enriched Gene Sets Between Original Database and Reconstructed Database for ORA")
#print(ora_just_counts)

########################################## PLOTS AFTER FRACTIONAL AMOUNTS #####################################################3

############## ADDING FRACTIONAL AMOUNTS -- > COMPARING to union amount
#summary_df = summary_df %>% 
#    mutate(int_union_frac_wsc = shared_wsc/union_wsc,
#           int_union_frac_ap = shared_ap/union_ap,
#           int_union_frac_fdr = shared_fdr/union_fdr,
#           added_union_frac_wsc = added_wsc/union_wsc,
#           added_union_frac_ap = added_ap/union_ap,
#           added_union_frac_fdr = added_fdr/union_fdr,
#           removed_union_frac_wsc = removed_wsc/union_wsc,
#           removed_union_frac_ap = removed_ap/union_ap,
#           removed_union_frac_fdr = removed_fdr/union_fdr,
#           )
#
#union_frac_df = summary_df %>% 
#  select(c(1:4, grep("union_frac", names(summary_df))))
#
#union_frac_df_long = union_frac_df %>% 
#  pivot_longer(cols = 5:13,
#               names_to = "metric",
#               values_to = "value")
#
#
#union_frac_df_long = union_frac_df_long %>% 
#  mutate(redundancy_reduction = str_extract(metric, pattern = "(?<=_)[a-z]+$"),
#         comparison = str_extract(metric, pattern = "[a-z]+(?=_)"),
#         comparison = case_when(comparison == "int" ~ "shared",
#                                T ~ comparison)
#         )
#
#union_frac_df_long$redundancy_reduction = factor(union_frac_df_long$redundancy_reduction,
#                                                 levels = c("fdr", "ap", "wsc"))
#### Now creating plots from union fractional amount 
#
### GSEA OF COMPARISON/UNION
#gsea_union_frac_counts = union_frac_df_long %>% 
#  filter(method == "GSEA") %>% 
#  ggplot(aes(x = redundancy_reduction, y = value, fill = comparison)) +
#  geom_boxplot() +
#  #geom_jitter() +
#  facet_wrap(~geo) +
#  scale_fill_viridis_d() +
#  scale_y_continuous(limits = c(0,1)) +
#  labs(x = "Redundancy Reduction",
#       y = "Enriched Gene Set Count / Union Gene Set Count",
#       fill = "Comparison",
#       title = "Comparing Enriched Gene Sets Between Original GO Annotation and LLM Reconstructed Annotation for GSEA")
#print(gsea_union_frac_counts)
#
#ggsave(plot = gsea_union_frac_counts,
#       filename = "plots/gsea_gs_comparison_union_frac.png", 
#       dpi = 600,
#       width = 8, 
#       height = 6, 
#       units = "cm",
#       scale = 3.5
#)


