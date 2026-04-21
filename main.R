source("functions.R")
source("functions_analysis.R")

# experiment name: experiment_group (DIA1-3, DDA3-4), MS, specification (DIA: ups2 or ups3, DDA: diann or mq)
# DIA1_popeye, DIA2_olive, DIA2_oliveslice, DIA2_fozzie, DIA2_sam
# DIA3_olive_ups2, DIA3_olive_ups3, DIA3_oliveslice_ups2, DIA3_oliveslice_ups3
# DIA3_fozzie_ups2, DIA3_fozzie_ups3, DIA3_sam_ups2, DIA3_sam_ups3
# DDA3_fozzie_diann, DDA3_sam_diann, DDA3_fozzie_mq, DDA3_sam_mq
# DDA4_fozzie_diann, DDA4_sam_diann, DDA4_fozzie_mq, DDA4_sam_mq

{
experiment <- "DIA3_fozzie_ups2"
  
source("load_data.R")
source("filter_overview.R")
# save stuff function 1
print(attr(data, "anz_plot"))
ggsave_png()
save_table(attr(data, "anz"), "anz")
save_table(attr(data, "anz_condition"), "anz_condition")
save_table(attr(data, "ups_protein_info"), "ups_protein_info")
}

score_name <- "iBAQ"
# score_name <- "topN"
# score_name <- "mean_int"
{
data <- data %>% 
  { if (standard == "UPS2") {
      ups2_score(., mean_score = "no")       # apply function to the current data frame
    } else {.}} %>%    # otherwise return df without doing anything
  { if (standard == "UPS3") {
      ups3_score(.)       # apply function to the current data frame
  } else {.}} %>% 
  u2os_score() %>%  # function 3
  u2os_absolut()   # function 4
 
# safe results 
# function 2
{ if (standard == "UPS2") {
  save_table(attr(data, "ups2_score"), paste0("ups2_", score_name, "_scores"))
  save_table(attr(data, "ups2_score_info"), paste0("ups2_", score_name, "_slope_info"))
  print(attr(data, "ups2_score_plot"))   
  ggsave_png()
} else if (standard == "UPS3") {
  save_table(attr(data, "ups3_score"), paste0("ups3_", score_name))
  save_table(attr(data, "ups3_score_info"), paste0("ups3_", score_name, "_slope_info"))
  print(attr(data, "ups3_score_plot"))  
  ggsave_png()
  }}
# function 3
print(attr(data, "u2os_score_plot"))
ggsave_png()
# function 4
print(attr(data, "u2os_absolut_fmol_plot"))
ggsave_png()
save_table(attr(data, "u2os_iqr_fmol"), paste0("u2os_iqr_fmol_", score_name))
save_table(attr(data, "histone_iqr_fmol"), paste0("histone_iqr_fmol_", score_name))
save_table(attr(data, "anz_condition_filtered"), paste0("anz_condition_filtered_", score_name))
print(attr(data, "u2os_absolut_pg_plot"))
ggsave_png()
save_table(attr(data, "u2os_iqr_pg"), paste0("u2os_iqr_pg_", score_name))
save_table(attr(data, "histone_iqr_pg"), paste0("histone_iqr_pg_", score_name))
}
# copynumber_estimate <- attr(data, "u2os_absolut_pg")
# save_table(attr(data, "u2os_absolut_pg"), paste0("u2os_absolut_pg_", score_name))


# plots for the BM
data <- data %>% 
  ups2_peptides() %>%   # additional BM plots
  ups2_score() %>%   # part 2
  u2os_score()   # part3

print(attr(data, "peptides_fc_plot"))
ggsave_png()
print(attr(data, "peptides_int_plot"))
ggsave_png()
print(attr(data, "ups2_proteins_plot"))
ggsave_png()
save_table(attr(data, "ups2_score"), "ups2_score")
print(attr(data, "ups2_score_plot"))
ggsave_png()
# function 3
print(attr(data, "u2os_score_plot"))
ggsave_png()



# function 2b
data_pg <- load_data(str_c(path, "/txt/proteinGroups.txt")) %>% 
  filter(!grepl("CON_", `Majority protein IDs`)) %>% 
  select(-(starts_with("Intensity") | starts_with("Peptides") | starts_with("Razor") | starts_with("Unique peptides") | starts_with("iBAQ") | starts_with("Top3") | starts_with("LFQ") | starts_with("MS/MS count")| starts_with("Sequence coverage")) | contains("UPS3")) %>%   # use only columns with UPS3
  ups3_MQ_score(score = "Top3")
# ups3_MQ_score(score = "iBAQ")
print(attr(data_pg, "ups3_MQ_score_plot"))
ggsave_png()

ups3_factor_pmol <- if (any(grepl("fozzie", path))) {95.91} else if (any(grepl("sam", path))) {239.77} else {NULL}
ups_factor <- if (any(grepl("fozzie", path))) {0.09591} else if (any(grepl("sam", path))) {0.23977} else if (any(grepl("olive", path))) {0.23977} else {NULL}
histone_theoretical <- if (any(grepl("fozzie", path))) {333} else if (any(grepl("sam", path))) {133} else if (any(grepl("olive", path))) {133} else {NULL}

