
source("functions.R")
source("functions_analysis.R")

# experiment name: experiment_group (DIA1-3, DDA3-4), MS, specification (DIA: ups2 or ups3, DDA: diann or mq)
# DIA1_popeye, DIA2_olive, DIA2_oliveslice, DIA2_fozzie, DIA2_sam
# DIA3_olive_ups2, DIA3_olive_ups3, DIA3_oliveslice_ups2, DIA3_oliveslice_ups3
# DIA3_fozzie_ups2, DIA3_fozzie_ups3, DIA3_sam_ups2, DIA3_sam_ups3
# DDA3_fozzie_diann, DDA3_sam_diann, DDA3_fozzie_mq, DDA3_sam_mq
# DDA4_fozzie_diann, DDA4_sam_diann, DDA4_fozzie_mq, DDA4_sam_mq

{
experiment <- "DIA3_olive_ups2"
  
source("load_data.R")
source("filter_overview.R")
# save results function 1
print(attr(data, "anz_plot"))
ggsave_png()
save_table(attr(data, "anz"), "anz")
save_table(attr(data, "anz_condition"), "anz_condition")
save_table(attr(data, "ups_protein_info"), "ups_protein_info")
}

score_name <- "iBAQ"

{
data <- data %>% 
  { if (standard == "UPS2") {
      ups2_score(., mean_score = "no")       # apply function to the current data frame
    } else {.}} %>%    # otherwise return df without doing anything
  { if (standard == "UPS3") {
      ups3_score(., mean_score = "no")       # apply function to the current data frame
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
save_table(attr(data, "core_histone_sum_fmol"), paste0("core_histone_sum_fmol_", score_name))
save_table(attr(data, "core_histone_sum_pg"), paste0("core_histone_sum_pg_", score_name))
}
# copynumber_estimate <- attr(data, "u2os_absolut_pg")
# save_table(attr(data, "u2os_absolut_pg"), paste0("u2os_absolut_pg_", score_name))

# plots for the sample BM
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
