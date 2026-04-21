# load experiment path corresponding to selected experiment (and list of conditions and label for final plots)
if (experiment == "DIA1_popeye") {
  path_report <- "E:/MaxQuant/ACB/DIA1_popeye/20250619_combinedlibtxt_frag_v182/report.tsv" 
  conditions <- c("5ng_UPS2", "10ng_UPS2", "20ng_UPS2", "200ng_HSH", "HSH_55ng_UPS2", "HSH_110ng_UPS2", "HSH_220ng_UPS2")
  condition_labels <- c("5ng_UPS2", "10ng_UPS2", "20ng_UPS2", "U2OS", "UPS2:U2OS = 1:18", "UPS2:U2OS = 1:9", "UPS2:U2OS = 1:4.5")
  remove_conditions <- c("5ng_UPS2", "10ng_UPS2", "20ng_UPS2")
  condition_filter = "HSH_220ng_UPS2"
  
} else if (experiment == "DIA2_olive"| experiment == "DIA2_olive_ups2") {
  path_report <- "E:/MaxQuant/ACB/DIA2_olive/20251210_olive_dia_s11_v182/report.tsv" 
  conditions <- c("dia_HSH_025ng_UPS2", "dia_HSH_1ng_UPS2", "dia_HSH_5ng_UPS2", "dia_HSH_10ng_UPS2", "dia_HSH", "dia_UPS2")
  condition_labels <- c("0.25 ng U2OS+UPS2", "1 ng U2OS+UPS2", "5 ng U2OS+UPS2", "10 ng U2OS+UPS2", "U2OS", "UPS2")
  remove_conditions <- c("dia_HSH", "dia_UPS2")
  condition_filter = "dia_HSH_10ng_UPS2"
  ups_factor <- 6883
  histone_theoretical <- 5.57
} else if (experiment == "DIA2_oliveslice"| experiment == "DIA2_oliveslice_ups2") {
  path_report <- "E:/MaxQuant/ACB/DIA2_olive/20251210_olive_slice_s11_v182/report.tsv" 
  conditions <- c("slice_HSH_025ng_UPS2", "slice_HSH_1ng_UPS2", "slice_HSH_5ng_UPS2", "slice_HSH_10ng_UPS2", "slice_HSH", "slice_UPS2")
  condition_labels <- c("0.25 ng U2OS+UPS2", "1 ng U2OS+UPS2", "5 ng U2OS+UPS2", "10 ng U2OS+UPS2", "U2OS", "UPS2")
  remove_conditions <- c("slice_HSH", "slice_UPS2")
  condition_filter = "slice_HSH_10ng_UPS2"
  ups_factor <- 6883
  histone_theoretical <- 5.57
} else if (experiment == "DIA2_fozzie"| experiment == "DIA2_fozzie_ups2") {
  path_report <- "E:/MaxQuant/ACB/DIA2_fozzie/20251210_fozzie_s11_v20/report.parquet"
  conditions <- c("HSH_025ng_UPS2", "HSH_1ng_UPS2", "HSH_5ng_UPS2", "HSH_25ng_UPS2", "HSH_50ng_UPS2", "HSH_100ng_UPS2", "HSH", "UPS2")
  condition_labels <- c("0.25 ng U2OS+UPS2", "1 ng U2OS+UPS2", "5 ng U2OS+UPS2", "25 ng U2OS+UPS2", "50 ng U2OS+UPS2", "100 ng U2OS+UPS2", "U2OS", "UPS2")
  remove_conditions <- c("HSH", "UPS2")
  condition_filter = "HSH_100ng_UPS2"
  ups_factor <- 689
  histone_theoretical <- 55.72
  } else if (experiment == "DIA2_sam"| experiment == "DIA2_sam_ups2") {
  path_report <- "E:/MaxQuant/ACB/DIA2_sam/20251210_sam_s11_v20/report.parquet"
  conditions <- c("HSH_025ng_UPS2", "HSH_1ng_UPS2", "HSH_5ng_UPS2", "HSH_25ng_UPS2", "HSH_50ng_UPS2", "HSH_100ng_UPS2", "HSH", "UPS2")
  condition_labels <- c("0.25 ng U2OS+UPS2", "1 ng U2OS+UPS2", "5 ng U2OS+UPS2", "25 ng U2OS+UPS2", "50 ng U2OS+UPS2", "100 ng U2OS+UPS2", "U2OS", "UPS2")
  remove_conditions <- c("HSH", "UPS2")
  condition_filter = "HSH_100ng_UPS2"
  ups_factor <- 689
  histone_theoretical <- 55.72
  
  } else if (experiment == "DIA3_olive_ups2"| experiment == "DIA3_olive") {
  path_report <- "E:/MaxQuant/ACB/DIA3_olive/20251124_olive_v18_ups2/report.tsv"
  conditions <- c("dia_10ng_U2OS_UPS2", "dia_50ng_U2OS_UPS2")
  condition_labels <- c("10 ng U2OS+UPS2", "50 ng U2OS+UPS2")
  condition_filter = "dia_10ng_U2OS_UPS2"
  ups_factor <- 5782
  histone_theoretical <- 5.39
} else if (experiment == "DIA3_oliveslice_ups2"| experiment == "DIA3_oliveslice") {
  path_report <- "E:/MaxQuant/ACB/DIA3_olive_sliced/20251124_olive_v18_ups2_slice/report.tsv" 
  conditions <- c("slice_10ng_U2OS_UPS2")
  condition_labels <- c("10 ng U2OS+UPS2")
  condition_filter = "slice_10ng_U2OS_UPS2"
  ups_factor <- 5782
  histone_theoretical <- 5.39
} else if (experiment == "DIA3_fozzie_ups2"| experiment == "DIA3_fozzie") {
  path_report <- "E:/MaxQuant/ACB/DIA3_fozzie/20251124_fozzie_v20_ups2/report.parquet"
  conditions <- c("100ng_U2OS_UPS2", "200ng_U2OS_UPS2")
  condition_labels <- c("100 ng U2OS+UPS2", "200 ng U2OS+UPS2")
  # remove_conditions <- c("100ng_U2OS_UPS2")
  condition_filter = "200ng_U2OS_UPS2"
  ups_factor <- 289
  histone_theoretical <- 107.77
} else if (experiment == "DIA3_sam_ups2"| experiment == "DIA3_sam") {
  path_report <- "E:/MaxQuant/ACB/DIA3_sam/20251124_sam_v20_ups2/report.parquet"
  conditions <- c("U2OS_UPS2")
  condition_labels <- c("100 ng U2OS+UPS2")
  condition_filter = "U2OS_UPS2"
  ups_factor <- 578
  histone_theoretical <- 53.89
  
} else if (experiment == "DIA3_olive_ups3") {
  path_report <- "E:/MaxQuant/ACB/DIA3_olive/20260324_olive_v18_ups3/report.tsv" 
  conditions <- c("dia_10ng_U2OS_UPS3", "dia_50ng_U2OS_UPS3")
  condition_labels <- c("10 ng U2OS+UPS3", "50 ng U2OS+UPS3")
  condition_filter = "dia_10ng_U2OS_UPS3"
  ups_factor <- 4.8
  histone_theoretical <- 5.49
  } else if (experiment == "DIA3_oliveslice_ups3") {
  path_report <- "E:/MaxQuant/ACB/DIA3_olive_sliced/20251124_olive_v18_ups3_slice/report.tsv" 
  conditions <- c("slice_10ng_U2OS_UPS3")
  condition_labels <- c("10 ng U2OS+UPS3")
  condition_filter = "slice_10ng_U2OS_UPS3"
  ups_factor <- 4.8
  histone_theoretical <- 5.49
  } else if (experiment == "DIA3_fozzie_ups3") {
  path_report <- "E:/MaxQuant/ACB/DIA3_fozzie/20251124_fozzie_v20_ups3/report.parquet"
  conditions <- c("200ng_U2OS_UPS3")
  condition_labels <- c("200 ng U2OS+UPS3")
  condition_filter = "200ng_U2OS_UPS3"
  ups_factor <- 0.24
  histone_theoretical <- 109.78
  } else if (experiment == "DIA3_sam_ups3") {
  path_report <- "E:/MaxQuant/ACB/DIA3_sam/20251124_sam_v20_ups3/report.parquet"
  conditions <- c("U2OS_UPS3")
  condition_labels <- c("100 ng U2OS+UPS3")
  condition_filter = "U2OS_UPS3"
  ups_factor <- 0.48
  histone_theoretical <- 54.89
  
  
  } else if (experiment == "DDA3_fozzie_diann") {
    path_report <- "E:/MaxQuant/ACB/DDA/test3_fozzie/diann/combined_report.tsv"
    conditions <- c("DDA_500ng_U2OS_UPS2")
    condition_labels <- c("500 ng U2OS+UPS2 DiaNN")
    ups_factor <- 116
    histone_theoretical <- 263.47
  } else if (experiment == "DDA3_sam_diann") {
    path_report <- "E:/MaxQuant/ACB/DDA/test3_sam/diann/combined_report.tsv"
    conditions <- c("160ngHSH_40ngUPS2", "196ngHSH_4ngUPS2")
    condition_labels <- c("200 ng U2OS+UPS2 DiaNN", "200 ng U2OS+UPS3 DiaNN (1:50)")
    remove_conditions <- c("196ngHSH_4ngUPS2")
    ups_factor <- 265
    histone_theoretical <- 105.39
  } else if (experiment == "DDA3_fozzie_mq") {
    path_report <- "E:/MaxQuant/ACB/DDA/test3_fozzie/MQ/combined/txt/peptides.txt"
    conditions <- c("DDA_500ng_U2OS_UPS2")
    condition_labels <- c("500 ng U2OS+UPS2 MQ")
    ups_factor <- 116
    histone_theoretical <- 263.47
  } else if (experiment == "DDA3_sam_mq") {
    path_report <- "E:/MaxQuant/ACB/DDA/test3_sam/MQ/combined/txt/peptides.txt"
    conditions <- c("160ngHSH_40ngUPS2")
    condition_labels <- c("200 ng U2OS+UPS2 MQ")
    ups_factor <- 265
    histone_theoretical <- 105.39
  
} else if (experiment == "DDA4_fozzie_diann") {
  path_report <- "E:/MaxQuant/ACB/DDA/test4_fozzie/20251016_fozzie_diann_U2OS_UPS3/report.parquet"
  conditions <- c("500ng_U2OS_UPS3")
  condition_labels <- c("500 ng U2OS+UPS3 DiaNN")
    ups_factor <- 0.096
  histone_theoretical <- 274.44
  } else if (experiment == "DDA4_sam_diann") {
  path_report <- "E:/MaxQuant/ACB/DDA/test4_sam/20251016_sam_diann_U2OS_UPS3/report.parquet"
  conditions <- c("200ng_U2OS_UPS3")
  condition_labels <- c("200 ng U2OS+UPS3 DiaNN")
  ups_factor <- 0.24
  histone_theoretical <- 109.78
  } else if (experiment == "DDA4_fozzie_mq") {
  path_report <- "E:/MaxQuant/ACB/DDA/test4_fozzie/20251017_fozzie_MQ/txt/peptides.txt"
  conditions <- c("500ng_U2OS_UPS3")
  condition_labels <- c("500 ng U2OS+UPS3 MQ")
  ups_factor <- 0.096
  histone_theoretical <- 274.44
  } else if (experiment == "DDA4_sam_mq") {
  path_report <- "E:/MaxQuant/ACB/DDA/test4_sam/20251016_sam_MQ/txt/peptides.txt"
  conditions <- c("UPS3")
  condition_labels <- c("200 ng U2OS+UPS3 MQ")
  ups_factor <- 0.24
  histone_theoretical <- 109.78
  } else {
  stop("no data for this experiment (", experiment, ")")
}


# get metadata for selected experiment
data <- load_data(path_report)
experiment_group <- str_split_fixed(experiment, "_", 3)[,1]
MS <- str_split_fixed(experiment, "_", 3)[,2]
specification <- str_split_fixed(experiment, "_", 3)[,3]
path <- dirname(path_report)
rm(path_report)


# folder for results (it doesn't exist)
# if (!dir.exists(str_c(path, "/results_thesis"))) {
#   dir.create(str_c(path, "/results_thesis"))
# }
results_path <- str_c("E:/MaxQuant/ACB/results_thesis/", experiment, "_", format(Sys.Date(), "%Y%m%d"))
if (!dir.exists(results_path)) {
  dir.create(results_path)
}

# load refined speclib and theoretical peptides for iBAQ (and MS name for captions)
if (MS == "popeye") {
  u2os_speclib <- load_data("E:/MaxQuant/ACB/test/speclib_olive/20250428/speclib_tsv/lib.tsv")
  MS_name <- "timsTOF_HT"
} else if (MS == "olive") {
  u2os_speclib <- load_data("E:/MaxQuant/ACB/speclib/speclib_olive/as_txt/lib.tsv")
  MS_name <- "timsTOF_Ultra2"
} else if (MS == "oliveslice") {
  u2os_speclib <- load_data("E:/MaxQuant/ACB/speclib/speclib_olive_slice/as_txt/lib.tsv")
  MS_name <- "timsTOF_Ultra2_slicePASEF"
} else if (MS == "fozzie") {
  u2os_speclib <- load_data("E:/MaxQuant/ACB/speclib/speclib_fozzie/as_txt/lib.tsv")
  MS_name <- "Orbitrap_Exploris"
} else if (MS == "sam") {
  u2os_speclib <- load_data("E:/MaxQuant/ACB/speclib/speclib_sam/as_txt/lib.tsv")
  MS_name <- "Orbitrap_Astral"
} else {
  stop("no refined speclib for this experiment (", MS, ")")
}

# u2os_tp_1MC <- theoretical_peptides("E:/MaxQuant/ACB/databases/uniprotkb_proteome_UP000005640_AND_revi_2025_08_12.fasta") %>%
  # mutate(Protein.Group_short = str_extract(Protein.Group, "(?<=\\|)[^|]+(?=\\|)"))
# write_tsv(u2os_tp_1MC, "E:/MaxQuant/ACB/databases/u2os_tp_1MC.tsv")
# u2os_tp_2MC <- theoretical_peptides("E:/MaxQuant/ACB/databases/uniprotkb_proteome_UP000005640_AND_revi_2025_08_12.fasta", missed_cleavages = 2) %>%
#   mutate(Protein.Group_short = str_extract(Protein.Group, "(?<=\\|)[^|]+(?=\\|)"))
# write_tsv(u2os_tp_2MC, "E:/MaxQuant/ACB/databases/u2os_tp_2MC.tsv")
u2os_tp <- load_data("E:/MaxQuant/ACB/databases/u2os_tp_2MC.tsv")    # takes too much time to calculate every time, choose 1 or 2 missed cleavages


# u2os_mw <- readAAStringSet("E:/MaxQuant/ACB/databases/uniprotkb_proteome_UP000005640_AND_revi_2025_08_12.fasta") %>%
#   {
#     tibble(protein_id = names(.), sequence = as.character(.))
#   } %>%
# mutate(clean_sequence = gsub("[^ACDEFGHIKLMNPQRSTVWY]", "", sequence),
#        seq_length = nchar(clean_sequence),
#        protein_mw_Da = Peptides::mw(clean_sequence),
#        Protein.Group_short = str_extract(protein_id, "(?<=\\|)[^|]+(?=\\|)")) %>%
# mutate(protein_mw = protein_mw_Da * (1.66*10^(-12)))
# write_tsv(u2os_mw, "E:/MaxQuant/ACB/databases/u2os_mw.tsv")
u2os_mw <- load_data("E:/MaxQuant/ACB/databases/u2os_mw.tsv")    # takes too much time to calculate every time
u2os_mw_longest <- u2os_mw %>%
  group_by(Protein.Group_short) %>%
  slice_max(seq_length, n = 1, with_ties = FALSE) %>%   # bei mehreren gleich langen wird nur eine (erste) behalten
  ungroup()

# load standard information and speclib
if (any(grepl("UPS2", conditions))) {
  ups2_standard <- load_data("E:/MaxQuant/ACB/databases/ups2_standard.txt") 
  ups2_speclib <- load_data("E:/MaxQuant/ACB/speclib/predicted_speclib/UPS2/as_txt/lib.tsv")
  ups2_tp <- theoretical_peptides("E:/MaxQuant/ACB/databases/ups2.fasta") %>% 
    mutate(Protein.Group_short = str_extract(Protein.Group, "^[^|]+"))
  standard = "UPS2"}

if (any(grepl("UPS3", conditions))) {
  ups3_standard <- load_data("E:/MaxQuant/ACB/databases/ups3_standard.txt")
  ups3_speclib <- load_data("E:/MaxQuant/ACB/speclib/predicted_speclib/UPS3/as_txt/lib.tsv")
  ups3_tp <- theoretical_peptides("E:/MaxQuant/ACB/databases/ups3.fasta") %>% 
    mutate(Protein.Group_short = str_extract(Protein.Group, "^[^|]+"))
  standard = "UPS3"}

if (any(grepl("UPS2", conditions)) & any(grepl("UPS3", conditions))) {
  standard = "both standards in the data"}
