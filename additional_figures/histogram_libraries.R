# speclib histograms
source("functions.R")
results_path <- "E:/MaxQuant/ACB/results_thesis/speclib_histograms/"

# load predicted and refined speclibs
u2os_speclib_predicted <- load_data("E:/MaxQuant/ACB/speclib/speclib_sam/as_txt/lib.tsv") %>% 
  filter(!str_detect(ProteinGroup, ";"))
u2os_speclib_popeye <- load_data("E:/MaxQuant/ACB/test/speclib_olive/20250428/speclib_tsv/lib.tsv") %>% 
  filter(!str_detect(ProteinGroup, ";"))
u2os_speclib_olive <- load_data("E:/MaxQuant/ACB/speclib/speclib_olive/as_txt/lib.tsv") %>% 
  filter(!str_detect(ProteinGroup, ";"))
u2os_speclib_oliveslice <- load_data("E:/MaxQuant/ACB/speclib/speclib_olive_slice/as_txt/lib.tsv") %>% 
  filter(!str_detect(ProteinGroup, ";"))
u2os_speclib_fozzie <- load_data("E:/MaxQuant/ACB/speclib/speclib_fozzie/as_txt/lib.tsv") %>% 
  filter(!str_detect(ProteinGroup, ";"))
u2os_speclib_sam <- load_data("E:/MaxQuant/ACB/speclib/speclib_sam/as_txt/lib.tsv") %>% 
  filter(!str_detect(ProteinGroup, ";"))


speclibs_list <- list(
  predicted_library = u2os_speclib_predicted,
  DIA1_library = u2os_speclib_popeye,
  timsTOF_Ultra2 = u2os_speclib_olive,
  timsTOF_Ultra2_slicePASEF = u2os_speclib_oliveslice,
  Orbitrap_Exploris = u2os_speclib_fozzie,
  Orbitrap_Astral = u2os_speclib_sam)

# Count unique values per data frame
unique_peptideSequence <- map_dbl(speclibs_list, ~ n_distinct(.x$PeptideSequence))
unique_ProteinGroup <- map_dbl(speclibs_list, ~ n_distinct(.x$ProteinGroup))

# Turn into a data frame for plotting
plot_df <- tibble::tibble(
  dataset = names(unique_peptideSequence),
  PeptideSequence = unique_peptideSequence,
  ProteinGroup = unique_ProteinGroup)
plot_df$dataset <- factor(plot_df$dataset, levels = plot_df$dataset)


ggplot(plot_df, aes(dataset, PeptideSequence)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_small() +
  # scale_y_log10() +
  labs(title = "Number of unique peptide sequences in refined spectral libraries", 
       x = "",
       y = "Number of peptide sequences",
       caption = "spectral_libraries_peptides")
ggsave_png()

ggplot(plot_df, aes(dataset, ProteinGroup)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_small() +
  # scale_y_log10() +
  labs(title = "Number of unique protein groups in refined spectral libraries", 
       x = "",
       y = "Number of protein groups",
       caption = "spectral_libraries_proteins")
ggsave_png()

save_table(plot_df, "anz_speclibs")
