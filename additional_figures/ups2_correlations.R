#run main script before
ups2_score_H <- data %>% 
  filter(grepl("ups", Protein.Group), (Channel == "H" | Channel == "H, H"), intensity > 0) %>% 
  get_score(score_name = score_name, ups2_tp) %>%
  rename(int_score_H = int_score) %>% 
  # condition_name() %>%
  mutate(protein = sub("ups.*", "", Protein.Group))

ups2_score <- attr(data, "ups2_score") %>% 
  left_join(ups2_score_H %>% select(protein, Run, `int_score_H`), by = c("protein", "Run")) %>% 
  mutate(ups2_PTMs = ups2_standard$`Potential PTMs`[match(protein, ups2_standard$`UniProt Accession Number`)]) %>% 
  mutate(ups2_MW = ups2_standard$`MW (Da)`[match(protein, ups2_standard$`UniProt Accession Number`)]) %>% 
  mutate(ups2_source = ups2_standard$`Source or recombinant`[match(protein, ups2_standard$`UniProt Accession Number`)])

ups2_score_info <- attr(data, "ups2_score_info") 


# adjust color, subtitle and caption
ggplot(ups2_score, aes(log10(ups2_amount), log10(int_score),  color = log10(int_score_H))) +
  geom_point() +
  facet_wrap(~ factor(condition_label, levels = condition_labels)) +
  labs(title = paste0("correlation ", score_name, " with ups2 amount"), 
       subtitle = "colored by iBAQ score of the H protein",
       caption = paste0("ups2_", score_name, "_", MS_name, "_score_H"), 
       x = "log10 ups2 in standard [fmol in vial]", y = paste0("log10 ", score_name, " score")) +
  theme_large() +
  scale_fill_gradientn(colours = rev(c("#f7fbff", "#6baed6", "#08306b"))) +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 10)) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
  geom_text(data = ups2_score_info, aes(x = Inf, y = Inf, label = cor_score_2),
            hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) 

ggsave_png()





acetylated_proteins <- ups2_score %>% 
  filter(ups2_PTMs == "Acetylation") 
acetylated_proteins_list <- unique(acetylated_proteins$Protein.Group)

acetylated_proteins_peptides <- data %>% 
  filter(Protein.Group %in% acetylated_proteins_list, (Channel == "L" | Channel == "L, L"), intensity > 0) %>% 
  group_by(condition, Run, Protein.Group, Sequence) %>% 
  summarise(intensity_sequence = sum(intensity)) %>% 
  ungroup() %>% 
  # group_by(condition, Protein.Group, Sequence) %>% 
  # summarise(intensity_sequence_median = median(intensity_sequence)) %>% 
  # ungroup() %>% 
  label_conditions() %>% 
  mutate(
    seq_class = case_when(
      !grepl("[KR]$", Sequence) ~ 1,   # does NOT end with K or R - possible end
      grepl("IKASFK$", Sequence) ~ 1,   # possible end for CAH2
      grepl("K", Sequence)      ~ 2,   # contains K anywhere
      TRUE                      ~ 3    # everything else
      ),
    seq_class = factor(seq_class, levels = c(1, 2, 3),
      labels = c("N-terminal peptide", "Contains K", "no possible Acetylation"))) %>% 
  mutate(protein_name = str_remove(Protein.Group, "^.*\\|") %>%   # remove everything before |
                 str_remove("_UPS$")) # and the _UPS in the end

  
ggplot(acetylated_proteins_peptides, aes(protein_name, log10(intensity_sequence),  color = seq_class)) +
  geom_point() +
  facet_wrap(~ factor(condition_label, levels = condition_labels))+
  labs(title = paste0("peptide intensity of poteintially acetylated UPS2 proteins"), 
       subtitle = "colored by the potential of acetylation",
       caption = paste0("ups2_peptides", MS_name, "_acetylation"), 
       x = "", y = paste0("median peptide intensity")) +
  theme_small() +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 10))

ggsave_png()
