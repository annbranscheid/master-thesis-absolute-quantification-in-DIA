source("functions.R")
source("functions_analysis.R")
library(patchwork)

experiment <- "DIA3_fozzie_ups2"
  
source("load_data.R")
data_UPS2 <- data %>% 
  #mutate(Channel = str_remove(str_extract(Precursor.Id, "-[LMH]"), "-")) %>%   # for olive
  filter(grepl("ups", Protein.Group), (Channel == "L" | Channel == "L, L"), Ms1.Area > 0)
source("filter_overview.R")
data_UPS2_filtered <- data %>% 
  #mutate(Channel = str_remove(str_extract(Precursor.Id, "-[LMH]"), "-")) %>%   # for olive
  filter(grepl("ups", Protein.Group), (Channel == "L" | Channel == "L, L"), intensity > 0)

detected_ups2 <- ups2_speclib %>% 
  distinct(PeptideSequence, .keep_all = TRUE) %>% 
  mutate(detected = PeptideSequence %in% data_UPS2$Stripped.Sequence, 
         processed = PeptideSequence %in% data_UPS2_filtered$Sequence) %>% 
  mutate(prot = (str_sub(ProteinGroup, end = str_locate(ProteinGroup, "ups") [,1] -1))) %>%   # basically UniProt ID
  mutate(ups2_amount = ups2_standard$`UPS2 Amount (fmol)`[match(prot, ups2_standard$`UniProt Accession Number`)],
         MW = ups2_standard$`MW (Da)`[match(prot, ups2_standard$`UniProt Accession Number`)],
         Source = ups2_standard$`Source or recombinant`[match(prot, ups2_standard$`UniProt Accession Number`)],
         Host = ups2_standard$Host[match(prot, ups2_standard$`UniProt Accession Number`)],
         Tag = ups2_standard$Tag[match(prot, ups2_standard$`UniProt Accession Number`)],
         PTMs = ups2_standard$`Potential PTMs`[match(prot, ups2_standard$`UniProt Accession Number`)])   # add the protein amount from standard


ggplot(detected_ups2, aes(Tr_recalibrated)) +
  geom_histogram(bins = 20, fill = "lightskyblue") +
  geom_histogram(data = subset(detected_ups2, (detected == TRUE)), bins = 20,
                 fill = "cornflowerblue", alpha = 0.6, color = NA) +   # Subset detected
  geom_histogram(data = subset(detected_ups2, (processed == TRUE)), bins = 20,
                 fill = "royalblue4", alpha = 0.6, color = NA) +   # Subset detected
  # geom_histogram(data = subset(detected_ups2, PTMs == "Acetylation"), bins = 20,
  #                fill = "cornflowerblue", alpha = 0.6, color = NA) +   # Subset detected
  # geom_histogram(data = subset(detected_ups2, (PTMs == "Acetylation" & detected == TRUE)), bins = 20,
  #                fill = "royalblue4", alpha = 0.6, color = NA) +   # Subset detected
  # geom_histogram(data = subset(detected_ups2, (PTMs == "Acetylation" & processed == TRUE)), bins = 20,
  #                fill = "black", alpha = 0.6, color = NA) +   # Subset detected
  labs(title = "Tr_recalibrated - Acetylation", x = "expected retention time", y = "Count", caption = paste0("Tr_recalibrated_", MS_name)) +
  theme_minimal()
ggsave_png()



pot_factors <- c("Tr_recalibrated", "IonMobility", "PrecursorMz")
detected_ups2_long <- detected_ups2 %>%
  pivot_longer(cols = all_of(pot_factors), 
               names_to = "variable", 
               values_to = "value")

ggplot(detected_ups2_long, aes(x = value)) +
  # All data (background)
  geom_histogram(bins = 20, fill = "lightskyblue") +
  # Detected subset
  geom_histogram(data = detected_ups2_long %>% filter(detected == TRUE),
    bins = 20, fill = "cornflowerblue", alpha = 0.6, color = NA) +
   # Processed subset
  geom_histogram(data = detected_ups2_long %>% filter(processed == TRUE),
    bins = 20, fill = "royalblue4", alpha = 0.6, color = NA) +
  facet_wrap(~variable, scales = "free_x") +
  labs(title = "Detected and processed peptides across variables", 
       subtitle = "Light: all peptides | Medium: detected | Dark: processed", x = NULL, 
       y = "Count", caption = "detected_peptides_ups2") +
  theme_large() 
ggsave_png()


ggsave_png <- function(name = last_plot()$labels$caption) {
  ggsave(
    # filename <- str_c(path, "/results_thesis/", format(Sys.Date(), "%Y%m%d"), "_", name, ".png"),   # to save in the same folder as data
    filename <- str_c(results_path, "/", format(Sys.Date(), "%Y%m%d"), "_", name, ".png"),   # to save results of all MS in one folder
    width = 17,
    height = 7.53,
    dpi = 300,
    units = "in",
    device = "png"
  )
}
