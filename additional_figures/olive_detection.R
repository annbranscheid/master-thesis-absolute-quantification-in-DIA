source("functions.R")
source("functions_analysis.R")
library(patchwork)

olive_ups2 <- load_data("E:/MaxQuant/ACB/speclib/speclib_olive/combined_lib_ups2.tsv")
olive_ups3 <- load_data("E:/MaxQuant/ACB/speclib/speclib_olive/combined_lib_ups3.tsv")
oliveslice_ups2 <- load_data("E:/MaxQuant/ACB/speclib/speclib_olive_slice/combined_lib_ups2.tsv")
oliveslice_ups3 <- load_data("E:/MaxQuant/ACB/speclib/speclib_olive_slice/combined_lib_ups3.tsv")

# DIA2 olive ----
experiment <- "DIA2_olive"
source("load_data.R")
source("filter_overview.R")
DIA2_olive_ups2 <- data %>% 
  filter(condition == "dia_HSH_10ng_UPS2") %>% 
  filter(!str_detect(Protein.Group, ";")) %>% 
  mutate(rt_prediction_accuracy = Predicted.RT - RT)

DIA2_olive_ups2_speclib <- olive_ups2 %>% 
  filter(!str_detect(ProteinGroup, ";")) %>% 
  distinct(PeptideSequence, .keep_all = TRUE) %>% 
  mutate(detected = PeptideSequence %in% DIA2_olive_ups2$Sequence)

ggplot(DIA2_olive_ups2, aes(rt_prediction_accuracy)) +
  geom_histogram(bins = 20, fill = "lightskyblue") +
  labs(title = "RT prediction accuracy - dia_HSH_10ng_UPS2", x = "RT - expected RT", y = "Count", caption = paste0("u2os_rt_prediction_accuracy", MS_name, "_dia_HSH_10ng_UPS2")) +
  theme_minimal()
ggsave_png()

ggplot(DIA2_olive_ups2_speclib, aes(Tr_recalibrated)) +
  geom_histogram(bins = 20, fill = "lightskyblue") +
  geom_histogram(data = subset(DIA2_olive_ups2_speclib, (detected == TRUE)), bins = 20,
                 fill = "cornflowerblue", alpha = 0.6, color = NA) +   # Subset detected
  labs(title = "Tr_recalibrated - dia_HSH_10ng_UPS2", x = "expected retention time", y = "Count", caption = paste0("u2os_Tr_recalibrated", MS_name, "_dia_HSH_10ng_UPS2")) +
  theme_minimal()
ggsave_png()

# DIA3 olive UPS2 ----
experiment <- "DIA3_olive_ups2"
source("load_data.R")
source("filter_overview.R")
DIA3_olive_ups2 <- data %>% 
  filter(condition == "dia_10ng_U2OS_UPS2") %>% 
  filter(!str_detect(Protein.Group, ";")) %>% 
  mutate(rt_prediction_accuracy = Predicted.RT - RT)

DIA3_olive_ups2_speclib <- olive_ups2 %>% 
  filter(!str_detect(ProteinGroup, ";")) %>% 
  distinct(PeptideSequence, .keep_all = TRUE) %>% 
  mutate(detected = PeptideSequence %in% DIA3_olive_ups2$Sequence)

ggplot(DIA3_olive_ups2, aes(rt_prediction_accuracy)) +
  geom_histogram(bins = 20, fill = "lightskyblue") +
  labs(title = "RT prediction accuracy - dia_10ng_U2OS_UPS2", x = "RT - expected RT", y = "Count", caption = paste0("u2os_rt_prediction_accuracy", MS_name, "_dia_10ng_U2OS_UPS2")) +
  theme_minimal()
ggsave_png()

ggplot(DIA3_olive_ups2_speclib, aes(Tr_recalibrated)) +
  geom_histogram(bins = 20, fill = "lightskyblue") +
  geom_histogram(data = subset(DIA3_olive_ups2_speclib, (detected == TRUE)), bins = 20,
                 fill = "cornflowerblue", alpha = 0.6, color = NA) +   # Subset detected
  labs(title = "Tr_recalibrated - dia_10ng_U2OS_UPS2", x = "expected retention time", y = "Count", caption = paste0("u2os_Tr_recalibrated", MS_name, "_dia_10ng_U2OS_UPS2")) +
  theme_minimal()
ggsave_png()

# DIA3 olive UPS3 ----
experiment <- "DIA3_olive_ups3"
source("load_data.R")
source("filter_overview.R")
DIA3_olive_ups3 <- data %>% 
  filter(condition == "dia_10ng_U2OS_UPS3") %>% 
  filter(!str_detect(Protein.Group, ";")) %>% 
  mutate(rt_prediction_accuracy = Predicted.RT - RT)

DIA3_olive_ups3_speclib <- olive_ups3 %>% 
  filter(!str_detect(ProteinGroup, ";")) %>% 
  distinct(PeptideSequence, .keep_all = TRUE) %>% 
  mutate(detected = PeptideSequence %in% DIA3_olive_ups3$Sequence)

ggplot(DIA3_olive_ups3, aes(rt_prediction_accuracy)) +
  geom_histogram(bins = 20, fill = "lightskyblue") +
  labs(title = "RT prediction accuracy - dia_10ng_U2OS_UPS3", x = "RT - expected RT", y = "Count", caption = paste0("u2os_rt_prediction_accuracy", MS_name, "_dia_10ng_U2OS_UPS3")) +
  theme_minimal()
ggsave_png()

ggplot(DIA3_olive_ups3_speclib, aes(Tr_recalibrated)) +
  geom_histogram(bins = 20, fill = "lightskyblue") +
  geom_histogram(data = subset(DIA3_olive_ups3_speclib, (detected == TRUE)), bins = 20,
                 fill = "cornflowerblue", alpha = 0.6, color = NA) +   # Subset detected
  labs(title = "Tr_recalibrated - dia_10ng_U2OS_UPS3", x = "expected retention time", y = "Count", caption = paste0("u2os_Tr_recalibrated", MS_name, "_dia_10ng_U2OS_UPS3")) +
  theme_minimal()
ggsave_png()

# DIA2 olive slice UPS2 ----
experiment <- "DIA2_oliveslice"
source("load_data.R")
source("filter_overview.R")
DIA2_oliveslice <- data %>% 
  filter(condition == "slice_HSH_10ng_UPS2") %>% 
  filter(!str_detect(Protein.Group, ";")) %>% 
  mutate(rt_prediction_accuracy = Predicted.RT - RT)

DIA2_oliveslice_speclib <- oliveslice_ups2 %>% 
  filter(!str_detect(ProteinGroup, ";")) %>% 
  distinct(PeptideSequence, .keep_all = TRUE) %>% 
  mutate(detected = PeptideSequence %in% DIA2_oliveslice$Sequence)

ggplot(DIA2_oliveslice, aes(rt_prediction_accuracy)) +
  geom_histogram(bins = 20, fill = "lightskyblue") +
  labs(title = "RT prediction accuracy - slice_HSH_10ng_UPS2", x = "RT - expected RT", y = "Count", caption = paste0("u2os_rt_prediction_accuracy", MS_name, "_slice_HSH_10ng_UPS23")) +
  theme_minimal()
ggsave_png()

ggplot(DIA2_oliveslice_speclib, aes(Tr_recalibrated)) +
  geom_histogram(bins = 20, fill = "lightskyblue") +
  geom_histogram(data = subset(DIA2_oliveslice_speclib, (detected == TRUE)), bins = 20,
                 fill = "cornflowerblue", alpha = 0.6, color = NA) +   # Subset detected
  labs(title = "Tr_recalibrated - slice_HSH_10ng_UPS2", x = "expected retention time", y = "Count", caption = paste0("u2os_Tr_recalibrated", MS_name, "_slice_HSH_10ng_UPS2")) +
  theme_minimal()
ggsave_png()

# DIA3 olive slice UPS2 ----
experiment <- "DIA3_oliveslice_ups2"
source("load_data.R")
source("filter_overview.R")
DIA3_oliveslice_ups2 <- data %>% 
  filter(condition == "slice_10ng_U2OS_UPS2") %>% 
  filter(!str_detect(Protein.Group, ";")) %>% 
  mutate(rt_prediction_accuracy = Predicted.RT - RT)

DIA3_oliveslice_ups2_speclib <- oliveslice_ups2 %>% 
  filter(!str_detect(ProteinGroup, ";")) %>% 
  distinct(PeptideSequence, .keep_all = TRUE) %>% 
  mutate(detected = PeptideSequence %in% DIA3_oliveslice_ups2$Sequence)

ggplot(DIA3_oliveslice_ups2, aes(rt_prediction_accuracy)) +
  geom_histogram(bins = 20, fill = "lightskyblue") +
  labs(title = "RT prediction accuracy - slice_10ng_U2OS_UPS2", x = "RT - expected RT", y = "Count", caption = paste0("u2os_rt_prediction_accuracy", MS_name, "_slice_10ng_U2OS_UPS2")) +
  theme_minimal()
ggsave_png()

ggplot(DIA3_oliveslice_ups2_speclib, aes(Tr_recalibrated)) +
  geom_histogram(bins = 20, fill = "lightskyblue") +
  geom_histogram(data = subset(DIA3_oliveslice_ups2_speclib, (detected == TRUE)), bins = 20,
                 fill = "cornflowerblue", alpha = 0.6, color = NA) +   # Subset detected
  labs(title = "Tr_recalibrated - slice_10ng_U2OS_UPS2", x = "expected retention time", y = "Count", caption = paste0("u2os_Tr_recalibrated", MS_name, "_slice_10ng_U2OS_UPS2")) +
  theme_minimal()
ggsave_png()

# DIA3 olive slice UPS3 ----
experiment <- "DIA3_oliveslice_ups3"
source("load_data.R")
source("filter_overview.R")
DIA3_oliveslice_ups3 <- data %>% 
  filter(condition == "slice_10ng_U2OS_UPS3") %>% 
  filter(!str_detect(Protein.Group, ";")) %>% 
  mutate(rt_prediction_accuracy = Predicted.RT - RT)

DIA3_oliveslice_ups3_speclib <- oliveslice_ups3 %>% 
  filter(!str_detect(ProteinGroup, ";")) %>% 
  distinct(PeptideSequence, .keep_all = TRUE) %>% 
  mutate(detected = PeptideSequence %in% DIA3_oliveslice_ups3$Sequence)

ggplot(DIA3_oliveslice_ups3, aes(rt_prediction_accuracy)) +
  geom_histogram(bins = 20, fill = "lightskyblue") +
  labs(title = "RT prediction accuracy - slice_10ng_U2OS_UPS3", x = "RT - expected RT", y = "Count", caption = paste0("u2os_rt_prediction_accuracy", MS_name, "_slice_10ng_U2OS_UPS3")) +
  theme_minimal()
ggsave_png()

ggplot(DIA3_oliveslice_ups3_speclib, aes(Tr_recalibrated)) +
  geom_histogram(bins = 20, fill = "lightskyblue") +
  geom_histogram(data = subset(DIA3_oliveslice_ups3_speclib, (detected == TRUE)), bins = 20,
                 fill = "cornflowerblue", alpha = 0.6, color = NA) +   # Subset detected
  labs(title = "Tr_recalibrated - slice_10ng_U2OS_UPS3", x = "expected retention time", y = "Count", caption = paste0("u2os_Tr_recalibrated", MS_name, "_slice_10ng_U2OS_UPS3")) +
  theme_minimal()
ggsave_png()



# note ----
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

