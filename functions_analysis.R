# part 1 ----
prepare_data <- function(data) {
  data <- data %>%
    { if (any(grepl("DIA", path))) rename(.,intensity = Ms1.Area, Sequence = Stripped.Sequence) else . } %>%
    { if (any(grepl("MQ", path))) rename(.,Protein.Group = `Leading razor protein`) else . } %>%   # for DDA MQ
    { if (any(grepl("diann", path))) rename(.,intensity = Ms1.Area, Sequence = Stripped.Sequence) else . } %>%   # for DDA diann
    condition_name() %>%
    filter(!condition %in% remove_conditions) %>%
    filter(!str_detect(Protein.Group, ";"), intensity > 0) %>%
    normalisation() %>%
    mutate(factor_run = attr(., "scaling_factor")$factor_run[match(Run, attr(., "scaling_factor")$Run)]) %>%
    mutate(intensity = intensity * factor_run)
return(data)}

# DiaNN ups2
anz_ups2 <- function(data) {
  anz_u2os <- data %>%
    filter(!grepl("ups", Protein.Group)) %>%
    distinct(Run, Protein.Group, Sequence) %>%  # unique peptide-protein-run combinations
    group_by(Run) %>%
    summarise(
      u2os_proteins = n_distinct(Protein.Group),   # number of proteins
      u2os_peptides  = n_distinct(Sequence),  # number of peptides
      .groups = "drop")
  anz_ups2 <- data %>%
    filter(grepl("ups", Protein.Group), (Channel == "L" | Channel == "L, L")) %>%
    distinct(Run, Protein.Group, Sequence) %>%  # unique peptide-protein-run combinations
    group_by(Run) %>%
    summarise(
      ups2_proteins = n_distinct(Protein.Group),   # number of proteins
      ups2_peptides  = n_distinct(Sequence),  # number of peptides
      .groups = "drop")
  anz <- left_join(anz_u2os, anz_ups2, by = "Run") %>%     # combine u2os and ups2 & return
    condition_name() %>%
    sort_by_run()
  
  attr(data, "anz") <- anz
  attr(data, "anz_plot") <- anz_hist_ups2(anz = anz)
  
  anz_condition_u2os <- data %>%
    filter(!grepl("ups", Protein.Group)) %>%
    condition_name() %>%
    distinct(condition, Protein.Group, Sequence) %>%  # unique peptide-protein-condition combinations
    group_by(condition) %>%
    summarise(
      u2os_proteins = n_distinct(Protein.Group),   # number of proteins
      u2os_peptides  = n_distinct(Sequence),  # number of peptides
      .groups = "drop")
  anz_condition_ups2 <- data %>%
    filter(grepl("ups", Protein.Group), (Channel == "L" | Channel == "L, L")) %>%
    distinct(condition, Protein.Group, Sequence) %>%  # unique peptide-protein-condition combinations
    group_by(condition) %>%
    summarise(
      ups2_proteins = n_distinct(Protein.Group),   # number of proteins
      ups2_peptides  = n_distinct(Sequence),  # number of peptides
      .groups = "drop")
  anz_condition <- left_join(anz_condition_u2os, anz_condition_ups2, by = "condition") %>%     # combine u2os and ups2 & return
    mutate(condition = factor(condition, levels = conditions)) %>%
    arrange(condition)
  
  attr(data, "anz_condition") <- anz_condition
  
  data_ups2 <- data %>% 
    filter(grepl("ups", Protein.Group), (Channel == "L" | Channel == "L, L"))
  run <- unique(data_ups2$Run)
  prot <- unique(data_ups2$Protein.Group)
  ups2_protein_info <- prot %>%
    map_dfr(~ {data_ups2 %>%
        filter(grepl(.x, Protein.Group)) %>%
        distinct(Run, Sequence) %>%
        count(Run, name = "anz_peptides") %>%
        mutate(prot = .x)}) %>% 
    pivot_wider(names_from  = "Run", values_from  = "anz_peptides")
  
  attr(data, "ups_protein_info") <- ups2_protein_info
  
  return(data)}


anz_hist_ups2 <- function(anz) {
  anz_long <- anz %>%
    pivot_longer(cols = c(u2os_proteins, u2os_peptides, ups2_proteins, ups2_peptides), names_to = "plot", values_to = "Count")
  
  anz_long$Run <- factor(anz_long$Run, levels = unique(anz_long$Run))
 
   anz_plot <- ggplot(anz_long, aes(Run, Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    facet_wrap(~ plot, scales = "free") +
    theme_small() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "protein and peptide Counts", caption = paste0("anz_", MS_name), x = "", y = "count")
}



# DiaNN ups3
anz_ups3 <- function(data) {
  anz_u2os <- data %>%
    filter(!grepl("ups3", Protein.Group)) %>% 
    distinct(Run, Protein.Group, Sequence) %>%  # unique peptide-protein-run combinations
    group_by(Run) %>%
    summarise(
      u2os_proteins = n_distinct(Protein.Group),   # number of proteins
      u2os_peptides  = n_distinct(Sequence),  # number of peptides
      .groups = "drop")
  anz_ups3 <- data %>%
    filter(grepl("ups3", Protein.Group)) %>% 
    distinct(Run, Protein.Group, Sequence) %>%  # unique peptide-protein-run combinations
    group_by(Run) %>%
    summarise(
      ups3_proteins = n_distinct(Protein.Group),   # number of proteins
      ups3_peptides  = n_distinct(Sequence),  # number of peptides
      .groups = "drop")
  anz <- left_join(anz_u2os, anz_ups3, by = "Run") %>%    # combine u2os and ups3 & return
    condition_name() %>%
    sort_by_run()
  
  attr(data, "anz") <- anz
  attr(data, "anz_plot") <- anz_hist_ups3(anz = anz)
  
  anz_condition_u2os <- data %>%
    filter(!grepl("ups3", Protein.Group)) %>%
    condition_name() %>%
    distinct(condition, Protein.Group, Sequence) %>%  # unique peptide-protein-condition combinations
    group_by(condition) %>%
    summarise(
      u2os_proteins = n_distinct(Protein.Group),   # number of proteins
      u2os_peptides  = n_distinct(Sequence),  # number of peptides
      .groups = "drop")
  anz_condition_ups2 <- data %>%
    filter(grepl("ups3", Protein.Group)) %>%
    distinct(condition, Protein.Group, Sequence) %>%  # unique peptide-protein-condition combinations
    group_by(condition) %>%
    summarise(
      ups2_proteins = n_distinct(Protein.Group),   # number of proteins
      ups2_peptides  = n_distinct(Sequence),  # number of peptides
      .groups = "drop")
  anz_condition <- left_join(anz_condition_u2os, anz_condition_ups2, by = "condition") %>%     # combine u2os and ups2 & return
    mutate(condition = factor(condition, levels = conditions)) %>%
    arrange(condition)
  
  attr(data, "anz_condition") <- anz_condition
  
  data_ups3 <- data %>% 
    filter(grepl("ups3", Protein.Group)) 
  run <- unique(data_ups3$Run)
  prot <- unique(data_ups3$Protein.Group)
  ups3_protein_info <- prot %>%
    map_dfr(~ {data_ups3 %>%
        filter(grepl(.x, Protein.Group)) %>%
        distinct(Run, Sequence) %>%
        count(Run, name = "anz_peptides") %>%
        mutate(prot = .x)}) %>% 
    pivot_wider(names_from  = "Run", values_from  = "anz_peptides")
  
  attr(data, "ups_protein_info") <- ups3_protein_info
  
  return(data)}


anz_hist_ups3 <- function(anz) {
  anz_long <- anz %>%
    pivot_longer(cols = c(u2os_proteins, u2os_peptides, ups3_proteins, ups3_peptides), names_to = "plot", values_to = "Count")
  
  anz_long$Run <- factor(anz_long$Run, levels = unique(anz_long$Run))
  
  anz_plot <- ggplot(anz_long, aes(Run, Count)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    facet_wrap(~ plot, scales = "free") +
    theme_small() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "protein and peptide Counts", caption = paste0("anz_", MS_name), x = "", y = "count")
}

# MQ ups2
anz_MQ <- function(data) {
  anz_u2os <- data %>%
    filter(Channel == "H") %>%
    distinct(Run, Protein.Group, Sequence) %>%  # unique peptide-protein-run combinations
    group_by(Run) %>%
    summarise(
      u2os_proteins = n_distinct(Protein.Group),   # number of proteins
      u2os_peptides  = n_distinct(Sequence),  # number of peptides
      .groups = "drop")
  anz_ups2 <- data %>%
    filter(grepl("ups", Protein.Group), (Channel == "L")) %>%
    distinct(Run, Protein.Group, Sequence) %>%  # unique peptide-protein-run combinations
    group_by(Run) %>%
    summarise(
      ups2_proteins = n_distinct(Protein.Group),   # number of proteins
      ups2_peptides  = n_distinct(Sequence),  # number of peptides
      .groups = "drop")
  anz <- left_join(anz_u2os, anz_ups2, by = "Run") %>%     # combine u2os and ups2 & return
    condition_name() %>%
    sort_by_run()
  
  attr(data, "anz") <- anz
  attr(data, "anz_plot") <- anz_hist_ups2(anz = anz)
  
  anz_condition_u2os <- data %>%
    filter(Channel == "H") %>%
    condition_name() %>%
    distinct(condition, Protein.Group, Sequence) %>%  # unique peptide-protein-condition combinations
    group_by(condition) %>%
    summarise(
      u2os_proteins = n_distinct(Protein.Group),   # number of proteins
      u2os_peptides  = n_distinct(Sequence),  # number of peptides
      .groups = "drop")
  anz_condition_ups2 <- data %>%
    filter(grepl("ups", Protein.Group), (Channel == "L")) %>%
    distinct(condition, Protein.Group, Sequence) %>%  # unique peptide-protein-condition combinations
    group_by(condition) %>%
    summarise(
      ups2_proteins = n_distinct(Protein.Group),   # number of proteins
      ups2_peptides  = n_distinct(Sequence),  # number of peptides
      .groups = "drop")
  anz_condition <- left_join(anz_condition_u2os, anz_condition_ups2, by = "condition") %>%     # combine u2os and ups2 & return
    mutate(condition = factor(condition, levels = conditions)) %>%
    arrange(condition)
  
  attr(data, "anz_condition") <- anz_condition
  
  data_ups2 <- data %>% 
    filter(grepl("ups", Protein.Group), (Channel == "L"))
  run <- unique(data_ups2$Run)
  prot <- unique(data_ups2$Protein.Group)
  ups2_protein_info <- prot %>%
    map_dfr(~ {data_ups2 %>%
        filter(grepl(.x, Protein.Group)) %>%
        distinct(Run, Sequence) %>%
        count(Run, name = "anz_peptides") %>%
        mutate(prot = .x)}) %>% 
    pivot_wider(names_from  = "Run", values_from  = "anz_peptides")
  
  attr(data, "ups_protein_info") <- ups2_protein_info
  
  return(data)}


# part 2 ----
ups2_score <- function(data, protein_tp = ups2_tp, ups2_standard = get("ups2_standard", envir = .GlobalEnv), score_name = get("score_name", envir = parent.frame()), mean_score = "no") {
  ups2_score <- data %>% 
    filter(grepl("ups", Protein.Group), (Channel == "L" | Channel == "L, L"), intensity > 0) %>% 
    # group_by(Run, Protein.Group) %>%
    # filter(abs(intensity - median(intensity)) / median(intensity) <= 10) %>%   # filter peptides whose intensity lies outside the 100% range
    # ungroup() %>%
    get_score(score_name = score_name, protein_tp) %>%
    condition_name() %>%
    {if(mean_score == "yes"){   # mean score for each protein in each condition
      group_by(., condition, Protein.Group) %>%
      summarise(int_score = mean(int_score), .groups = "drop")
    } else {
      .}} %>%
    mutate(protein = sub("ups.*", "", Protein.Group)) %>% 
    mutate(ups2_amount = ups2_standard$`UPS2 Amount (fmol)`[match(protein, ups2_standard$`UniProt Accession Number`)]) %>% 
    label_conditions()

  attr(data, "ups2_score") <- ups2_score

  ups2_score <- ups2_score %>%
    slopes_score(data_log10 = "no") %>%
    correlation_score(data_log10 = "no")

  ups2_score_info <- cbind(attr(ups2_score, "slopes"), attr(ups2_score, "correlation_score")[,-1 ]) %>% 
    label_conditions()

  attr(data, "ups2_score_info") <- ups2_score_info

  attr(data, "ups2_score_plot") <- ggplot(ups2_score, aes(log10(ups2_amount), log10(int_score),  color = protein)) +
    geom_point() +
    facet_wrap(~ factor(condition_label, levels = condition_labels)) +
    # scale_x_log10() +
    # ylim((6), (12)) +
    # scale_y_log10() +
    labs(title = paste0("correlation ", score_name, " with ups2 amount"), caption = paste0("ups2_", score_name, "_", MS_name), x = "log10 ups2 in standard [fmol in vial]", y = paste0("log10 ", score_name, " score")) +
    theme_large() +
    theme(legend.position = "none") +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
    geom_text(data = ups2_score_info, aes(x = Inf, y = Inf, label = cor_score_2),
              hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) # +
    # geom_text(data = ups2_score_info, aes(x = Inf, y = Inf, label = slope),
    #           inherit.aes = FALSE, hjust = 1.1, vjust = 2.5, size = 3, color = "black") +
    # geom_text(data = ups2_score_info, aes(x = Inf, y = Inf, label = intercept),
    #           inherit.aes = FALSE, hjust = 1.1, vjust = 3.9, size = 3, color = "black")
  
  return(data)}


ups3_score <- function(data, protein_tp = ups3_tp, ups3_standard = get("ups3_standard", envir = .GlobalEnv), score_name = get("score_name", envir = parent.frame())) {
  ups3_score <- data %>% 
    filter(grepl("ups3", Protein.Group), intensity > 0) %>% 
    # group_by(Run, Protein.Group) %>%
    # filter(abs(intensity - median(intensity)) / median(intensity) <= 10) %>%   # filter peptides whose intensity lies outside the 100% range
    # ungroup() %>% 
    get_score(score_name = score_name, protein_tp) %>% 
    mutate(protein = sub("ups3.*", "", Protein.Group)) %>% 
    condition_name() %>%
    mutate(ups3_amount = ups3_standard$`amount[pmol]`[match(protein, ups3_standard$protein)]) %>% 
    label_conditions()
  
  
  attr(data, "ups3_score") <- ups3_score
  
  ups3_score <- ups3_score %>% 
    slopes_score(data_log10 = "no", ups3 = "yes") %>% 
    correlation_score(data_log10 = "no", ups3 = "yes")
  
  ups3_score_info <- cbind(attr(ups3_score, "slopes"), attr(ups3_score, "correlation_score")[,-1 ])
  
  attr(data, "ups3_score_info") <- ups3_score_info
  
  attr(data, "ups3_score_plot") <- ggplot(ups3_score, aes(log10(ups3_amount), log10(int_score),  color = protein)) +
    geom_point() +
    facet_wrap(~ factor(condition_label, levels = condition_labels)) +
    # scale_x_log10() +
    #ylim((6), (12)) +
    # scale_y_log10() +
    theme_large() +
    labs(title = paste0("correlation ", score_name, " with ups3 amount"), caption = paste0("ups3_", score_name, "_", MS_name), x = "log10 ups3 in standard [pmol in vial]", y = paste0("log10 ", score_name, " score")) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
    geom_text(data = ups3_score_info, aes(x = Inf, y = Inf, label = cor_score_2),
              hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) # +
    # geom_text(data = ups3_score_info, aes(x = Inf, y = Inf, label = slope),
    #           inherit.aes = FALSE, hjust = 1.1, vjust = 2.5, size = 3, color = "black") +
    # geom_text(data = ups3_score_info, aes(x = Inf, y = Inf, label = intercept),
    #           inherit.aes = FALSE, hjust = 1.1, vjust = 3.9, size = 3, color = "black")
  
  return(data)}


# part 2b
ups3_MQ_score <- function(data_pg, score = "iBAQ", ups3_standard = get("ups3_standard", envir = .GlobalEnv)) {
  ups3_MQ_score <- data_pg %>% 
    filter(grepl("ups3", `Majority protein IDs`)) %>% 
    pivot_longer(cols = starts_with(score), names_to = "Run", values_to = "int_score") %>% 
    condition_name() %>% 
    mutate(protein = sub("ups3.*", "", `Majority protein IDs`)) %>% 
    mutate(ups3_amount = ups3_standard$`amount[pmol]`[match(protein, ups3_standard$protein)]) %>% 
    label_conditions()
  
  attr(data_pg, "ups3_MQ_score") <- ups3_MQ_score
  
  
  ups3_MQ_score <- ups3_MQ_score %>% 
    slopes_score(data_log10 = "no", ups3 = "yes") %>% 
    correlation_score(data_log10 = "no", ups3 = "yes")
  
  ups3_MQ_score_info <- cbind(attr(ups3_MQ_score, "slopes"), attr(ups3_MQ_score, "correlation_score")[,-1 ])
  
  attr(data_pg, "ups3_MQ_score_info") <- ups3_MQ_score_info
  
  attr(data_pg, "ups3_MQ_score_plot") <- ggplot(ups3_MQ_score, aes(ups3_amount, log10(int_score),  color = protein)) +
    geom_point() +
    facet_wrap(~ factor(condition, levels = conditions)) +
    scale_x_log10() +
    ylim((6),(12)) +
    # scale_y_log10() +
    theme_large() +
    labs(title = paste0("correlation ", score, " with ups3 amount"), caption = paste0("ups3_MQ_", score, "_", MS_name), x = "ups3 in standard [pmol in vial]", y = paste0(score, " score (log10)")) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 1) +
    geom_text(data = ups3_MQ_score_info, aes(x = Inf, y = Inf, label = cor_score_2),
              hjust = 1.1, vjust = 1.1, size = 3, inherit.aes = FALSE) +
    geom_text(data = ups3_MQ_score_info, aes(x = Inf, y = Inf, label = slope),
              inherit.aes = FALSE, hjust = 1.1, vjust = 2.5, size = 3, color = "black") +
    geom_text(data = ups3_MQ_score_info, aes(x = Inf, y = Inf, label = intercept),
              inherit.aes = FALSE, hjust = 1.1, vjust = 3.9, size = 3, color = "black")
  
  return(data_pg)} 




# part 3 ----
u2os_score <- function(data, protein_tp = u2os_tp, standard = get("standard", envir = parent.frame()), condition_filter = get("condition_filter", envir = parent.frame()), score_name = get("score_name", envir = parent.frame())) {

  u2os_prot_groups <- data %>%
    { if (standard == "UPS2") filter(., !grepl("ups", Protein.Group), !(Channel %in% c("L", "L, L"))) else . } %>%
    { if (standard == "UPS3") filter(., !grepl("ups3", Protein.Group)) else . } %>%
    filter(!str_detect(Protein.Group, ";"), intensity > 0)  %>%
    find_ribosome_histone_proteasome()

    ribosomes_df <- u2os_prot_groups$ribosomes_df %>%
      mutate(ribosome_type = case_when(
        str_detect(description, regex("mitochondrial", ignore_case = TRUE)) ~ "Mitochondrial",
        str_detect(external_gene_name, "^RPS") ~ "40S",
        str_detect(external_gene_name, "^RPL") ~ "60S",
        TRUE ~ "Unknown")) %>%
      mutate(ribosome_class = case_when(
        ribosome_type == "40S" ~ 1,
        ribosome_type == "60S" ~ 2,
        ribosome_type == "Mitochondrial" ~ 3,
        TRUE ~ 0))
    histones_df <- u2os_prot_groups$histones_df %>%
      mutate(histone_type = case_when(
        str_detect(external_gene_name, "^H2A") ~ "H2A",
        str_detect(external_gene_name, "^H2B") ~ "H2B",
        str_detect(external_gene_name, "^H3") ~ "H3",
        str_detect(external_gene_name, "^H4") ~ "H4",
        TRUE ~ "Unknown")) %>%
      mutate(histone_class = case_when(
        histone_type == "H2A" ~ 1,
        histone_type == "H2B" ~ 2,
        histone_type == "H3" ~ 3,
        histone_type == "H4" ~ 4,
        TRUE ~ 0))
    proteasome_df <- u2os_prot_groups$proteasome_df %>%
      mutate(proteasome_type = case_when(
        str_detect(description, regex("mitochondrial", ignore_case = TRUE)) ~ "Mitochondrial",   # mitochondrial proteasome-related
        str_detect(external_gene_name, "^PSMA") | str_detect(external_gene_name, "^PSMB") ~ "20S core",   # 20S core particle (PSMA / PSMB)
        str_detect(external_gene_name, "^PSMC") | str_detect(external_gene_name, "^PSMD") ~ "19S regulatory",   # 19S regulatory particle (PSMC / PSMD)
        TRUE ~ "Unknown")) %>%
      mutate(proteasome_class = case_when(
        proteasome_type == "20S core" ~ 1,
        proteasome_type == "19S regulatory" ~ 2,
        proteasome_type == "Mitochondrial" ~ 3,
        TRUE ~ 0))
    all_annotated <- u2os_prot_groups$annotations_df
    
  # attr(data, "histones_df") <- histones_df

  u2os_score <- all_annotated %>%
    filter(!str_detect(Protein.Group, ";")) %>% 
    get_score(score_name = score_name, protein_tp) %>% 
    filter(!is.na(int_score)) %>% 
    condition_name() %>%
    group_by(condition, Protein.Group) %>%
    summarise(median_int_score = median(int_score)) %>%   # changed to median (also in plot and on y axis and in all boxplots if i show them)
    ungroup() %>%
    label_conditions() %>%
    mutate(highlight = case_when(
      Protein.Group %in% ribosomes_df$Protein.Group ~ 1,
      Protein.Group %in% histones_df$Protein.Group ~ 2,
      Protein.Group %in% proteasome_df$Protein.Group ~ 3,
      TRUE ~ 0)) %>%
    mutate(ribosome_class = ribosomes_df$ribosome_class[match(Protein.Group, ribosomes_df$Protein.Group)]) %>%
    mutate(histone_class = histones_df$histone_class[match(Protein.Group, histones_df$Protein.Group)]) %>% 
    mutate(proteasome_class = proteasome_df$proteasome_class[match(Protein.Group, proteasome_df$Protein.Group)])

  attr(data, "u2os_score") <- u2os_score

  attr(data, "u2os_score_plot") <- ggplot(u2os_score, aes(condition_label, log10(median_int_score))) +
    # scale_x_discrete(limits = c("0.2ng U2OS", "0.8ng U2OS", "4.2ng U2OS", "21 ng U2OS", "42ng U2OS", "85ng U2OS")) +
    geom_violin(color = "black") +
    # geom_boxplot(width = 0.5) +
    # stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 2)), vjust = -0.5, size = 3) +
    # ylim((4), (12)) +
    # geom_boxplot(data = subset(u2os_score, (highlight == 1 & ribosome_class == 1)), color = "seagreen3", width = 0.05) +   # ribosomes 40S
    # geom_boxplot(data = subset(u2os_score, (highlight == 1 & ribosome_class == 2)), color = "darkgreen", width = 0.05) + 
    # geom_boxplot(data = subset(u2os_score, (highlight == 1 & ribosome_class == 3)), color = "palegreen1", width = 0.05) +   # ribosomes 60S
    # geom_boxplot(data = subset(u2os_score, (highlight == 2 & histone_class == 1)), color = "darkred", width = 0.05) +   # histones
    # geom_boxplot(data = subset(u2os_score, (highlight == 3 & proteasome_class == 1)), color = "dodgerblue3", width = 0.05) +   # proteasome 20S core
    # geom_boxplot(data = subset(u2os_score, (highlight == 3 & proteasome_class == 2)), color = "darkblue", width = 0.05) +   # proteasome 19S regulatory
    # geom_boxplot(data = subset(u2os_score, (highlight == 3 & proteasome_class == 3)), color = "lightblue2", width = 0.05) +   # proteasome Mitochondrial
    theme_large() +
    theme(axis.text.x = element_text(size = 11)) +   # x labels for sam are too big
    labs(x = "condition", y = paste0("log10 ", score_name),
         title = paste0(score_name, " of U2OS proteins (log10)"),
         # subtitle = "sore based on Ms1.Normalised",
         caption = paste0("log10(u2os_", score_name, ")", "_", MS_name))
  
  return(data)}





# part 4 ----
u2os_absolut <- function(data, standard = get("standard", envir = parent.frame()), condition_filter = get("condition_filter", envir = parent.frame()), ups_factor = get("ups_factor", envir = .GlobalEnv), histone_theoretical = get("histone_theoretical", envir = .GlobalEnv)) {

  # option for selecting best condition and use ups2 slope coeffitions for all runs
  #  if (standard == "UPS2") {
  #    ups2_score_info <- attr(data, "ups2_score_info")
  #    slope <- ups2_score_info[["estimate_ups_amount"]][ups2_score_info[["condition"]] == condition_filter]
  #    intercept <- ups2_score_info[["estimate_(Intercept)"]][ups2_score_info[["condition"]] == condition_filter]
  # } else if  (standard == "UPS3") {
  #   ups3_score_info <- attr(data, "ups3_score_info")
  #   slope <- ups3_score_info[["estimate_ups_amount"]][ups3_score_info[["condition"]] == condition_filter]
  #   intercept <- ups3_score_info[["estimate_(Intercept)"]][ups3_score_info[["condition"]] == condition_filter]
  # } else {print("choose one standard")}
  
  # option for using ups2 slope coefficients of each condition (with left_join before absolute amount calculation)
  if (standard == "UPS2") {
    ups_score_info <- attr(data, "ups2_score_info")
    ups_score <- attr(data, "ups2_score")
  } else if  (standard == "UPS3") {
    ups_score_info <- attr(data, "ups3_score_info")
    ups_score <- attr(data, "ups3_score")
  } else {print("choose one standard")}
  
  anz_condition_u2os <- attr(data, "u2os_score") %>%
    distinct(condition, Protein.Group) %>%  # unique u2os proteins 
    group_by(condition) %>%
    summarise(u2os_proteins = n_distinct(Protein.Group), .groups = "drop") %>% 
    ungroup()
  
  anz_condition_ups2 <- ups_score %>%
    distinct(condition, Protein.Group) %>%  # unique ups proteins 
    group_by(condition) %>%
    summarise(ups2_proteins = n_distinct(Protein.Group), .groups = "drop")%>% 
    ungroup()
  
  anz_condition_filtered <- left_join(anz_condition_u2os, anz_condition_ups2, by = "condition") %>%     # combine u2os and ups2 & return
    mutate(condition = factor(condition, levels = conditions)) %>%
    arrange(condition)
  
  attr(data, "anz_condition_filtered") <- anz_condition_filtered
  
  # option for calculating absolute amount from coeffitions calculated with log10 data
  # u2os_absolut_fmol <- attr(data, "u2os_score") %>%  
  #   condition_name() %>%
  #   label_conditions() %>% 
  #   mutate(log10_score = log10(median_int_score)) %>%   # based on median score between replicates 
  #   mutate(x = (log10_score - intercept) / slope) %>%   # log10 amount
  #   mutate(amount = 10^x) %>% 
  #   mutate(absolute = amount/ups_factor) %>% 
  #   mutate(absolute_2 = format(absolute, scientific = FALSE)) 
  

  # option for stayinf in log space
  u2os_absolut_fmol <- attr(data, "u2os_score") %>%
    left_join(ups_score_info %>%
                select(condition,
                       slope = estimate_ups_amount, intercept = `estimate_(Intercept)`),
              by = "condition") %>%
    mutate(log10_score = log10(median_int_score)) %>%   # based on median score between replicates
    mutate(x = (log10_score - intercept) / slope) %>%   # log10 amount
    mutate(amount = 10^x) %>%
    mutate(absolute = amount/ups_factor) %>%
    mutate(absolute_2 = format(absolute, scientific = FALSE))
  
  # option for staying in lin space
  # u2os_absolut_fmol <- attr(data, "u2os_score") %>%
  #   left_join(ups_score_info %>%
  #               select(condition,
  #                      slope = estimate_ups_amount, intercept = `estimate_(Intercept)`),
  #             by = "condition") %>%
  #   mutate(amount = (median_int_score - intercept) / slope) %>%  # changed to median of replicates
  #   mutate(absolute = amount/ups_factor) %>%
  #   mutate(absolute_2 = format(absolute, scientific = FALSE))
  # 

  attr(data, "u2os_absolut_fmol") <- u2os_absolut_fmol
  
  # get one intensity for each core histone (summed isoformes)
  core_histone_sum_fmol <- u2os_absolut_fmol %>%
    filter(histone_class > 0) %>% 
    # filter(histone_class == 1 | histone_class == 2 | histone_class == 3) %>% 
    group_by(condition_label, histone_class) %>%   # Run, if i look at replicates
    summarise(absolute = sum(absolute, na.rm = TRUE)) %>%
    ungroup() %>% 
    mutate(histone_label = case_when(histone_class == 1 ~ "H2A", histone_class == 2 ~ "H2B", histone_class == 3 ~ "H3", histone_class == 4 ~ "H4", TRUE ~ NA_character_))
  
  attr(data, "u2os_absolut_fmol_plot") <- ggplot(u2os_absolut_fmol, aes(condition_label, log10(absolute))) +
    # scale_x_discrete(limits = c("0.2ng U2OS", "0.8ng U2OS", "4.2ng U2OS", "21 ng U2OS", "42ng U2OS", "85ng U2OS")) +
    geom_violin(color = "black") +
    geom_hline(yintercept = log10(histone_theoretical)) +
    # ylim((-3), (4)) +
    #theme(title = element_text(size = 14), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
    geom_boxplot(data = subset(u2os_absolut_fmol, (highlight == 1 & ribosome_class == 1)), color = "seagreen3", width = 0.05) +   # ribosomes 40S
    geom_boxplot(data = subset(u2os_absolut_fmol, (highlight == 1 & ribosome_class == 2)), color = "darkgreen", width = 0.05) +
    geom_boxplot(data = subset(u2os_absolut_fmol, (highlight == 1 & ribosome_class == 3)), color = "palegreen1", width = 0.05) +   # ribosomes 60S
    # geom_boxplot(data = subset(u2os_absolut_fmol, (highlight == 2 & histone_class == 1)), color = "darkred", width = 0.05) +   # histones
    geom_boxplot(data = subset(u2os_absolut_fmol, (highlight == 3 & proteasome_class == 1)), color = "dodgerblue3", width = 0.05) +   # proteasome 20S core
    geom_boxplot(data = subset(u2os_absolut_fmol, (highlight == 3 & proteasome_class == 2)), color = "darkblue", width = 0.05) +   # proteasome 19S regulatory
    geom_boxplot(data = subset(u2os_absolut_fmol, (highlight == 3 & proteasome_class == 3)), color = "lightblue2", width = 0.05) +   # proteasome Mitochondrial
    # geom_text_repel(data = subset(u2os_absolut_fmol, highlight == 2 & histone_class == 1), aes(label = Protein.Group), vjust = -1.5, hjust = 1, color = "darkred", size = 2, max.overlaps = Inf) +   # max.overlaps = Inf allows many labels if needed
    geom_point(data = core_histone_sum_fmol, aes(condition_label, log10(absolute)), color = "darkred", size = 3) +
    geom_text_repel(data = core_histone_sum_fmol, aes(condition_label, log10(absolute), label = histone_label), color = "darkred", size = 4) +
    theme_large() +
    labs(x = "condition", y = paste0("log10 protein amount [fmol]"),
         title = paste0("absolute amount of U2OS proteins (log10(fmol))"), subtitle = "red - histones, green - ribosomes, blue - proteasome",
        caption = paste0("log10(u2os_amount_fmol)_", score_name, "_", MS_name))

  u2os_iqr_fmol <- u2os_absolut_fmol %>%
    group_by(condition, highlight, ribosome_class, histone_class, proteasome_class) %>%
    summarise(n = n(),
      Q1 = quantile(log10(absolute), 0.25, na.rm = TRUE),
      Median = median(log10(absolute), na.rm = TRUE),
      Q3 = quantile(log10(absolute), 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      .groups = "drop")
  attr(data, "u2os_iqr_fmol") <- u2os_iqr_fmol
  
  histone_iqr_fmol <- core_histone_sum_fmol %>%
    group_by(condition_label) %>% 
    summarise(
      n = n(),
      Q1 = quantile(log10(absolute), 0.25, na.rm = TRUE),
      Median = median(log10(absolute), na.rm = TRUE),
      Q3 = quantile(log10(absolute), 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      .groups = "drop")
  attr(data, "histone_iqr_fmol") <- histone_iqr_fmol
  attr(data, "core_histone_sum_fmol") <- core_histone_sum_fmol
  
  
  
  # option for using the slope
  u2os_absolut_pg <- attr(data, "u2os_score") %>%
    left_join(ups_score_info %>%
                select(condition,
                       slope = estimate_ups_amount, intercept = `estimate_(Intercept)`),
              by = "condition") %>%
    mutate(log10_score = log10(median_int_score)) %>%   # based on median score between replicates
    mutate(x = (log10_score - intercept) / slope) %>%   # log10 amount
    mutate(tCN = 10^x) %>%
    # mutate(protein_mw = u2os_mw$protein_mw[match(Protein.Group, u2os_mw$Protein.Group_short)]) %>%   # new column with mw
    left_join(u2os_mw_longest %>% 
        select(Protein.Group_short, protein_mw), by = c("Protein.Group" = "Protein.Group_short")) %>% 
    mutate(tPM = tCN * protein_mw) %>% 
    group_by(condition) %>%   # changed to condition when i use median
    mutate(tPM_total = sum(tPM)) %>%   # changed to median 
    ungroup() %>%
    mutate(tPM_proportion = tPM / tPM_total) %>%
    mutate(absolute = 300 * tPM_proportion) %>% 
    mutate(copy_number = absolute / protein_mw) 
  
  attr(data, "u2os_absolut_pg") <-u2os_absolut_pg
  
  # MW matching check  
  # Wie viele Protein-Gruppen hatten mehrere Isoformen?
  u2os_mw %>%
    count(Protein.Group_short) %>%
    filter(n > 1)
  # Wurden alle Proteine gematcht?
  n_missing_mw <- sum(is.na(u2os_absolut_pg$protein_mw))
  if (n_missing_mw != 0) {
    warning(sprintf(n_missing_mw, " Proteins don't hava a assigned MW"), call. = FALSE)}
  
  # u2os_absolut_pg <- attr(data, "u2os_score") %>%
  #   group_by(condition) %>%   # changed to condition when i use median
  #   mutate(score_total = sum(int_score)) %>%   # changed to median 
  #   ungroup() %>% 
  #   mutate(score_proportion = int_score / score_total) %>%    # changed to median
  #   mutate(absolute = 300 * score_proportion) %>%    # absolute amount in pg per cell
  #   mutate(protein_mw = u2os_mw$protein_mw[match(Protein.Group, u2os_mw$Protein.Group_short)]) %>%   # new column with mw
  #   mutate(copy_number = absolute / protein_mw)   # protein copy number = pg per cell / mw
  
  attr(data, "u2os_absolut_pg") <- u2os_absolut_pg
  
  # get one intensity for each core histone (summed isoformes)
  core_histone_sum_pg <- u2os_absolut_pg %>%
    filter(histone_class > 0) %>% 
    # filter(histone_class == 1 | histone_class == 2 | histone_class == 3) %>% 
    group_by(condition_label, histone_class) %>%   # Run, if i look at replicates
    summarise(absolute = sum(absolute, na.rm = TRUE)) %>%
    ungroup() %>% 
    mutate(histone_label = case_when(histone_class == 1 ~ "H2A", histone_class == 2 ~ "H2B", histone_class == 3 ~ "H3", histone_class == 4 ~ "H4", TRUE ~ NA_character_))
  
  
  attr(data, "u2os_absolut_pg_plot") <- ggplot(u2os_absolut_pg, aes(condition_label, log10(absolute))) +
    # scale_x_discrete(limits = c("0.2ng U2OS", "0.8ng U2OS", "4.2ng U2OS", "21 ng U2OS", "42ng U2OS", "85ng U2OS")) +
    geom_violin(color = "black") +
    # ylim((-3), (4)) +
    #theme(title = element_text(size = 14), axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10), axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16)) +
    geom_boxplot(data = subset(u2os_absolut_pg, (highlight == 1 & ribosome_class == 1)), color = "seagreen3", width = 0.05) +   # ribosomes 40S
    geom_boxplot(data = subset(u2os_absolut_pg, (highlight == 1 & ribosome_class == 2)), color = "darkgreen", width = 0.05) +
    geom_boxplot(data = subset(u2os_absolut_pg, (highlight == 1 & ribosome_class == 3)), color = "palegreen1", width = 0.05) +   # ribosomes 60S
    # geom_boxplot(data = subset(u2os_absolut_pg, (highlight == 2 & histone_class == 1)), color = "darkred", width = 0.05) +   # histones
    geom_boxplot(data = subset(u2os_absolut_pg, (highlight == 3 & proteasome_class == 1)), color = "dodgerblue3", width = 0.05) +   # proteasome 20S core
    geom_boxplot(data = subset(u2os_absolut_pg, (highlight == 3 & proteasome_class == 2)), color = "darkblue", width = 0.05) +   # proteasome 19S regulatory
    geom_boxplot(data = subset(u2os_absolut_pg, (highlight == 3 & proteasome_class == 3)), color = "lightblue2", width = 0.05) +   # proteasome Mitochondrial
    # geom_text_repel(data = subset(u2os_absolut_pg, highlight == 2 & histone_class == 1), aes(label = Protein.Group), vjust = -1.5, hjust = 1, color = "darkred", size = 2, max.overlaps = Inf) +   # max.overlaps = Inf allows many labels if needed
    geom_point(data = core_histone_sum_pg, aes(condition_label, log10(absolute)), color = "darkred", size = 3) +
    geom_text_repel(data = core_histone_sum_pg, aes(condition_label, log10(absolute), label = histone_label), color = "darkred", size = 4) +
    theme_large() +
    labs(x = "conditon", y = paste0("log10 protein amount [pg]"),
         title = paste0("absolute amount of U2OS proteins (log10(pg))"), subtitle = "red - histones, green - ribosomes, blue - proteasome",
         caption = paste0("log10(u2os_amount_pg)_", score_name, "_", MS_name))
  
  
  u2os_iqr_pg <- u2os_absolut_pg %>%
    group_by(condition, highlight, ribosome_class, histone_class, proteasome_class) %>%
    summarise(n = n(),
              Q1 = quantile(log10(absolute), 0.25, na.rm = TRUE),
              Median = median(log10(absolute), na.rm = TRUE),
              Q3 = quantile(log10(absolute), 0.75, na.rm = TRUE),
              IQR = Q3 - Q1,
              .groups = "drop")
  attr(data, "u2os_iqr_pg") <- u2os_iqr_pg
  
  histone_iqr_pg <- core_histone_sum_pg %>%
    group_by(condition_label) %>% 
    summarise(
      n = n(),
      Q1 = quantile(log10(absolute), 0.25, na.rm = TRUE),
      Median = median(log10(absolute), na.rm = TRUE),
      Q3 = quantile(log10(absolute), 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      .groups = "drop")
  attr(data, "histone_iqr_pg") <- histone_iqr_pg
  attr(data, "core_histone_sum_pg") <- core_histone_sum_pg
  
  return(data)}

# BM ----

ups2_peptides <- function(data) {
  peptides_fc <- data %>% 
    mutate(ups2 = Sequence %in% ups2_speclib$PeptideSequence) %>%   # look for all ups2 sequences (H and L)
    filter(ups2 == TRUE) %>% 
    label_conditions() %>% 
    
    # find H and L pairs for peptides with ups2 sequence
    rowwise() %>%
    filter(length(unique(Channel)) == 1) %>%   # only single labelled here
    ungroup() %>% 
    group_split(Run) %>%  # Split into list of data frames by 'Run'
    map(~ {
      peptides_h <- filter(.x, grepl("H", Channel))
      peptides_l <- filter(.x, grepl("L", Channel)) %>% 
        left_join(peptides_h,
                  by = "Stripped.Sequence.Charge",   # H and L pairs based on stripped sequence charge
                  suffix = c(".l", ".h"))
    }) %>%
    bind_rows() %>% 
    
    # calculate log2 FC
    mutate(log2_MS1_L = log2(intensity.l), log2_MS1_H = log2(intensity.h)) %>%   # intensities in log2
    # mutate(Ms1.Translated.h = if_else(is.na(Ms1.Translated.h), Channel.H.l, Ms1.Translated.h)) %>%   # if only light found - use Channel H value
    # mutate(log2_MS1_L = log2(Ms1.Translated.l), log2_MS1_H = log2(Ms1.Translated.h)) %>% 
    mutate(log2_MS1_Ratio = log2_MS1_L - log2_MS1_H) %>%   # calculate ratios
    mutate(across(c(log2_MS1_Ratio, log2_MS1_L), ~ ifelse(is.finite(.), ., NA))) %>%
    
    # some extra information on peptide pairs
    mutate(protein = sub("ups.*", "", Protein.Group.l)) %>% 
    mutate(ups2_amount = ups2_standard$`UPS2 Amount (fmol)`[match(protein, ups2_standard$`UniProt Accession Number`)],
           MW = ups2_standard$`MW (Da)`[match(protein, ups2_standard$`UniProt Accession Number`)],
           Source = ups2_standard$`Source or recombinant`[match(protein, ups2_standard$`UniProt Accession Number`)],
           Host = ups2_standard$Host[match(protein, ups2_standard$`UniProt Accession Number`)],
           Tag = ups2_standard$Tag[match(protein, ups2_standard$`UniProt Accession Number`)],
           PTMs = ups2_standard$`Potential PTMs`[match(protein, ups2_standard$`UniProt Accession Number`)]) %>%   # add the protein amount from standard
    mutate(prot_amount = str_c(protein, " (", ups2_amount, ")")) %>% 
    group_by(Protein.Group.l) %>%   # start counting for every protein new
    mutate(color_prot_pep = dense_rank(Sequence.l))   # for each protein: same sequence - same color (over all runs)
  
  attr(data, "peptides_fc") <- peptides_fc
  
  # FC plot
  attr(data, "peptides_fc_plot") <- ggplot(peptides_fc, aes(factor(condition_label.l, levels = condition_labels), 
                                                            log2_MS1_Ratio)) +
    geom_boxplot() +
    geom_hline(yintercept = 0) +
    stat_summary(fun = median, geom = "text", aes(label = round(after_stat(y), 2)), size = 4, vjust = -0.6, color = "red") +   # show mean
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +   # angle for runs on x axis 
    theme_small() +
    labs(x = "Condition", y = "log2(L/H)",
         title = "Peptide ratios (UPS2/U2OS) in each condition", caption = "pep_ratio")
  
  
  # get median ratio (here MS1) for each peptide in one condition
  peptides_fc_condition <- peptides_fc %>% 
    group_by(condition_label.l, Sequence.l) %>%
    summarise(peptide_ratio = median(c(log2_MS1_Ratio), na.rm = TRUE),
              log2_MS1_L_median = median(log2_MS1_L, na.rm = TRUE),
              log2_MS1_H_median = median(log2_MS1_H, na.rm = TRUE),
              across(c(protein, ups2_amount), first),
              .groups = "drop")
  
  attr(data, "peptides_fc_condition") <- peptides_fc_condition
  
  # intensity over ups2 amount plot
  attr(data, "peptides_int_plot") <- ggplot(peptides_fc_condition, aes(peptide_ratio, log2_MS1_L_median, color = factor(ups2_amount))) +
    geom_point() +
    facet_wrap(~ factor(condition_label.l, levels = condition_labels)) + # plot for each condition
    theme_minimal() +
    scale_color_manual(values = c(
      "0.5" = "lightblue",
      "5" = "lightskyblue",
      "50" = "royalblue2",
      "500" = "royalblue4",
      "5000" = "navyblue",
      "50000" = "black"
    ), name = "UPS2 amount [fmol]") +
    theme_small() +
    labs(x = "log2(L/H)", y = "log2 UPS2 int", title = "Peptide intensity over ratio", 
         subtitle = "peptide intensity: median of replicates", caption = "pep_int_ratio_amount")
  
  # pep_prot_ratio_cond (UPS2 Proteine)
  attr(data, "ups2_proteins_plot") <- ggplot(peptides_fc, aes(factor(condition_label.l, levels = condition_labels), log2_MS1_Ratio,
                             color = factor(color_prot_pep))) +
      geom_point(size = 1.2) +
      facet_wrap(~prot_amount) + # plot for each protein
      theme_minimal() +
      scale_x_discrete(labels = c("200ng_HSH" = "HSH", "HSH_55ng_UPS2" = "55 ng", "HSH_110ng_UPS2" = "110 ng", "HSH_220ng_UPS2" = "220 ng")) +   # rename conditions
    theme_small() +  
    labs(x = "UPS2 Concentraions", y = "median log2(L/H) (PT and MS1)",
           title = "ratios over UPS2 concentrations (peptide level)", subtitle = "colored by same sequence over all runs", caption = "pep_ratio_protein") +
      scale_color_manual(values = my_colors) +   # color palette
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(legend.position = "none")


  

  
  return(data)}

