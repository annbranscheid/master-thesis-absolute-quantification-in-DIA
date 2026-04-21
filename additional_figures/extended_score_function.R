# load the function files, experiment, load-data and data-overview before

get_score <- function(df, score_name = "topN", protein_tp){
  if (score_name == "topN") {
    df <- df %>%
      group_by(Run, Sequence, Protein.Group) %>%
      summarise(topN_intensity.sequence = sum(intensity), anz = n(), .groups = "drop") %>%
      ungroup() %>%
      group_by(Run, Protein.Group) %>%
      filter(n_distinct(Sequence) >= 2) %>%   # rn top 2-3
      slice_max(order_by = topN_intensity.sequence, n = 3, with_ties = FALSE) %>%
      summarise(topN_intensity = sum(topN_intensity.sequence), anz = n(), .groups = "drop") %>%
      mutate(int_score = topN_intensity / anz)
  } else if (score_name == "top3") {
    df <- df %>%
      group_by(Run, Sequence, Protein.Group) %>%
      summarise(topN_intensity.sequence = sum(intensity), anz = n(), .groups = "drop") %>%
      ungroup() %>%
      group_by(Run, Protein.Group) %>%
      filter(n_distinct(Sequence) >= 3) %>% 
      slice_max(order_by = topN_intensity.sequence, n = 3, with_ties = FALSE) %>%
      summarise(top3_intensity = sum(topN_intensity.sequence), anz = n(), .groups = "drop") %>%
      mutate(int_score = top3_intensity / anz)
  } else if (score_name == "iBAQ") {
    df <- df %>%
      group_by(Run, Sequence, Protein.Group) %>%
      summarise(iBAQ_intensity.sequence = sum(intensity), anz = n(), .groups = "drop") %>%
      ungroup() %>%
      group_by(Run, Protein.Group) %>%
      filter(n_distinct(Sequence) >= 1) %>%
      summarise(iBAQ_intensity = sum(iBAQ_intensity.sequence), anz  = n(), .groups = "keep") %>%   # using sum of all precursor rn
      ungroup() %>% 
      mutate(Protein.Group_short = str_extract(Protein.Group, "^[^|]+")) %>%
      left_join(select(protein_tp, Protein.Group_short, theoretical_peptides),
                by = join_by(Protein.Group_short)) %>%
      mutate(int_score = iBAQ_intensity / theoretical_peptides)
  } else if (score_name == "riBAQ") {
    df <- df %>%
      group_by(Run, Sequence, Protein.Group) %>%
      summarise(iBAQ_intensity.sequence = sum(intensity), anz = n(), .groups = "drop") %>%
      ungroup() %>%
      group_by(Run, Protein.Group) %>%
      filter(n_distinct(Sequence) >= 1) %>%
      summarise(iBAQ_intensity = sum(iBAQ_intensity.sequence), anz  = n(), .groups = "keep") %>%   # using sum of all precursor rn
      ungroup() %>%
      mutate(Protein.Group_short = str_extract(Protein.Group, "^[^|]+")) %>%
      left_join(select(protein_tp, Protein.Group_short, theoretical_peptides),
                by = join_by(Protein.Group_short)) %>%
      mutate(iBAQ = iBAQ_intensity / theoretical_peptides) %>% 
      group_by(Run) %>%
      mutate(sum_iBAQ = sum(iBAQ, na.rm = TRUE),  # sum of iBAQ within each Run
        int_score = iBAQ / sum_iBAQ) %>%
      ungroup()
  } else if (score_name == "iBAQ_noise") {
    df <- df %>%
      {noise <<- min(.$intensity[.$intensity != 0])
      .} %>%
      group_by(Run, Sequence, Protein.Group) %>%
      summarise(iBAQ_intensity.sequence = sum(intensity), anz = n(), .groups = "drop") %>%
      ungroup() %>%
      group_by(Run, Protein.Group) %>%
      filter(n_distinct(Sequence) >= 1) %>%
      summarise(iBAQ_intensity = sum(iBAQ_intensity.sequence), anz  = n(), .groups = "keep") %>%   # using sum of all precursor rn
      ungroup() %>% 
      mutate(Protein.Group_short = str_extract(Protein.Group, "^[^|]+")) %>%
      left_join(select(protein_tp, Protein.Group_short, theoretical_peptides),
                by = join_by(Protein.Group_short)) %>%
      mutate(int_score =  (iBAQ_intensity + (noise * (theoretical_peptides-anz)))/ theoretical_peptides)
  } else if (score_name == "mean_int") {
    df <- df %>%
      group_by(Run, Protein.Group) %>%
      filter(n_distinct(Sequence) >= 2) %>%
      summarise(sum_intensity = sum(intensity), anz  = n(), .groups = "drop") %>%   # using sum of all precursor rn
      mutate(int_score = sum_intensity / anz)
  } else {
    stop("score is not in the list")
  }
}



score_name <- "iBAQ"
score_name <- "topN"
score_name <- "top3"
score_name <- "mean_int"
score_name <- "riBAQ"
score_name <- "iBAQ_noise"

{data <- data %>% 
  ups2_score(., mean_score = "no")  
 
save_table(attr(data, "ups2_score"), paste0("ups2_", score_name, "_scores"))
save_table(attr(data, "ups2_score_info"), paste0("ups2_", score_name, "_slope_info"))
print(attr(data, "ups2_score_plot"))   
ggsave_png()
}
