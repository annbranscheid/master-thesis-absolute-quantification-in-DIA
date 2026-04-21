library(tidyverse)
library(usethis)
#library(devtools)
library(conflicted)
# library(Polychrome)   # for color palettes
library(broom)   # for slopes
library(arrow)   # to read parquet files
library(biomaRt)   # annotation
library(ggrepel)   # prevent labels from overlapping
library(Biostrings)   # for reading FASTA 
library(cleaver)   # for in-silico digestion
library(Peptides)   # for MW function
conflicted::conflicts_prefer(
  dplyr::filter,
  dplyr::select,
  dplyr::rename,
  dplyr::slice,
  dplyr::summarise,
  dplyr::mutate,
  dplyr::arrange,
  dplyr::desc,
  dplyr::first
)

# look at the packages and verions used
# library(renv)
# renv::dependencies("functions.R")
# renv::init()
# renv::snapshot()
# deps <- renv::dependencies("functions.R")$Package
# lock <- jsonlite::fromJSON("renv.lock")
# pkg_names <- names(lock$Packages)
# data.frame(Package = deps,
#            Version = sapply(deps, function(p) {
#              if (p %in% pkg_names) {
#                lock$Packages[[p]]$Version
#                } else {NA_character_}}))


remove_conditions <- c()

theme_small <- function(base_size = 11) {
  theme_grey(base_size = base_size) +
    theme(
      plot.title    = element_text(size = 1.2 * base_size, face = "bold"),
      plot.subtitle = element_text(size = 1 * base_size),
      axis.title    = element_text(size = 0.9 * base_size),
      axis.text.x   = element_text(size = 0.7 * base_size),
      axis.text.y   = element_text(size = 0.8 * base_size),
      plot.caption  = element_text(size = 8, hjust = 1),
      # Facet strips (if i want to change the facet title background)
      # strip.background = element_rect(fill  = "grey85",colour = "grey50"),
      # strip.text = element_text(face = "bold", size = 0.9 * base_size)
    )
}

theme_large <- function(base_size = 16) {
  theme_grey(base_size = base_size) +
    theme(
      plot.title    = element_text(size = 1.2 * base_size, face = "bold"),
      plot.subtitle = element_text(size = 1 * base_size),
      axis.title    = element_text(size = 0.9 * base_size),
      axis.text.x   = element_text(size = 0.7 * base_size),
      axis.text.y   = element_text(size = 0.8 * base_size),
      plot.caption  = element_text(size = 8, hjust = 1),
      # Facet strips (if i want to change the facet title background)
      # strip.background = element_rect(fill  = "grey85",colour = "grey50"),
      # strip.text = element_text(face = "bold", size = 0.9 * base_size)
    )
}

my_colors <- c(
  "dodgerblue", "firebrick", "purple", "salmon", "slateblue", "darkred",
  "steelblue", "rosybrown", "mediumorchid", "forestgreen", "mediumseagreen", "deeppink",
  "saddlebrown", "cadetblue", "orangered", "seagreen", "peru", "tomato",
  "paleturquoise", "gray40", "lightcoral", "maroon", "blue", "olivedrab",
  "indianred", "navy", "lightsteelblue", "mediumpurple", "coral", "green",
  "plum", "slategray", "darkcyan", "orange", "mediumturquoise", "darkgreen",
  "thistle", "darkorange", "sienna", "chocolate", "midnightblue", "mediumvioletred",
  "turquoise", "red", "royalblue", "darkviolet", "brown", "mediumblue"
)

load_data <- function(file_path) {
  ext <- tools::file_ext(file_path)
  if (ext == "tsv") {
    df <- read_tsv(file_path, show_col_types = FALSE)
  } else if (ext %in% c("xls", "xlsx")) {
    df <- readxl::read_excel(file_path, show_col_types = FALSE)
  } else if (ext == "parquet") {
    df <- read_parquet(file_path, show_col_types = FALSE)
  } else if (ext == "txt") {
    df <- read_tsv(file_path, show_col_types = FALSE)
  } else {
    stop("Unsupported file type: ", ext)
  }
  
  return(df)
}

ggsave_png <- function(name = last_plot()$labels$caption) {
  ggsave(
    # filename <- str_c(path, "/results_thesis/", format(Sys.Date(), "%Y%m%d"), "_", name, ".png"),   # to save in the same folder as data
    filename <- str_c(results_path, "/", format(Sys.Date(), "%Y%m%d"), "_", name, ".png"),   # to save results of all MS in one folder
    width = 8.5,
    height = 7.53,
    dpi = 300,
    units = "in",
    device = "png"
  )
}

save_table <- function(df, name = NULL) {
  if(is.null(name)) {
    name <- deparse(substitute(df))
    # Clean up any special characters for filenames
    name <- gsub("[^A-Za-z0-9_]", "_", name)
  }
  write_tsv(df,
            # filename <- str_c(path, "/results_thesis/", format(Sys.Date(), "%Y%m%d"), "_", name, ".tsv")   # to save in the same folder as data
            filename <- str_c(results_path, "/", format(Sys.Date(), "%Y%m%d"), "_", name, ".tsv")   # to save results of all MS in one folder
            )
}

# run name function
run_name <- function(df, col = "Run", start_from = 1, end_from_end = 0) {
  df %>%
    mutate({{col}} := sapply(str_split(.data[[col]], "_"), function(parts) {
        n <- length(parts)
        end_pos <- if (end_from_end == 0) n else n - end_from_end   # if 0, don't remove anything
        if(n >= start_from - end_pos) {   # string long enough?
          str_c(parts[start_from:end_pos], collapse = "_")
        } else {NA}}))}

condition_name <- function(df, col = "Run", MS = get("MS", envir = parent.frame())) {
  cut_from_end <-
    if (MS == "olive" || MS == "oliveslice" || MS == "popeye") {2} else if (MS == "fozzie" || MS == "sam") {1} else {0}
  df %>%
    mutate(condition = sapply(str_split(.data[[col]], "_"), function(parts) {
      n <- length(parts)# position to keep up to (from the left)
      end_pos <- n - cut_from_end
      if (end_pos >= 1) {str_c(parts[1:end_pos], collapse = "_")} else {NA}}))}

slopes_score <- function(df, data_log10 = "yes", ups3 = "no") {
  slopes <- df  %>%
    { if (ups3 == "no") {mutate(., ups_amount = ups2_amount)} else . } %>%
    { if (ups3 == "yes") {mutate(., ups_amount = ups3_amount)} else . } %>%
    { if (data_log10 == "no") {mutate(., int_score = log10(int_score), ups_amount = log10(ups_amount))} else . } %>%
    group_by(condition) %>%
    group_modify(~ {
      df2 <- drop_na(.x, int_score, ups_amount)
      if (nrow(df2) == 0) return(tibble())  # returns empty group summary
      tidy(lm(int_score ~ ups_amount, data = df2))}) %>%
    #filter(term == "ups_amount") %>%
    select(c(condition, term, estimate, std.error)) %>%
    pivot_wider(names_from = term, values_from = c(estimate, std.error)) %>%
    mutate(slope = paste0("slope = ", round(estimate_ups_amount, 2)),
           intercept = paste0("intercept = ", round(`estimate_(Intercept)`, 2)),
           std_ups = paste0("std = ", round(estimate_ups_amount, 2)))
  attr(df, "slopes") <- slopes
  return(df)
}

correlation_score <- function(df, score = int_score, data_log10 = "yes", ups3 = "no"){
  correlation_score <- df %>% #new df with correlations as labels for the plot (pearson - linear relationship between MS1 and PT expected, or spearman)
    { if (!"int_score" %in% names(.)) mutate(., int_score = {{score}}) else . } %>% 
    { if (ups3 == "no") {mutate(., ups_amount = ups2_amount)} else . } %>%
    { if (ups3 == "yes") {mutate(., ups_amount = ups3_amount)} else . } %>%
    { if (data_log10 == "no") {mutate(., int_score = log10(int_score), ups_amount = log10(ups_amount))} else . } %>%
    group_by(condition) %>%
    summarize(cor_score = cor(int_score, ups_amount, use = "complete.obs", method = c("pearson"))) %>%
    mutate(cor_score_1 = paste0("r = ", round(cor_score, 2)), cor_score_2 = paste0("r^2 = ", round((cor_score**2), 2))) 
  attr(df, "correlation_score") <- correlation_score
  return(df)
}

theoretical_peptides <- function(fasta_path, protease = "trypsin", missed_cleavages = 1, min_length = 7, max_length = 30) {
  proteins <- readAAStringSet(fasta_path)   # read FASTA
  
  # Digest all proteins
  peptide_tbl <- tibble(
    Protein.Group = names(proteins),
    Sequence  = as.character(proteins)) %>%   # tibble for proteins and sequence
    mutate(peptides = map(Sequence, ~ {
      p <- cleave(.x, enzym = protease, missedCleavages = missed_cleavages)
      # Flatten in case cleave returns a list of length 1
      if (is.list(p) && length(p) == 1) p[[1]] else p})) %>%
    unnest(peptides) %>%
    mutate(peptide_length = nchar(peptides)) %>%
    filter(peptide_length >= min_length, peptide_length <= max_length)
  
  # Count peptides per protein
  peptide_tbl %>%
    group_by(Protein.Group) %>%
    summarise(theoretical_peptides = n(), .groups = "drop")
}

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

find_ribosome_histone_proteasome <- function(protein_df, id_column = "Protein.Group") {
  # Extract UniProt IDs, removes duplicates, NAs, or empty strings
  flat_ids <- protein_df %>%
    select(all_of(id_column)) %>%   # take uniprot name column
    mutate(ProteinID = strsplit(as.character(.data[[id_column]]), ";")) %>%   # shouldn't be any in there
    unnest(ProteinID) %>%   # split multiple IDs per row into separate rows
    mutate(ProteinID = str_trim(ProteinID)) %>%
    filter(ProteinID != "") %>% 
    distinct(ProteinID)   # remove dublicates from replicate runs
  
  unique_ids <- flat_ids$ProteinID   # IDs as a list
  mart <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")   # Connect to Ensembl biomart, remove mirror if ensemble works again   , mirror = "asia"
  annotations <- getBM(attributes = c("uniprotswissprot", "external_gene_name", "description"),
                       filters = "uniprotswissprot", values = unique_ids, mart = mart)   # get annotations from BioMart
  
  # Ensure one row per ProteinID
  annotations <- annotations %>%
    distinct(uniprotswissprot, .keep_all = TRUE)
  
  # Merge annotations with protein list from data
  annotated_ids <- flat_ids %>%
    left_join(annotations, by = c("ProteinID" = "uniprotswissprot"))
  
  # Join back to original data (replicates preserved)
  annotated_df <- protein_df %>%
    #distinct(ProteinID, .keep_all = TRUE) %>%
    left_join(annotated_ids, by = setNames("ProteinID", id_column))
  
  stopifnot(nrow(annotated_df) == nrow(protein_df))
  # annotated_df <- protein_df %>%
  #   mutate(
  #     match_idx = match(.data[[id_column]], annotated_ids$ProteinID),
  #     annotation = annotated_ids$annotation[match_idx]
  #   ) %>%
  #   select(-match_idx)
  
  # Filter for ribosomal proteins
  ribosomes_df <- annotated_df %>%
    filter(str_detect(description, regex("ribosomal", ignore_case = TRUE)))
  histones_df <- annotated_df %>%
    filter(str_detect(description, regex("histone", ignore_case = TRUE)))
  proteasome_df <- annotated_df %>%
    filter(str_detect(description,regex("proteasom|ubiquitin", ignore_case = TRUE)))
  
  # Return a list of results
  list(annotations_df = annotated_df, ribosomes_df = ribosomes_df, histones_df = histones_df, proteasome_df = proteasome_df)
}


label_conditions <- function(data, conditions = get("conditions", envir = .GlobalEnv), condition_labels = get("condition_labels", envir = .GlobalEnv)) {
  stopifnot(length(conditions) == length(condition_labels), "condition" %in% names(data))
  
  # Create new column with relabeled factor
  data$condition_label <- factor(
    data$condition,
    levels = conditions,
    labels = condition_labels
  )
  
  return(data)
}

sort_by_run <- function(df, cond = conditions) {
  df <- df %>%
    mutate(Run_number = as.numeric(str_match(Run, "_(\\d+)$")[,2]),
           Run_order = match(condition, cond) * 10 + Run_number) %>% 
    arrange(Run_order) %>%
    select(- c(Run_order, Run_number))
}

normalisation <- function(df) {
  scaling_factor <- df %>% 
    group_by(condition, Run) %>% 
    summarise(median_int = median(intensity)) %>%   # mean int of each run
    ungroup()
  
  condition_ref <- scaling_factor %>% 
    group_by(condition) %>% 
    summarise(ref_int = median(median_int))   #mean run as reference for the other two
  
  scaling_factor <- scaling_factor %>% 
    mutate(ref_int = condition_ref$ref_int[match(condition, condition_ref$condition)]) %>%
    mutate(factor_run = ref_int/median_int)
  
  attr(df, "scaling_factor") <- scaling_factor
  return(df)
}
