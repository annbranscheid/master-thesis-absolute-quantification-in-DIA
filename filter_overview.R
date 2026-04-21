# filter data (contaminants and Q-vlaue) and get information about identified peptide and protein amount

if (experiment_group == "DIA1") {
  # BM: DiaNN 1.8 used with SILAC
  data <- data %>% 
    filter(!str_detect(Protein.Group, "Cont")) %>% 
    mutate(Stripped.Sequence.Charge = str_remove_all(Precursor.Id, "\\(SILAC-[KR]-[LH]\\)"), 
           Channel = str_extract_all(Precursor.Id, "(?<=-)[LH]")) %>%
    mutate(Run =  str_replace_all(Run, "pg", "ng")) %>% 
    filter(Q.Value < 0.01,  PG.Q.Value < 0.05) %>% 
    filter((Channel == "H" | Channel == "H, H" | ((Channel == "L" | Channel == "L, L") & Channel.Q.Value < 0.03))) %>% 
    # filter(!grepl("L", Channel) | (grepl("L", Channel) & Channel.Q.Value < 0.03)) %>%   # filter channel q value only for standard peptides
    run_name(start_from = 6, end_from_end = 1) %>% 
    prepare_data() %>% 
    anz_ups2()
  
} else if ((experiment_group == "DIA2" | experiment_group == "DIA3") & (MS == "olive" | MS == "oliveslice") & (specification == "" | specification == "ups2")) {
  # DiaNN 1.8 used with SILAC
  data <- data %>% 
    filter(!str_detect(Protein.Group, "Cont")) %>% 
    mutate(Stripped.Sequence.Charge = str_remove_all(Precursor.Id, "\\(SILAC-[KR]-[LH]\\)"), 
           Channel = str_extract_all(Precursor.Id, "(?<=-)[LH]")) %>%
    filter(Q.Value < 0.01,  PG.Q.Value < 0.05) %>% 
    filter((Channel == "H" | Channel == "H, H" | ((Channel == "L" | Channel == "L, L") & Channel.Q.Value < 0.03))) %>%
    run_name(start_from = 6, end_from_end = 1) %>% 
    prepare_data() %>% 
    anz_ups2()
  
} else if (experiment_group == "DIA3" & (MS == "olive" | MS == "oliveslice") & specification == "ups3") {
  # DiaNN 1.8 used with only light proteins
  data <- data %>% 
    filter(!str_detect(Protein.Group, "Cont")) %>% 
    mutate(Stripped.Sequence.Charge = str_remove_all(Precursor.Id, "\\(SILAC-[KR]-[LH]\\)"), 
           Channel = str_extract_all(Precursor.Id, "(?<=-)[LH]")) %>%
    filter(Q.Value < 0.01, PG.Q.Value < 0.05) %>% 
    run_name(start_from = 6, end_from_end = 1) %>% 
    prepare_data() %>% 
    anz_ups3()
  
} else if ((experiment_group == "DIA2" | experiment_group == "DIA3") & (MS == "fozzie" | MS == "sam") & (specification == "" | specification == "ups2")) {
  # DiaNN 2.0 used with SILAC
  data <- data %>% 
    filter(!str_detect(Protein.Group, "Cont")) %>% 
    filter(Q.Value < 0.01,PG.Q.Value < 0.05) %>% 
    filter((Channel == "H" | Channel == "H, H" | ((Channel == "L" | Channel == "L, L") & Channel.Q.Value < 0.03))) %>% 
    # filter(!grepl("L", Channel) | (grepl("L", Channel) & Channel.Q.Value < 0.03)) %>%
    run_name(start_from = 6) %>% 
    prepare_data() %>% 
    anz_ups2()
  
} else if (experiment_group == "DIA3" & (MS == "fozzie" | MS == "sam") & specification == "ups3") {
  # DiaNN 2.0 used with only light proteins
  data <- data %>%
    filter(!str_detect(Protein.Group, "Cont")) %>%
    filter(Q.Value < 0.01, PG.Q.Value < 0.05) %>%
    run_name(start_from = 6) %>%
    prepare_data() %>%
    anz_ups3()
  
} else if (experiment_group == "DDA3" & MS == "fozzie" & specification == "diann") { 
  # DiaNN 2.3 used with only light proteins
  standard <- "UPS2"
  data <- data %>% 
    filter(Precursor.Charge > 1, !str_detect(Protein.Group, "cRAP")) %>% 
    filter(Q.Value < 0.01, PG.Q.Value < 0.05) %>% 
    mutate(Run =  str_replace_all(Run, "UPS3", "UPS2")) %>%
    run_name(start_from = 6) %>%
    prepare_data() %>% 
    anz_ups2()
} else if (experiment_group == "DDA3" & MS == "sam" & specification == "diann") { 
  # DiaNN 2.3 used with only light proteins
  data <- data %>% 
    filter(Precursor.Charge > 1, !str_detect(Protein.Group, "cRAP")) %>% 
    filter(Q.Value < 0.01, PG.Q.Value < 0.05) %>% 
    run_name(start_from = 7) %>%
    prepare_data() %>% 
    anz_ups2()

} else if (experiment_group == "DDA3" & specification == "mq") {
  # MQ used with only light proteins
  data <- data %>% 
    # select(-(starts_with("Experiment") | starts_with("Intensity")) | contains("UPS3")) %>%   # use only columns with UPS3
    filter(!grepl("\\+", `Potential contaminant`)) %>% 
    pivot_longer(cols = starts_with("Intensity"), names_to = "Run", values_to = "intensity") %>% 
    filter(!intensity == 0) %>% 
    mutate(Run = sub("^Intensity ", "", Run)) %>%
    filter(!Run %in% c("L", "H")) %>% 
    mutate(Channel = case_when(grepl("L", Run) ~ "L", grepl("H", Run) ~ "H", TRUE ~ NA_character_)) %>% 
    mutate(Run = substr(Run, 3, nchar(Run))) %>% 
    mutate(Run = str_c(conditions, "_", Run)) %>% 
    prepare_data() %>% 
    anz_MQ()
  
} else if (experiment_group == "DDA4" & MS == "fozzie" & specification == "diann") { 
  # DiaNN 2.3 used with only light proteins
  data <- data %>% 
    filter(Precursor.Charge > 1, !str_detect(Protein.Group, "cRAP")) %>% 
    filter(Q.Value < 0.01, PG.Q.Value < 0.05) %>% 
    run_name(start_from = 6) %>%
    prepare_data() %>% 
    anz_ups3()
  
} else if (experiment_group == "DDA4" & MS == "sam" & specification == "diann") { 
  # DiaNN 2.3 used with only light proteins
  data <- data %>% 
    filter(Precursor.Charge > 1, !str_detect(Protein.Group, "cRAP")) %>% 
    filter(Q.Value < 0.01, PG.Q.Value < 0.05) %>% 
    run_name(start_from = 7) %>%
    prepare_data() %>% 
    anz_ups3()
  
} else if (experiment_group == "DDA4" & specification == "mq") {
  # MQ used with only light proteins
  data <- data %>% 
    select(-(starts_with("Experiment") | starts_with("Intensity")) | contains("UPS3")) %>%   # use only columns with UPS3
    filter(!grepl("\\+", `Potential contaminant`)) %>% 
    pivot_longer(cols = starts_with("Intensity"), names_to = "Run", values_to = "intensity") %>% 
    filter(!intensity == 0) %>% 
    mutate(Run = sub("^Intensity ", "", Run)) %>%
    prepare_data() %>% 
    anz_ups3()
  
} else {
  stop("no filtering defined for this experiment (", experiment_group, ")")
}

