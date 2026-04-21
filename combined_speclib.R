# adjust speclib path (if changing ups speclib also change output name)
# load data----
path_lib <- "E:/MaxQuant/ACB/speclib/speclib_olive_slice/"
lib_ref <- read_tsv(str_c(path_lib, "as_txt/lib.tsv")) %>%
  mutate(fragment_info = (str_sub(transition_name, str_locate(transition_name, "_") [,1] + 1)))

# test protein number
protein_number <- lib_ref %>%
  filter(!str_detect(ProteinGroup, ";")) %>%
  summarise(n = n_distinct(ProteinGroup)) %>%
  pull(n)

ups <- read_tsv("E:/MaxQuant/ACB/speclib/predicted_speclib/UPS2/as_txt/lib.tsv") %>%
  mutate(fragment_info = (str_sub(transition_name, str_locate(transition_name, "_") [,1] + 1)))

# rows that should be the same but have different masses----
ups_in_lib_ref <- lib_ref %>%
  inner_join(ups, by = c("transition_name")) 

# combine libraries----
{lib_ref_wo_ups <- lib_ref %>%
  anti_join(ups_in_lib_ref, by = c("transition_name"))

if (nrow(lib_ref_wo_ups)==(nrow(lib_ref)-nrow(ups_in_lib_ref))) {
  combined_lib <- bind_rows(lib_ref_wo_ups, ups)
  rm(lib_ref_wo_ups)
  rm(ups_in_lib_ref)
} else {
  print("there is something wrong with the lib")
}}

write_tsv(combined_lib, str_c(path_lib, "combined_lib_ups2.tsv"))
