source("functions.R")

report_u2os <- load_data("E:/MaxQuant/ACB/DDA/test3_fozzie/diann/20251127_fozzie_u2os_v23/report.parquet") %>% 
  mutate(Channel = "H")
report_ups2 <- load_data("E:/MaxQuant/ACB/DDA/test3_fozzie/diann/20251126_fozzie_ups2_v23/report.parquet") %>% 
  mutate(Channel = "L")


combined_reprort <- bind_rows(report_u2os, report_ups2)

write_tsv(combined_reprort, str_c("E:/MaxQuant/ACB/DDA/test3_fozzie/diann/", "combined_report.tsv"))



source("functions.R")

report_u2os <- load_data("E:/MaxQuant/ACB/DDA/test3_sam/diann/20250926_sam_u2os_v23/report.parquet") %>% 
  mutate(Channel = "H")
report_ups2 <- load_data("E:/MaxQuant/ACB/DDA/test3_sam/diann/20250926_sam_ups2b_v23/report.parquet") %>% 
  mutate(Channel = "L")


combined_reprort <- bind_rows(report_u2os, report_ups2)

write_tsv(combined_reprort, str_c("E:/MaxQuant/ACB/DDA/test3_sam/diann/", "combined_report.tsv"))
