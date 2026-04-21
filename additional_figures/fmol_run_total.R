a <- attr(data, "u2os_absolut_fmol") %>%
  # filter(condition == "200ng_U2OS_UPS2") %>%
  left_join(u2os_mw_longest %>% 
              select(Protein.Group_short, protein_mw_Da), by = c("Protein.Group" = "Protein.Group_short")) %>% 
  mutate(prot_mass_pg = absolute * protein_mw_Da * 10^(-3))

protein_mass_pg <- sum(a$prot_mass_pg)
# b is total prot mass in pg

protein_mass_ng <- protein_mass_pg/1000
# c is total prot amount in ng
