source("scripts/00_pkg_functions.R")

# Run phenesse

# Ischnura verticalis

iv <- read.csv("downloaded_data/Ischnura verticalis.csv", stringsAsFactors = FALSE) %>% 
  dplyr::select(-X)

iv_100km <- plt_summary(cell_size = 100000, dat = iv, n_per_cell = 10)

duration_iv <- run_phenesse(df = iv_100km$dat_to_use, minimum_obs = 10, earliest_year = 2015, 
                            last_year = 2019,n_item = 5000, onset_perct = 0.05, 
                            num_cores = 30, offset_perct = 0.95)

