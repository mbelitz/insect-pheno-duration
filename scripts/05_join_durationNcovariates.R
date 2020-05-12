library(dplyr)
library(sf)
library(ggplot2)

total_duration <- read.csv("phenesse_outputs/total/total_duration_V2_withCoordinates.csv", stringsAsFactors = F)

td <- total_duration %>% 
  mutate(XX = round(X, digits = 0),
         YY = round(Y, digits = 0))

# read in climate data

clim <- read.csv("data/gridded/clim_grid.csv", stringsAsFactors = FALSE)

pop <- read.csv("data/gridded/pop_grid.csv", stringsAsFactors = FALSE)

cd <- clim %>% 
  mutate(XX = round(X, digits = 0),
         YY = round(Y, digits = 0))

pd <- pop %>% 
  mutate(XX = round(X, digits = 0),
         YY = round(Y, digits = 0))

## join to make a super table

dur_clim <- left_join(td, cd, by = c("XX", "YY"))
dur_vars <- left_join(dur_clim, pd, by = c("XX", "YY")) %>% 
  distinct()

total_dur_vars <- dur_vars %>% 
  mutate(prec = case_when(year == 2015 ~ prcp_2015,
                          year == 2016 ~ prcp_2016,
                          year == 2017 ~ prcp_2017,
                          year == 2018 ~ prcp_2018,
                          year == 2019 ~ prcp_2019)) %>% 
  mutate(temp = case_when(year == 2015 ~ tmp_2015,
                          year == 2016 ~ tmp_2016,
                          year == 2017 ~ tmp_2017,
                          year == 2018 ~ tmp_2018,
                          year == 2019 ~ tmp_2019)) %>% 
  mutate(bio4 = case_when(year == 2015 ~ bio4_2015,
                          year == 2016 ~ bio4_2016,
                          year == 2017 ~ bio4_2017,
                          year == 2018 ~ bio4_2018,
                          year == 2019 ~ bio4_2019)) %>% 
  mutate(bio15 = case_when(year == 2015 ~ bio15_2015,
                           year == 2016 ~ bio15_2016,
                           year == 2017 ~ bio15_2017,
                           year == 2018 ~ bio15_2018,
                           year == 2019 ~ bio15_2019))

total_dur_vars <- total_dur_vars %>% 
  dplyr::select(onset, onset_low, onset_high, rowname, year, offset, offset_low, offset_high, 
                duration, scientificName, Order, lon, lat, OS, prec, temp,
                bio4, bio15, pop)

write.csv(total_dur_vars, file = "data/model_dfs/duration_climate_population_data.csv", row.names = FALSE)
