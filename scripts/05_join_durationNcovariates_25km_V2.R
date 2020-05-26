library(raster)
library(dplyr)
library(sf)
library(ggplot2)
library(raster)

## Model Temp and Precip for estimates
dur <- read.csv("phenesse_outputs/total_25km/total_duration_25km_V2.csv", stringsAsFactors = FALSE)
dur_sf <- st_as_sf(dur, coords = c("X", "Y"))

# read in climate data
clim <- read.csv("data/gridded/clim_grid.csv", stringsAsFactors = FALSE)
clim_sf <- st_as_sf(clim, coords = c("X", "Y"))

pop <- read.csv("data/gridded/pop_grid.csv", stringsAsFactors = FALSE)
pop_sf <- st_as_sf(pop, coords = c("X", "Y")) 

near_pop <- st_nearest_feature(dur_sf, pop_sf)
near_clim <- st_nearest_feature(dur_sf, clim_sf)

pop <- pop %>% 
  tibble::rownames_to_column() %>% 
  mutate(nearest_pop = as.integer(rowname))
clim <- clim %>% 
  tibble::rownames_to_column() %>% 
  mutate(nearest_clim = as.integer(rowname))

dur <- dur %>% 
  mutate(nearest_pop = near_pop) %>% 
  mutate(nearest_clim = near_clim)

## join to make a super table
dur_clim <- left_join(dur, clim, by = "nearest_clim")
dur_clim_pop <- left_join(dur_clim, pop, by = "nearest_pop")

total_dur_clim_pop <- dur_clim_pop %>% 
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

total_dur_clim_pop <- total_dur_clim_pop %>% 
  dplyr::select(onset, year, offset, 
                duration, scientificName, Order, lon, lat, OS, 
                prec, temp, bio4, bio15, pop)

ggplot(total_dur_clim_pop) + 
  geom_tile(mapping = aes(x = lon, y = lat, fill = log10(pop))) +
  geom_point(mapping = aes(x = lon, y = lat, color = duration), size = 2) + 
  scale_fill_viridis_c() +
  coord_equal()

write.csv(total_dur_clim_pop, file = "data/model_dfs/duration_climate_population_data_25km_V2.csv", row.names = FALSE)
