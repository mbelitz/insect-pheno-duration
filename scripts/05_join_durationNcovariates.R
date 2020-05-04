library(raster)
library(dplyr)
library(sf)
library(ggplot2)
library(raster)

## Model Temp and Precip for estimates

dur <- read.csv("phenesse_outputs/total/total_duration_id_cellFIX.csv", stringsAsFactors = FALSE)


# make grid cells
na <-  rnaturalearth::ne_countries(country = c("United States of America", "Mexico", "Canada"),
                                   returnclass = "sp")
na_map <-  rnaturalearth::ne_countries(country = c("United States of America", "Mexico", "Canada"),
                                       returnclass = "sf")

# make grid over na
na_map <- st_transform(na_map, crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")
grids <-  st_make_grid(na_map, cellsize = c(100000, 100000))
grids_sf <-  mutate(st_sf(geometry = grids), id_cells = 1:n())

coords <- st_centroid(grids) %>% 
  st_coordinates() %>% 
  data.frame()

grid_coords <- bind_cols(grids_sf, coords) %>% 
  st_set_geometry(NULL)

total_duration <- left_join(dur, grid_coords) %>% 
  mutate(lon = X, lat = Y) %>% 
  mutate(OS = paste(Order, scientificName, sep = ":"))

# read in climate data

clim <- read.csv("data/gridded/clim_grid.csv", stringsAsFactors = FALSE)

pop <- read.csv("data/gridded/pop_grid.csv", stringsAsFactors = FALSE)

## join to make a super table

dur_clim <- left_join(total_duration, clim, by = "id_cells")

total_dur_clim <- dur_clim %>% 
  mutate(prec = case_when(year == 2015 ~ prcp_2015,
                          year == 2016 ~ prcp_2016,
                          year == 2017 ~ prcp_2017,
                          year == 2018 ~ prcp_2018,
                          year == 2019 ~ prcp_2019)) %>% 
  mutate(temp = case_when(year == 2015 ~ tmp_2015,
                          year == 2016 ~ tmp_2016,
                          year == 2017 ~ tmp_2017,
                          year == 2018 ~ tmp_2018,
                          year == 2019 ~ tmp_2019))

total_dur_clim <- total_dur_clim %>% 
  dplyr::select(onset, onset_low, onset_high, rowname, year, offset, offset_low, offset_high, 
                duration, scientificName, Order, id_cells, lon, lat, OS, prec, temp)

total_dur_clim_pop <- left_join(total_dur_clim, pop) %>% 
  dplyr::select(-X, -Y)

write.csv(total_dur_clim_pop, file = "data/model_dfs/duration_climate_population_data.csv", row.names = FALSE)
