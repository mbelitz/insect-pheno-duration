library(raster)
library(dplyr)
library(sf)
library(ggplot2)
library(raster)

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

##

pop <- raster("data/population/gpw_v4_population_density_rev11_2020_15_min.tif")

na_pop <- crop(pop, na)

na_pop_proj <- projectRaster(from = na_pop, crs = crs(grids_sf))
napopdf <- as.data.frame(na_pop_proj, xy = TRUE)
pop_sf <- st_as_sf(napopdf, coords = c(x = 'x', y = 'y'),
                   crs = st_crs(na_map)) %>% 
  na.omit() %>% 
  rename(pop = 1)

pop_join <- st_join(x = pop_sf, y = grids_sf, st_intersects)
pop_joindf <- st_drop_geometry(pop_join) %>% 
  as.data.frame()

pop_df <- pop_joindf %>% 
  group_by(id_cells) %>% 
  summarise(pop = mean(pop))

na_grid <- st_join(grids_sf, na_map, st_intersects) %>% 
  filter(!is.na(sovereignt))

nacoords <- st_centroid(na_grid) %>% 
  st_coordinates() %>% 
  data.frame()

na_grid_coords <- bind_cols(na_grid, nacoords) %>% 
  st_set_geometry(NULL) %>% 
  dplyr::select(id_cells, X, Y)

pop_grid <- left_join(na_grid_coords, pop_df)

ggplot() +
  geom_tile(pop_grid, mapping = aes(x = X, y = Y, fill = pop)) +
  scale_fill_viridis_c()

write.csv(pop_grid, file = "data/gridded/pop_grid.csv", row.names = FALSE)