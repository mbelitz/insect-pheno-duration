library(raster)
library(dplyr)
library(sf)
library(ggplot2)
library(raster)


## Model Temp and Precip for estimates

dur <- read.csv("phenesse_outputs/total/total_duration_allspp.csv", stringsAsFactors = FALSE)


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


### raster stack all of our rasters together

p1 <- raster('data/daymet/daymet_v3_prcp_annttl_2015_na.tif')
p2 <- raster('data/daymet/daymet_v3_prcp_annttl_2016_na.tif')
p3 <- raster('data/daymet/daymet_v3_prcp_annttl_2017_na.tif')
p4 <- raster('data/daymet/daymet_v3_prcp_annttl_2018_na.tif')
p5 <- raster('data/daymet/daymet_v3_prcp_annttl_2019_na.tif')

t1 <- raster("data/daymet/daymet_v3_tmax_annavg_2015_na.tif")
t2 <- raster("data/daymet/daymet_v3_tmax_annavg_2016_na.tif")
t3 <- raster("data/daymet/daymet_v3_tmax_annavg_2017_na.tif")
t4 <- raster("data/daymet/daymet_v3_tmax_annavg_2018_na.tif")
t5 <- raster("data/daymet/daymet_v3_tmax_annavg_2019_na.tif")

rstack <- stack(p1, p2, p3, p4, p5, t1, t2, t3, t4, t5)
rstack_proj <- projectRaster(from = rstack, crs = crs(grids_sf), res = 25000)
rstackdf <- as.data.frame(rstack_proj, xy = TRUE)
rstackdf_sf <- st_as_sf(rstackdf, coords = c(x = 'x', y = 'y'),
                    crs = st_crs(rstack_proj)) %>% 
  na.omit()

rstack_join <- st_join(x = rstackdf_sf, y = grids_sf, st_intersects)
rstack_joindf <- st_drop_geometry(rstack_join) %>% 
  rename(prcp_2015 = 1, prcp_2016 = 2, prcp_2017 = 3, prcp_2018 = 4, prcp_2019 = 5,
         tmp_2015 = 6, tmp_2016 = 7, tmp_2017 = 8, tmp_2018 = 9, tmp_2019 = 10, id_cells = 11) %>% 
  as.data.frame()

clim_df <- rstack_joindf %>% 
  group_by(id_cells) %>% 
  summarise(prcp_2015 = mean(prcp_2015),
            prcp_2016 = mean(prcp_2016),
            prcp_2017 = mean(prcp_2017),
            prcp_2018 = mean(prcp_2018),
            prcp_2019 = mean(prcp_2019),
            tmp_2015 = mean(tmp_2015),
            tmp_2016 = mean(tmp_2016),
            tmp_2017 = mean(tmp_2017),
            tmp_2018 = mean(tmp_2018),
            tmp_2019 = mean(tmp_2019))

clim_grid <- left_join(grid_coords, clim_df)

write.csv(clim_grid, file = "data/gridded/clim_grid.csv", row.names = FALSE)
