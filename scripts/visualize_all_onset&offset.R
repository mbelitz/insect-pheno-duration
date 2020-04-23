library(dplyr)
library(sf)
library(ggplot2)
library(raster)
library(ggforce)

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

ggplot(filter(total_duration, year == 2019)) + 
  geom_point(aes(x = onset, y = lat, shape = as.character(year), color = Order), alpha = 0.25) +
  geom_point(aes(x = offset, y = lat, shape = as.character(year), color = Order), alpha = 0.25) +
  geom_smooth(aes(x = onset, y = lat, color = Order), method = "lm", se = FALSE) +
  geom_smooth(aes(x = offset, y = lat, color = Order),  method = "lm", se = FALSE) +
  scale_shape_discrete(name = "Year") +
  labs(x = "Onset & Offset", y = "Latitude (m)") +
  scale_color_viridis_d(option = "plasma") +
  theme_bw() 

ggsave(filename = "outputs_to_share/models/all_spp_together.png", width = 8, height = 7)


###### Now to plot species individually ############

p <- ggplot(total_duration) + 
  geom_point(aes(x = onset, y = lat, shape = as.character(year)), color = "purple", alpha = 0.25) +
  geom_point(aes(x = offset, y = lat, shape = as.character(year)),color = "black", alpha = 0.25) +
  geom_smooth(aes(x = onset, y = lat), method = "lm",color = "purple", se = FALSE) +
  geom_smooth(aes(x = offset, y = lat),  method = "lm", color = "black", se = FALSE) +
  scale_shape_discrete(name = "Year") +
  labs(x = "Onset & Offset", y = "Latitude (m)") +
  scale_color_viridis_d(option = "plasma") +
  theme_bw()  +
  facet_wrap_paginate(~OS, scales = "free",
                      nrow = 1, ncol =1)


required_n_pages <- n_pages(p)

  
 g <-  ggplot(total_duration) + 
    geom_point(aes(x = onset, y = lat, shape = as.character(year)), color = "purple", alpha = 0.25) +
    geom_point(aes(x = offset, y = lat, shape = as.character(year)),color = "black", alpha = 0.25) +
    geom_smooth(aes(x = onset, y = lat), method = "lm",color = "purple", se = FALSE) +
    geom_smooth(aes(x = offset, y = lat),  method = "lm", color = "black", se = FALSE) +
    scale_shape_discrete(name = "Year") +
    labs(x = "Onset & Offset", y = "Latitude (m)") +
    scale_color_viridis_d(option = "plasma") +
    theme_bw()  +
    facet_wrap_paginate(~OS, scales = "free",
                        nrow = 5, ncol =5, page = 6)
 

# Visualize for temperate 
  
ggplot(total_dur_clim) + 
   geom_point(aes(y = onset, x = temp, shape = as.character(year)), color = "purple", alpha = 0.25) +
   geom_point(aes(y = offset, x = temp, shape = as.character(year)),color = "black", alpha = 0.25) +
   geom_smooth(aes(y = onset, x = temp), method = "lm",color = "purple", se = FALSE) +
   geom_smooth(aes(y = offset, x = temp),  method = "lm", color = "black", se = FALSE) +
   scale_shape_discrete(name = "Year") +
   labs(x = "Temp", y = "Onset & Offset") +
   scale_color_viridis_d(option = "plasma") +
   theme_bw()  +
   facet_wrap_paginate(~OS, scales = "free",
                       nrow = 5, ncol =5, page = 1)
 
onof <- ggplot(total_dur_clim) + 
  geom_point(aes(y = onset, x = temp, shape = as.character(year)), color = "purple", alpha = 0.25) +
  geom_point(aes(y = offset, x = temp, shape = as.character(year)),color = "black", alpha = 0.25) +
  geom_smooth(aes(y = onset, x = temp), method = "lm",color = "purple", se = FALSE) +
  geom_smooth(aes(y = offset, x = temp),  method = "lm", color = "black", se = FALSE) +
  scale_shape_discrete(name = "Year") +
  labs(x = "Temp", y = "Onset & Offset") +
  scale_color_viridis_d(option = "plasma") +
  theme_bw()  +
  facet_wrap(~OS, scales = "free", ncol = 6)

ggsave('outputs_to_share/models/onset_temp.pdf', plot = onof, width = 16, height = 36)


onof_precip <- ggplot(total_dur_clim) + 
  geom_point(aes(y = onset, x = prec, shape = as.character(year)), color = "purple", alpha = 0.25) +
  geom_point(aes(y = offset, x = prec, shape = as.character(year)),color = "black", alpha = 0.25) +
  geom_smooth(aes(y = onset, x = prec), method = "lm",color = "purple", se = FALSE) +
  geom_smooth(aes(y = offset, x = prec),  method = "lm", color = "black", se = FALSE) +
  scale_shape_discrete(name = "Year") +
  labs(x = "Precipitation", y = "Onset & Offset") +
  scale_color_viridis_d(option = "plasma") +
  theme_bw()  +
  facet_wrap(~OS, scales = "free", ncol = 6)

ggsave('outputs_to_share/models/onset_precip.pdf', plot = onof_precip, width = 16, height = 36)
 
duration_temp <- ggplot(total_dur_clim) + 
  geom_point(aes(y = duration, x = temp, shape = as.character(year)), alpha = 0.5) +
  geom_smooth(aes(y = duration, x = temp), method = "lm") +
  scale_shape_discrete(name = "Year") +
  labs(x = "Temp", y = "Duration") +
  theme_bw()  +
  facet_wrap(~OS, scales = "free", ncol = 6)

ggsave('outputs_to_share/models/duration_temp.pdf', plot = duration_temp, width = 16, height = 36)
