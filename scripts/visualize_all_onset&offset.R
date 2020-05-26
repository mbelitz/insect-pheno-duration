library(dplyr)
library(sf)
library(ggplot2)
library(raster)
library(ggforce)

## Model Temp and Precip for estimates

dur <- read.csv("data/model_dfs/duration_climate_population_data.csv", stringsAsFactors = FALSE)

ggplot(filter(dur, year == 2019)) + 
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

p <- ggplot(dur) + 
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
 

# Visualize for temperate and stuff
 
dur <- read.csv("data/model_dfs/duration_climate_population_data.csv", stringsAsFactors = FALSE)
  
ggplot(dur) + 
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
 
onof <- ggplot(dur) + 
  geom_point(aes(y = onset, x = temp, shape = as.character(year)), color = "purple", alpha = 0.25) +
  geom_point(aes(y = offset, x = temp, shape = as.character(year)),color = "black", alpha = 0.25) +
  geom_smooth(aes(y = onset, x = temp), method = "lm",color = "purple", se = FALSE) +
  geom_smooth(aes(y = offset, x = temp),  method = "lm", color = "black", se = FALSE) +
  scale_shape_discrete(name = "Year") +
  labs(x = "Temp", y = "Onset & Offset") +
  theme_bw()  +
  facet_wrap(~OS, scales = "free", ncol = 6)

ggsave('outputs_to_share/models/onset_temp.pdf', plot = onof, width = 16, 
       height = 65, limitsize = FALSE)

onof_precip <- ggplot(dur) + 
  geom_point(aes(y = onset, x = prec, shape = as.character(year)), color = "purple", alpha = 0.25) +
  geom_point(aes(y = offset, x = prec, shape = as.character(year)),color = "black", alpha = 0.25) +
  geom_smooth(aes(y = onset, x = prec), method = "lm",color = "purple", se = FALSE) +
  geom_smooth(aes(y = offset, x = prec),  method = "lm", color = "black", se = FALSE) +
  scale_shape_discrete(name = "Year") +
  labs(x = "Precipitation", y = "Onset & Offset") +
  theme_bw()  +
  facet_wrap(~OS, scales = "free", ncol = 6)

ggsave('outputs_to_share/models/onset_precip.pdf', plot = onof_precip, width = 16, 
       height = 65, limitsize = FALSE)

dur <- dur %>% 
  filter(year >= 2015)

duration_temp <- ggplot(dur) + 
  geom_point(aes(y = duration, x = temp, shape = as.character(year)), alpha = 0.5) +
  geom_smooth(aes(y = duration, x = temp), method = "lm") +
  scale_shape_discrete(name = "Year") +
  labs(x = "Temp", y = "Duration") +
  theme_bw()  +
  facet_wrap(~OS, scales = "free", ncol = 6)

ggsave('outputs_to_share/models/duration_temp.pdf', plot = duration_temp, width = 16, height = 49)

onof_pop <- ggplot(dur) + 
  geom_point(aes(y = onset, x = pop, shape = as.character(year)), color = "purple", alpha = 0.25) +
  geom_point(aes(y = offset, x = pop, shape = as.character(year)),color = "black", alpha = 0.25) +
  geom_smooth(aes(y = onset, x = pop), method = "lm",color = "purple", se = FALSE) +
  geom_smooth(aes(y = offset, x = pop),  method = "lm", color = "black", se = FALSE) +
  scale_shape_discrete(name = "Year") +
  labs(x = "Log10 Population Density", y = "Onset & Offset") +
  scale_x_log10() +
  theme_bw()  +
  facet_wrap(~OS, scales = "free", ncol = 6)

ggsave('outputs_to_share/models/onset_pop.pdf', plot = onof_pop, width = 16, height = 36)
