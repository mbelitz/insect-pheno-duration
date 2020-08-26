library(tidyverse)
library(sf)
library(rnaturalearth)

# study location map

pglmm_data <- read.csv("PGLMM_data/model_df_5cells_treespp.csv")

spp_year_cell <- pglmm_data %>% 
  group_by(lat, lon) %>% 
  summarise(Combinations = n(), years = length(unique(year)), 
            spp = length(unique(scientificName)))


na2 <- maps::map(regions=c("usa", "mexico", "canada"), plot = FALSE, fill = TRUE) %>%
  st_as_sf()
na_equalArea <- st_transform(na2, crs="+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")


study_area <- ggplot() + 
  geom_sf(na_equalArea, mapping = aes(), fill = "grey96") + 
  geom_tile(spp_year_cell, 
             mapping = aes(x = lon, y = lat, fill = Combinations)) +
  coord_sf(xlim = c(-2303077, 2471923), ylim = c(-2475719,1449281)) +
  scale_fill_viridis_c(option = "plasma", trans = "log") + 
  labs(x = "", y = "", fill = "Species x Year Combinations")+
  theme_minimal()

ggsave(filename = "Tables&Figures/studyarea.png", dpi = 400)  
  
  
