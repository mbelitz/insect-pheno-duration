library(rgbif)
library(dplyr)
library(lubridate)
library(ggplot2)
library(rnaturalearth)
library(sf)

# function to map gbif occurrence data

map_occ_data <- function(binomial){
  dat <- occ_search(scientificName = binomial,
                    hasCoordinate = TRUE, country = "US",
                    limit = 200000)
  
  sp_df <- dat$data %>% 
    mutate(doy = lubridate::day(eventDate)) %>% 
    filter(!is.na(year)) %>% 
    filter(!is.na(doy))
  
  
  
  us <- rnaturalearth::ne_countries(country = "United States of America", 
                                 returnclass = "sf")
  
  ggplot() + 
    geom_sf(us, mapping = aes()) + 
    geom_point(sp_df, mapping = aes(x = decimalLongitude, y = decimalLatitude, color = basisOfRecord)) + 
    coord_sf( xlim = c(-125, -65), ylim = c(20, 60))+
    ggtitle(label = binomial) +
    facet_wrap(~year)
  
  
  
}

# testing
# map_occ_data("Popillia japonica")

