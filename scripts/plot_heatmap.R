library(ridigbio)
library(rgbif)
library(dplyr)
library(ggplot2)
library(hexbin)


# ridigbio
io_idig <- idig_search_records(limit = 10000, rq = list(scientificname = "Automeris io"), 
            fields = c("data.dwc:decimalLongitude", "data.dwc:decimalLatitude",
                       "data.dwc:year", "data.dwc:month", "data.dwc:eventDate",
                       "data.dwc:scientificName", "data.dwc:basisOfRecord")) 

names(io_idig) <- c("decimalLongitude", "decimalLatitude", "year", "month", "eventDate",
                    "scientificName", "basisOfRecord")

#rgbif
io_gbif <- occ_search(scientificName = "Automeris io",
                      hasCoordinate = TRUE, country = "US",
                      limit = 200000)

sp_df <- io_gbif$data %>% 
  filter(basisOfRecord == "HUMAN_OBSERVATION") %>% 
  select(decimalLongitude, decimalLatitude, year, month, eventDate, 
         scientificName, basisOfRecord)

total <- rbind(sp_df, io_idig) %>% 
  filter(!is.na(year)) %>% 
  filter(!is.na(decimalLatitude)) %>% 
  filter(!is.na(decimalLongitude))

total$decimalLongitude <- as.numeric(total$decimalLongitude)
total$decimalLatitude <- as.numeric(total$decimalLatitude)
total$year <- as.numeric(total$year)

total <- total %>% 
  mutate(decade = floor(year/10)*10)
  

us <- rnaturalearth::ne_countries(country = "United States of America", 
                                  returnclass = "sf")

a <- ggplot() + 
  geom_sf(us, mapping = aes()) + 
  geom_point(total, mapping = aes(x = decimalLongitude, y = decimalLatitude, color = basisOfRecord)) +
  coord_sf( xlim = c(-125, -65), ylim = c(20, 60)) +
  facet_wrap(~decade)
a

ggplot() + 
  geom_sf(us, mapping = aes()) +
  geom_bin2d(total,mapping = aes(x = decimalLongitude, y = decimalLatitude), alpha = 0.9, binwidth =c(2,2)) + 
  coord_sf( xlim = c(-125, -65), ylim = c(20, 60)) +
  scale_fill_continuous(type = "viridis",trans = "log") + 
  facet_wrap(~decade)
