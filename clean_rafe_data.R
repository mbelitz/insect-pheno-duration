rafe_data <- read.csv("downloaded_data/insects_idigbio_no_geo_2020-03-10a.csv")

rafe_data <- rafe_data %>% 
  mutate(decimalLongitude = lng) %>% 
  mutate(decimalLatitude = lat) %>% 
  mutate(coordinateUncertaintyInMeters = uncert) %>% 
  rename(eventDate = dwc.eventDate, year = dwc.year,
         month = dwc.month, scientificName = dwc.scientificName,
         basisOfRecord = dwc.basisOfRecord, institutionCode = dwc.institutionCode)

  
rafe_data2 <- rafe_data %>% 
  filter(!is.na(decimalLongitude)) %>% 
  filter(!is.na(decimalLatitude)) %>% 
  filter(!is.na(eventDate)) %>% 
  mutate(decade = floor(year/10)*10) %>% 
  mutate(lon_bin = floor(decimalLongitude/1)*1) %>% 
  mutate(lat_bin = floor(decimalLatitude/1)*1)

  
rafe_data2 <- rafe_data2 %>% 
  dplyr::select(decimalLongitude, decimalLatitude, year, month, eventDate, 
                scientificName, basisOfRecord, coordinateUncertaintyInMeters, 
                institutionCode, decade, lat_bin, lon_bin)

write.csv(rafe_data2, "downloaded_data/cleaned_georeferenced_idigbio.csv", row.names = FALSE)
