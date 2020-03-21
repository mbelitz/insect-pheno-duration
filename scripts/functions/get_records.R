library(rgbif)
library(ridigbio)
library(dplyr)

get_records <- function(binomial){

  idig <- idig_search_records(limit = 20000, rq = list(scientificname = binomial), 
                            fields = c("data.dwc:decimalLongitude", "data.dwc:decimalLatitude",
                                       "data.dwc:year", "data.dwc:month", "data.dwc:eventDate",
                                       "data.dwc:scientificName", "data.dwc:basisOfRecord", 
                                       "data.dwc:coordinateUncertaintyInMeters", "data.dwc:institutionCode" )) 

  names(idig) <- c("decimalLongitude", "decimalLatitude", "year", "month", "eventDate",
                   "scientificName", "basisOfRecord", "coordinateUncertaintyInMeters", "institutionCode")
  
  #rgbif
  gbif <- occ_search(scientificName = binomial,
                     hasCoordinate = TRUE, country = "US",
                     limit = 200000)
  
  gbif_df <- gbif$data %>% 
    dplyr::filter(institutionCode == "iNaturalist") %>% 
    dplyr::select(decimalLongitude, decimalLatitude, year, month, eventDate, 
           scientificName, basisOfRecord, coordinateUncertaintyInMeters, institutionCode)
  
  total <- rbind(gbif_df, idig) %>% 
    filter(!is.na(year)) %>% 
    filter(!is.na(decimalLatitude)) %>% 
    filter(!is.na(decimalLongitude)) %>% 
    filter(!is.na(eventDate)) %>% 
    filter(decimalLatitude != 0 | decimalLongitude != 0)
  
  total$decimalLongitude <- as.numeric(total$decimalLongitude)
  total$decimalLatitude <- as.numeric(total$decimalLatitude)
  total$year <- as.numeric(total$year)
  
  total <- total %>% 
    mutate(decade = floor(year/10)*10) %>% 
    mutate(lon_bin = floor(decimalLongitude/1)*1) %>% 
    mutate(lat_bin = floor(decimalLatitude/1*1))
  
  return(total)
}

# example
#
#en <- get_records("Epicauta normalis")
#head(en)
#write.csv(en, "downloaded_data/Epicauta_normalis.csv", row.names = FALSE)
