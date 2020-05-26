library(dplyr)

## read in spp list from google drive 
spp_list <- read.csv('data/traits/spp_list.csv', stringsAsFactors = FALSE) %>% 
  rename(scientificName = scientific_name)

## read in phenesse outputs
model_df <- read.csv(file = "data/model_dfs/duration_climate_population_data.csv",
                     stringsAsFactors = FALSE) 

id_cells <- model_df %>% 
  group_by(lon, lat) %>% 
  summarise(count = n()) %>% 
  tibble::rownames_to_column() %>% 
  rename(id_cells = rowname)

model_df <- left_join(model_df, id_cells)

model_df2 <- model_df %>% 
  na.omit() %>% 
  mutate(temp = scale(temp),
         prec = scale(prec),
         pop = scale(log10(pop)),
         prec_seas = scale(bio15),
         temp_seas = scale(bio4))

datadens <- model_df2 %>% 
  group_by(scientificName, Order) %>% 
  summarise(count = n())

has_10_cells <- filter(datadens, count >= 10) %>% 
  filter (scientificName != "Apis mellifera") # 145

model_df2 <- filter(model_df2, scientificName %in% has_10_cells$scientificName)

### combine model_df2 w/ traits
model_df3 <- filter(spp_list, scientificName %in% has_10_cells$scientificName)

# save csv

write.csv(model_df3, "data/traits/spp_traits.csv", row.names = FALSE)
