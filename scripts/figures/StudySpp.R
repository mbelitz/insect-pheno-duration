library(tidyverse)

# How many species are we losing because we don't have trait data?

model_df <- read.csv(file = "data/model_dfs/duration_climate_population_data_25km.csv",
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

has_5_cells <- filter(datadens, count >= 5) %>% 
  filter (scientificName != "Apis mellifera") # 185

model_df_5cells <- filter(model_df2, scientificName %in% has_5_cells$scientificName)

na <- rnaturalearth::ne_countries(continent = "North America", returnclass = "sf") %>% 
  sf::st_transform(crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" )

has_pheno_data <- filter(model_df_5cells) %>% 
  distinct(scientificName, .keep_all = TRUE) %>% 
  group_by(Order) %>% 
  summarise(pheno_count = n())

## remove migratory species
# make season trait

spp_traits <- read.csv('data/traits/insect_traits_july10.csv', stringsAsFactors = FALSE) %>% 
  na_if("") %>% 
  dplyr::select(1:7)

spp_seas <- model_df2 %>% 
  group_by(scientificName) %>% 
  summarise(ave_on = mean(onset))

spp_seas2 <- spp_seas %>% 
  mutate(seas = case_when(ave_on <= 125 ~ "Spring",
                          ave_on > 125 & ave_on <= 175 ~ "Summer",
                          ave_on > 175 ~ "Fall"))

# combine traits and results
model_df3_5cells <- left_join(model_df_5cells, spp_traits)
traits <- left_join(model_df3_5cells, spp_seas2) %>% 
  select(-ave_on) %>% 
  na.omit()

# unique(model_df3_5cells$scientificName) %>% length() 129 if Migratory could be included

no_traits <- model_df_5cells %>%
  filter(!scientificName %in% traits$scientificName)

unique(no_traits$scientificName) %>% length

no_traits_bp <- no_traits %>% 
  distinct(scientificName, .keep_all = TRUE) %>% 
  group_by(Order) %>% 
  summarise(count = n()) %>% 
  mutate(Included = "No - Missing Traits")

## do not overwinter
noOverwinter_spp_bp <- model_df3_5cells %>% 
  filter(diapause.stage == "None")

noOverwinter_spp_bp <- noOverwinter_spp_bp %>% 
  distinct(scientificName, .keep_all = TRUE) %>% 
  group_by(Order) %>% 
  summarise(count = n()) %>% 
  mutate(Included = "No - NoOverwinter")

  
## Mig spp

# migratory spp

mig_spp <- c("Anax junius", "Pantala flavenscens", "Pantala hymenaea", "Tramea lacerata", 
             "Sympetrum corruptum", "Sympetrum vicinum", "Libellula pulchella", "Libellula vibrans",
             "Tramea lacerata", "Tramea onusta", "Pantala flavescens", "Libellula quadrimaculata",
             "Ertythrodiplax umbrata", "Epiaeschna heros", "Tramea carolina",
             "Libellula semifasciata", "Pantala hymenaea", 
             "Spoladea recurvalis", "Ponoquina ocola", "Plutella xylostella",
             "Chrysodeixis includens", "Phoebis sennae", "Abaeis nicippe",
             "Libytheana carinenta", "Agraulis vanillae", "Junonia coenia",
             "Danaus plexippus", "Vanessa virginiensis", "Vanessa cardui",
             "Vanessa atalanta", "Danaus gilippus", "Nymphalis antiopa", 
             "Polygonia interrogationis", "Lerema accius")

mig_spp_bp <- model_df_5cells %>% 
  filter(scientificName %in% mig_spp)

mig_spp_bp <- mig_spp_bp %>% 
  distinct(scientificName, .keep_all = TRUE) %>% 
  group_by(Order) %>% 
  summarise(count = n()) %>% 
  mutate(Included = "No - Migratory")

mig_spp_bp2 <- full_join(mig_spp_bp, noOverwinter_spp_bp)

mig_spp_bp3 <- data.frame(
  Order = c("Lepidoptera", "Odonata" , "Diptera"),
  count = c(13, 11, 1),
  Included = c("No - Migratory/No Diapuase", "No - Migratory/No Diapuase", "No - Migratory/No Diapuase")
)

# Species included

final_spp <- model_df_5cells %>% 
  distinct(scientificName, .keep_all = TRUE) %>% 
  group_by(Order) %>% 
  summarise(count = n()) %>% 
  mutate(Included = "Yes")

## Total SPP df

t_spp_df <- rbind(final_spp, mig_spp_bp3, no_traits_bp)


spp_fig <- ggplot(t_spp_df) + 
  geom_bar(stat = "identity", aes (x = Order, y = count, fill = Included)) +
  labs(y = "Number of Species") +
  scale_fill_brewer(palette = 3) +
  scale_y_continuous(expand = c(0,0)) +
  theme_classic()

ggsave(plot = spp_fig, filename = "Tables&Figures/StudySpecies.png", dpi = 400,
       width = 8, height = 5)




