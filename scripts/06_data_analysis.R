library(lme4)
library(lmerTest)
library(dplyr)
library(effects)
library(ggplot2)
library(car)

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
         temp_seas = scale(bio4),
         onsetci = onset_high - onset_low,
         offsetci = offset_high - offset_low)

datadens <- model_df2 %>% 
  group_by(scientificName, Order) %>% 
  summarise(count = n())

has_5_cells <- filter(datadens, count >= 5) %>% 
  filter (scientificName != "Apis mellifera") #186 species 

ggplot(has_5_cells, aes(x = Order)) + 
  geom_bar(stat = "count", aes (fill = Order)) +
  labs(y = "Number of Species") +
  ggtitle("At Least 5 Comparable Cells - 186 total spp") +
  theme_bw()


has_10_cells <- filter(datadens, count >= 10) %>% 
  filter (scientificName != "Apis mellifera") # 105

ggplot(has_10_cells, aes(x = Order)) + 
  geom_bar(stat = "count", aes (fill = Order)) +
  labs(y = "Number of Species") +
  ggtitle("At Least 10 Comparable Cells - 105 total spp") +
  theme_bw()

model_df2 <- filter(model_df2, scientificName %in% has_10_cells$scientificName)

########## Onset Models ###########

m1on <- lmer(onset ~ temp + pop + prec + prec_seas + temp_seas +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName) + 
                (0 + prec_seas | scientificName),
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m1on) 

m2on <- lmer(onset ~ temp + pop + prec + prec_seas + temp_seas + temp:pop +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + pop | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + prec_seas | scientificName) + 
               (0 + temp:pop | scientificName),
             data = model_df2, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(m2on) 

m3on <- lmer(onset ~ temp + pop + prec + prec_seas + temp_seas + temp:prec +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + pop | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + prec_seas | scientificName) + 
               (0 + temp:prec | scientificName),
             data = model_df2, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(m3on) 


m4on <- lmer(onset ~ temp + pop + prec + temp_seas + temp:pop +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + pop | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + temp:pop | scientificName),
             data = model_df2, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(m4on) 

m5on <- lmer(onset ~ temp + pop + prec + temp_seas + temp:prec +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + pop | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + temp:prec | scientificName),
             data = model_df2, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(m5on) 


m6on <- lmer(onset ~ temp + pop + prec + temp_seas + 
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + pop | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName),
             data = model_df2, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(m6on) 

# model rank top models 
AIC(m1on, m2on, m3on, m4on, m5on, m6on)
MuMIn::Weights(AIC(m1on, m2on, m3on, m4on, m5on, m6on))

summary(m3on)
plot(effects::effect("temp:prec", m3on), multiline = T)
MuMIn::r.squaredGLMM(m3on)
car::vif(m3on)
class(m3on) = "lmerMod"
performance::check_collinearity(m3on)
performance::model_performance(m3on)




######### OFFSET MODELS

m1off <- lmer(offset ~ temp + pop + prec + prec_seas + temp_seas +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + pop | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + prec_seas | scientificName),
             data = model_df2, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(m1off) 

m2off <- lmer(offset ~ temp + pop + prec + prec_seas + temp_seas + temp:pop +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName) + 
                (0 + prec_seas | scientificName) +
                (0 + temp:pop | scientificName),
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m2off) 

m1off <- lmer(offset ~ temp + pop + prec + prec_seas + temp_seas +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName) + 
                (0 + prec_seas | scientificName),
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m1off) 