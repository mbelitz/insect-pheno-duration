library(dplyr)

## read in spp list from google drive 
spp_traits <- read.csv('data/traits/spp_traits.csv', stringsAsFactors = FALSE)

## read in phenesse outputs
model_df <- read.csv(file = "data/model_dfs/duration_climate_population_data_25km_V2.csv",
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

# make season trait

spp_seas <- model_df2 %>% 
  group_by(scientificName) %>% 
  summarise(ave_on = mean(onset))

hist(spp_seas$ave_on)

spp_seas2 <- spp_seas %>% 
  mutate(seas = case_when(ave_on <= 125 ~ "Spring",
                          ave_on > 125 & ave_on <= 175 ~ "Summer",
                          ave_on > 175 ~ "Fall"))

# combine traits and results

model_df3 <- left_join(model_df2, spp_traits)
model_df3 <- left_join(model_df3, spp_seas2) %>% 
  select(-ave_on)

# Onset model

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

#step(m3on)
# Model found:
#onset ~ temp + pop + prec + temp_seas + temp:prec +
#  (1 | id_cells) + (1 | scientificName) + 
#  (0 + temp | scientificName) + 
#  (0 + pop |scientificName) + 
#  (0 + prec | scientificName) + 
#  (0 + temp_seas | scientificName) +
#  (0 + temp:prec | sceintificName)

## final model

m_on <- lmer(onset ~ temp + pop + prec +  temp_seas + temp:prec +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + pop | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + temp:prec | scientificName),
             data = model_df3, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(m_on)
# step(m_onset2) # confirmed that no interactions needed
car::vif(m_on)
suppressWarnings(MuMIn::r.squaredGLMM(m_on))

### Extract Random term estimates

on_coef <- coef(m_on)$scientificName %>% 
  rownames_to_column("scientificName") %>% 
  rename(intercept_ave_onset = "(Intercept)") %>% 
  left_join(spp_traits, by = "scientificName") %>% 
  left_join(spp_seas2, by = "scientificName")

pivot_longer(on_coef, cols = c("higher_taxon", "larval.diet", "gen.time", 
                               "immature.habitat", "development",
                               "diapause.stage", "Migrate.", "seas")) %>% 
  ggplot(aes(x = value, y = intercept_ave_onset)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") 


pivot_longer(on_coef, cols = c("higher_taxon", "larval.diet", "gen.time", 
                               "immature.habitat", "development",
                               "diapause.stage", "Migrate.", "seas")) %>% 
  ggplot(aes(x = value, y = temp)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") +
  labs(y = "Sensentivity to temp")


pivot_longer(on_coef, cols = c("higher_taxon", "larval.diet", "gen.time", 
                               "immature.habitat", "development",
                               "diapause.stage", "Migrate.", "seas")) %>% 
  ggplot(aes(x = value, y = pop)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") +
  labs(y = "Sensentivity to pop")

pivot_longer(on_coef, cols = c("higher_taxon", "larval.diet", "gen.time", 
                               "immature.habitat", "development",
                               "diapause.stage", "Migrate.", "seas")) %>% 
  ggplot(aes(x = value, y = prec)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") +
  labs(y = "Sensentivity to prec")


pivot_longer(on_coef, cols = c("higher_taxon", "larval.diet", "gen.time", 
                               "immature.habitat", "development",
                               "diapause.stage", "Migrate.", "seas")) %>% 
  ggplot(aes(x = value, y = temp_seas)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") +
  labs(y = "Sensentivity to temp_seas")

######## Add traits to model 

m_on_traits <- lmer(onset ~ temp + pop + prec +  temp_seas + temp:prec +
               development + temp:development + pop:development + prec:development + 
               diapause.stage + temp:diapause.stage + pop:diapause.stage + prec:diapause.stage +
               gen.time + temp:gen.time + pop:gen.time + prec:gen.time + 
               immature.habitat + temp:immature.habitat + pop:immature.habitat + prec:immature.habitat +
               larval.diet + temp:larval.diet + pop:larval.diet + prec:larval.diet +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + pop | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + temp:prec | scientificName),
             data = model_df3, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

m_on_traits_s = step(m_on_traits, reduce.random = T)

summary(get_model(m_on_traits_s))
car::Anova(get_model(m_on_traits_s))
MuMIn::r.squaredGLMM(get_model(m_on_traits_s))
plot(effects::effect("temp:prec", get_model(m_on_traits_s)), lines=list(multiline=TRUE))

plot(effects::predictorEffects(get_model(m_on_traits_s), ~ temp), multiline = T)



## Explore traits for offset data

m_off <- lmer(offset ~ temp + pop + prec +  temp_seas + temp:prec +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + pop | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + temp:prec | scientificName),
             data = model_df3, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(m_off)

m_off_s <- step(m_off)

m_off <- lmer(offset ~ prec +  temp_seas +
                (1|id_cells) + (1|scientificName) +
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName),
              data = model_df3, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

# step(m_onset2) # confirmed that no interactions needed
car::vif(m_off)
suppressWarnings(MuMIn::r.squaredGLMM(m_off))

### Extract Random term estimates

off_coef <- coef(m_off)$scientificName %>% 
  rownames_to_column("scientificName") %>% 
  rename(intercept_ave_offset = "(Intercept)") %>% 
  left_join(spp_traits, by = "scientificName") %>% 
  left_join(spp_seas2, by = "scientificName")

pivot_longer(off_coef, cols = c("higher_taxon", "larval.diet", "gen.time", 
                               "immature.habitat", "development",
                               "diapause.stage", "Migrate.", "seas")) %>% 
  ggplot(aes(x = value, y = intercept_ave_offset)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") 


pivot_longer(off_coef, cols = c("higher_taxon", "larval.diet", "gen.time", 
                               "immature.habitat", "development",
                               "diapause.stage", "Migrate.", "seas")) %>% 
  ggplot(aes(x = value, y = prec)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") +
  labs(y = "Sensentivity to prec")


pivot_longer(off_coef, cols = c("higher_taxon", "larval.diet", "gen.time", 
                               "immature.habitat", "development",
                               "diapause.stage", "Migrate.", "seas")) %>% 
  ggplot(aes(x = value, y = temp_seas)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") +
  labs(y = "Sensentivity to temp_seas")

######## Add traits to model 

m_off_traits <- lmer(offset ~ temp + pop + prec +  temp_seas + temp:prec +
                      development + temp:development + pop:development + prec:development + 
                      diapause.stage + temp:diapause.stage + pop:diapause.stage + prec:diapause.stage +
                      gen.time + temp:gen.time + pop:gen.time + prec:gen.time + 
                      immature.habitat + temp:immature.habitat + pop:immature.habitat + prec:immature.habitat +
                      larval.diet + temp:larval.diet + pop:larval.diet + prec:larval.diet +
                      (1|id_cells) + (1|scientificName) +
                      (0 + temp | scientificName) + 
                      (0 + pop | scientificName) + 
                      (0 + prec | scientificName) +
                      (0 + temp_seas | scientificName) + 
                      (0 + temp:prec | scientificName),
                    data = model_df3, REML = FALSE, 
                    lmerControl(optimizer = "bobyqa"))

m_on_traits_s = step(m_on_traits, reduce.random = T)

summary(get_model(m_on_traits_s))
car::Anova(get_model(m_on_traits_s))
MuMIn::r.squaredGLMM(get_model(m_on_traits_s))
plot(effects::effect("temp:prec", get_model(m_on_traits_s)), lines=list(multiline=TRUE))

plot(effects::predictorEffects(get_model(m_on_traits_s), ~ temp), multiline = T)
