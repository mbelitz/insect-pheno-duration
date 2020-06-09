library(tidyverse)
library(lme4)
library(lmerTest)
library(sjPlot)

## read in spp list from google drive 
spp_traits <- read.csv('data/traits/insect_traits.csv', stringsAsFactors = FALSE) %>% 
  na_if("")

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
mon <- lmer(onset ~ temp + prec + prec_seas + temp_seas + temp:prec +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + prec_seas | scientificName) + 
               (0 + temp:prec | scientificName),
             data = model_df3, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

mon2 <- lmer(onset ~ temp + prec + prec_seas + temp_seas +
              (1|id_cells) + (1|scientificName) +
              (0 + temp | scientificName) + 
              (0 + prec | scientificName) +
              (0 + temp_seas | scientificName) + 
              (0 + prec_seas | scientificName) + 
              (0 + temp:prec | scientificName),
            data = model_df3, REML = FALSE, 
            lmerControl(optimizer = "bobyqa"))

car::vif(mon)
car::vif(mon2)
MuMIn::Weights(AIC(mon, mon2))

#step(mon)
#Model found:
#  onset ~ temp + prec + temp_seas + (1 | id_cells) + (1 | scientificName) + 
#  (0 + temp | scientificName) + (0 + prec | scientificName) + 
#  (0 + temp_seas | scientificName) + (0 + prec_seas | scientificName) + 
#  (0 + temp:prec | scientificName) + temp:prec

## final model

m_on_final <- lmer(onset ~ temp + prec +  temp_seas + temp:prec +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp:prec | scientificName),
             data = model_df3, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(m_on_final)
car::vif(m_on_final)
suppressWarnings(MuMIn::r.squaredGLMM(m_on_final))
plot_model(m_on_final, type = "pred", terms = c("temp", "prec"))

### Extract Random term estimates
on_coef <- coef(m_on_final)$scientificName %>% 
  tibble::rownames_to_column("scientificName") %>% 
  rename(intercept_ave_onset = "(Intercept)") %>% 
  left_join(spp_traits, by = "scientificName") %>% 
  left_join(spp_seas2, by = "scientificName")

# clean on_coef
on_coef <- on_coef %>% 
  na.omit %>% 
  filter(diapause.stage != "None")

pivot_longer(on_coef, cols = c("higher_taxon", "larval.diet", "flights", 
                               "immature.habitat", "development",
                               "diapause.stage", "seas")) %>% 
  ggplot(aes(x = value, y = intercept_ave_onset)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") 


pivot_longer(on_coef, cols = c("higher_taxon", "larval.diet", "flights", 
                               "immature.habitat", "development",
                               "diapause.stage", "seas")) %>% 
  ggplot(aes(x = value, y = temp)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") +
  labs(y = "Sensentivity to temp")


pivot_longer(on_coef, cols = c("higher_taxon", "larval.diet", "flights", 
                               "immature.habitat", "development",
                               "diapause.stage", "seas")) %>% 
  ggplot(aes(x = value, y = prec)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") +
  labs(y = "Sensentivity to prec")


pivot_longer(on_coef, cols = c("higher_taxon", "larval.diet", "flights", 
                               "immature.habitat", "development",
                               "diapause.stage", "seas")) %>% 
  ggplot(aes(x = value, y = temp_seas)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") +
  labs(y = "Sensentivity to temp_seas")

######## Add traits to model 
model_df4 <- na.omit(model_df3) %>% 
  filter(diapause.stage != "None")

m_on_traits <- lmer(onset ~ temp + prec +  temp_seas + temp:prec +
                      development + temp:development +  prec:development + temp_seas:development +
                      diapause.stage + temp:diapause.stage + prec:diapause.stage + temp_seas:diapause.stage +
                      flights + temp:flights + prec:flights + temp_seas:flights +
                      immature.habitat + temp:immature.habitat + prec:immature.habitat + temp_seas:immature.habitat +
                      larval.diet + temp:larval.diet + prec:larval.diet + temp_seas:larval.diet +
                      (1|id_cells) + (1|scientificName) +
                      (0 + temp | scientificName) + 
                      (0 + prec | scientificName) +
                      (0 + temp_seas | scientificName) + 
                      (0 + temp:prec | scientificName),
                    data = model_df4, REML = FALSE, 
                    lmerControl(optimizer = "bobyqa"))

# m_on_traits_s = step(m_on_traits, reduce.random = T)

# m_on_traits_s
# summary(get_model(m_on_traits_s))

m_on_traits2 <- lmer(onset ~ temp + prec + temp_seas + temp:prec +
                       diapause.stage + immature.habitat +
                       temp:diapause.stage + prec:diapause.stage +
                       (1 | id_cells) + (1 | scientificName) + 
                       (0 + temp | scientificName) + 
                       (0 + prec | scientificName) +
                       (0 + temp_seas | scientificName) + 
                       (0 + temp:prec | scientificName),
                     data = model_df4, REML = FALSE, 
                     lmerControl(optimizer = "bobyqa"))

# make sure we have the top model & no more interactions are needed
# m_on_traits_s2 = step(m_on_traits2, reduce.random = T)
# m_on_traits_s2

car::Anova(m_on_traits2)
MuMIn::r.squaredGLMM(m_on_traits2)
plot_model(m_on_traits2, type = "pred", terms = c("temp", "diapause.stage"), ci.lvl = NA)
plot_model(m_on_traits2, type = "pred", terms = c("prec", "diapause.stage"), ci.lvl = NA)

## Explore traits for offset data

m_off <- lmer(offset ~ temp + prec +  temp_seas + prec_seas + temp:prec +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + prec | scientificName) +
                (0 + prec_seas | scientificName) +
                (0 + temp_seas | scientificName) + 
                (0 + temp:prec | scientificName),
              data = model_df3, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))


summary(m_off)

m_off_s <- step(m_off)
m_off_s

# Final model without traits!
m_off_final <-lmer(offset ~ prec +  temp_seas +
              (1|id_cells) + (1|scientificName) +
              (0 + prec | scientificName) +
              (0 + temp_seas | scientificName),
               data = model_df3, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

car::vif(m_off_final)
suppressWarnings(MuMIn::r.squaredGLMM(m_off_final))

### Extract Random term estimates
off_coef <- coef(m_off_final)$scientificName %>% 
  rownames_to_column("scientificName") %>% 
  rename(intercept_ave_offset = "(Intercept)") %>% 
  left_join(spp_traits, by = "scientificName") %>% 
  left_join(spp_seas2, by = "scientificName")

pivot_longer(off_coef, cols = c("higher_taxon", "larval.diet", "flights", 
                                "immature.habitat", "development",
                                "diapause.stage", "seas")) %>% 
  ggplot(aes(x = value, y = intercept_ave_offset)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") 


pivot_longer(off_coef, cols = c("higher_taxon", "larval.diet", "flights", 
                                "immature.habitat", "development",
                                "diapause.stage", "seas")) %>% 
  ggplot(aes(x = value, y = prec)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") +
  labs(y = "Sensentivity to prec")


pivot_longer(off_coef, cols = c("higher_taxon", "larval.diet", "flights", 
                                "immature.habitat", "development",
                                "diapause.stage", "seas")) %>% 
  ggplot(aes(x = value, y = temp_seas)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") +
  labs(y = "Sensentivity to temp_seas")

######## Add traits to model 
m_off_traits <- lmer(offset ~ prec +  temp_seas +
                      development + prec:development + temp_seas:development +
                      diapause.stage + prec:diapause.stage + temp_seas:diapause.stage +
                      flights + prec:flights + temp_seas:flights +
                      immature.habitat + prec:immature.habitat + temp_seas:immature.habitat +
                      larval.diet +  prec:larval.diet + temp_seas:larval.diet +
                      (1|id_cells) + (1|scientificName) +
                      (0 + prec | scientificName) +
                      (0 + temp_seas | scientificName),
                    data = model_df4, REML = FALSE, 
                    lmerControl(optimizer = "bobyqa"))


m_off_traits_s = step(m_off_traits, reduce.random = T)

m_off_traits_final <- lmer(offset ~ prec +  temp_seas +
                    development + diapause.stage + flights + immature.habitat + larval.diet +
                    (1|id_cells) + (1|scientificName) +
                    (0 + temp_seas | scientificName) + 
                    (0 + prec | scientificName) + 
                    temp_seas:development + prec:diapause.stage + temp_seas:flights +
                      prec:immature.habitat,
                    data = model_df4, REML = FALSE, 
                    lmerControl(optimizer = "bobyqa"))

summary(m_off_traits_final)
car::Anova(m_off_traits_final)
MuMIn::r.squaredGLMM(m_off_traits_final)
plot_model(m_off_traits_final, type = "pred", terms = c("temp_seas", "development"), ci.lvl = NA)
plot_model(m_off_traits_final, type = "pred", terms = c("temp_seas", "flights"), ci.lvl = NA)
plot_model(m_off_traits_final, type = "pred", terms = c("prec", "diapause.stage"),ci.lvl = NA)
plot_model(m_off_traits_final, type = "pred", terms = c("prec", "immature.habitat"), ci.lvl = NA)


############ Duration

m_dur <- lmer(duration ~ temp + prec + temp_seas + prec_seas + temp:prec +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + prec | scientificName) +
                (0 + prec_seas | scientificName) +
                (0 + temp_seas | scientificName) + 
                (0 + temp:prec | scientificName),
              data = model_df4, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

m_dur_s <- step(m_dur)
m_dur_s

## Final Duration Model No Traits ##
m_dur_final <- lmer(duration ~ temp + prec + temp_seas + temp:prec +
                      (1|id_cells) + (1|scientificName) +
                      (0 + temp | scientificName) + 
                      (0 + prec | scientificName) +
                      (0 + temp_seas | scientificName) + 
                      (0 + temp:prec | scientificName),
                    data = model_df4, REML = FALSE, 
                    lmerControl(optimizer = "bobyqa"))

summary(m_dur_final)
car::vif(m_dur_final)
suppressWarnings(MuMIn::r.squaredGLMM(m_dur_final))

### Extract Random term estimates
dur_coef <- coef(m_dur_final)$scientificName %>% 
  rownames_to_column("scientificName") %>% 
  rename(intercept_ave_duration = "(Intercept)") %>% 
  left_join(spp_traits, by = "scientificName") %>% 
  left_join(spp_seas2, by = "scientificName")

pivot_longer(dur_coef, cols = c("higher_taxon", "larval.diet", "flights", 
                                "immature.habitat", "development",
                                "diapause.stage", "seas")) %>% 
  ggplot(aes(x = value, y = intercept_ave_duration)) +
  geom_boxplot() + geom_jitter() +
  facet_wrap(~name, scales = "free") 

## Add traits to duration model
m_dur_traits <- lmer(duration ~ temp + prec +  temp_seas + temp:prec +
                       development + prec:development + temp_seas:development + temp:development +
                       diapause.stage + prec:diapause.stage + temp_seas:diapause.stage + temp:diapause.stage +
                       flights + prec:flights + temp_seas:flights + temp:flights +
                       immature.habitat + prec:immature.habitat + temp_seas:immature.habitat + temp:immature.habitat +
                       larval.diet +  prec:larval.diet + temp_seas:larval.diet + temp:larval.diet +
                       (1|id_cells) + (1|scientificName) +
                       (0 + prec | scientificName) +
                       (0 + temp | scientificName) + 
                       (0 + prec_seas | scientificName) +
                       (0 + temp_seas | scientificName) +
                       (0 + temp:prec | scientificName),
                     data = model_df4, REML = FALSE, 
                     lmerControl(optimizer = "bobyqa"))

m_dur_traits_s = step(m_dur_traits, reduce.random = T)
m_dur_traits_s

m_dur_traits_final <- lmer(duration ~ temp + prec +  temp_seas + 
        development + diapause.stage + flights + immature.habitat + 
          (1 | id_cells) + (1 | scientificName) +
          (0 + prec | scientificName) +
          (0 + temp | scientificName) + 
          (0 + prec_seas | scientificName) +
          (0 + temp_seas | scientificName) +
          (0 + temp:prec | scientificName) + 
          temp:prec +
          temp:development +
          temp_seas:diapause.stage +
          temp:diapause.stage + 
          temp:flights,
          data = model_df4, REML = FALSE, 
          lmerControl(optimizer = "bobyqa"))

summary(m_dur_traits_final)
car::Anova(m_dur_traits_final)
MuMIn::r.squaredGLMM(m_dur_traits_final)
plot_model(m_dur_traits_final, type = "pred", terms = c("temp", "prec"), ci.lvl = NA)
plot_model(m_dur_traits_final, type = "pred", terms = c("temp_seas", "diapause.stage"), ci.lvl = NA)
plot_model(m_dur_traits_final, type = "pred", terms = c("temp", "diapause.stage"),ci.lvl = NA)
plot_model(m_dur_traits_final, type = "pred", terms = c("temp", "development"), ci.lvl = NA)
plot_model(m_dur_traits_final, type = "pred", terms = c("temp", "flights"), ci.lvl = NA)