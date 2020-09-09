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
unique(model_df2$scientificName) %>% length() # 217 spp

## remove migratory species
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

model_df2 <- filter(model_df2, !scientificName %in% mig_spp)
unique(model_df2$scientificName) %>% length() # left with 194 spp


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
              (0 + prec_seas | scientificName),
            data = model_df3, REML = FALSE, 
            lmerControl(optimizer = "bobyqa"))

car::vif(mon)
car::vif(mon2)
MuMIn::Weights(AIC(mon, mon2))

step(mon)
#Model found:
#  onset ~ temp + prec + temp_seas + (1 | id_cells) + (1 | scientificName) + 
#  (0 + temp | scientificName) + (0 + prec | scientificName) + 
#  (0 + temp_seas | scientificName) + (0 + prec_seas | scientificName) + 
#  (0 + temp:prec | scientificName) + temp:prec

## final model

m_on_final <- lmer(onset ~ temp + prec +  temp_seas + 
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp:prec | scientificName) +
                temp:prec,
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

m_on_traits_s = step(m_on_traits, reduce.random = T)
m_on_traits_s
# Model found:
#   onset ~ temp + prec + temp_seas + diapause.stage + immature.habitat + 
#   (1 | id_cells) + (1 | scientificName) + (0 + temp | scientificName) + 
#   (0 + prec | scientificName) + (0 + temp_seas | scientificName) + 
#   (0 + temp:prec | scientificName) + temp:prec + temp_seas:diapause.stage + 
#   temp:immature.habitat + prec:immature.habitat

m_on_traits_final <- lmer(onset ~ temp + prec + temp_seas + temp:prec +
                       diapause.stage + immature.habitat +
                       (1 | id_cells) + (1 | scientificName) + 
                       (0 + temp | scientificName) + 
                       (0 + prec | scientificName) +
                       (0 + temp_seas | scientificName) + 
                       (0 + temp:prec | scientificName) + 
                       temp_seas:diapause.stage + 
                       temp:flights +
                       temp:immature.habitat + 
                       prec:immature.habitat,
                     data = model_df4, REML = FALSE, 
                     lmerControl(optimizer = "bobyqa"))

# First check vifs to see if final model will stand
car::vif(m_on_traits_final) ## Generalized Collinearity Diagnostics (Fox and Monette 1992)
# looks good so now make sure we have the top model & no more interactions are needed
m_on_traits_s2 = step(m_on_traits_final, reduce.random = T)
m_on_traits_s2 # looks like temp:flights is not needed interaction remove and test again

m_on_traits_final <- lmer(onset ~ temp + prec + temp_seas + temp:prec +
                            diapause.stage + immature.habitat +
                            (1 | id_cells) + (1 | scientificName) + 
                            (0 + temp | scientificName) + 
                            (0 + prec | scientificName) +
                            (0 + temp_seas | scientificName) + 
                            (0 + temp:prec | scientificName) + 
                            temp_seas:diapause.stage + 
                            temp:immature.habitat + 
                            prec:immature.habitat,
                          data = model_df4, REML = FALSE, 
                          lmerControl(optimizer = "bobyqa"))

m_on_traits_s3 = step(m_on_traits_final, reduce.random = T)
m_on_traits_s3 # looks like we found the top model

summary(m_on_traits_final)
car::Anova(m_on_traits_final) # temp_seas:diapuase stage & temp:diapause stage are sig interactions
car::vif(m_on_traits_final) ## Generalized Collinearity Diagnostics (Fox and Monette 1992)
MuMIn::r.squaredGLMM(m_on_traits_final)
plot_model(m_on_traits_final, type = "pred", terms = c("temp_seas", "diapause.stage"), ci.lvl = NA)
plot_model(m_on_traits_final, type = "pred", terms = c("temp", "immature.habitat"), ci.lvl = NA)
plot_model(m_on_traits_final, type = "pred", terms = c("prec", "immature.habitat"), ci.lvl = NA)

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

#Model found:
#  offset ~ prec + temp_seas + (1 | id_cells) + (1 | scientificName) + 
#  (0 + temp | scientificName) + (0 + prec_seas | scientificName) + 
#  (0 + temp_seas | scientificName) + (0 + temp:prec | scientificName)

# Final model without traits!
m_off_final <-lmer(offset ~ prec +  temp_seas +
              (1|id_cells) + (1|scientificName) +
              (0 + prec | scientificName) +
              (0 + temp_seas | scientificName),
               data = model_df3, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(m_off_final)
car::vif(m_off_final)
suppressWarnings(MuMIn::r.squaredGLMM(m_off_final))

### Extract Random term estimates
off_coef <- coef(m_off_final)$scientificName %>% 
  rownames_to_column("scientificName") %>% 
  rename(intercept_ave_offset = "(Intercept)") %>% 
  left_join(spp_traits, by = "scientificName") %>% 
  left_join(spp_seas2, by = "scientificName")

# clean off_coef
off_coef <- off_coef %>% 
  na.omit %>% 
  filter(diapause.stage != "None")


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
                      seas + prec:seas + temp_seas:seas +
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

# stepwise regression to select model
m_off_traits_s = step(m_off_traits, reduce.random = T)
m_off_traits_s

# Model found:
#  offset ~ prec + temp_seas + development + diapause.stage + flights + 
#  immature.habitat + larval.diet + (1 | id_cells) + (1 | scientificName) + 
#  (0 + prec | scientificName) + (0 + temp_seas | scientificName) + 
#  temp_seas:development + prec:diapause.stage + prec:immature.habitat + 
#  temp_seas:immature.habitat + prec:larval.diet

# Final offset model with traits
m_off_traits_final <- lmer(offset ~ prec +  temp_seas + 
                    seas + development +  flights + immature.habitat +
                    (1|id_cells) + (1|scientificName) +
                    (0 + temp_seas | scientificName) + 
                    (0 + prec | scientificName) + 
                    prec:seas +
                    temp_seas:seas +
                    temp_seas:development +
                    prec:immature.habitat + 
                    temp_seas:immature.habitat,
                    data = model_df4, REML = FALSE, 
                    lmerControl(optimizer = "bobyqa"))
#first see if VIFs are okay
car::vif(m_off_traits_final) # they are not. temp_seas:development is inflated, so removed

# Final offset model with traits
m_off_traits_final <- lmer(offset ~ prec +  temp_seas + 
                             seas + development +  flights + immature.habitat +
                             (1|id_cells) + (1|scientificName) +
                             (0 + temp_seas | scientificName) + 
                             (0 + prec | scientificName) + 
                             prec:seas +
                             temp_seas:seas +
                             prec:immature.habitat + 
                             temp_seas:immature.habitat,
                           data = model_df4, REML = FALSE, 
                           lmerControl(optimizer = "bobyqa"))
# check vifs again
car::vif(m_off_traits_final) # we're good. Double check best model

# stepwise regression to select model
m_off_traits_s2 = step(m_off_traits_final, reduce.random = T)
m_off_traits_s2

#Model found:
#  offset ~ prec + temp_seas + seas + flights + immature.habitat + 
#  (1 | id_cells) + (1 | scientificName) + 
# (0 + temp_seas | scientificName) + (0 + prec | scientificName) + prec:seas + 
#  temp_seas:seas + temp_seas:immature.habitat

m_off_traits_final <- lmer(offset ~ prec + temp_seas +
                             seas + flights + immature.habitat + 
                             (1|id_cells) + (1|scientificName) +
                             (0 + temp_seas | scientificName) + 
                             (0 + prec | scientificName) + 
                             prec:seas + 
                             temp_seas:seas + 
                             temp_seas:immature.habitat,
                           data = model_df4, REML = FALSE, 
                           lmerControl(optimizer = "bobyqa"))

# check vifs again
car::vif(m_off_traits_final)

summary(m_off_traits_final)
car::vif(m_off_traits_final)
car::Anova(m_off_traits_final)
MuMIn::r.squaredGLMM(m_off_traits_final)
plot_model(m_off_traits_final, type = "pred", terms = c("temp_seas", "seas"), ci.lvl = NA)
plot_model(m_off_traits_final, type = "pred", terms = c("temp_seas", "immature.habitat"), ci.lvl = NA)
plot_model(m_off_traits_final, type = "pred", terms = c("prec", "seas"),ci.lvl = NA)

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

# Model found:
#   duration ~ temp + prec + temp_seas + (1 | id_cells) + (1 | scientificName) + 
#   (0 + temp | scientificName) + (0 + prec_seas | scientificName) + 
#   (0 + temp_seas | scientificName) + (0 + temp:prec | scientificName) + 
#   temp:prec

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
                       seas + prec:seas + temp_seas:seas + temp:seas + 
                       development + prec:development + temp_seas:development + temp:development +
                       diapause.stage + prec:diapause.stage + temp_seas:diapause.stage + temp:diapause.stage +
                       flights + prec:flights + temp_seas:flights + temp:flights +
                       immature.habitat + prec:immature.habitat + temp_seas:immature.habitat + temp:immature.habitat +
                       larval.diet +  prec:larval.diet + temp_seas:larval.diet + temp:larval.diet +
                       (1|id_cells) + (1|scientificName) +
                       (0 + prec | scientificName) +
                       (0 + temp | scientificName) + 
                       (0 + temp_seas | scientificName) +
                       (0 + temp:prec | scientificName),
                     data = model_df4, REML = FALSE, 
                     lmerControl(optimizer = "bobyqa"))

m_dur_traits_s = step(m_dur_traits, reduce.random = T)
m_dur_traits_s

m_dur_traits_final <- lmer(duration ~ temp + prec +  temp_seas + 
        seas + development + diapause.stage + flights + immature.habitat + larval.diet +
          (1 | id_cells) + (1 | scientificName) +
          (0 + prec | scientificName) +
          (0 + temp | scientificName) +
          (0 + temp_seas | scientificName) +
          (0 + temp:prec | scientificName) + 
          temp_seas:seas + temp:development + 
          prec:diapause.stage + temp:flights +
          temp_seas:immature.habitat + temp:immature.habitat + 
          prec:larval.diet,
          data = model_df4, REML = FALSE, 
          lmerControl(optimizer = "bobyqa"))
#check vifs
car::vif(m_dur_traits_final) # temp:development is inflated remove

m_dur_traits_final <- lmer(duration ~ temp + prec +  temp_seas + 
                             seas + development + diapause.stage + flights + immature.habitat + larval.diet +
                             (1 | id_cells) + (1 | scientificName) +
                             (0 + prec | scientificName) +
                             (0 + temp | scientificName) +
                             (0 + temp_seas | scientificName) +
                             (0 + temp:prec | scientificName) + 
                             temp_seas:seas + 
                             prec:diapause.stage + temp:flights +
                             temp_seas:immature.habitat + temp:immature.habitat + 
                             prec:larval.diet,
                           data = model_df4, REML = FALSE, 
                           lmerControl(optimizer = "bobyqa"))

#check vifs again
car::vif(m_dur_traits_final) #  development is inflated remove

m_dur_traits_final <- lmer(duration ~ temp + prec +  temp_seas + 
                             seas +  diapause.stage + flights + immature.habitat + larval.diet +
                             (1 | id_cells) + (1 | scientificName) +
                             (0 + prec | scientificName) +
                             (0 + temp | scientificName) +
                             (0 + temp_seas | scientificName) +
                             temp_seas:seas + 
                             prec:diapause.stage + temp:flights +
                             temp_seas:immature.habitat + temp:immature.habitat + 
                             prec:larval.diet,
                           data = model_df4, REML = FALSE, 
                           lmerControl(optimizer = "bobyqa"))

#check vifs again
car::vif(m_dur_traits_final) # good now 

## see if final model remains the same
m_dur_traits_s2 = step(m_dur_traits_final, reduce.random = T)
m_dur_traits_s2

#Model found:
#  duration ~ temp + prec + temp_seas + seas + diapause.stage + 
#  flights + immature.habitat + larval.diet + 
#  (1 | id_cells) + (1 | scientificName) + 
#  (0 + prec | scientificName) + (0 + temp | scientificName) + (0 + temp_seas | scientificName) + 
#  temp_seas:seas + prec:diapause.stage + temp:flights + temp_seas:immature.habitat + 
#  temp:immature.habitat + prec:larval.diet

m_dur_traits_final <- lmer(duration ~ temp + prec +  temp_seas + 
                             seas +  diapause.stage + flights + immature.habitat + larval.diet +
                             (1 | id_cells) + (1 | scientificName) +
                             (0 + temp | scientificName) +
                             (0 + temp_seas | scientificName) +
                             (0 + prec | scientificName) + 
                             temp_seas:seas + temp_seas:immature.habitat +
                             temp:flights + temp:immature.habitat + 
                             prec:larval.diet + prec:diapause.stage,
                           data = model_df4, REML = FALSE, 
                           lmerControl(optimizer = "bobyqa"))

summary(m_dur_traits_final)
car::vif(m_dur_traits_final)
car::Anova(m_dur_traits_final)
MuMIn::r.squaredGLMM(m_dur_traits_final)
plot_model(m_dur_traits_final, type = "pred", terms = c("temp_seas", "seas"), ci.lvl = NA)
plot_model(m_dur_traits_final, type = "pred", terms = c("temp_seas", "immature.habitat"), ci.lvl = NA)
plot_model(m_dur_traits_final, type = "pred", terms = c("temp", "flights"),ci.lvl = NA)
plot_model(m_dur_traits_final, type = "pred", terms = c("temp", "immature.habitat"),ci.lvl = NA)
plot_model(m_dur_traits_final, type = "pred", terms = c("prec", "diapause.stage"), ci.lvl = NA)
plot_model(m_dur_traits_final, type = "pred", terms = c("prec", "larval.diet"),ci.lvl = NA)

### NOW ADD PHYLOGENY ########

## try the ROTL package
library(rotl)
nmz <- tnrs_match_names(names = unique(model_df4$scientificName), context_name = "Animals")
nmz2 <- filter(nmz, ott_id != 7146697)
insect_tree <- tol_induced_subtree(ott_ids = nmz2$ott_id)
insect_tree$tip.label <- word(insect_tree$tip.label, start = 1, end = 2, sep = "_") %>% 
  sub(pattern = "_", replacement = " ")

plot(insect_tree, type = "fan")
plot(insect_tree)


insect_tree_bl <- ape::compute.brlen(insect_tree)

insect_tree_bl

on_coef
physig_col = function(inte){
  col_to = names(inte)[-1]
  out = vector("list", length = length(col_to))
  names(out) = col_to
  for(i in col_to){
    x = inte[[i]]
    names(x) = inte$scientificName
    xp1 = phytools::phylosig(insect_tree_bl, x, method = "K", test = T)
    xp2 = phytools::phylosig(insect_tree_bl, x, method = "lambda", test = T)
    out[[i]] = tibble::tibble(statistic = c(xp1$K, xp2$lambda),
                              P = c(ifelse(xp1$P < 0.5, xp1$P, 1 - xp1$P), xp2$P),
                              test = c("K", "Lambda"))
  }
  bind_rows(out, .id = "terms")
}

inte = left_join(dplyr::select(on_coef, scientificName, intercept_ave_onset),
                 dplyr::select(off_coef, scientificName, intercept_ave_offset)) %>% 
  left_join(dplyr::select(dur_coef, scientificName, intercept_ave_duration))

physig_intercept = physig_col(inte) %>% 
  arrange(test) %>% 
  mutate(statistic = round(statistic, 4),
         P = round(P, 3)) 
physig_intercept

physig_slopes = left_join(dplyr::select(on_coef, scientificName, temp_onset = temp, prec_onset = prec, 
                                        temp_seas_onset = temp_seas, tempprecint_onset = temp:prec),
                          dplyr::select(off_coef, scientificName, prec_offset = prec, temp_seas_offset = temp_seas,
                          )) %>% 
  left_join(dplyr::select(dur_coef, scientificName, temp_dur = temp, prec_dur = prec, 
                          temp_seas_dur = temp_seas, tempprecint_dur = temp:prec,
  )) %>% 
  physig_col() %>% 
  arrange(test) %>% 
  mutate(statistic = round(statistic, 4),
         P = round(P, 3)) 
as.data.frame(physig_slopes)



####### PGLMM
# onset 
library(phyr)
library(INLA)

# Burnsius communis not in phylogeny, remove
model_df5 <- filter(model_df4, scientificName %in% insect_tree_bl$tip.label)

a <- insect_tree_bl$tip.label %in% model_df5$scientificName
a[140]
a[32]

insect_tree_bl$tip.label[32]
insect_tree_bl$tip.label[140]


insect_tree_bl2 <- ape::drop.tip(insect_tree_bl, c("Enodia anthedon", "Aeshna cyanea")) # = (C:1,(D:1,E:1):1);


pm_onset <- pglmm(onset ~ temp + prec + temp_seas + temp:prec +
                    diapause.stage + immature.habitat +
                    (1 | id_cells) + (1 | scientificName__) + 
                    (0 + temp | scientificName__) + 
                    (0 + prec | scientificName__) +
                    (0 + temp_seas | scientificName__)  +    
                    temp_seas:diapause.stage + 
                    temp:immature.habitat + 
                    prec:immature.habitat,
                  data = model_df5, 
                  cov_ranef = list(scientificName = insect_tree_bl), 
                  bayes = TRUE)


ranef(pm_onset)
fixef(pm_onset) %>% knitr::kable()
inla_onset = pm_onset$inla.model
summary(inla_onset)
inla_onset$summary.fixed # kld are small, good
plot(inla_onset$marginals.fixed$`(Intercept)`)
plot(inla_onset$marginals.fixed$temp)
plot(inla_onset$marginals.fixed$`temp:prec`)
# 
names(inla_onset$marginals.fixed)
names(inla_onset$marginals.hyperpar)
length(inla_onset$summary.random)
inla_onset$marginals.random

install.packages("remotes")
remotes::install_github("julianfaraway/brinla")
library(brinla)

# bri.hyperpar.plot(inla_onset, F)
# bri.lmresid.plot(inla_onset)
inla_onset$marginals.fitted.values
invsqrt <- function(x) 1/sqrt(x)
invsqrt(inla_onset$summary.hyperpar[, -2]) # SD of random terms
# bri.hyperpar.summary(inla_onset)
inla_onset$marginals.random # species level random term
# bri.random.plot(inla_onset)
# 

# PGLMM OFFSET #
pm_offset <- pglmm(offset ~ prec + temp_seas + seas + flights + immature.habitat +
                    (1 | id_cells) + (1 | scientificName__) + 
                    (0 + prec | scientificName__) +
                    (0 + temp_seas | scientificName__)  +    
                     prec:seas + 
                     temp_seas:seas + 
                     temp_seas:immature.habitat,
                  data = model_df5, 
                  cov_ranef = list(scientificName = insect_tree_bl), 
                  bayes = TRUE)

ranef(pm_offset)
fixef(pm_offset) %>% knitr::kable()


# PGLMM Duration #
pm_dur <- pglmm(duration ~ temp + prec + temp_seas + 
                     seas + diapause.stage + flights + immature.habitat + larval.diet +
                     (1 | id_cells) + (1 | scientificName__) + 
                     (0 + temp | scientificName__) + 
                     (0 + prec | scientificName__) +
                     (0 + temp_seas | scientificName__)  +    
                     temp_seas:seas + 
                     temp_seas:immature.habitat + 
                     temp:flights + 
                     temp:immature.habitat + 
                     prec:larval.diet + 
                     prec:diapause.stage,
                   data = model_df5, 
                   cov_ranef = list(scientificName = insect_tree_bl), 
                   bayes = TRUE)

ranef(pm_dur)
fixef(pm_dur) %>% knitr::kable()
