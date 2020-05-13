library(lme4)
library(lmerTest)
library(dplyr)
library(effects)
library(ggplot2)
library(car)

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
         temp_seas = scale(bio4),
         onsetci = onset_high - onset_low,
         offsetci = offset_high - offset_low)

## currently how many unique species are there?
unique(model_df2$scientificName) %>% length() # 369

datadens <- model_df2 %>% 
  group_by(scientificName, Order) %>% 
  summarise(count = n())

has_5_cells <- filter(datadens, count >= 5) %>% 
  filter (scientificName != "Apis mellifera")#300 species 

ggplot(has_5_cells, aes(x = Order)) + 
  geom_bar(stat = "count", aes (fill = Order)) +
  labs(y = "Number of Species") +
  ggtitle("At Least 5 Comparable Cells - 300 total spp") +
  theme_bw()

has_10_cells <- filter(datadens, count >= 10) %>% 
  filter (scientificName != "Apis mellifera") # 217

ggplot(has_10_cells, aes(x = Order)) + 
  geom_bar(stat = "count", aes (fill = Order)) +
  labs(y = "Number of Species") +
  ggtitle("At Least 10 Comparable Cells - 217 total spp") +
  theme_bw()

model_df2 <- model_df2 %>% 
  filter(scientificName %in% has_10_cells$scientificName)

unique(model_df2$scientificName) %>% length()

# onset models

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

m2on <- lmer(onset ~ temp + prec + temp_seas + 
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) +
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) ,
             data = model_df2, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(m2on) 

m3on <- lmer(onset ~ temp + prec + temp_seas + temp:prec +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) +
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) +
               (0 + temp:prec | scientificName),
             data = model_df2, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(m3on) 

m4on <- lmer(onset ~ temp + prec_seas + (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + prec_seas | scientificName),
             REML = F, 
             lmerControl(optimizer = "bobyqa"),
             data = model_df2)

summary(m4on) 

m5on <- lmer(onset ~ temp + temp_seas +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + temp_seas),
             REML = F, 
             lmerControl(optimizer = "bobyqa"),
             data = model_df2)

summary(m5on) 

m6on <- lmer(onset ~ temp + prec + 
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + prec | scientificName), 
             REML = F,
             lmerControl(optimizer = "bobyqa"),
             data = model_df2)

summary(m6on) 


m7on <- lmer(onset ~ temp + temp_seas + temp:temp_seas + 
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + temp_seas | scientificName) + 
               (0 + temp:temp_seas),
             REML = F, 
             lmerControl(optimizer = "bobyqa"),
             data = model_df2)

summary(m7on) 

m8on <- lmer(onset ~ temp + prec + prec_seas + temp_seas + 
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + prec_seas | scientificName),
             data = model_df2, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(m8on) 

m9on <- lmer(onset ~ temp + pop + prec + prec_seas + temp_seas + temp:pop +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + pop | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + prec_seas | scientificName) + 
               (0 + temp:pop | scientificName),
             data = model_df2, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(m9on) 

# model rank top models 
MuMIn::AICc(m1on, m2on, m3on, m4on, m5on, m6on, m7on, m8on, m9on)
MuMIn::Weights(AIC(m1on, m2on, m3on, m4on, m5on, m6on, m7on, m8on, m9on))

summary(m9on)
MuMIn::r.squaredGLMM(m9on)
car::vif(m9on)
class(m9on) = "lmerMod"
performance::check_collinearity(m9on)
performance::model_performance(m9on)


## Now try offset

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

m2off <- lmer(offset ~ temp + pop + prec + temp_seas +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName),
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m2off) 

m3off <- lmer(offset ~ pop + prec + temp_seas +
                (1|id_cells) + (1|scientificName) +
                (0 + pop | scientificName) + 
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName) + ,
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m3off) 

m4off <- lmer(offset ~ temp + pop + prec + temp_seas + temp:pop +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName) + 
                (0 + temp:pop | scientificName),
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m4off) 

m5off <- lmer(offset ~ temp + pop + prec + temp_seas + temp_seas:pop +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName) +
                (0 + temp_seas:pop | scientificName),
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m5off) 

m6off <- lmer(offset ~ temp + pop + prec + temp_seas + temp_seas:prec +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName) + +
                (0 + temp_seas:prec | scientificName),
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m6off) 

m7off <- lmer(offset ~ temp + pop + prec + prec_seas + temp_seas + temp:pop +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName) + 
                (0 + prec_seas | scientificName) +
                (0 + temp:pop | scientificName),
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m7off) 

# model rank top models 
MuMIn::AICc(m1off, m2off, m3off, m4off, m5off, m6off, m7off)
MuMIn::Weights(AIC(m1off, m2off, m3off, m4off, m5off, m6off, m7off))

summary(m7off)
MuMIn::r.squaredGLMM(m7off)
car::vif(m7off)
class(m7off) = "lmerMod"
performance::check_collinearity(m7off)
performance::model_performance(m7off)


### Duration Models ### Giddy Up ###

m1dur <- lmer(duration ~ temp + pop + prec + prec_seas + temp_seas +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName) + 
                (0 + prec_seas | scientificName),
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m1dur) 

m2dur <- lmer(duration ~ temp + pop + prec +  temp_seas +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                
                (0 + temp_seas | scientificName) + 
                (0 + prec_seas | scientificName),
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m2dur) 

m3dur <- lmer(duration ~ temp + pop + temp_seas +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                
                (0 + temp_seas | scientificName),
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m3dur) 

m4dur <- lmer(duration ~ temp + pop + temp_seas + temp:pop +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + temp:pop | scientificName) +
                (0 + temp_seas | scientificName),
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m4dur) 

m5dur <- lmer(duration ~ temp + pop + temp_seas + temp:pop + prec_seas + 
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) +
                (0 + prec_seas | scientificName) +
                (0 + temp:pop | scientificName) +
                (0 + temp_seas | scientificName),
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m5dur) 


m6dur <- lmer(duration ~ temp + pop + temp_seas + temp:pop + prec_seas + prec +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) +
                (0 + prec_seas | scientificName) +
                (0 + temp:pop | scientificName) +
                (0 + temp_seas | scientificName) +
                (0 + prec | scientificName),
              data = model_df2, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(m6dur) 
plot(effects::effect("temp:pop", m6dur), multiline = T)
# model rank top models 
MuMIn::AICc(m1dur, m2dur, m3dur, m4dur, m5dur, m6dur)
MuMIn::Weights(AIC(m1dur, m2dur, m3dur, m4dur, m5dur, m6dur))

summary(m6dur)
MuMIn::r.squaredGLMM(m6dur)
car::vif(m6dur)
class(m6dur) = "lmerMod"
performance::check_collinearity(m6dur)
performance::model_performance(m6dur)

### Check out full models for each species

cole <- model_df2 %>% 
  filter(Order == "Coleoptera")

dip <- model_df2 %>% 
  filter(Order == "Diptera")

hemi <- model_df2 %>% 
  filter(Order == "Hemiptera")

hym <- model_df2 %>% 
  filter(Order == "Hympenoptera")

leps <- model_df2 %>% 
  filter(Order == "Lepidoptera")

odes <- model_df2 %>% 
  filter(Order == "Odonata")


## beetles first
cole_on <- lmer(onset ~ temp + pop + prec +
                   (1|id_cells) + (1|scientificName),
                 data = model_df2, REML = FALSE, 
                 lmerControl(optimizer = "bobyqa"))

summary(cole_on) 

cole_offset <- lmer(offset ~ temp + pop + prec +
                   (1|id_cells) + (1|scientificName),
                 data = cole, REML = FALSE, 
                 lmerControl(optimizer = "bobyqa"))

summary(cole_off) 

cole_dur <- lmer(duration ~ temp + pop + prec +
                (1|id_cells) + (1|scientificName),
              data = cole, REML = FALSE, 
              lmerControl(optimizer = "bobyqa"))

summary(cole_dur) 

# model rank top models 
MuMIn::r.squaredGLMM(cole_dur)
MuMIn::r.squaredGLMM(cole_on)
MuMIn::r.squaredGLMM(cole_off) 

## flies next
dip_on <- lmer(onset ~ temp + pop + prec +
                   (1|id_cells) + (1|scientificName),
                 data = dip, REML = FALSE, 
                 lmerControl(optimizer = "bobyqa"))

summary(dip_on) 

dip_off <- lmer(offset ~ temp + pop + prec +
                   (1|id_cells) + (1|scientificName),
                 data = dip, REML = FALSE, 
                 lmerControl(optimizer = "bobyqa"))

summary(dip_off) 

dip_dur <- lmer(duration ~ temp + pop + prec +
                   (1|id_cells) + (1|scientificName),
                 data = dip, REML = FALSE, 
                 lmerControl(optimizer = "bobyqa"))

summary(dip_dur) 
# model rank top models 
MuMIn::r.squaredGLMM(dip_dur)
MuMIn::r.squaredGLMM(dip_on)
MuMIn::r.squaredGLMM(dip_off) 



## bees next
hym_on <- lmer(onset ~ temp + pop + prec +
                 (1|id_cells) + (1|scientificName),
               data = hym, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(hym_on) 

hym_off <- lmer(offset ~ temp + pop + prec +
                  (1|id_cells) + (1|scientificName),
                data = hym, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(hym_off) 

hym_dur <- lmer(duration ~ temp + pop + prec +
                  (1|id_cells) + (1|scientificName),
                data = hym, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(hym_dur) 
# model rank top models 
MuMIn::r.squaredGLMM(hym_dur)
MuMIn::r.squaredGLMM(hym_on)
MuMIn::r.squaredGLMM(hym_off) 



## leps 
leps_on <- lmer(onset ~ temp + pop + prec +
                 (1|id_cells) + (1|scientificName),
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(leps_on) 

leps_off <- lmer(offset ~ temp + pop + prec +
                  (1|id_cells) + (1|scientificName),
                data = leps, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(leps_off) 

leps_dur <- lmer(duration ~ temp + pop + prec +
                  (1|id_cells) + (1|scientificName),
                data = leps, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(leps_dur) 
# model rank top models 
MuMIn::r.squaredGLMM(leps_dur)
MuMIn::r.squaredGLMM(leps_on)
MuMIn::r.squaredGLMM(leps_off) 



## DRAGONflies next
odes_on <- lmer(onset ~ temp + pop + prec +
                 (1|id_cells) + (1|scientificName),
               data = odes, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(odes_on) 

odes_off <- lmer(offset ~ temp + pop + prec +
                  (1|id_cells) + (1|scientificName),
                data = odes, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(odes_off) 

odes_dur <- lmer(duration ~ temp + pop + prec +
                  (1|id_cells) + (1|scientificName),
                data = odes, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(odes_dur) 
# model rank top models 
MuMIn::r.squaredGLMM(odes_dur)
MuMIn::r.squaredGLMM(odes_on)
MuMIn::r.squaredGLMM(odes_off) 