library(lme4)
library(lmerTest)
library(dplyr)
library(effects)
library(ggplot2)
library(car)

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

has_5_cells <- filter(datadens, count >= 5) %>% 
  filter (scientificName != "Apis mellifera") #234 species 

ggplot(has_5_cells, aes(x = Order)) + 
  geom_bar(stat = "count", aes (fill = Order)) +
  labs(y = "Number of Species") +
  ggtitle("At Least 5 Comparable Cells - 186 total spp") +
  theme_bw()


has_10_cells <- filter(datadens, count >= 10) %>% 
  filter (scientificName != "Apis mellifera") # 145

ggplot(has_10_cells, aes(x = Order)) + 
  geom_bar(stat = "count", aes (fill = Order)) +
  labs(y = "Number of Species") +
  ggtitle("At Least 10 Comparable Cells - 145 total spp") +
  theme_bw()

model_df2 <- filter(model_df2, scientificName %in% has_10_cells$scientificName)

odes <- filter(model_df2, Order == "Odonata")
leps <- filter(model_df2, Order == "Lepidoptera")
dips <- filter(model_df2, Order == "Diptera")
hyms <- filter(model_df2, Order == "Hymenoptera")
coles <- filter(model_df2, Order == "Coleoptera")

########## Onset Models ###########

ode1on <- lmer(onset ~ temp + pop + prec + prec_seas + temp_seas +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + pop | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + prec_seas | scientificName),
             data = odes, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(ode1on) 

ode2on <- lmer(onset ~ temp + pop + prec + prec_seas + temp_seas + temp:pop +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + pop | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + prec_seas | scientificName) + 
               (0 + temp:pop | scientificName),
             data = odes, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(ode2on) 


ode3on <- lmer(onset ~ temp + pop + prec_seas + temp_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + temp_seas | scientificName) + 
                 (0 + prec_seas | scientificName) ,
               data = odes, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(ode3on) 


ode4on <- lmer(onset ~ temp + pop + prec + temp_seas + temp:pop +
               (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + pop | scientificName) + 
               (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + 
               (0 + temp:pop | scientificName),
             data = odes, REML = FALSE, 
             lmerControl(optimizer = "bobyqa"))

summary(ode4on) 


ode5on <- lmer(onset ~ temp + pop + prec_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec_seas | scientificName) ,
               data = odes, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(ode5on) 


ode6on <- lmer(onset ~ temp + prec_seas + 
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + prec_seas | scientificName) ,
               data = odes, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(ode6on) 


# model rank top models 
AIC(ode1on, ode2on, ode3on, ode4on, ode5on, ode6on)
MuMIn::Weights(AIC(ode1on, ode2on, ode3on, ode4on, ode5on, ode6on))

summary(ode3on)
MuMIn::r.squaredGLMM(ode3on)
car::vif(ode3on)
performance::check_collinearity(ode3on)
performance::model_performance(ode3on)

# Leps onset

lep1on <- lmer(onset ~ temp + pop + prec + prec_seas + temp_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec | scientificName) +
                 (0 + temp_seas | scientificName) + 
                 (0 + prec_seas | scientificName),
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep1on) 

lep1.5on <- lmer(onset ~ temp + pop + prec + prec_seas + temp_seas + temp:prec +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec | scientificName) +
                 (0 + temp_seas | scientificName) + 
                 (0 + prec_seas | scientificName) + 
                 (0 + temp:prec | scientificName),
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep1.5on) 


lep2on <- lmer(onset ~ temp + pop + prec + prec_seas + temp_seas + temp:pop +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec | scientificName) +
                 (0 + temp_seas | scientificName) + 
                 (0 + prec_seas | scientificName) + 
                 (0 + temp:pop | scientificName),
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep2on) 


lep3on <- lmer(onset ~ temp + pop + prec_seas + temp_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + temp_seas | scientificName) + 
                 (0 + prec_seas | scientificName) ,
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep3on) 


lep4on <- lmer(onset ~ temp + pop + prec + temp_seas + temp:pop +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec | scientificName) +
                 (0 + temp_seas | scientificName) + 
                 (0 + temp:pop | scientificName),
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep4on) 


lep5on <- lmer(onset ~ temp + pop + prec_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec_seas | scientificName) ,
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep5on) 


lep6on <- lmer(onset ~ temp + prec_seas + 
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + prec_seas | scientificName) ,
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep6on) 


# model rank top models 
AIC(lep1on, lep1.5on, lep2on, lep3on, lep4on, lep5on, lep6on)
MuMIn::Weights(AIC(lep1on, lep1.5on, lep2on, lep3on, lep4on, lep5on, lep6on))

summary(lep1.5on)
plot(effects::effect("temp:prec", lep1.5on), multiline = T)
MuMIn::r.squaredGLMM(lep1.5on)
car::vif(lep1.5on)
performance::check_collinearity(lep1.5on)


#### Hymenoptera

#hym1on <- lmer(onset ~ temp + pop + prec + prec_seas + temp_seas +
#                 (1|id_cells) + (1|scientificName) +
#                 (0 + temp | scientificName) + 
#                 (0 + pop | scientificName) + 
#                 (0 + prec | scientificName) +
#                 (0 + temp_seas | scientificName) + 
#                 (0 + prec_seas | scientificName),
#               data = hyms, REML = FALSE, 
#               lmerControl(optimizer = "bobyqa"))
#
#summary(hym1on) #isSingular

#hym1on <- lmer(onset ~ temp + pop + prec + temp_seas +
#                 (1|id_cells) + (1|scientificName) +
#                 (0 + temp | scientificName) + 
#                 (0 + pop | scientificName) + 
#                 (0 + prec | scientificName) +
#                 (0 + temp_seas | scientificName),
#               data = hyms, REML = FALSE, 
#               lmerControl(optimizer = "bobyqa"))
#
#summary(hym1on) # is singluar

hym1on <- lmer(onset ~ temp + pop +  temp_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + temp_seas | scientificName),
               data = hyms, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(hym1on) 

hym2on <- lmer(onset ~ temp + pop +  temp_seas + temp:pop +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + temp_seas | scientificName) + 
                 (0 + temp:pop | scientificName),
               data = hyms, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(hym2on) 

hym3on <- lmer(onset ~pop +  temp_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + pop | scientificName) + 
                 (0 + temp_seas | scientificName),
               data = hyms, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(hym3on) 


# Model rank top models 
AIC(hym1on, hym2on, hym3on)
MuMIn::Weights(AIC(hym1on, hym2on, hym3on))

summary(hym1on)
MuMIn::r.squaredGLMM(hym1on)
car::vif(hym1on)
performance::check_collinearity(hym1on)
performance::model_performance(hym1on)

# Coleoptera

# cole1on <- lmer(onset ~ temp + pop + prec + prec_seas + temp_seas +
#                  (1|id_cells) + (1|scientificName) +
#                  (0 + temp | scientificName) + 
#                  (0 + pop | scientificName) + 
#                  (0 + prec | scientificName) +
#                  (0 + temp_seas | scientificName) + 
#                  (0 + prec_seas | scientificName),
#                data = coles, REML = FALSE, 
#                lmerControl(optimizer = "bobyqa"))
# 
# summary(cole1on) # isSingluar

#cole1on <- lmer(onset ~ temp + pop + prec + prec_seas + 
#                (1|id_cells) + (1|scientificName) +
#                (0 + temp | scientificName) + 
#                (0 + pop | scientificName) + 
#                (0 + prec | scientificName) +
#                (0 + prec_seas | scientificName),
#              data = coles, REML = FALSE, 
#              lmerControl(optimizer = "bobyqa"))
#
#summary(cole1on) # isSingluar

# cole1on <- lmer(onset ~ temp + pop + prec_seas + 
#                   (1|id_cells) + (1|scientificName) +
#                   (0 + temp | scientificName) + 
#                   (0 + pop | scientificName) + 
#                   (0 + prec_seas | scientificName),
#                 data = coles, REML = FALSE, 
#                 lmerControl(optimizer = "bobyqa"))
# 
# summary(cole1on) # isSingluar

cole1on <- lmer(onset ~ temp + temp:pop + 
                  (1|id_cells) + (1|scientificName) +
                  (0 + temp | scientificName) + 
                  (0 + temp:pop | scientificName),
                data = coles, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(cole1on) 

cole2on <- lmer(onset ~ temp + 
                  (1|id_cells) + (1|scientificName) +
                  (0 + temp | scientificName) ,
                data = coles, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(cole2on) 



# mcolel rank top mcolels 
AIC(cole1on, cole2on)
MuMIn::Weights(AIC(cole1on, cole2on))

summary(cole2on)
MuMIn::r.squaredGLMM(cole2on)
performance::model_performance(cole2on)


# Diptera 

# dip1on <- lmer(onset ~ temp + pop + prec + prec_seas + temp_seas +
#                  (1|id_cells) + (1|scientificName) +
#                  (0 + temp | scientificName) + 
#                  (0 + pop | scientificName) + 
#                  (0 + prec | scientificName) +
#                  (0 + temp_seas | scientificName) + 
#                  (0 + prec_seas | scientificName),
#                data = dips, REML = FALSE, 
#                lmerControl(optimizer = "bobyqa"))
# 
# summary(dip1on) # isSingular

# dip1on <- lmer(onset ~ temp + prec + prec_seas + temp_seas +
#                  (1|id_cells) + (1|scientificName) +
#                  (0 + temp | scientificName) + 
#                  (0 + prec | scientificName) +
#                  (0 + temp_seas | scientificName) + 
#                  (0 + prec_seas | scientificName),
#                data = dips, REML = FALSE, 
#                lmerControl(optimizer = "bobyqa"))
# 
# summary(dip1on) # isSummary

dip1on <- lmer(onset ~ temp + pop + prec + prec_seas + 
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec | scientificName) + 
                 (0 + prec_seas | scientificName),
               data = dips, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(dip1on) 

dip2on <- lmer(onset ~ temp + prec + prec_seas + 
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + prec | scientificName) + 
                 (0 + prec_seas | scientificName),
               data = dips, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(dip2on) 

dip3on <- lmer(onset ~ temp + prec +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + prec | scientificName),
               data = dips, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(dip3on) 

dip4on <- lmer(onset ~ temp + prec + temp:prec + 
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + prec | scientificName),
               data = dips, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(dip4on) 

# mcolel rank top mcolels 
AIC(dip1on, dip2on, dip3on, dip4on)
MuMIn::Weights(AIC(dip1on, dip2on, dip3on, dip4on))

summary(dip3on)
MuMIn::r.squaredGLMM(dip3on)
performance::model_performance(dip3on)


##### Offset

########## Offset Models ###########

# ode1off <- lmer(offset ~ temp + pop + prec + prec_seas + temp_seas +
#                  (1|id_cells) + (1|scientificName) +
#                  (0 + temp | scientificName) + 
#                  (0 + pop | scientificName) + 
#                  (0 + prec | scientificName) +
#                  (0 + temp_seas | scientificName) + 
#                  (0 + prec_seas | scientificName),
#                data = odes, REML = FALSE, 
#                lmerControl(optimizer = "bobyqa"))
# 
# summary(ode1off) # is Singular

ode1off <- lmer(offset ~ pop + prec + prec_seas + temp_seas + 
                 (1|id_cells) + (1|scientificName) +
                 (0 + pop | scientificName) + 
                 (0 + prec | scientificName) +
                 (0 + temp_seas | scientificName) + 
                 (0 + prec_seas | scientificName),
               data = odes, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(ode1off) 

ode2off <- lmer(offset ~prec + prec_seas + temp_seas + 
                  (1|id_cells) + (1|scientificName) +
                  (0 + prec | scientificName) +
                  (0 + temp_seas | scientificName) + 
                  (0 + prec_seas | scientificName),
                data = odes, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(ode2off) 

ode3off <- lmer(offset ~ prec + temp_seas + 
                  (1|id_cells) + (1|scientificName) +
                  (0 + prec | scientificName) +
                  (0 + temp_seas | scientificName),
                data = odes, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(ode3off) 

ode4off <- lmer(offset ~ temp_seas + 
                  (1|id_cells) + (1|scientificName)  +
                  (0 + temp_seas | scientificName),
                data = odes, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(ode4off) 

# model rank top models 
AIC(ode1off, ode2off, ode3off, ode4off)
MuMIn::Weights(AIC(ode1off, ode2off, ode3off, ode4off))

summary(ode1off)
MuMIn::r.squaredGLMM(ode1off)
car::vif(ode1off)
performance::check_collinearity(ode1off)
performance::model_performance(ode1off)

# Leps offset

lep1off <- lmer(offset ~ temp + pop + prec + prec_seas + temp_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec | scientificName) +
                 (0 + temp_seas | scientificName) + 
                 (0 + prec_seas | scientificName),
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep1off) 

lep1.5off <- lmer(offset ~ temp + pop + prec + prec_seas + temp_seas + temp:prec +
                   (1|id_cells) + (1|scientificName) +
                   (0 + temp | scientificName) + 
                   (0 + pop | scientificName) + 
                   (0 + prec | scientificName) +
                   (0 + temp_seas | scientificName) + 
                   (0 + prec_seas | scientificName) + 
                   (0 + temp:prec | scientificName),
                 data = leps, REML = FALSE, 
                 lmerControl(optimizer = "bobyqa"))

summary(lep1.5off) 


lep2off <- lmer(offset ~ temp + pop + prec + prec_seas + temp_seas + temp:pop +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec | scientificName) +
                 (0 + temp_seas | scientificName) + 
                 (0 + prec_seas | scientificName) + 
                 (0 + temp:pop | scientificName),
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep2off) 


lep3off <- lmer(offset ~ temp + prec + prec_seas + temp_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + prec | scientificName) + 
                 (0 + temp_seas | scientificName) + 
                 (0 + prec_seas | scientificName) ,
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep3off) 


lep4off <- lmer(offset ~ temp + prec + prec_seas + temp_seas + temp:prec +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + prec | scientificName) +
                 (0 + prec_seas | scientificName) +
                 (0 + temp_seas | scientificName) + 
                 (0 + temp:prec | scientificName),
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep4off) 

lep5off <- lmer(offset ~ temp + prec + temp_seas + temp:prec +
                  (1|id_cells) + (1|scientificName) +
                  (0 + temp | scientificName) + 
                  (0 + prec | scientificName) +
                  (0 + temp_seas | scientificName) + 
                  (0 + temp:prec | scientificName),
                data = leps, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(lep5off) 

lep6off <- lmer(offset ~ temp + prec + temp_seas +
                  (1|id_cells) + (1|scientificName) +
                  (0 + temp | scientificName) + 
                  (0 + prec | scientificName) +
                  (0 + temp_seas | scientificName),
                data = leps, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(lep6off) 

# model rank top models 
AIC(lep1off, lep1.5off, lep2off, lep3off, lep4off, lep5off, lep6off)
MuMIn::Weights(AIC(lep1off, lep1.5off, lep2off, lep3off, lep4off, lep5off, lep6off))

summary(lep1.5off)
plot(effects::effect("temp:prec", lep1.5off), multiline = T)
MuMIn::r.squaredGLMM(lep1.5off)
car::vif(lep1.5off)
performance::check_collinearity(lep1.5off)


#### Hymenoptera

#hym1off <- lmer(offset ~ temp + pop + prec + prec_seas + temp_seas +
#                 (1|id_cells) + (1|scientificName) +
#                 (0 + temp | scientificName) + 
#                 (0 + pop | scientificName) + 
#                 (0 + prec | scientificName) +
#                 (0 + temp_seas | scientificName) + 
#                 (0 + prec_seas | scientificName),
#               data = hyms, REML = FALSE, 
#               lmerControl(optimizer = "bobyqa"))
#
#summary(hym1off) #isSingular

#hym1off <- lmer(offset ~ temp + pop + prec_seas + temp_seas +
#                (1|id_cells) + (1|scientificName) +
#                (0 + temp | scientificName) + 
#                (0 + pop | scientificName) + 
#                (0 + prec_seas | scientificName) +
#                (0 + temp_seas | scientificName),
#              data = hyms, REML = FALSE, 
#              lmerControl(optimizer = "bobyqa"))
#
#summary(hym1off) # is singluar

# hym1off <- lmer(offset ~ temp + prec_seas +  temp_seas +
#                  (1|id_cells) + (1|scientificName) +
#                  (0 + temp | scientificName) + 
#                  (0 + prec_seas | scientificName) + 
#                  (0 + temp_seas | scientificName),
#                data = hyms, REML = FALSE, 
#                lmerControl(optimizer = "bobyqa"))
# 
# summary(hym1off) 

# hym1off <- lmer(offset ~ temp +  temp_seas + 
#                  (1|id_cells) + (1|scientificName) +
#                  (0 + temp | scientificName) + 
#                  (0 + temp_seas | scientificName) ,
#                data = hyms, REML = FALSE, 
#                lmerControl(optimizer = "bobyqa"))
# 
# summary(hym1off) # isSingular

hym1off <- lmer(offset ~ temp +  
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) ,
               data = hyms, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(hym1off) 


# Model rank top models 
AIC(hym1off)
MuMIn::Weights(AIC(hym1off))

summary(hym1off)
MuMIn::r.squaredGLMM(hym1off)
car::vif(hym1off)
performance::check_collinearity(hym1off)
performance::model_performance(hym1off)

# Coleoptera

# cole1off <- lmer(offset ~ temp + pop + prec + prec_seas + temp_seas +
#                  (1|id_cells) + (1|scientificName) +
#                  (0 + temp | scientificName) + 
#                  (0 + pop | scientificName) + 
#                  (0 + prec | scientificName) +
#                  (0 + temp_seas | scientificName) + 
#                  (0 + prec_seas | scientificName),
#                data = coles, REML = FALSE, 
#                lmerControl(optimizer = "bobyqa"))
# 
# summary(cole1off) # isSingluar

#cole1off <- lmer(offset ~ temp + pop + prec + prec_seas + 
#                (1|id_cells) + (1|scientificName) +
#                (0 + temp | scientificName) + 
#                (0 + pop | scientificName) + 
#                (0 + prec | scientificName) +
#                (0 + prec_seas | scientificName),
#              data = coles, REML = FALSE, 
#              lmerControl(optimizer = "bobyqa"))
#
#summary(cole1off) # isSingluar

#cole1off <- lmer(offset ~ temp + pop + prec_seas + 
#                   (1|id_cells) + (1|scientificName) +
#                   (0 + temp | scientificName) + 
#                   (0 + pop | scientificName) + 
#                   (0 + prec_seas | scientificName),
#                 data = coles, REML = FALSE, 
#                 lmerControl(optimizer = "bobyqa"))
# 
# summary(cole1off) # isSingluar

#cole1off <- lmer(offset ~ temp + temp:pop + 
#                  (1|id_cells) + (1|scientificName) +
#                  (0 + temp | scientificName) + 
#                  (0 + temp:pop | scientificName),
#                data = coles, REML = FALSE, 
#                lmerControl(optimizer = "bobyqa"))
#
#summary(cole1off) # is singular

cole1off <- lmer(offset ~ temp + 
                  (1|scientificName) +
                  (0 + temp | scientificName) ,
                data = coles, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(cole1off) 



# mcolel rank top mcolels 
AIC(cole1off)
MuMIn::Weights(AIC(cole1off))

summary(cole1off)
MuMIn::r.squaredGLMM(cole1off)
performance::model_performance(cole1off)


# Diptera 

# dip1off <- lmer(offset ~ temp + pop + prec + prec_seas + temp_seas +
#                  (1|id_cells) + (1|scientificName) +
#                  (0 + temp | scientificName) + 
#                  (0 + pop | scientificName) + 
#                  (0 + prec | scientificName) +
#                  (0 + temp_seas | scientificName) + 
#                  (0 + prec_seas | scientificName),
#                data = dips, REML = FALSE, 
#                lmerControl(optimizer = "bobyqa"))
# 
# summary(dip1off) # isSingular

#dip1off <- lmer(offset ~ temp + prec + prec_seas + temp_seas +
#                 (1|id_cells) + (1|scientificName) +
#                 (0 + temp | scientificName) + 
#                 (0 + prec | scientificName) +
#                 (0 + temp_seas | scientificName) + 
#                 (0 + prec_seas | scientificName),
#               data = dips, REML = FALSE, 
#               lmerControl(optimizer = "bobyqa"))
#
#summary(dip1off) # isSummary

#dip1off <- lmer(offset ~ temp + pop + prec + prec_seas + 
#                 (1|id_cells) + (1|scientificName) +
#                 (0 + temp | scientificName) + 
#                 (0 + pop | scientificName) + 
#                 (0 + prec | scientificName) + 
#                 (0 + prec_seas | scientificName),
#               data = dips, REML = FALSE, 
#               lmerControl(optimizer = "bobyqa"))
#
#summary(dip1off) 
#
#dip1off <- lmer(offset ~ temp + prec + prec_seas + 
#                 (1|id_cells) + (1|scientificName) +
#                 (0 + temp | scientificName) + 
#                 (0 + prec | scientificName) + 
#                 (0 + prec_seas | scientificName),
#               data = dips, REML = FALSE, 
#               lmerControl(optimizer = "bobyqa"))
#
#summary(dip1off) 

dip1off <- lmer(offset ~ prec_seas +
                 (1|scientificName) ,
               data = dips, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(dip1off) 

dip4off <- lmer(offset ~ prec_seas + 
                 
                 (0 + prec_seas | scientificName),
               data = dips, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(dip4off) 

# mcolel rank top mcolels 
AIC(dip1off, dip4off)
MuMIn::Weights(AIC(dip1off,dip4off))

summary(dip1off)
MuMIn::r.squaredGLMM(dip1off)
performance::model_performance(dip1off)


########## Duration Models ###########

ode1dur <- lmer(duration ~ temp + pop + prec + prec_seas + temp_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec | scientificName) +
                 (0 + temp_seas | scientificName) + 
                 (0 + prec_seas | scientificName),
               data = odes, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(ode1dur) 

ode2dur <- lmer(duration ~ temp + temp_seas + temp:prec +
                  (1|id_cells) + (1|scientificName) +
                  (0 + temp | scientificName) + 
                  (0 + temp_seas | scientificName) + 
                  (0 + temp:prec | scientificName),
                data = odes, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(ode2dur) 


ode3dur <- lmer(duration ~ temp + pop + prec_seas + temp_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + temp_seas | scientificName) + 
                 (0 + prec_seas | scientificName) ,
               data = odes, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(ode3dur) 


ode4dur <- lmer(duration ~ temp + pop + prec + temp_seas + temp:pop +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec | scientificName) +
                 (0 + temp_seas | scientificName) + 
                 (0 + temp:pop | scientificName),
               data = odes, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(ode4dur) 


ode5dur <- lmer(duration ~ temp + pop + prec_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec_seas | scientificName) ,
               data = odes, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(ode5dur) 


ode6dur <- lmer(duration ~ temp + prec_seas + 
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + prec_seas | scientificName) ,
               data = odes, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(ode6dur) 


# model rank top models 
AIC(ode1dur, ode2dur, ode3dur, ode4dur, ode5dur, ode6dur)
MuMIn::Weights(AIC(ode1dur, ode2dur, ode3dur, ode4dur, ode5dur, ode6dur))

summary(ode1dur)
MuMIn::r.squaredGLMM(ode1dur)
car::vif(ode1dur)
performance::check_collinearity(ode1dur)
performance::model_performance(ode1dur)

# Leps duration

lep1dur <- lmer(duration ~ temp + pop + prec + prec_seas + temp_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec | scientificName) +
                 (0 + temp_seas | scientificName) + 
                 (0 + prec_seas | scientificName),
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep1dur) 

lep1.5dur <- lmer(duration ~ temp + pop + prec + prec_seas + temp_seas + temp:prec +
                   (1|id_cells) + (1|scientificName) +
                   (0 + temp | scientificName) + 
                   (0 + pop | scientificName) + 
                   (0 + prec | scientificName) +
                   (0 + temp_seas | scientificName) + 
                   (0 + prec_seas | scientificName) + 
                   (0 + temp:prec | scientificName),
                 data = leps, REML = FALSE, 
                 lmerControl(optimizer = "bobyqa"))

summary(lep1.5dur) 

lep2dur <- lmer(duration ~ temp + pop + prec + prec_seas + temp_seas + temp:pop +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec | scientificName) +
                 (0 + temp_seas | scientificName) + 
                 (0 + prec_seas | scientificName) + 
                 (0 + temp:pop | scientificName),
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep2dur) 

lep3dur <- lmer(duration ~ temp + pop + prec_seas + temp_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + temp_seas | scientificName) + 
                 (0 + prec_seas | scientificName) ,
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep3dur) 

lep4dur <- lmer(duration ~ temp + pop + prec + temp_seas + temp:pop +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec | scientificName) +
                 (0 + temp_seas | scientificName) + 
                 (0 + temp:pop | scientificName),
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep4dur) 

lep5dur <- lmer(duration ~ temp + pop + prec_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + prec_seas | scientificName) ,
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep5dur) 

lep6dur <- lmer(duration ~ temp + prec_seas + 
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + prec_seas | scientificName) ,
               data = leps, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(lep6dur) 

# model rank top models 
AIC(lep1dur, lep1.5dur, lep2dur, lep3dur, lep4dur, lep5dur, lep6dur)
MuMIn::Weights(AIC(lep1dur, lep1.5dur, lep2dur, lep3dur, lep4dur, lep5dur, lep6dur))

summary(lep1.5dur)
plot(effects::effect("temp:prec", lep1.5dur), multiline = T)
MuMIn::r.squaredGLMM(lep1.5dur)
car::vif(lep1.5dur)
performance::check_collinearity(lep1.5dur)


#### Hymenoptera

# hym1dur <- lmer(duration ~ temp + pop + prec + prec_seas + temp_seas +
#                  (1|id_cells) + (1|scientificName) +
#                  (0 + temp | scientificName) + 
#                  (0 + pop | scientificName) + 
#                  (0 + prec | scientificName) +
#                  (0 + temp_seas | scientificName) + 
#                  (0 + prec_seas | scientificName),
#                data = hyms, REML = FALSE, 
#                lmerControl(optimizer = "bobyqa"))
# 
# summary(hym1dur) #isSingular

# hym1dur <- lmer(duration ~ temp + pop + prec_seas + temp_seas +
#                  (1|id_cells) + (1|scientificName) +
#                  (0 + temp | scientificName) + 
#                  (0 + pop | scientificName) + 
#                  (0 + prec_seas | scientificName) +
#                  (0 + temp_seas | scientificName),
#                data = hyms, REML = FALSE, 
#                lmerControl(optimizer = "bobyqa"))
# 
# summary(hym1dur) # is singluar

hym1dur <- lmer(duration ~ temp + pop +  temp_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName) + 
                 (0 + temp_seas | scientificName),
               data = hyms, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(hym1dur) 

hym2dur <- lmer(duration ~ temp + pop +
                 (1|id_cells) + (1|scientificName) +
                 (0 + temp | scientificName) + 
                 (0 + pop | scientificName),
               data = hyms, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(hym2dur) 

hym3dur <- lmer(duration ~pop +  temp_seas +
                 (1|id_cells) + (1|scientificName) +
                 (0 + pop | scientificName) + 
                 (0 + temp_seas | scientificName),
               data = hyms, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(hym3dur) 

hym4dur <- lmer(duration ~ temp +  temp:pop +
                  (1|id_cells) + (1|scientificName) +
                  (0 + temp | scientificName) + 
                  (0 + temp:pop | scientificName),
                data = hyms, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(hym4dur) 

# Model rank top models 
AIC(hym1dur, hym2dur, hym3dur, hym4dur)
MuMIn::Weights(AIC(hym1dur, hym2dur, hym3dur, hym4dur))

summary(hym2dur)
MuMIn::r.squaredGLMM(hym2dur)
car::vif(hym2dur)
performance::check_collinearity(hym2dur)
performance::model_performance(hym2dur)

# Coleoptera

# cole1dur <- lmer(duration ~ temp + pop + prec + prec_seas + temp_seas +
#                  (1|id_cells) + (1|scientificName) +
#                  (0 + temp | scientificName) + 
#                  (0 + pop | scientificName) + 
#                  (0 + prec | scientificName) +
#                  (0 + temp_seas | scientificName) + 
#                  (0 + prec_seas | scientificName),
#                data = coles, REML = FALSE, 
#                lmerControl(optimizer = "bobyqa"))
# 
# summary(cole1dur) # isSingluar

# cole1dur <- lmer(duration ~ temp_seas + pop + prec + prec_seas + 
#                 (1|id_cells) + (1|scientificName) +
#                 (0 + temp_seas | scientificName) + 
#                 (0 + pop | scientificName) + 
#                 (0 + prec | scientificName) +
#                 (0 + prec_seas | scientificName),
#               data = coles, REML = FALSE, 
#               lmerControl(optimizer = "bobyqa"))
# 
# summary(cole1dur) # isSingluar

# cole1dur <- lmer(duration ~ temp_seas + pop + prec_seas + 
#                   (1|id_cells) + (1|scientificName) +
#                   (0 + temp_seas | scientificName) + 
#                   (0 + pop | scientificName) + 
#                   (0 + prec_seas | scientificName),
#                 data = coles, REML = FALSE, 
#                 lmerControl(optimizer = "bobyqa"))
# 
# summary(cole1dur) # isSingluar

cole1dur <- lmer(duration ~ temp + temp:pop + 
                  (1|id_cells) + (1|scientificName) +
                  (0 + temp | scientificName) + 
                  (0 + temp:pop | scientificName),
                data = coles, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(cole1dur) 

cole2dur <- lmer(duration ~ temp + 
                  (1|id_cells) + (1|scientificName) +
                  (0 + temp | scientificName) ,
                data = coles, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(cole2dur) 

cole3dur <- lmer(duration ~ prec_seas  + 
                   (1|id_cells) + (1|scientificName) +
                   (0 + prec_seas | scientificName),
                 data = coles, REML = FALSE, 
                 lmerControl(optimizer = "bobyqa"))

summary(cole3dur) 


# mcolel rank top mcolels 
AIC(cole1dur, cole2dur, cole3dur)
MuMIn::Weights(AIC(cole1dur, cole2dur, cole3dur))

summary(cole2dur)
MuMIn::r.squaredGLMM(cole2dur)
performance::model_performance(cole2dur)


# Diptera 

#dip1dur <- lmer(duration ~ temp + pop + prec + prec_seas + temp_seas +
#                 (1|id_cells) + (1|scientificName) +
#                 (0 + temp | scientificName) + 
#                 (0 + pop | scientificName) + 
#                 (0 + prec | scientificName) +
#                 (0 + temp_seas | scientificName) + 
#                 (0 + prec_seas | scientificName),
#               data = dips, REML = FALSE, 
#               lmerControl(optimizer = "bobyqa"))
#
#summary(dip1dur) # isSingular

# dip1dur <- lmer(duration ~ temp + prec + prec_seas + temp_seas +
#                  (1|id_cells) + (1|scientificName) +
#                  (0 + temp | scientificName) + 
#                  (0 + prec | scientificName) +
#                  (0 + temp_seas | scientificName) + 
#                  (0 + prec_seas | scientificName),
#                data = dips, REML = FALSE, 
#                lmerControl(optimizer = "bobyqa"))
# 
# summary(dip1dur) # isSummary

# dip1dur <- lmer(duration ~ temp + prec + temp_seas + 
#                  (1|id_cells) + (1|scientificName) +
#                  (0 + temp | scientificName) + 
#                  (0 + pop | scientificName) + 
#                  (0 + prec | scientificName) + 
#                  (0 + temp_seas | scientificName),
#                data = dips, REML = FALSE, 
#                lmerControl(optimizer = "bobyqa"))
# 
# summary(dip1dur) 

# dip1dur <- lmer(duration ~ temp + prec + 
#                   (1|id_cells) + (1|scientificName) +
#                   (0 + temp | scientificName) + 
#                   (0 + prec | scientificName),
#                 data = dips, REML = FALSE, 
#                 lmerControl(optimizer = "bobyqa"))
# 
# summary(dip1dur) 

dip1dur <- lmer(duration ~ pop  +
                  (1|scientificName) +
                 (0 + pop | scientificName),
               data = dips, REML = FALSE, 
               lmerControl(optimizer = "bobyqa"))

summary(dip1dur) 

dip2dur <- lmer(duration ~ prec  +
                  (0 + prec | scientificName),
                data = dips, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(dip2dur) 

dip3dur <- lmer(duration ~ temp  +
                  (0 + temp | scientificName),
                data = dips, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(dip3dur) 

dip4dur <- lmer(duration ~ temp_seas  +
                  (0 + temp_seas | scientificName),
                data = dips, REML = FALSE, 
                lmerControl(optimizer = "bobyqa"))

summary(dip4dur) 


# mcolel rank top mcolels 
AIC(dip1dur, dip2dur, dip3dur, dip4dur)
MuMIn::Weights(AIC(dip1dur, dip2dur, dip3dur, dip4dur))

summary(dip3dur)
MuMIn::r.squaredGLMM(dip3dur)
performance::model_performance(dip3dur)




