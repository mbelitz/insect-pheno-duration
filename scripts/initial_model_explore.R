library(lme4)
library(lmerTest)
library(dplyr)
library(effects)
library(ggplot2)
library(car)

model_df <- read.csv(file = "data/model_dfs/duration_climate_population_data.csv",
                     stringsAsFactors = FALSE)

model_df2 <- model_df %>% 
  na.omit() %>% 
  mutate(temp = scale(temp),
         prec = scale(prec),
         pop = scale(pop),
         prec_seas = scale(bio15),
         temp_seas = scale(bio4),
         onsetci = onset_high - onset_low,
         offsetci = offset_high - offset_low)

m1on <- lmer(onset ~ temp + pop + prec + prec_seas + temp_seas + (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + (0 + pop | scientificName) + (0 + prec | scientificName) +
               (0 + temp_seas | scientificName) + (0 + prec_seas | scientificName),
             data = model_df2, REML = FALSE)

summary(m1on) ## model failed to converge but usful in showing temp and prec_seas as most important predictors

m2on <- lmer(onset ~ temp:prec + prec + prec_seas + pop + (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + (0 + pop | scientificName) + (0 + prec | scientificName) +
               (0 + prec_seas | scientificName) +
               (0 + temp:prec | scientificName), REML = F,
             data = model_df2)

summary(m2on)

m3on <- lmer(onset ~ temp:pop + prec:temp + temp + prec + pop + (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + (0 + pop | scientificName) + (0 + prec | scientificName) +
               (0 + temp:pop | scientificName) + (0 + prec:temp | scientificName), REML = F,
             data = model_df2)

summary(m3on) ## Failed to converge 

m3on <- lmer(onset ~ temp + prec_seas + pop + (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + (0 + pop | scientificName) + (0 + prec_seas | scientificName)
              , REML = F,
             data = model_df2)

summary(m3on) 

m4on <- lmer(onset ~ temp + prec_seas + (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + prec_seas | scientificName),
             REML = F,
             data = model_df2)

summary(m4on) 

m5on <- lmer(onset ~ temp + prec_seas + temp:prec_seas + (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + 
               (0 + prec_seas | scientificName) + (0 + temp:prec_seas),
             REML = F,
             data = model_df2)

summary(m5on) 

# model rank top models 
MuMIn::AICc(m2on, m3on, m4on, m5on)
MuMIn::Weights(AIC(m2on, m3on, m4on, m5on))

summary(m3on)
MuMIn::r.squaredGLMM(m3on)
car::vif(m3on)
class(m3on) = "lmerMod"
performance::check_collinearity(m3on)
performance::model_performance(m3on)
performance::check_normality(m3on)

### Now try offset

m1off <- lmer(offset ~ temp + pop + prec + prec_seas + temp_seas + (1|id_cells) + (1|scientificName) +
               (0 + temp | scientificName) + (0 + pop | scientificName) + (0 + prec | scientificName) +
                (0 + prec_seas | scientificName) + (0 + temp_seas | scientificName), REML = FALSE, 
             data = model_df2)

summary(m1off) ## Did not converged using to check correlation and examine what variables are most important
# removing prec_seas

m1off <- lmer(offset ~ temp + pop + prec + temp_seas + (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + (0 + pop | scientificName) + (0 + prec | scientificName) +
                (0 + temp_seas | scientificName), REML = FALSE, 
              data = model_df2)

summary(m1off) ## still no convergence

m1off <- lmer(offset ~ temp + pop + temp_seas + (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + (0 + pop | scientificName) +
                (0 + temp_seas | scientificName), REML = FALSE, 
              data = model_df2)

summary(m1off)

m1.5off <- lmer(offset ~ temp + pop + temp_seas + scale(onset) + (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + (0 + pop | scientificName) +
                (0 + temp_seas | scientificName), REML = FALSE, 
              data = model_df2)

summary(m1.5off)

m2off <- lmer(offset ~ temp + pop + prec + (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + (0 + pop | scientificName) +
                (0 + prec | scientificName), REML = FALSE, 
              data = model_df2)

summary(m2off)

m3off <- lmer(offset ~ temp + pop + temp:pop + (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + (0 + pop | scientificName) +
                (0 + temp:pop | scientificName), REML = FALSE, 
              data = model_df2)

summary(m3off)

m4off <- lmer(offset ~ temp + pop + (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + (0 + pop | scientificName), 
              REML = FALSE, 
              data = model_df2)

summary(m4off)

m5off <- lmer(offset ~ temp + temp_seas + 
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + temp_seas | scientificName), REML = FALSE, 
              data = model_df2)

summary(m5off)

# model rank top models 
AIC(m1off, m1.5off, m2off, m3off, m4off, m5off)
MuMIn::Weights(AIC(m1off, m1.5off, m2off, m3off, m4off, m5off))

summary(m1.5off)
MuMIn::r.squaredGLMM(m1.5off)
car::vif(m1off)
class(m1off) = "lmerMod"
performance::check_collinearity(m1off)
performance::model_performance(m1off)
performance::check_normality(m1off)

## Duration 

### Now try offset

m1dur <- lmer(duration ~ temp + pop + prec + temp_seas + prec_seas + scale(onset) +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName) +
                (0 + prec_seas), REML = FALSE, 
              data = model_df2)
summary(m1dur)

m2dur <- lmer(duration ~ temp + pop + prec + temp_seas + scale(onset) +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + prec | scientificName) +
                (0 + temp_seas | scientificName), REML = FALSE, 
              data = model_df2)
summary(m2dur) ## Failed to converge

m2dur <- lmer(duration ~ temp + pop + temp_seas + scale(onset) +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + temp_seas | scientificName), REML = FALSE, 
              data = model_df2)
summary(m2dur)

m3dur <- lmer(duration ~ temp + pop + temp:pop + temp_seas + scale(onset) +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + temp_seas | scientificName)+
                (0 + temp:pop | scientificName), REML = FALSE, 
              data = model_df2)
summary(m3dur) ## Failed to converge

m3dur <- lmer(duration ~ temp + pop + temp:pop + scale(onset) +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName) + 
                (0 + temp:pop | scientificName), REML = FALSE, 
              data = model_df2)
summary(m3dur) 

m4dur <- lmer(duration ~ temp + pop + scale(onset) +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName) + 
                (0 + pop | scientificName), REML = FALSE, 
              data = model_df2)
summary(m4dur) 

m5dur <- lmer(duration ~ temp + scale(onset) +
                (1|id_cells) + (1|scientificName) +
                (0 + temp | scientificName), REML = FALSE, 
              data = model_df2)
summary(m5dur) 

m6dur <- lmer(duration ~ pop  + scale(onset) +
                (1|id_cells) + (1|scientificName) +
                (0 + pop | scientificName), REML = FALSE, 
              data = model_df2)
summary(m6dur) 

# model rank top models 
AIC(m1dur, m2dur, m3dur, m4dur, m5dur, m6dur)
MuMIn::Weights(AIC(m1dur, m2dur, m3dur, m4dur, m5dur, m6dur))

summary(m2dur)
MuMIn::r.squaredGLMM(m2dur)
car::vif(m2dur)
class(m2dur) = "lmerMod"
performance::check_collinearity(m2dur)
performance::model_performance(m2dur)
performance::check_normality(m2dur)
