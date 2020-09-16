library(tidyverse)
library(gridExtra)

# Interaction Response Curves
mdf <- readr::read_csv("PGLMM_data/model_df_5cells_treespp.csv")

# Onset interaction response curves

tmp.x <- mdf$temp

## pop interactions
pop_high <- 0.95
pop <- -0.02
pop_low <- -1

onset.low = 117.4678473 + 
  (-15.3536044 * tmp.x) + 
  (-2.5310345 * pop_low) + 
  (-3.0480112 * pop_low * tmp.x)

onset = 117.4678473 + 
  (-15.3536044 * tmp.x) + 
  (-2.5310345 * pop) + 
  (-3.0480112 * pop * tmp.x)

onset.high = 117.4678473 + 
  (-15.3536044 * tmp.x) + 
  (-2.5310345 * pop_high) + 
  (-3.0480112 * pop_high * tmp.x)

tmpz <- c(tmp.x, tmp.x, tmp.x)
onset_pop <- c(onset.high, onset, onset.low)
popz <- c(rep(.95, 3171), rep(-0.02, 3171), rep(-1, 3171))

pop_temp <- data.frame(temp = tmpz, onset = onset_pop, Pop = as.factor(popz)) 

tpop <- ggplot() + 
  geom_line(pop_temp, mapping = aes(x = temp, y = onset, color = Pop), size = 1.05) +
  labs(x = "Temp", y = "Onset") + 
  scale_color_manual(values = c("turquoise4", "grey5", "violetred2")) +
  theme_bw()
   
tpop

## Precipitation

## prec interactions
prec_high <- 0.99
prec <- 0.09
prec_low <- -0.81

onset.low.prec = 117.4678473 + 
  (-15.3536044 * tmp.x) + 
  (0.7572964 * prec_low) + 
  (-2.5246943 * prec_low * tmp.x)

onset.prec = 117.4678473 + 
  (-15.3536044 * tmp.x) + 
  (0.7572964 * prec) + 
  (-2.5246943 * prec * tmp.x)

onset.high.prec = 117.4678473 + 
  (-15.3536044 * tmp.x) + 
  (0.7572964 * prec_high) + 
  (-2.5246943 * prec_high * tmp.x)

tmpz <- c(tmp.x, tmp.x, tmp.x)
onset_prec <- c(onset.high.prec, onset.prec, onset.low.prec)
precz <- c(rep(.99, 3171), rep(0.09, 3171), rep(-0.81, 3171))

prec_temp <- data.frame(temp = tmpz, onset = onset_prec, Prec = as.factor(precz)) 


tprec <- ggplot() + 
  geom_line(prec_temp, mapping = aes(x = temp, y = onset_prec, color = Prec), size = 1.05) +
  labs(x = "Temp", y = "Onset") + 
  scale_color_manual(values = c("brown", "cyan", "Blue")) +
  theme_bw()
   

tprec

## pop and immature habitat

## immature habitat interactions
Fresh.Water <- 2.262
Nest <- -6.611
Underground <- -3.045

pop.x2 <- mdf$pop

onset_fw <- 117.4678473 + 
  (-2.531 * pop.x2) + 
  (-0.206) + 
  (Fresh.Water * pop.x2)

onset_nest <- 117.4678473 + 
  (-2.531 * pop.x2) + 
  (9.645) + 
  (Nest * pop.x2)

onset_underground <- 117.4678473 + 
  (-2.531 * pop.x2) + 
  (8.818) + 
  (Underground * pop.x2)

onset_above.ground.vegetation <- 117.4678473 + 
  (-2.531 * pop.x2) 

popz2 <- c(pop.x2, pop.x2, pop.x2, pop.x2)
onset_ih <- c(onset_above.ground.vegetation, onset_fw, onset_nest, onset_underground)
Immature.Habitat <- c(rep("Above Ground Vegetation", 3171), rep("Freshwater", 3171), 
                      rep("Nest", 3171), rep("Underground", 3171))

pop_immature.habitat <- data.frame(popz2 = popz2, onset = onset_ih, Immature.Habitat = as.factor(Immature.Habitat)) 


pih <- ggplot() + 
  geom_line(pop_immature.habitat, mapping = aes(x = popz2, y = onset, color = Immature.Habitat), size = 1.05) +
  labs(x = "Pop", y = "Onset") + 
  scale_color_viridis_d() +
  theme_bw()

onset_interactions <- egg::ggarrange(tpop, tprec, pih, ncol = 1,
                                     labels = c("A"," B", "C"))

ggsave(plot = onset_interactions, filename = "Tables&Figures/onset_interactions.png",
       dpi = 400, width = 8, height = 8)

### Offset

## temp Seas and immature habitat

## immature habitat interactions & tempseas
Fresh.Water.ts <- -4.571766
Nest.ts <- 8.441451
Underground.ts <- 11.794879

tmpseas <- mdf$temp_seas

offset_fw_ts <- 281.049541 + 
  (-9.336082 * tmpseas) + 
  (14.61098) + 
  (Fresh.Water.ts * tmpseas)

offset_nest_ts <- 281.049541 + 
  (-9.336082 * tmpseas) + 
  (3.404545) + 
  (Nest.ts * tmpseas)

offset_underground_ts <- 281.049541 + 
  (-9.336082 * tmpseas) + 
  (-17.517044) + 
  (Underground.ts * tmpseas)

offset_above.ground.vegetation_ts <- 281.049541 + 
  (-9.336082 * tmpseas) 

tmpseas.x <- c(tmpseas, tmpseas,tmpseas,tmpseas)
offset_ih_ts <- c(offset_above.ground.vegetation_ts, offset_fw_ts, offset_nest_ts, offset_underground_ts)
Immature.Habitat <- c(rep("Above Ground Vegetation", 3171), rep("Freshwater", 3171), 
                      rep("Nest", 3171), rep("Underground", 3171))

tempseas_immature.habitat <- data.frame(tempseaz = tmpseas.x, offset = offset_ih_ts, 
                                        Immature.Habitat = as.factor(Immature.Habitat)) 


tsih <- ggplot() + 
  geom_line(tempseas_immature.habitat, mapping = aes(x = tempseaz, y = offset, color = Immature.Habitat), size = 1.05) +
  labs(x = "Temp Seas", y = "Offset") + 
  scale_color_viridis_d() +
  theme_bw()
   

tsih


## immature habitat interactions x prec
Fresh.Water.prec <- -2.920
Nest.prec <- -3.798
Underground.prec <- -10.138

pcip <- mdf$prec

offset_fw_prec <- 281.049541 + 
  (3.657138 * pcip) + 
  (14.61098) + 
  (Fresh.Water.prec * pcip)

offset_nest_prec <- 281.049541 + 
  (3.657138 * pcip) + 
  (3.404545) + 
  (Nest.prec * pcip)

offset_underground_prec <- 281.049541 + 
  (3.657138 * pcip) + 
  (-17.517044) + 
  (Underground.prec * pcip)

offset_above.ground.vegetation_prec <- 281.049541 + 
  (3.657138 * pcip) 

prec.x <- c(pcip, pcip, pcip, pcip)
offset_ih_prec <- c(offset_above.ground.vegetation_prec, offset_fw_prec, offset_nest_prec, offset_underground_prec)
Immature.Habitat <- c(rep("Above Ground Vegetation", 3171), rep("Freshwater", 3171), 
                      rep("Nest", 3171), rep("Underground", 3171))

prec_immature.habitat <- data.frame(precz = prec.x, offset = offset_ih_prec, 
                                    Immature.Habitat = as.factor(Immature.Habitat)) 


prec_ih <- ggplot() + 
  geom_line(prec_immature.habitat, mapping = aes(x = precz, y = offset, color = Immature.Habitat), size = 1.05) +
  labs(x = "Prec", y = "Offset") + 
  scale_color_viridis_d() +
  theme_bw()
   

prec_ih

offset_interactions <- egg::ggarrange(tsih, prec_ih, ncol = 1,
                                      labels = c("A"," B"))

ggsave(plot = offset_interactions, filename = "Tables&Figures/offset_interactions.png",
       dpi = 400, width = 8, height = 6)


## Duration

## Temp:pop

tmp.x <- mdf$temp

## pop interactions
pop_high <- 0.95
pop <- -0.02
pop_low <- -1

duration.low.pop = 139.1467457 + 
  (16.1542697 * tmp.x) + 
  (2.3289379 * pop_low) + 
  (3.3671402 * pop_low * tmp.x)

duration.pop = 139.1467457 + 
  (16.1542697 * tmp.x) + 
  (2.3289379 * pop) + 
  (3.3671402 * pop * tmp.x)

duration.high.pop = 139.1467457 + 
  (16.1542697 * tmp.x) + 
  (2.3289379 * pop_high) + 
  (3.3671402 * pop_high * tmp.x)

tmpz <- c(tmp.x, tmp.x, tmp.x)
duration_pop <- c(duration.high.pop, duration.pop, duration.low.pop)
popz <- c(rep(.95, 3171), rep(-0.02, 3171), rep(-1, 3171))

pop_temp_dur <- data.frame(temp = tmpz, onset = duration_pop, Pop = as.factor(popz)) 


dur_pop <- ggplot() + 
  geom_line(pop_temp_dur, mapping = aes(x = temp, y = onset, color = Pop), size = 1.05) +
  labs(x = "Temp", y = "Duration") + 
  scale_color_manual(values = c("turquoise4", "grey5", "violetred2")) +
  theme_bw()
   
dur_pop

# Temp Seas : Activity Seas
Spring.ts <- -12.1436149
Summer.ts <- -3.1365185

tmpseas <- mdf$temp_seas

dur_Spring.ts <- 139.1467457 +  
  (-4.2899232 * tmpseas) + 
  (10.8594987) + 
  (Spring.ts * tmpseas)

dur_Summer_ts <- 139.1467457 + 
  (-4.2899232 * tmpseas) +
  (4.1971951) + 
  (Summer.ts * tmpseas)

dur_Fall_ts <- 139.1467457 + 
  (-4.2899232 * tmpseas)

tmpseas.x.seas <- c(tmpseas, tmpseas,tmpseas)
dur_seas_ts <- c(dur_Spring.ts, dur_Summer_ts, dur_Fall_ts)
Seasonality <- c(rep("Spring", 3171), rep("Summer", 3171), 
                 rep("Fall", 3171))

tempseas_seas_dur <- data.frame(tempseaz = tmpseas.x.seas, duration = dur_seas_ts, 
                                Seasonality = as.factor(Seasonality)) 

tmpseas_seas <- ggplot() + 
  geom_line(tempseas_seas_dur, mapping = aes(x = tempseaz, y = duration, color = Seasonality), size = 1.05) +
  labs(x = "Temp Seas", y = "Duration") + 
  theme_bw()
   
tmpseas_seas

# Temp:Voltinism
univoltine <- -10.2772051

tmp.x <- mdf$temp

dur_univoltine_temp <- 139.1467457 +  
  (16.1542697 * tmp.x) + 
  (-11.8102269) + 
  (univoltine * tmp.x)

dur_not_univoltine <- 139.1467457 + 
  (16.15542697 * tmp.x)

tmp_volt.x <- c(tmp.x, tmp.x)
dur_tmp_vol <- c(dur_univoltine_temp, dur_not_univoltine)
Voltinism <- c(rep("Univoltine", 3171), rep("Not Univoltine", 3171))

temp_voltine_dur <- data.frame(Temp = tmp_volt.x, duration = dur_tmp_vol, 
                               Voltinism = as.factor(Voltinism)) 

tvd <- ggplot() + 
  geom_line(temp_voltine_dur, mapping = aes(x = Temp, y = duration, color = Voltinism), size = 1.05) +
  labs(x = "Temp", y = "Duration") + 
  scale_color_manual(values = c("maroon", "gold")) +
  theme_bw()

tvd

## Temp Seas : Voltinism

univoltine.ts <- -9.5115716

ts.x <- mdf$temp_seas

dur_univoltine_ts <- 139.1467457 +  
  (-4.2899232 * ts.x) + 
  (-11.8102269) + 
  (univoltine.ts * ts.x)

dur_not_univoltine_ts <- 139.1467457 + 
  (-4.2899232 * ts.x)

ts_volt.x <- c(ts.x, ts.x)
dur_tmpseas_vol <- c(dur_univoltine_ts, dur_not_univoltine_ts)
Voltinism <- c(rep("Univoltine", 3171), rep("Not Univoltine", 3171))

ts_voltine_dur <- data.frame(Temp_Seas = ts_volt.x, duration = dur_tmpseas_vol, 
                             Voltinism = as.factor(Voltinism)) 

tsvd <- ggplot() + 
  geom_line(ts_voltine_dur, mapping = aes(x = Temp_Seas, y = duration, color = Voltinism), size = 1.05) +
  labs(x = "Temp Seas", y = "Duration") + 
  scale_color_manual(values = c("maroon", "gold")) +
  theme_bw()
   

tsvd

## pop and immature habitat for duration

## immature habitat interactions
Fresh.Water_dur <- 3.019070
Nest_dur <- 9.0204673   
Underground_dur <- 11.4712092

pop.x2 <- mdf$pop

dur_fw <- 139.1467457 + 
  (2.3289379 * pop.x2) + 
  (4.2725748) + 
  (Fresh.Water_dur * pop.x2)

dur_nest <- 139.1467457 + 
  (2.3289379 * pop.x2) + 
  (-2.6734721) + 
  (Nest_dur * pop.x2)

dur_underground <-  139.1467457 + 
  (2.3289379 * pop.x2) + 
  (-20.0519249) + 
  (Underground_dur * pop.x2)

dur_above.ground.vegetation <- 139.1467457 + 
  (2.3289379 * pop.x2) 

popz2 <- c(pop.x2, pop.x2, pop.x2, pop.x2)
dur_ih <- c(dur_above.ground.vegetation, dur_fw, dur_nest, dur_underground)
Immature.Habitat <- c(rep("Above Ground Vegetation", 3171), rep("Freshwater", 3171), 
                      rep("Nest", 3171), rep("Underground", 3171))

dur_pop_immature.habitat <- data.frame(popz2 = popz2, Duration = dur_ih, Immature.Habitat = as.factor(Immature.Habitat)) 


dur_pih <- ggplot() + 
  geom_line(dur_pop_immature.habitat, mapping = aes(x = popz2, y = Duration, color = Immature.Habitat), size = 1.05) +
  labs(x = "Pop", y = "Duration") + 
  scale_color_viridis_d() +
  theme_bw()
   

dur_pih


# temp:larval diet
detritivorous <- 25.7738613
live.plants <- 0.4720430   

tmp.x2 <- mdf$temp

dur_detriv <- 139.1467457 + 
  (16.1542697 * tmp.x2) + 
  (30.1744343) + 
  (detritivorous * tmp.x2)

dur_live.plants <- 139.1467457 + 
  (16.1542697 * tmp.x2) +
  (-10.2299286) + 
  (live.plants * pop.x2)

dur_carnivorous <-  139.1467457 + 
  (16.1542697 * tmp.x2) 

tmpz2 <- c(tmp.x2, tmp.x2, tmp.x2)
dur_diet <- c(dur_detriv, dur_live.plants, dur_carnivorous)
Immature.Diet <- c(rep("Detritivorous", 3171), rep("Herbivorous", 3171), 
                   rep("Carnivorous", 3171))

dur_tmp_larval.diet <- data.frame(Temp = tmpz2, Duration = dur_diet, Immature.Diet = as.factor(Immature.Diet)) 


dur_ld <- ggplot() + 
  geom_line(dur_tmp_larval.diet, mapping = aes(x = Temp, y = Duration, color = Immature.Diet), size = 1.05) +
  labs(x = "Temp", y = "Duration") + 
  scale_color_viridis_d(option = "plasma") +
  theme_bw()
   
dur_ld

dur_interactions <- egg::ggarrange(dur_pop, tmpseas_seas, tvd, 
                                   tsvd, dur_pih, dur_ld, 
                                   ncol = 2,
                                   labels = c("A"," B", "C", "D", "E", "F")
                                  )

ggsave(plot = dur_interactions, filename = "Tables&Figures/duration_interactions.png",
       dpi = 400, width = 10, height = 6)


## Grid arrange over a 6x3 matrix
##lay <- rbind(c(1,2,3),
#             c(4,5,6),
#             c(7, NA, 8),
#             c(NA, NA, 9),
#             c(NA, NA, 10),
#             c(NA, NA, 11))
#
#
#
##ga <- grid.arrange(arrangeGrob(tpop, top = "Onset"), arrangeGrob(tsih, top = "Offset"), 
#             arrangeGrob(dur_pop, top = "Duration"), tprec, prec_ih,
#             tmpseas_seas, pih, tvd, tsvd, dur_pih, dur_ld, layout_matrix = lay)
#
#ggsave(plot = ga, filename = "Tables&Figures/InteractiveEffects.png", dpi = 400, 
#       width = 12, height = 10)

