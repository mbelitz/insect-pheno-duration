library(tidyverse)
library(ggeffects)

# Interaction Response Curves
mdf <- readr::read_csv("PGLMM_data/model_df_5cells_treespp.csv")

common_spp <- mdf %>% 
  group_by(scientificName, Order) %>% 
  summarise(count = n())

# most common spp
get_temp_coef_onset <- function(binomial){
  df <- dplyr::filter(mdf, scientificName == binomial)
  mod <- lm(onset ~ temp, data = df)
  coef <- unname(mod$coefficients['temp'])
  return(coef)
}

get_temp_coef_se_onset <- function(binomial){
  df <- dplyr::filter(mdf, scientificName == binomial)
  mod <- lm(onset ~ temp, data = df)
  se <- coef(summary(mod))[,"Std. Error"]['temp'] 
  return(se)
}

spp_list <- c("Pachydiplax longipennis",
              "Bombus impatiens",
              "Papilio glaucus",
              "Diabrotica undecimpunctata",
              "Delphinia picta",
              "Neotibicen superbus")

model_coef <- lapply(X = spp_list, FUN = get_temp_coef_onset)
model_se <- lapply(X = spp_list, FUN = get_temp_coef_se_onset)

temp_coef_df_onset <- data.frame(
  scientificName = c("Pachydiplax longipennis",
                     "Bombus impatiens",
                     "Papilio glaucus",
                     "Diabrotica undecimpunctata",
                     "Delphinia picta",
                     "Neotibicen superbus"),
  
  coef = c(model_coef[[1]],
           model_coef[[2]],
           model_coef[[3]],
           model_coef[[4]],
           model_coef[[5]],
           model_coef[[6]]),
  
  high_se =c(model_coef[[1]] + (2 * model_se[[1]]),
             model_coef[[2]] + (2 * model_se[[2]]),
             model_coef[[3]] + (2 * model_se[[3]]),
             model_coef[[4]] + (2 * model_se[[4]]),
             model_coef[[5]] + (2 * model_se[[5]]),
             model_coef[[6]] + (2 * model_se[[6]])),
  
   low_se =c(model_coef[[1]] - (2 * model_se[[1]]),
             model_coef[[2]] - (2 * model_se[[2]]),
             model_coef[[3]] - (2 * model_se[[3]]),
             model_coef[[4]] - (2 * model_se[[4]]),
             model_coef[[5]] - (2 * model_se[[5]]),
             model_coef[[6]] - (2 * model_se[[6]]))
) %>% 
  mutate(Phenometric = "Onset")

get_temp_coef_offset <- function(binomial){
  df <- dplyr::filter(mdf, scientificName == binomial)
  mod <- lm(offset ~ temp, data = df)
  coef <- unname(mod$coefficients['temp'])
  return(coef)
}

get_temp_coef_se_offset <- function(binomial){
  df <- dplyr::filter(mdf, scientificName == binomial)
  mod <- lm(offset ~ temp, data = df)
  se <- coef(summary(mod))[,"Std. Error"]['temp'] 
  return(se)
}

spp_list <- c("Pachydiplax longipennis",
              "Bombus impatiens",
              "Papilio glaucus",
              "Diabrotica undecimpunctata",
              "Delphinia picta",
              "Neotibicen superbus")

model_coef_off <- lapply(X = spp_list, FUN = get_temp_coef_offset)
model_se_off <- lapply(X = spp_list, FUN = get_temp_coef_se_offset)


temp_coef_df_offset <- data.frame(
  scientificName = c("Pachydiplax longipennis",
                     "Bombus impatiens",
                     "Papilio glaucus",
                     "Diabrotica undecimpunctata",
                     "Delphinia picta",
                     "Neotibicen superbus"),
  coef = c(model_coef_off[[1]],
           model_coef_off[[2]],
           model_coef_off[[3]],
           model_coef_off[[4]],
           model_coef_off[[5]],
           model_coef_off[[6]]),
  
  high_se =c(model_coef_off[[1]] + (2 * model_se_off[[1]]),
             model_coef_off[[2]] + (2 * model_se_off[[2]]),
             model_coef_off[[3]] + (2 * model_se_off[[3]]),
             model_coef_off[[4]] + (2 * model_se_off[[4]]),
             model_coef_off[[5]] + (2 * model_se_off[[5]]),
             model_coef_off[[6]] + (2 * model_se_off[[6]])),
  
  low_se =c(model_coef_off[[1]] - (2 * model_se_off[[1]]),
            model_coef_off[[2]] - (2 * model_se_off[[2]]),
            model_coef_off[[3]] - (2 * model_se_off[[3]]),
            model_coef_off[[4]] - (2 * model_se_off[[4]]),
            model_coef_off[[5]] - (2 * model_se_off[[5]]),
            model_coef_off[[6]] - (2 * model_se_off[[6]]))
  
) %>% 
  mutate(Phenometric = "Offset")


get_temp_coef_duration <- function(binomial){
  df <- dplyr::filter(mdf, scientificName == binomial)
  bimod <- lm(duration ~ temp, data = df)
  coef <- unname(bimod$coefficients['temp'])
  return(coef)
}

get_temp_coef_se_duration <- function(binomial){
  df <- dplyr::filter(mdf, scientificName == binomial)
  mod <- lm(duration ~ temp, data = df)
  se <- coef(summary(mod))[,"Std. Error"]['temp'] 
  return(se)
}

spp_list <- c("Pachydiplax longipennis",
              "Bombus impatiens",
              "Papilio glaucus",
              "Diabrotica undecimpunctata",
              "Delphinia picta",
              "Neotibicen superbus")

model_coef_dur <- lapply(X = spp_list, FUN = get_temp_coef_duration)
model_se_dur <- lapply(X = spp_list, FUN = get_temp_coef_se_duration)

temp_coef_df_duration <- data.frame(
  scientificName = c("Pachydiplax longipennis",
                     "Bombus impatiens",
                     "Papilio glaucus",
                     "Diabrotica undecimpunctata",
                     "Delphinia picta",
                     "Neotibicen superbus"),
  coef = c(model_coef_dur[[1]],
           model_coef_dur[[2]],
           model_coef_dur[[3]],
           model_coef_dur[[4]],
           model_coef_dur[[5]],
           model_coef_dur[[6]]),
  
  high_se =c(model_coef_dur[[1]] + (2 * model_se_dur[[1]]),
             model_coef_dur[[2]] + (2 * model_se_dur[[2]]),
             model_coef_dur[[3]] + (2 * model_se_dur[[3]]),
             model_coef_dur[[4]] + (2 * model_se_dur[[4]]),
             model_coef_dur[[5]] + (2 * model_se_dur[[5]]),
             model_coef_dur[[6]] + (2 * model_se_dur[[6]])),
  
  low_se =c(model_coef_dur[[1]] - (2 * model_se_dur[[1]]),
            model_coef_dur[[2]] - (2 * model_se_dur[[2]]),
            model_coef_dur[[3]] - (2 * model_se_dur[[3]]),
            model_coef_dur[[4]] - (2 * model_se_dur[[4]]),
            model_coef_dur[[5]] - (2 * model_se_dur[[5]]),
            model_coef_dur[[6]] - (2 * model_se_dur[[6]]))
) %>% 
  mutate(Phenometric = "Duration")


tot_coef_df <- rbind(temp_coef_df_onset, temp_coef_df_offset, temp_coef_df_duration)
tot_coef_df$Phenometric <- factor(tot_coef_df$Phenometric, levels = c("Onset", "Offset", "Duration"))


p <- ggplot(tot_coef_df, mapping = aes(x = coef, y = scientificName, color = Phenometric)) + 
  geom_vline(xintercept = 0, alpha = 0.5, linetype = "dashed") + 
  geom_point() +  
  geom_errorbarh(mapping = aes(xmin = low_se, xmax = high_se, y = scientificName),
                 height = 0.5) +
  scale_color_manual(values = c("slateblue2", "tomato1", "springgreen2")) + 
  labs(x = "Temp Coef", y = "") +
  theme_minimal()
p

pg <- "https://upload.wikimedia.org/wikipedia/commons/7/7e/Eastern_Tiger_Swallowtail_Papilio_glaucus_Female_2838px.jpg"
pl <- "https://upload.wikimedia.org/wikipedia/commons/8/8d/Pachydiplax_longipennis_Blue_Dasher_1500px.jpg"
ns <- "https://upload.wikimedia.org/wikipedia/commons/7/73/Neotibicen_superbus_P1500620a.jpg"
du <- "https://upload.wikimedia.org/wikipedia/commons/b/b2/Spotted_Cucumber_Beetle_-_Diabrotica_undecimpunctata%2C_near_Buffalo%2C_Texas.jpg"
dp <- "https://upload.wikimedia.org/wikipedia/commons/6/67/Delphinia_picta_P1320004b.jpg"
bi <- "https://upload.wikimedia.org/wikipedia/commons/c/ce/Common_Eastern_Bumble_Bee_%28Bombus_impatiens%29_-_Kitchener%2C_Ontario_01.jpg"

## add in images
library(cowplot)
pimage <- axis_canvas(p, axis = 'y') + 
  draw_image(pg, y = 5.5, scale = 1) +
  draw_image(pl, y = 4.5, scale = 1) +
  draw_image(ns, y = 3.5, scale = 1) +
  draw_image(du, y = 2.5, scale = 1) +
  draw_image(dp, y = 1.5, scale = 1) +
  draw_image(bi, y = 0.5, scale = 1) 
  
ggdraw(insert_yaxis_grob(p, pimage, position = "left"))

ggsave(filename = "Tables&Figures/spp_temp_sensitivity.png", width = 10, height = 7)
