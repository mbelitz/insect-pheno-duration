
d_onset_25km = filter(all_envi_25km, pheno == "onset_date", type == "normal_base") %>% 
  st_drop_geometry() %>% na.omit()
d_onset_25km[, c("ave_tmean", "pop_25km_log")] = scale(d_onset_25km[, c("ave_tmean", "pop_25km_log")])

m_onset_25km = lmer(julian_day ~ ave_tmean * pop_25km_log + (1|id_cells) + (1|scientific_name) +
                      (0 + ave_tmean | scientific_name) + (0 + pop_25km_log | scientific_name) +
                      (0 + ave_tmean : pop_25km_log | scientific_name), 
                    data = d_onset_25km) 
summary(m_onset_25km)

m_onset_25km = lmer(julian_day ~ ave_tmean + pop_25km_log + (1|id_cells) + (1|scientific_name) +
                      (0 + ave_tmean | scientific_name) + (0 + pop_25km_log | scientific_name), 
                    data = d_onset_25km) 
summary(m_onset_25km)

m_onset_25km = lmer(julian_day ~ ave_tmean + (1|id_cells) + (1|scientific_name) +
                      (0 + ave_tmean | scientific_name), 
                    data = d_onset_25km) 
summary(m_onset_25km)


d_offset_25km = filter(all_envi_25km, pheno == "offset_date", type == "normal_base") %>% 
  st_drop_geometry() %>% na.omit()
d_offset_25km[, c("ave_tmean", "pop_25km_log")] = scale(d_offset_25km[, c("ave_tmean", "pop_25km_log")])
cor.test(d_offset_25km$ave_tmean, d_offset_25km$pop_25km_log)

m_offset_25km = lmer(julian_day ~ ave_tmean * pop_25km_log + (1|id_cells) + (1|scientific_name) +
                       (0 + ave_tmean | scientific_name) + (0 + pop_25km_log | scientific_name) +
                      (0 + ave_tmean : pop_25km_log | scientific_name), 
                    data = d_offset_25km) 
summary(m_offset_25km)

m_offset_25km = lmer(julian_day ~ ave_tmean + pop_25km_log + (1|id_cells) + (1|scientific_name) +
                      (0 + ave_tmean | scientific_name) + (0 + pop_25km_log | scientific_name), 
                    data = d_offset_25km) 
summary(m_offset_25km)


d_dur_25km = filter(all_envi_25km, pheno == "obs_date_dur", type == "normal_base") %>% 
  st_drop_geometry() %>% na.omit()
d_dur_25km[, c("ave_tmean", "pop_25km_log")] = scale(d_dur_25km[, c("ave_tmean", "pop_25km_log")])
cor.test(d_dur_25km$ave_tmean, d_dur_25km$pop_25km_log)
n_distinct(d_dur_25km$id_cells)

ggplot(d_dur_25km, aes(x = ave_tmean, y = pop_25km_log)) +
  geom_point(alpha = 0.7) + 
  geom_smooth(method = "lm") +
  facet_wrap(~scientific_name, scales = "free")
d_dur_25km %>% 
  dplyr::select(scientific_name, ave_tmean, pop_25km_log) %>% 
  na.omit() %>% 
  group_by(scientific_name) %>% 
  do(broom::tidy(cor.test(~ ave_tmean + pop_25km_log, data = .))) %>% 
  View()

m_dur_25km_inter = lmer(dur ~ ave_tmean * pop_25km_log + (1|id_cells) + (1|scientific_name) +
                    (0 + ave_tmean | scientific_name) + (0 + pop_25km_log | scientific_name) +
                       (0 + ave_tmean : pop_25km_log | scientific_name), 
                     data = d_dur_25km, REML = F) 
getME(m_dur_25km_inter, "theta")
getME(m_dur_25km_inter, "lower")
xm = lmerTest::step(m_dur_25km_inter, reduce.fixed = T)
get_model(xm)

m_dur_25km_inter1 = lmer(dur ~ ave_tmean * pop_25km_log + (1|id_cells) + (1|scientific_name) +
                          (0 + ave_tmean | scientific_name) + #(0 + pop_25km_log | scientific_name) +
                          (0 + ave_tmean : pop_25km_log | scientific_name), 
                        data = d_dur_25km, REML = F) 
m_dur_25km_inter2 = lmer(dur ~ ave_tmean * pop_25km_log + (1|id_cells) + (1|scientific_name) +
                          # (0 + ave_tmean | scientific_name) + 
                           (0 + pop_25km_log | scientific_name) +
                          (0 + ave_tmean : pop_25km_log | scientific_name), 
                        data = d_dur_25km, REML = F) 
m_dur_25km_inter3 = lmer(dur ~ ave_tmean * pop_25km_log + (1|id_cells) + (1|scientific_name) +
                           # (0 + ave_tmean | scientific_name) + 
                           # (0 + pop_25km_log | scientific_name) +
                           (0 + ave_tmean : pop_25km_log | scientific_name), 
                         data = d_dur_25km, REML = F) 
MuMIn::Weights(AIC(m_dur_25km_inter0, m_dur_25km_inter1, m_dur_25km_inter2, m_dur_25km_inter3))
summary(m_dur_25km_inter3)


m_dur_25km_temp_pop = lmer(dur ~ ave_tmean + pop_25km_log + (1|id_cells) + (1|scientific_name) +
                       (0 + ave_tmean | scientific_name) + (0 + pop_25km_log | scientific_name), 
                     data = d_dur_25km, REML = F) 
summary(m_dur_25km_temp_pop)
xm = lmerTest::step(m_dur_25km_temp_pop, reduce.fixed = F)
get_model(xm)

# 27 sp 
# Estimate Std. Error     df t value Pr(>|t|)    
# (Intercept)    67.495      7.356 22.590   9.175 4.45e-09 ***
#   ave_tmean      23.485      7.856 25.766   2.989  0.00608 ** 
#   pop_25km_log    1.960      2.240 20.895   0.875  0.39141   

m_dur_25km_temp = lmer(dur ~ ave_tmean + (1|id_cells) + (1|scientific_name) +
                             (0 + ave_tmean | scientific_name), 
                           data = d_dur_25km, REML = F, 
                       control = lmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 2e5))) 
summary(m_dur_25km_temp)

getME(m_dur_25km_temp, "theta")
getME(m_dur_25km_temp, "lower")

m_dur_25km_pop = lmer(dur ~ pop_25km_log + (1|id_cells) + (1|scientific_name) +
                            (0 + pop_25km_log | scientific_name), 
                           data = d_dur_25km, REML = F) 
summary(m_dur_25km_pop)

MuMIn::Weights(AIC(m_dur_25km_inter, m_dur_25km_temp_pop, m_dur_25km_temp, m_dur_25km_pop))
# 0.12, 0.88, 0, 0

m_dur_25km_temp2 = lmer(dur ~ ave_tmean * bloomer + (1|id_cells) + (1|scientific_name) +
                         (0 + ave_tmean | scientific_name), 
                       data = d_dur_25km, REML = F, 
                       control = lmerControl(optimizer = "bobyqa",
                                             optCtrl = list(maxfun = 2e5))) 
summary(m_dur_25km_temp2)

library(effects)
plot(allEffects(m_dur_25km_temp2), multiline = T)
