
d_onset_25km = filter(all_envi_25km, pheno == "onset_date", type == "normal_base") %>% 
  st_drop_geometry() %>% na.omit()
d_onset_25km[, c("ave_tmean", "pop_25km_log")] = scale(d_onset_25km[, c("ave_tmean", "pop_25km_log")])

m_onset_25km = lmer(julian_day ~ ave_tmean * pop_25km_log + (1|id_cells) + (1|sp) +
                      (0 + ave_tmean : pop_25km_log | sp), 
                    data = d_onset_25km) 
summary(m_onset_25km)

m_onset_25km = lmer(julian_day ~ ave_tmean + pop_25km_log + (1|id_cells) + (1|sp) +
                      (0 + ave_tmean | sp) + (0 + pop_25km_log | sp), 
                    data = d_onset_25km) 
summary(m_onset_25km)

m_onset_25km = lmer(julian_day ~ ave_tmean + (1|id_cells) + (1|sp) +
                      (0 + ave_tmean | sp), 
                    data = d_onset_25km) 
summary(m_onset_25km)


d_offset_25km = filter(all_envi_25km, pheno == "offset_date", type == "normal_base") %>% 
  st_drop_geometry() %>% na.omit()
d_offset_25km[, c("ave_tmean", "pop_25km_log")] = scale(d_offset_25km[, c("ave_tmean", "pop_25km_log")])
cor.test(d_offset_25km$ave_tmean, d_offset_25km$pop_25km_log)

m_offset_25km = lmer(julian_day ~ ave_tmean * pop_25km_log + (1|id_cells) + (1|sp) +
                       (0 + ave_tmean | sp) + (0 + pop_25km_log | sp) +
                      (0 + ave_tmean : pop_25km_log | sp), 
                    data = d_offset_25km) 
summary(m_offset_25km)

m_offset_25km = lmer(julian_day ~ ave_tmean + pop_25km_log + (1|id_cells) + (1|sp) +
                      (0 + ave_tmean | sp) + (0 + pop_25km_log | sp), 
                    data = d_offset_25km) 
summary(m_offset_25km)


d_dur = d_dur_25km %>% 
  mutate(sp = as.character(sp)) %>% 
  # filter(sp != "Fragaria vesca") %>% 
  # st_drop_geometry() %>%
  na.omit() %>% as.data.frame()
d_dur[, c("pop_25km_log", "bio_1", "bio_4", "bio_12", "bio_15")] = 
  scale(d_dur[, c("pop_25km_log", "bio_1", "bio_4", "bio_12", "bio_15")])
n_distinct(d_dur$id_cells)
names(d_dur) = gsub("_", "", names(d_dur))
d_dur = mutate(d_dur, bio1pop25kmlog = bio1 * pop25kmlog, bio4pop25kmlog = bio4 * pop25kmlog)
cor(unique(d_dur[, c("pop25kmlog", "bio1", "bio4", "bio12", "bio15")])) %>% 
corrplot::corrplot(method = "color", type = "upper", 
                   addCoef.col = "gray50", diag = FALSE)
m_dur_lat = lmer(
  dur ~ bio1 + bio4 + bio12 + bio15 + pop25kmlog + (1 | idcells) + (1 | sp) + 
    (0 + bio1 | sp) + (0 + bio4 | sp) + (0 + bio12 | sp) + (0 + bio15 | sp) + (0 + pop25kmlog | sp),
  data = d_dur,
  REML = T,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))) 
summary(m_dur_lat)

anova(
  lmer(dur ~ lat + bio1 + bio4 + (1 | idcells) + (1 | sp) + 
         (0 + lat | sp) + (0 + bio1 | sp) + (0 + bio4 | sp), data = d_dur, REML = F),
  lmer(dur ~ lat + bio4 + (1 | idcells) + (1 | sp) + 
         (0 + lat | sp) + (0 + bio4 | sp), data = d_dur, REML = F)
)
summary(lmer(dur ~ lat + (1 | idcells) + (1 | sp) + 
               (0 + lat | sp), data = d_dur, REML = F))
summary(lmer(dur ~ bio1 + (1 | idcells) + (1 | sp) + 
               (0 + bio1 | sp), data = d_dur, REML = F))

anova(
  lmer(dur ~ bio1 + bio4 + (1 | idcells) + (1 | sp) + 
         (0 + bio1 | sp) + (0 + bio4 | sp), data = d_dur, REML = F),
  lmer(dur ~ bio4 + (1 | idcells) + (1 | sp) + 
         (0 + bio4 | sp), data = d_dur, REML = F)
)
# bio4 is stronger
anova(
  lmer(dur ~ bio1 + (1 | idcells) + (1 | sp) + 
         (0 + bio1 | sp), data = d_dur, REML = F),
  lmer(dur ~ bio4 + (1 | idcells) + (1 | sp) + 
         (0 + bio4 | sp), data = d_dur, REML = F)
)

m_dur = lmer(
  dur ~ bio1 + pop25kmlog + bio1pop25kmlog + (1 | idcells) + (1 | sp) + 
    (0 + bio1 | sp) + (0 + pop25kmlog | sp), # +
   # (0 + bio1pop25kmlog | sp),
  data = d_dur,
  REML = F,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))) 
summary(m_dur)
plot(effects::effect("pop25kmlog", m_dur))
step(m_dur)

m_dur = lmer(
  dur ~ bio1 + bio4 + pop25kmlog + bio1pop25kmlog + bio4pop25kmlog + (1 | idcells) + (1 | sp) + 
    (0 + bio1 | sp) + (0 + bio4 | sp) + (0 + pop25kmlog | sp) +
   (0 + bio1pop25kmlog | sp) + (0 + bio4pop25kmlog | sp),
  data = d_dur,
  REML = F,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))) 

m_dur = lme4::lmer(
  dur ~ 
    bio1 +
    bio4 +
    pop25kmlog + 
    bio1:pop25kmlog +
    bio4:pop25kmlog +
    (1 | idcells) + (1 | sp) + 
    # (0 + bio1*pop25kmlog | sp) + (0 + bio4*pop25kmlog | sp),
    (0 + bio1 | sp) +
    (0 + bio4 | sp) +
    (0 + pop25kmlog | sp) +
    (0 + bio1:pop25kmlog | sp) +
    (0 + bio4:pop25kmlog | sp),
  data = d_dur,
  REML = F,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))) 
summary(m_dur)
plot(effects::effect("bio1:pop25kmlog", m_dur), multiline = T, x.var = "pop25kmlog")
plot(effects::effect("bio4:pop25kmlog", m_dur), multiline = T, x.var = "pop25kmlog")
MuMIn::r.squaredGLMM(m_dur)
car::vif(m_dur)
class(m_dur) = "lmerMod"
performance::check_collinearity(m_dur)
performance::model_performance(m_dur)
performance::check_normality(m_dur)

m_dur2 = lmer(
  dur ~ 
    bio1 +
    bio4 +
    pop25kmlog + 
    bio1:pop25kmlog +
    bio4:pop25kmlog +
    (1 | idcells) + (1 | sp) +
    # (0 + bio1*pop25kmlog | sp) + (0 + bio4*pop25kmlog | sp),
    (0 + bio1 | sp) +
    (0 + bio4 | sp) +
    (0 + bio1:pop25kmlog | sp), # +
    # (0 + bio4:pop25kmlog | sp) +
    # (0 + pop25kmlog | sp),
  data = d_dur,
  REML = F,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))) 

summary(m_dur2)
anova(m_dur, m_dur2)

anova(m_dur, lmer(formula = dur ~ bio1 + bio4 + pop25kmlog + bio1pop25kmlog + 
               bio4pop25kmlog + (1 | idcells) + (1 | sp) + (0 + bio1 | sp) + 
               (0 + bio4 | sp) + (0 + pop25kmlog | sp), data = d_dur, REML = F, 
             control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05))))
m_dur2 = lmer(
  dur ~ bio1 + bio4 + pop25kmlog + bio1:pop25kmlog + (1 | idcells) + (1 | sp) + 
    (0 + bio1 | sp) + (0 + bio4 | sp) + (0 + pop25kmlog | sp) +
    (0 + bio1:pop25kmlog | sp),
  data = d_dur,
  REML = F,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))) 
summary(m_dur2)
anova(m_dur, m_dur2)
m_dur_f = lmerTest::step(m_dur)

# no 3-way interactions
m_dur = lmer(
  dur ~ bio1 * bio4 * pop25kmlog + (1 | id_cells) + (1 | sp) + 
    (0 + bio1 | sp) + (0 + bio4 | sp) + (0 + pop25kmlog | sp) +
    (0 + bio1:pop25kmlog | sp) + (0 + bio4:pop25kmlog | sp) +
    (0 + bio1:bio4:pop25kmlog | sp),
  data = d_dur,
  REML = T,
  control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))) 
summary(m_dur)
plot(effects::effect("bio_1:pop_25km_log", m_dur), multiline = T)

m_dur_25km_inter = lmer(dur ~ ave_tmean * pop_25km_log + (1|id_cells) + (1|sp) +
                    (0 + ave_tmean | sp) + (0 + pop_25km_log | sp) +
                       (0 + ave_tmean : pop_25km_log | sp), 
                     data = d_dur, REML = F) 
m_dur_25km_inter1 = lmer(dur ~ ave_tmean * pop_25km_log + (1|id_cells) + (1|sp) +
                          (0 + ave_tmean | sp) + #(0 + pop_25km_log | sp) +
                          (0 + ave_tmean : pop_25km_log | sp), 
                        data = d_dur, REML = F) 
m_dur_25km_inter2 = lmer(dur ~ ave_tmean * pop_25km_log + (1|id_cells) + (1|sp) +
                          # (0 + ave_tmean | sp) + 
                           (0 + pop_25km_log | sp) +
                          (0 + ave_tmean : pop_25km_log | sp), 
                        data = d_dur, REML = F) 
m_dur_25km_inter3 = lmer(dur ~ ave_tmean * pop_25km_log + (1|id_cells) + (1|sp) +
                           # (0 + ave_tmean | sp) + 
                           # (0 + pop_25km_log | sp) +
                           (0 + ave_tmean : pop_25km_log | sp), 
                         data = d_dur, REML = F) 
MuMIn::Weights(AIC(m_dur_25km_inter0, m_dur_25km_inter1, m_dur_25km_inter2, m_dur_25km_inter3))
summary(m_dur_25km_inter3)


m_dur_25km_temp_pop = lmer(dur ~ ave_tmean + pop_25km_log + (1|id_cells) + (1|sp) +
                       (0 + ave_tmean | sp) + (0 + pop_25km_log | sp), 
                     data = d_dur, REML = F) 
summary(m_dur_25km_temp_pop)

# 27 sp 
# Estimate Std. Error     df t value Pr(>|t|)    
# (Intercept)    67.495      7.356 22.590   9.175 4.45e-09 ***
#   ave_tmean      23.485      7.856 25.766   2.989  0.00608 ** 
#   pop_25km_log    1.960      2.240 20.895   0.875  0.39141   

m_dur_25km_temp = lmer(dur ~ ave_tmean + (1|id_cells) + (1|sp) +
                             (0 + ave_tmean | sp), 
                           data = d_dur_25km, REML = F) 
summary(m_dur_25km_temp)

m_dur_25km_pop = lmer(dur ~ pop_25km_log + (1|id_cells) + (1|sp) +
                            (0 + pop_25km_log | sp), 
                           data = d_dur_25km, REML = F) 
summary(m_dur_25km_pop)

MuMIn::Weights(AIC(m_dur_25km_inter, m_dur_25km_temp_pop, m_dur_25km_temp, m_dur_25km_pop))
# 0.12, 0.88, 0, 0