
# explore data ----
all_envi_25km = readRDS("data/all_envi_25km.rds")
all_envi_25km = rbind(readRDS("data/all_envi_25km.rds"),
                      readRDS("data/all_envi_25km_12sp.rds"))
saveRDS(all_envi_25km, "data_output/all_envi_25km_39sp.rds")
n_distinct(all_envi_25km$scientific_name)
table(all_envi_25km$scientific_name)

# onset
d_onset_25km = filter(all_envi_25km, pheno == "onset_date")
d_onset_order = d_onset_25km %>% 
  st_drop_geometry() %>% 
  group_by(scientific_name) %>% 
  summarise(earliest_onset = min(julian_day, na.rm = T)) %>% 
  mutate(bloomer = ifelse(earliest_onset < 60, "early",
                          ifelse(earliest_onset < 100, "mid", "late")))

all_envi_25km = all_envi_25km %>% 
  left_join(d_onset_order) %>% 
  mutate(scientific_name = forcats::fct_reorder(scientific_name, earliest_onset))

d_onset_25km = filter(all_envi_25km, pheno == "onset_date")

ggplot(d_onset_25km, aes(x = ave_tmean, y = julian_day)) +
  geom_point(aes(color = yr)) + geom_smooth(method = "lm") +
  facet_wrap(~scientific_name, scales = "free")
ggsave("figures/onset_temp_25km.pdf", width = 16, height = 16)

p_onset_normal = ggplot(filter(d_onset_25km, type == "normal_base", julian_day < 250), 
       aes(x = ave_tmean, y = julian_day)) +
  geom_point(aes(color = yr)) + geom_smooth(method = "lm") +
  facet_wrap(~scientific_name, scales = "free")
ggsave("figures/onset_temp_25km_normal.pdf", plot = p_onset_normal, width = 16, height = 16)

# # offset
# d_offset_25km = filter(all_envi_25km, pheno == "offset_date")
# 
# ggplot(d_offset_25km, aes(x = ave_tmean, y = julian_day)) +
#   geom_point(aes(color = yr)) + geom_smooth(method = "lm") +
#   facet_wrap(~scientific_name, scales = "free")
# ggsave("figures/offset_temp_25km.pdf", width = 16, height = 16)
# 
# ggplot(filter(d_offset_25km, type == "normal_base"), 
#        aes(x = ave_tmean, y = julian_day)) +
#   geom_point(aes(color = yr)) + geom_smooth(method = "lm") +
#   facet_wrap(~scientific_name, scales = "free")
# ggsave("figures/offset_temp_25km_normal.pdf", width = 16, height = 16)

# onset & offset
d_onset_offset_25km = filter(all_envi_25km, pheno != "obs_date_dur")
d_onset_offset_25km$type[d_onset_offset_25km$scientific_name == "Houstonia pusilla" &
                           d_onset_offset_25km$julian_day > 300] = "shifted_base"
d_onset_offset_25km$type[d_onset_offset_25km$scientific_name == "Salvia lyrata" &
                           d_onset_offset_25km$julian_day > 300] = "shifted_base"
d_onset_offset_25km$type[d_onset_offset_25km$scientific_name == "Anemone berlandieri" &
                           d_onset_offset_25km$julian_day > 220] = "shifted_base"
d_onset_offset_25km$type[d_onset_offset_25km$scientific_name == "Castilleja exserta" &
                           d_onset_offset_25km$julian_day > 300] = "shifted_base"
d_onset_offset_25km$type[d_onset_offset_25km$scientific_name == "Sanguinaria canadensis" &
                           d_onset_offset_25km$julian_day > 190] = "shifted_base"
d_onset_offset_25km$type[d_onset_offset_25km$scientific_name == "Cephalanthus occidentalis" &
                           d_onset_offset_25km$julian_day > 340] = "shifted_base"

d_onset_offset_25km$julian_day[d_onset_offset_25km$scientific_name == "Acmispon glaber" &
                                 d_onset_offset_25km$pheno == "onset_date" &
                           d_onset_offset_25km$julian_day > 340] = 
  d_onset_offset_25km$julian_day[d_onset_offset_25km$scientific_name == "Acmispon glaber" &
                                   d_onset_offset_25km$pheno == "onset_date" &
                                   d_onset_offset_25km$julian_day > 340] - 365
d_onset_offset_25km$type[d_onset_offset_25km$scientific_name == "Acmispon glaber" &
                           d_onset_offset_25km$pheno == "offset_date" &
                           d_onset_offset_25km$julian_day < 90] = "shifted_base"
d_onset_offset_25km$type[d_onset_offset_25km$scientific_name == "Sambucus cerulea" &
                           d_onset_offset_25km$pheno == "offset_date" &
                           d_onset_offset_25km$julian_day < 90] = "shifted_base"
d_onset_offset_25km$type[d_onset_offset_25km$scientific_name == "Sisyrinchium bellum" &
                           d_onset_offset_25km$pheno == "offset_date" &
                           (d_onset_offset_25km$julian_day < 90 |
                              d_onset_offset_25km$julian_day > 300)] = "shifted_base"
d_onset_offset_25km$type[d_onset_offset_25km$scientific_name == "Salvia lyrata" &
                           d_onset_offset_25km$pheno == "offset_date" &
                           d_onset_offset_25km$julian_day > 250] = "shifted_base"
d_onset_offset_25km$type[d_onset_offset_25km$scientific_name == "Cephalanthus occidentalis" &
                           d_onset_offset_25km$pheno == "offset_date" &
                           d_onset_offset_25km$julian_day < 180] = "shifted_base"
d_onset_offset_25km$type[d_onset_offset_25km$scientific_name == "Cypripedium acaule" &
                           d_onset_offset_25km$pheno == "onset_date" &
                           d_onset_offset_25km$julian_day < 50] = "shifted_base"


ggplot(d_onset_offset_25km, aes(x = ave_tmean, y = julian_day, color = pheno)) +
  geom_point(aes(shape = yr)) + geom_smooth(method = "lm") +
  facet_wrap(~scientific_name, scales = "free")
ggsave("figures/onset_offset_temp_25km.pdf", width = 16, height = 16)

p_on_off_normal = ggplot(filter(d_onset_offset_25km, type == "normal_base"), 
       aes(x = ave_tmean, y = julian_day, color = pheno)) +
  geom_point(aes(shape = yr)) + geom_smooth(method = "lm") +
  facet_wrap(~scientific_name, scales = "free")
ggsave("figures/onset_offset_temp_25km_normal.pdf", plot = p_on_off_normal, width = 16, height = 16)

ggplot(filter(d_onset_offset_25km, type == "normal_base"), 
       aes(x = ave_tmean, y = julian_day, color = pheno)) +
  geom_point(aes(shape = yr)) + geom_smooth(method = "lm") +
  facet_wrap(~scientific_name) +
  coord_cartesian(ylim = c(0, 365))
ggsave("figures/onset_offset_temp_25km_normal_fixed_scale.pdf", width = 16, height = 16)

ggplot(d_onset_offset_25km, aes(x = pop_25km_log, y = julian_day, color = pheno)) +
  geom_point(aes(shape = yr)) + geom_smooth(method = "lm") +
  facet_wrap(~scientific_name, scales = "free")
ggsave("figures/onset_offset_pop_25km.pdf", width = 16, height = 16)

p_on_off_normal_pop = ggplot(filter(d_onset_offset_25km, type == "normal_base"), 
       aes(x = pop_25km_log, y = julian_day, color = pheno)) +
  geom_point(aes(shape = yr)) + geom_smooth(method = "lm") +
  facet_wrap(~scientific_name, scales = "free")
ggsave("figures/onset_offset_pop_25km_normal.pdf", plot = p_on_off_normal_pop, width = 16, height = 16)


# duration
all_envi_25km$type[all_envi_25km$scientific_name == "Houstonia pusilla" &
                           all_envi_25km$julian_day > 300] = "shifted_base"
all_envi_25km$type[all_envi_25km$scientific_name == "Salvia lyrata" &
                           all_envi_25km$julian_day > 300] = "shifted_base"
all_envi_25km$type[all_envi_25km$scientific_name == "Anemone berlandieri" &
                           all_envi_25km$julian_day > 220] = "shifted_base"
all_envi_25km$type[all_envi_25km$scientific_name == "Castilleja exserta" &
                           all_envi_25km$julian_day > 300] = "shifted_base"
all_envi_25km$type[all_envi_25km$scientific_name == "Sanguinaria canadensis" &
                           all_envi_25km$julian_day > 190] = "shifted_base"
all_envi_25km$type[all_envi_25km$scientific_name == "Cephalanthus occidentalis" &
                           all_envi_25km$julian_day > 340] = "shifted_base"

all_envi_25km$julian_day[all_envi_25km$scientific_name == "Acmispon glaber" &
                                 all_envi_25km$pheno == "onset_date" &
                                 all_envi_25km$julian_day > 340] = 
  all_envi_25km$julian_day[all_envi_25km$scientific_name == "Acmispon glaber" &
                                   all_envi_25km$pheno == "onset_date" &
                                   all_envi_25km$julian_day > 340] - 365
all_envi_25km$type[all_envi_25km$scientific_name == "Acmispon glaber" &
                           all_envi_25km$pheno == "offset_date" &
                           all_envi_25km$julian_day < 90] = "shifted_base"
all_envi_25km$type[all_envi_25km$scientific_name == "Sambucus cerulea" &
                           all_envi_25km$pheno == "offset_date" &
                           all_envi_25km$julian_day < 90] = "shifted_base"
all_envi_25km$type[all_envi_25km$scientific_name == "Sisyrinchium bellum" &
                           all_envi_25km$pheno == "offset_date" &
                           (all_envi_25km$julian_day < 90 |
                              all_envi_25km$julian_day > 300)] = "shifted_base"
all_envi_25km$type[all_envi_25km$scientific_name == "Salvia lyrata" &
                           all_envi_25km$pheno == "offset_date" &
                           all_envi_25km$julian_day > 250] = "shifted_base"
all_envi_25km$type[all_envi_25km$scientific_name == "Cephalanthus occidentalis" &
                           all_envi_25km$pheno == "offset_date" &
                           all_envi_25km$julian_day < 180] = "shifted_base"
all_envi_25km$type[all_envi_25km$scientific_name == "Cypripedium acaule" &
                           all_envi_25km$pheno == "onset_date" &
                           all_envi_25km$julian_day < 50] = "shifted_base"

d_dur_25km = filter(all_envi_25km, pheno == "obs_date_dur")

ggplot(d_dur_25km, aes(x = ave_tmean, y = dur)) +
  geom_point(aes(color = yr)) + geom_smooth(method = "lm") +
  facet_wrap(~scientific_name, scales = "free")
ggsave("figures/dur_temp_25km.pdf", width = 16, height = 16)

p_dur_normal_temp = ggplot(filter(d_dur_25km, type == "normal_base"), 
       aes(x = ave_tmean, y = dur)) +
  geom_point(aes(color = yr)) + geom_smooth(method = "lm") +
  facet_wrap(~scientific_name, scales = "free") +
  labs(x = "Average temperature", 
       y = "Flower durations",
       title = "Species are sorted by earliest onset date observed")
ggsave("figures/dur_temp_25km_normal.pdf", plot = p_dur_normal_temp, width = 16, height = 16)

p_dur_normal_pop = ggplot(filter(d_dur_25km, type == "normal_base"), 
                           aes(x = pop_25km_log, y = dur)) +
  geom_point(aes(color = yr)) + geom_smooth(method = "lm") +
  facet_wrap(~scientific_name, scales = "free") +
  labs(x = "Log human population density", 
       y = "Flower durations",
       title = "Species are sorted by earliest onset date observed")
ggsave("figures/dur_pop_25km_normal.pdf", plot = p_dur_normal_pop, width = 16, height = 16)
