d_claytonia_virginica = readr::read_csv("data/iNat_scores/Claytonia virginica_inat.csv",
                                        col_types = cols(
                                          id = col_character(),
                                          flowers = col_character(),
                                          scientific_name = col_character(),
                                          observed_on = col_date(format = ""),
                                          latitude = col_double(),
                                          longitude = col_double()
                                        )) %>% 
  dplyr::select(longitude, latitude, everything()) %>% 
  filter(flowers == "1") %>% 
  drop_na(longitude, latitude) %>% 
  rename(id_iNat = id)
cell_25k_claytonia_virginica = 
  plt_summary(cell_size = 25000, dat = d_claytonia_virginica, n_per_cell = 10)
grids_25km = cell_25k_claytonia_virginica$grids
saveRDS(grids_25km, "data_output/grids_25km.rds")
grids_25km_raster = raster(x = grids_25km, resolution = 25000)

all_25km_phenesse_df = bind_rows(
  readRDS("phenesse_outputs/all_25km_phenesse_df_12sp_2.rds"),
  readRDS("phenesse_outputs/all_25km_phenesse_df_12sp.rds"),
  readRDS("phenesse_outputs/all_25km_phenesse_df.rds")
)

# a mistake about base date
all_25km_phenesse_df$onset_date[all_25km_phenesse_df$scientific_name == "Epilobium canum" &
                                  all_25km_phenesse_df$type == "shifted_base"] = 
all_25km_phenesse_df$onset_date[all_25km_phenesse_df$scientific_name == "Epilobium canum" &
                                  all_25km_phenesse_df$type == "shifted_base"] + 365
all_25km_phenesse_df$offset_date[all_25km_phenesse_df$scientific_name == "Epilobium canum" &
                                  all_25km_phenesse_df$type == "shifted_base"] = 
  all_25km_phenesse_df$offset_date[all_25km_phenesse_df$scientific_name == "Epilobium canum" &
                                    all_25km_phenesse_df$type == "shifted_base"] + 365
all_25km_phenesse_df$base_date[all_25km_phenesse_df$scientific_name == "Epilobium canum" &
                                   all_25km_phenesse_df$type == "shifted_base"] = 
  all_25km_phenesse_df$base_date[all_25km_phenesse_df$scientific_name == "Epilobium canum" &
                                     all_25km_phenesse_df$type == "shifted_base"] + 365

est_25k <- all_25km_phenesse_df %>%
  dplyr::select(scientific_name, id_cells, observed_year, onset, offset, 
                dur = duration, onset_date, offset_date, type) %>%
  left_join(grids_25km, by = "id_cells") %>%
  unique() %>% 
  st_sf()

est_25k = mutate(est_25k, dur = ifelse(dur >= 365, 365, dur),
                 # offset = dur - onset,
                 offset_date = onset_date + dur) 

cent_25km = st_centroid(est_25k)
arrange(cent_25km, onset_date)


# extract climatic data ----

dys = cent_25km %>% 
  mutate(onset_date = as.character(onset_date),
         id_cells = as.character(id_cells),
         offset_date = as.character(offset_date)) %>% 
  dplyr::select(id_cells, onset_date, offset_date) %>% 
  unique() %>% 
  gather(key = "pheno", value = "dates", onset_date, offset_date) %>% 
  unique()
# some species have the same estimated onset/offset, thus the unique df is less

n_distinct(dys$dates)

# get all dates 90 days before onset/offset dates
dys2 = group_by(st_drop_geometry(dys), dates) %>% 
  do(try(get_date(unique(.$dates), ndays_ahead = 90))) %>% 
  ungroup() %>% 
  right_join(st_drop_geometry(dys), by = "dates") %>% 
  unique() %>% arrange(days_needed) 
table(dys2$pheno)
n_distinct(dys2$days_needed)
n_distinct(as.character(dys2$days_needed))

group_by(dys2, pheno, id_cells, dates) %>% tally() %>% pull(n) %>% unique()

# fixed window Nov 15 - Feb 28 for all id_cells
fixed_winter_dates = bind_rows(
  tibble(yr = "2017", days_needed = as.character(seq(as.Date("2016-11-15"), 
                                                    as.Date("2017-02-28"), by = 1))),
  tibble(yr = "2018", days_needed = as.character(seq(as.Date("2017-11-15"), 
                                                    as.Date("2018-02-28"), by = 1))),
  tibble(yr = "2019", days_needed = as.character(seq(as.Date("2018-11-15"), 
                                                    as.Date("2019-02-28"), by = 1)))
)
fixed_prism = dplyr::select(dys2, dates, id_cells) %>% 
  unique() %>% 
  mutate(yr = str_sub(dates, 1, 4)) %>% 
  dplyr::select(-dates) %>% 
  unique() %>% 
  mutate(yr = recode(yr, "2016" = "2017")) %>% 
  left_join(fixed_winter_dates, by = "yr") %>% 
  dplyr::select(-yr) %>% 
  unique()
dys_fixed = dplyr::select(dys2, dates, id_cells, pheno) %>% 
  unique() %>% 
  left_join(fixed_prism, by = "id_cells") %>% 
  unique()

# # duration: whole flower window + 90 days ahead
# dys3 = st_drop_geometry(cent_25km) %>% 
#   mutate(obs_date_dur = as.integer(offset_date - onset_date) + 90) %>% 
#   dplyr::select(id_cells, observed_year, offset_date, obs_date_dur) %>% 
#   unique() %>% 
#   rename(dates = offset_date) %>% 
#   mutate(dates = as.character(dates)) %>% 
#   group_by(id_cells, dates, obs_date_dur) %>% 
#   do(try(get_date(unique(.$dates), ndays_ahead = .$obs_date_dur))) %>% 
#   mutate(pheno = "obs_date_dur") %>% ungroup() 
# 
# n_distinct(dys3$days_needed)
# n_distinct(as.character(dys3$days_needed))
# filter(dys3, id_cells == 63)

# dys4 = bind_rows(mutate(dys2, id_cells = as.numeric(id_cells)), dys3) %>% 
#   left_join(unique(dplyr::select(cent_25km, id_cells))) %>% 
#   st_sf()
dys4 = bind_rows(mutate(dys2, days_needed = as.character(days_needed)), 
          dys_fixed) %>% 
  distinct() %>% 
  mutate(id_cells = as.numeric(id_cells)) %>% 
  left_join(unique(dplyr::select(cent_25km, id_cells))) %>% 
  st_sf()

table(dys4$pheno)
n_distinct(dys4$days_needed)
sort(unique(as.character(dys4$days_needed)))

dys4_list = split(dys4, dys4$days_needed)

# get_clim_per_date_flower_dur(loc_df = dys4_list[[1294]], grid_raster = grids_25km_raster)

if(!file.exists("data_output/all_prism_25km_df_all.rds")){
  # library(furrr)
  # future::plan(multisession, workers = 20)
  # 
  # all_prism_25km = future_map(dys4_list,
  #                             get_clim_per_date_flower_dur,
  #                             grid_raster = grids_25km_raster,
  #                             .progress = T)
  # 
  all_prism_25km = pbmclapply(dys4_list, get_clim_per_date_flower_dur, 
                              grid_raster = grids_25km_raster,
                              mc.cores = 20)
  map(all_prism_25km, class) %>% map_lgl(~"try-error" %in% .x) %>% sum()
  saveRDS(all_prism_25km, file = "data_output/all_prism_25km_all.rds")
  all_prism_25km_df = plyr::ldply(all_prism_25km) %>% as_tibble() %>%
    st_sf(crs = st_crs(grids_25km_raster))
  saveRDS(all_prism_25km_df, file = "data_output/all_prism_25km_df_all.rds")
  # 
  # all_prism_25km_df2 = na.omit(all_prism_25km_df)
  # all_prism_25km_df3 = anti_join(all_prism_25km_df, st_drop_geometry(all_prism_25km_df2))
  # 
  # all_prism_25km_df3_list2 = plyr::llply(all_prism_25km_df3_list, get_clim_per_date_flower_dur,
  #             grid_raster = grids_25km_raster,
  #             .progress = "text")
  
  
} else {
  all_prism_25km_df = readRDS("data_output/all_prism_25km_df_all.rds")
}

filter(all_prism_25km_df, is.na(ppt))
dim(na.omit(all_prism_25km_df))

# Nov 15 - Feb. 28
t_base = 5
prims_winter = st_drop_geometry(all_prism_25km_df) %>% 
  filter(days_needed %in% fixed_winter_dates$days_needed) %>% 
  group_by(id_cells, pheno, dates) %>% 
  arrange(days_needed) %>% # some dates say 2018-12-14 has more days
  slice(1:106) %>% 
  mutate(gdd = (tmax + tmin)/2 - t_base,
         gdd = ifelse(gdd < 0, 0, gdd)) %>% 
  summarise(n_chilling_days_winter = sum(tmean > 0 & tmean < 7.2, na.rm = T),
            tmean_winter = mean(tmean, na.rm = T),
            ppt_winter = mean(ppt, na.rm = T),
            GDD_winter = sum(gdd, na.rm = T)) %>% 
  ungroup()

all_prism_25km_df_90days = left_join(mutate(dys2, days_needed = as.character(days_needed),
                                            id_cells = as.numeric(id_cells)), 
                                     all_prism_25km_df)
# chilling days
# Chilling days are all days with average temperatures between 0 and 7.2 degrees C.

n_chill = all_prism_25km_df_90days %>% 
  group_by(id_cells, pheno, dates) %>% 
  summarise(n_chilling_days_90days = sum(tmean > 0 & tmean < 7.2, na.rm = T),
            nd = n()) %>% 
  ungroup()
unique(n_chill$nd)
if(unique(n_chill$nd) == 91){
  n_chill$nd = NULL
} else {
  error("Number of days not the same for all dates")
}

# GDD
t_base = 5
all_prism_25km_df_90days = group_by(all_prism_25km_df_90days, id_cells, pheno, dates) %>% 
  mutate(gdd = (tmax + tmin)/2 - t_base,
         gdd = ifelse(gdd < 0, 0, gdd)) %>% 
  ungroup()

all_prism_25km_df_90days_summary =
  group_by(all_prism_25km_df_90days, id_cells, pheno, dates) %>% 
  summarise(ave_ppt_90days = mean(ppt, na.rm = T),
            ave_tdmean_90days = mean(tdmean, na.rm = T),
            ave_tmax_90days = mean(tmax, na.rm = T),
            ave_tmean_90days = mean(tmean, na.rm = T),
            ave_tmin_90days = mean(tmin, na.rm = T),
            ave_vpdmax_90days = mean(vpdmax, na.rm = T),
            ave_vpdmin_90days = mean(vpdmin, na.rm = T),
            GDD_90days = sum(gdd, na.rm = T)
  ) %>% 
  ungroup()

# slopes of GDD and tmean
slope_tmean = group_by(all_prism_25km_df_90days, id_cells, pheno, dates) %>% 
  arrange(days_needed) %>% 
  mutate(nd = 1:n()) %>% 
  drop_na(tmean) %>% 
  do(broom::tidy(lm(tmean ~ nd, data = .))) %>% 
  ungroup() %>% 
  filter(term == "nd") %>% 
  dplyr::select(id_cells, pheno, dates, slope_tmean_90days = estimate)
slope_GDD = group_by(all_prism_25km_df_90days, id_cells, pheno, dates) %>% 
  arrange(days_needed) %>% 
  drop_na(gdd) %>% 
  mutate(nd = 1:n(), GDD = cumsum(gdd)) %>% 
  do(broom::tidy(lm(GDD ~ nd, data = .))) %>% 
  ungroup() %>% 
  filter(term == "nd") %>% 
  dplyr::select(id_cells, pheno, dates, slope_GDD_90days = estimate)
summary(slope_GDD$slope_GDD_90days)

all_prism_25km_df_90days_summary = left_join(all_prism_25km_df_90days_summary, 
                                             n_chill) %>% 
  left_join(slope_tmean) %>% 
  left_join(slope_GDD)
saveRDS(all_prism_25km_df_90days_summary, file = "data_output/all_prism_25km_df_90days_summary_all.rds")


ggplot() +
  geom_sf(data = grids_25km %>%
            mutate(to_color = id_cells %in% unique(filter(all_prism_25km_df_90days_summary, 
                                                          !is.finite(ave_tmean_90days))$id_cells)),
          aes(fill = to_color), size = 0.1, color = "gray70") +
  geom_sf(data = readRDS("data/usa_map.rds"), alpha = 0) +
  scale_fill_manual(values = c("gray70", "red")) 
# Mostly outside US

# annual prism data
cent_25km_cells = cent_25km %>% 
  dplyr::select(id_cells, observed_year) %>% 
  unique()

cent_25km_cells_2017 = filter(cent_25km_cells, observed_year == 2017)
cent_25km_cells_2018 = filter(cent_25km_cells, observed_year == 2018)
cent_25km_cells_2019 = filter(cent_25km_cells, observed_year == 2019)

# annual ppt and tmean
fnames = c("data/prism/PRISM_ppt_stable_4kmM3_2017_bil//PRISM_ppt_stable_4kmM3_2017_bil.bil" ,
  "data/prism/PRISM_tmean_stable_4kmM3_2017_bil//PRISM_tmean_stable_4kmM3_2017_bil.bil",
  "data/prism/PRISM_ppt_stable_4kmM3_2018_bil//PRISM_ppt_stable_4kmM3_2018_bil.bil" ,
  "data/prism/PRISM_tmean_stable_4kmM3_2018_bil//PRISM_tmean_stable_4kmM3_2018_bil.bil",
  "data/prism/PRISM_ppt_provisional_4kmM3_2019_bil/PRISM_ppt_provisional_4kmM3_2019_bil.bil",
  "data/prism/PRISM_tmean_provisional_4kmM3_2019_bil/PRISM_tmean_provisional_4kmM3_2019_bil.bil")

for(i in fnames){
  cat(i, "\n")
  xi = projectRaster(raster::raster(i), grids_25km_raster, method="bilinear")
  xname = names(xi)
  yyr = as.numeric(gsub("PRISM_.*_([0-9]{4}).*$", "\\1", xname))
  var = gsub("^PRISM_([a-z]*)_.*$", "\\1", xname)
  if(var == "ppt"){
    if(yyr == 2017){
      cent_25km_cells_2017 = mutate(cent_25km_cells_2017, 
                                    annual_ppt = raster::extract(xi, cent_25km_cells_2017))
    }
    if(yyr == 2018){
      cent_25km_cells_2018 = mutate(cent_25km_cells_2018, 
                                    annual_ppt = raster::extract(xi, cent_25km_cells_2018))
    }
    if(yyr == 2019){
      cent_25km_cells_2019 = mutate(cent_25km_cells_2019, 
                                    annual_ppt = raster::extract(xi, cent_25km_cells_2019))
    }
  }
  if(var == "tmean"){
    if(yyr == 2017){
      cent_25km_cells_2017 = mutate(cent_25km_cells_2017, 
                                    annual_tmean = raster::extract(xi, cent_25km_cells_2017))
    }
    if(yyr == 2018){
      cent_25km_cells_2018 = mutate(cent_25km_cells_2018, 
                                    annual_tmean = raster::extract(xi, cent_25km_cells_2018))
    }
    if(yyr == 2019){
      cent_25km_cells_2019 = mutate(cent_25km_cells_2019, 
                                    annual_tmean = raster::extract(xi, cent_25km_cells_2019))
    }
  }
}

# monthly ppt and tmean to calculate bio1, bio4, bio12, and bio15
f2017 = paste0("2017", c(paste0("0", 1:9), "10", "11", "12"))
f2018 = paste0("2018", c(paste0("0", 1:9), "10", "11", "12"))
# use 10-12 of 2018 for 2019
f2019 = c(paste0("2019", c(paste0("0", 1:9))), paste0("2018", c("10", "11", "12")))

fp_2017 = c(paste0("data/prism/PRISM_ppt_stable_4kmM3_", f2017, "_bil.bil"),
            paste0("data/prism/PRISM_tmean_stable_4kmM3_", f2017, "_bil.bil"))
for(i in fp_2017){
  cat(i, "\n")
  xi = projectRaster(raster::raster(i), grids_25km_raster, method="bilinear")
  mth = paste(sub(pattern = ".*PRISM_([a-z]*)_.*", "\\1", i), 
               str_sub(i, -10, -9), sep = "_")
  cent_25km_cells_2017 = mutate(cent_25km_cells_2017, 
                                vv = raster::extract(xi, cent_25km_cells_2017))
  names(cent_25km_cells_2017)[names(cent_25km_cells_2017) == "vv"] = mth
}

fp_2018 = c(paste0("data/prism/PRISM_ppt_stable_4kmM3_", f2018, "_bil.bil"),
            paste0("data/prism/PRISM_tmean_stable_4kmM3_", f2018, "_bil.bil"))
for(i in fp_2018){
  cat(i, "\n")
  xi = projectRaster(raster::raster(i), grids_25km_raster, method="bilinear")
  mth = paste(sub(pattern = ".*PRISM_([a-z]*)_.*", "\\1", i), 
              str_sub(i, -10, -9), sep = "_")
  cent_25km_cells_2018 = mutate(cent_25km_cells_2018, 
                                vv = raster::extract(xi, cent_25km_cells_2018))
  names(cent_25km_cells_2018)[names(cent_25km_cells_2018) == "vv"] = mth
}

fp_2019 = c(paste0("data/prism/PRISM_ppt_stable_4kmM3_", f2019, "_bil.bil"),
            paste0("data/prism/PRISM_tmean_stable_4kmM3_", f2019, "_bil.bil"))
for(i in fp_2019){
  cat(i, "\n")
  xi = projectRaster(raster::raster(i), grids_25km_raster, method="bilinear")
  mth = paste(sub(pattern = ".*PRISM_([a-z]*)_.*", "\\1", i), 
              str_sub(i, -10, -9), sep = "_")
  cent_25km_cells_2019 = mutate(cent_25km_cells_2019, 
                                vv = raster::extract(xi, cent_25km_cells_2019))
  names(cent_25km_cells_2019)[names(cent_25km_cells_2019) == "vv"] = mth
}

cells_25km_annual_prism = rbind(cent_25km_cells_2017, cent_25km_cells_2018, cent_25km_cells_2019) %>% 
  st_drop_geometry()

cells_25km_annual_prism = left_join(dplyr::select(cells_25km_annual_prism, id_cells, observed_year, 
                                                  bio_1 = annual_tmean, bio_12 = annual_ppt),
          left_join(
            dplyr::select(cells_25km_annual_prism, id_cells, observed_year, ppt_01:tmean_12) %>% 
              pivot_longer(cols = ppt_01:tmean_12) %>% 
              separate("name", c("var", "mth"), sep = "_") %>% 
              filter(var == "ppt") %>% 
              group_by(id_cells, observed_year, var) %>% 
              summarise(annual_ave = mean(value, na.rm = T),
                        annual_sd = sd(value, na.rm = T)) %>% 
              mutate(bio_15 = 100 * annual_sd / (annual_ave + 1)) %>% 
              dplyr::select(-starts_with("annual"), -var) %>% ungroup(),
            dplyr::select(cells_25km_annual_prism, id_cells, observed_year, ppt_01:tmean_12) %>% 
              pivot_longer(cols = ppt_01:tmean_12) %>% 
              separate("name", c("var", "mth"), sep = "_") %>% 
              filter(var == "tmean") %>% 
              group_by(id_cells, observed_year, var) %>% 
              summarise(annual_ave = mean(value, na.rm = T),
                        annual_sd = sd(value, na.rm = T)) %>% 
              mutate(bio_4 = 100 * annual_sd) %>% 
              dplyr::select(-starts_with("annual"), -var) %>% ungroup()
          )
)



# pop and elev
pop_25km = readRDS("data/popu_usa_m5_2010_25km.rds")
elev_25km = readRDS("data/elev_usa_25km.rds")

# estimated durations and dates; add pop. and elevation
est_25k_centroid = mutate(cent_25km, 
                          pop_25km = raster::extract(pop_25km, cent_25km),
                          elev = raster::extract(elev_25km, cent_25km)) %>% 
  st_drop_geometry() %>% 
  mutate(pop_25km_log = log10(pop_25km + 1),
         yr = as.factor(observed_year))

est_25k_centroid = mutate(est_25k_centroid, obs_date_dur = offset_date) %>% 
  pivot_longer(cols = c("onset_date", "offset_date", "obs_date_dur"),
               names_to = "pheno", values_to = "dates")

all_prism_25km_90days = readRDS("data_output/all_prism_25km_df_90days_summary_all.rds") 

all_envi_25km = right_join(all_prism_25km_90days, 
                          mutate(est_25k_centroid, dates = as.character(dates))) %>% 
  mutate(julian_day = yday(dates)) %>% 
  left_join(cells_25km_annual_prism) %>% 
  left_join(
    filter(grids_25km, id_cells %in% est_25k_centroid$id_cells)
  ) %>% 
  st_sf()
sum(is.na(filter(all_envi_25km, pheno == "onset_date")$ave_ppt_90days))

saveRDS(all_envi_25km, "data_output/all_envi_25km_all.rds")

