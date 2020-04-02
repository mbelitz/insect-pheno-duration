library(phenesse)
library(tidyverse)
library(sf)
library(lubridate)
library(pbmcapply)
library(raster)
library(rgbif)
library(ridigbio)
library(ggplot2)

#' Plot observation data and group by grid cells
#' 
#' @param na_map North American map.
#' @param cell_size Size (meters) of the square grid cells.
#' @param dat Data frame of observed information. The first column must be longitude;
#' the second column must be latitude.
#' @param n_per_cell Minimum records per cell to be used.
#' @param days_per_cell Minimum days of data per cell to be used.
#' @param sd_cutoff The cutoff of standard deviation of the number of records of 
#' all days for each year and cell combination; if higher than this cutoff, it
#' suggests that some days have way more records than the other days; we will thin
#' records for these days.
#' @param show_fig Plot the figure at the same time?
#' @param add_lakes_map Plot lakes?
#' @return A list of maps, a data frame to summarise number of records per cell,
#' and a data frame of the records that fall within cells with enough records.
#' 
plt_summary = function(na_map = rnaturalearth::ne_countries(country = c("United States of America", "Mexico", "Canada"),
                                                            returnclass = "sf"), 
                       cell_size = 100000, 
                       dat = tibble(long = runif(300, -110, -85),
                                    lat = runif(300, 26, 45), z = 1),
                       n_per_cell = 10,
                       show_fig = FALSE){
  # make grid over na
  na_map <- st_transform(na_map, crs = "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs ")
  grids = st_make_grid(na_map, cellsize = c(cell_size, cell_size))
  # add grid cell id
  grids = mutate(st_sf(geometry = grids), id_cells = 1:n())
  # convert lat/long to have the same crs
  gazette_idig <- 
  dat <- dat %>% 
    dplyr::rename(long = decimalLongitude) %>% 
    dplyr::rename(lat = decimalLatitude)
  
  if(!(grepl("^long", names(dat)[1], ignore.case = TRUE) &
       grepl("^lat", names(dat)[2], ignore.case = TRUE))){
    stop("The first two columns of dat must be longitude and latitude, respectively.")
  }

  dat = st_transform(st_as_sf(dat, coords = 1:2, crs = 4326), 
                     crs = st_crs(na_map)$proj4string)
  # which cell each point falls in?
  dat_cells = st_join(dat, grids)
  dat_cells_count = group_by(st_drop_geometry(dat_cells), id_cells) %>% 
    tally() %>% 
    arrange(desc(n)) %>% 
    mutate(enough_data = n >= n_per_cell)
  
  # cells with data
  cells_with_data = dplyr::filter(grids, id_cells %in% dat_cells_count$id_cells) %>% 
    left_join(dat_cells_count, by = "id_cells")
  # add centroid coords
  cells_with_data = bind_cols(cells_with_data, 
                              suppressWarnings(st_centroid(cells_with_data) %>% 
                                                 st_transform(4326) %>% 
                                                 st_coordinates() %>% 
                                                 as.data.frame() %>% 
                                                 rename(long_cell = X, lat_cell = Y)))
  
  # records fall within cells with >= n_per_cell records
  dat_to_use = filter(dat_cells, id_cells %in%
                        filter(cells_with_data, enough_data)$id_cells)
  
  plt_base = ggplot() +
    geom_sf(data = na_map) +
    geom_sf(data = grids, alpha = 0, size = 0.1, color = "gray") 
  
  plt = plt_base +
    geom_sf(data = dat, size = 0.5, alpha = 0.6) + 
    geom_sf(data = filter(cells_with_data, enough_data), alpha = 0, 
            size = 0.15, color = "red") +
    labs(title = paste(nrow(filter(cells_with_data, enough_data)),
                       "highlighted cells with records more than",
                       n_per_cell, 
                       "(Cell resolution:", cell_size/1000, "km by",
                       cell_size/1000, "km)",
                       collapse = " ")) 
  
  if(show_fig){
    print(plt)
  }
  
  cat(nrow(filter(cells_with_data, enough_data)), 
      "cells with records more than", n_per_cell, "\n")
  
  
  list(cells_with_data = cells_with_data, grids = grids,
       dat_to_use = dat_to_use, fig = plt, fig_base = plt_base)
}

#' Take the output from plt_summary function and run phenesse on the cell of
#' interest.
#' 
#' @param df The named output `dat_to_use` of the plt_summary function; after 
#' some examniation for potential outliers and decide a new base date if needed
#' to deal with flowering across calendar years.
#' @param minimum_obs Minimum records per cell to be used.
#' @param minimum_days Minimum days per cell to be used.
#' @param earliest_year Earliest year to be included in analysis
#' @param latest_year Latest year to be included in analysis
#' @param flowering_cutoff Day of year to filter observations by if next year's
#' flowering maybe occurring in December
#' @param n_item Number of iterations for phenesse.
#' @param onset_perct Which percentile to use for onset?
#' @param offset_perct Which percentile to use for offset?
#' @param num_cores Number of cores to use in calculation
#' @return A dataframe of the onset, offset, and duration
#'  calculation for the cells of interest.
#'  run_phenesse(cell_100k, 100, 2018, 2018, num_cores = 4)

run_phenesse <- function(df, minimum_obs = 10, 
                         earliest_year = 2015, last_year = 2019, n_item = 500,
                         onset_perct = 0, offset_perct = 1, num_cores){
  # make Daijiang function output into dataframe
  df <- df
  
  df <- df %>% 
    mutate(year = year(as_date(eventDate))) %>% 
    mutate(doy = yday(as_date(eventDate)))
  
  # filter data to only years of interest
  df <- df %>% 
    filter(year >= earliest_year & year <= last_year) 
  
  # filter to only include one observation per cell per year per day
  df2 <- df %>% 
    group_by(doy, id_cells, year) %>% 
    slice(1)
  
  # count number of records for each cell x year combination
  num_of_records <- df2 %>% 
    group_by(year, id_cells) %>% 
    summarise(count = n())
  
  not_enough_data <- num_of_records %>% 
    filter(count >= minimum_obs) %>% 
    as.data.frame()
  
  # remove cell, year combinations that do not have enough records
  enough_records <- left_join(df2, not_enough_data, by = c("year", "id_cells")) %>% 
    filter(!is.na(count)) %>% 
    filter(doy != 1) %>% 
    dplyr::select(-geometry.x, -geometry.y)
  
  # make list with all doy values in it for each cell x year combination
  species_cell_year <- split(enough_records, 
                             f = list(enough_records$year, enough_records$id_cells),
                             drop = TRUE)
  
  
  # lapply functions
  if(onset_perct == 0 | offset_perct == 1) {
    setestimator <- function(x, niter = n_item, perct = 0){
      tibble(est = weib_percentile(observations = x$doy, 
                                   iterations = niter, percentile = perct))
    } } else {
      setestimator <- function(x, niter = n_item, perct = 0){
        tibble(est = quantile_ci(observations = x$doy, 
                                 bootstraps = niter, percentile = perct))
      }}
  
  # Estimate onseet and offset
  if(num_cores > 1){
    onset <- pbmclapply(species_cell_year, setestimator, 
                        perct = onset_perct, mc.cores = num_cores)
    offset <- pbmclapply(species_cell_year, setestimator, 
                         perct = offset_perct, mc.cores = num_cores)
  } else{
    onset <- lapply(species_cell_year, setestimator, perct = onset_perct)
    offset <- lapply(species_cell_year, setestimator, perct = offset_perct)
  }
  
  # split outputs back to df
  onset_df <- as.data.frame(matrix(unlist(onset), nrow = length(onset), byrow = TRUE)) %>% 
    mutate(rowname = names(onset))
  onset_df <- onset_df %>% 
    mutate(year = substr(x = rowname, start = 0, stop = 4 )) %>% 
    mutate(id_cells = substr(x = rowname, start = 6, stop = nchar(rowname))) %>% 
    rename(onset = V1, onset_low = V2, onset_high = V3)
  onset_df$year <- as.numeric(onset_df$year)
  onset_df$id_cells <- as.numeric(onset_df$id_cells)
  offset_df <- as.data.frame(matrix(unlist(offset), nrow = length(offset), byrow = TRUE))%>% 
    mutate(rowname = names(onset))
  offset_df <- offset_df %>% 
    mutate(year = substr(x = rowname, start = 0, stop = 4 )) %>% 
    mutate(id_cells = substr(x = rowname, start = 6, stop = nchar(rowname))) %>% 
    rename(offset = V1, offset_low = V2, offset_high = V3)
  offset_df$year <- as.numeric(offset_df$year)
  offset_df$id_cells <- as.numeric(offset_df$id_cells)  
  
  # join estimates with original sf dataframe based on cell_ids and year
  cell_duration <- left_join(onset_df, offset_df) 
  cell_duration <- cell_duration %>% 
    mutate(duration = offset - onset)
  
  return(cell_duration)
}

#' Function to get gbif and idig bio records

get_records <- function(binomial){
  
  idig <- idig_search_records(limit = 20000, rq = list(scientificname = binomial), 
                              fields = c("data.dwc:decimalLongitude", "data.dwc:decimalLatitude",
                                         "data.dwc:year", "data.dwc:month", "data.dwc:eventDate",
                                         "data.dwc:scientificName", "data.dwc:basisOfRecord", 
                                         "data.dwc:coordinateUncertaintyInMeters", "data.dwc:institutionCode")) 
  
  names(idig) <- c("decimalLongitude", "decimalLatitude", "year", "month", "eventDate",
                   "scientificName", "basisOfRecord", "coordinateUncertaintyInMeters", "institutionCode")
  
  #rgbif
  gbif <- occ_search(scientificName = binomial,
                     hasCoordinate = TRUE, country = "US",
                     limit = 200000)
  
  gbif_df <- gbif$data %>% 
    dplyr::select(decimalLongitude, decimalLatitude, year, month, eventDate, 
                  scientificName, basisOfRecord, coordinateUncertaintyInMeters, institutionCode)
  
  total <- rbind(gbif_df, idig) %>% 
    filter(!is.na(year)) %>% 
    filter(!is.na(decimalLatitude)) %>% 
    filter(!is.na(decimalLongitude)) %>% 
    filter(!is.na(eventDate)) %>% 
    filter(decimalLatitude != 0 | decimalLongitude != 0)
  
  total$decimalLongitude <- as.numeric(total$decimalLongitude)
  total$decimalLatitude <- as.numeric(total$decimalLatitude)
  total$year <- as.numeric(total$year)
  
  total <- total %>% 
    mutate(decade = floor(year/10)*10) %>% 
    mutate(lon_bin = floor(decimalLongitude/1)*1) %>% 
    mutate(lat_bin = floor(decimalLatitude/1)*1)
  
  return(total)
}