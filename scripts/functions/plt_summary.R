library(rnaturalearth)
library(ggplot2)
library(sf)

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

## Test

#ten <- plt_summary(cell_size = 100000, dat = en, n_per_cell = 10)
#ten$fig
