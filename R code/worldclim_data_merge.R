# PURPOSE: merge past 6 months of temperature, precipitation, photoperiod data
# NOTES: b = breeding grounds, w = wintering grounds
# assuming birds decide fall departure date based on breeding ground conditions
# and spring departure date based on wintering ground conditions
# also, this is for the monthly worldclim data
# IMPORTANT: THIS HAS CODE FOR GETTING LAND COORDS FOR ERRONEOUS LLG OBSERVATIONS THAT WERE IN THE OCEAN

# load packages
library(raster)
library(dplyr)
library(lubridate)
library(sp)
library(geosphere)
library(ggplot2)
library(rnaturalearth)


# load llg data
llg <- read.csv('data/pheno_archival_raw_steph.csv')

# convert dates
llg$fall_dep_new
llg$fall_dep_new <- as.Date(llg$fall_dep_new, '%m/%d/%Y')
range(llg$fall_dep_new, na.rm = TRUE)

llg$spring_dep_new
llg$spring_dep_new <- as.Date(llg$spring_dep_new, '%m/%d/%Y')
range(llg$spring_dep_new, na.rm = TRUE)

llg$release_date
llg$release_date <- as.Date(llg$release_date, '%m/%d/%y')
range(llg$release_date, na.rm = TRUE)

# subset to necessary columns
llg_sub <- llg |> dplyr::select(release_date, spring_dep_new, release_GPS.N, 
                         release_GPS.W, fall_dep_new, wintering_lat, wintering_long) |> 
  filter(!is.na(spring_dep_new) | !is.na(fall_dep_new))

# prefixes for files
prec_name <- 'prec/wc2.1_cruts4.09_2.5m_prec_'
tmin_name <- 'tmin/wc2.1_cruts4.09_2.5m_tmin_'
tmax_name <- 'tmax/wc2.1_cruts4.09_2.5m_tmax_'

# path for data in flash
data_path <- '/Volumes/My Passport/SUMMER25/Delmore_Lab/data/'

# number of months to go back
n_mos <- 7

# to store loop data in for all birds and climate variables
clim_df <- as.data.frame(matrix(nrow = nrow(llg_sub), ncol = n_mos))
clim_data <- list(prec_b = clim_df, prec_w = clim_df, tmin_b = clim_df, 
                  tmin_w = clim_df, tmax_b = clim_df, tmax_w = clim_df)

# FUNCTIONS---------------------------------------------------------------------
# function to get nearest land (non-NA) pixel to given coordinates
nearest_land <- function(lat, lon, raster) {
  # using set buffer size
  max_dist <- 500000
  pt <- SpatialPoints(data.frame(x = lon, y = lat))
  
  # get non-NA pixels within buffer
  neighbors <- extract(raster, pt, buffer = max_dist, cellnumbers = TRUE, df = TRUE)
  colnames(neighbors) <- c('id', 'cell', 'value')
  neighbors <- as.data.frame(neighbors) |> filter(!is.na(value))
  
  # calculate distance from given coords
  neighbors_coords <- xyFromCell(raster, neighbors$cell)
  neighbors <- cbind(neighbors, neighbors_coords)
  neighbors$dist <- distm(neighbors_coords, c(lon, lat), fun = distHaversine)
  idx <- which(neighbors$dist == min(neighbors$dist))
  
  list(val = neighbors$value[idx], lat = neighbors$y[idx], lon = neighbors$x[idx])
}
# function to extract raster data from given location and date
extract_pt <- function(dt, lat, lon, path) {
  
  # list to store data
  clim_list <- list(prec = c(), tmin = c(), tmax = c(), land_lat = 0, land_lon = 0)
  
  # get location as sp object
  pt <- SpatialPoints(data.frame(x = lon, y = lat))
  
  # get dates for previous 6 months
  all_dt <- seq(dt, length = n_mos, by = '-31 days')
  
  # ocean flag 
  ocean_flag <- FALSE
  # get data for 6 previous months and store
  for (j in 1:length(all_dt)) {
    # get current date
    curr_dt <- all_dt[j]
    
    # load rasters
    prec_rast <- raster(paste0(path, prec_name, year(curr_dt), '-', format(curr_dt, '%m'), '.tif'))
    tmin_rast <-raster(paste0(path, tmin_name, year(curr_dt), '-', format(curr_dt, '%m'), '.tif'))
    tmax_rast <-raster(paste0(path, tmax_name, year(curr_dt), '-', format(curr_dt, '%m'), '.tif'))
    
    # get value at exact point
    clim_list$prec[j] <- raster::extract(prec_rast, pt)
    clim_list$tmin[j] <- raster::extract(tmin_rast, pt)
    clim_list$tmax[j] <- raster::extract(tmax_rast, pt)
    
    # break if ocean pixel
    if (is.na(clim_list$prec[j]) & is.na(clim_list$tmin[j]) & is.na(clim_list$tmax[j])) {
      ocean_flag <- TRUE
      break
    }
  }
  
  if (ocean_flag == TRUE) {
    land_lat <- lat
    land_lon <- lon
    for (j in 1:length(all_dt)) {
      # get current date
      curr_dt <- all_dt[j]
      
      # load rasters
      prec_rast <- raster(paste0(path, prec_name, year(curr_dt), '-', format(curr_dt, '%m'), '.tif'))
      tmin_rast <-raster(paste0(path, tmin_name, year(curr_dt), '-', format(curr_dt, '%m'), '.tif'))
      tmax_rast <-raster(paste0(path, tmax_name, year(curr_dt), '-', format(curr_dt, '%m'), '.tif'))
      
      # get lat + lon of nearest land pixel just first time
      if (j == 1) {
        land_info <- nearest_land(lat, lon, prec_rast)
        land_lat <- land_info$lat
        land_lon <- land_info$lon
        pt <- SpatialPoints(data.frame(x = land_lon, y = land_lat))
        clim_list$land_lat <- land_lat
        clim_list$land_lon <- land_lon
      }
      
      # get value at exact point
      clim_list$prec[j] <- raster::extract(prec_rast, pt)
      clim_list$tmin[j] <- raster::extract(tmin_rast, pt)
      clim_list$tmax[j] <- raster::extract(tmax_rast, pt)
    }
  }
  
  # return list of climate variables over past 6 months
  clim_list
}

photoperiod <- function(dt, lat) {
  lens <- c()
  all_dt <- seq(dt, length = n_mos, by = '-31 days')
  for (k in 1:length(all_dt)) {
    lens[k] <- daylength(lat, all_dt[k])
  }
  lens
}
#-------------------------------------------------------------------------------

# vectors to track any ocean points and the nearest land coordinates
land_coords <- list(lat_b = rep(NA, nrow(llg_sub)), lon_b = rep(NA, nrow(llg_sub)),
                    lat_w = rep(NA, nrow(llg_sub)), lon_w = rep(NA, nrow(llg_sub)))

# for loop to iterate through all birds
for (i in 1:nrow(llg_sub)) {
  # get current bird
  bird <- llg_sub[i, ]
  
  # extract datq
  dt_w <- bird$spring_dep_new
  dt_b <- bird$fall_dep_new
  lat_w <- bird$wintering_lat
  lat_b  <- bird$release_GPS.N
  lon_w <- bird$wintering_long
  lon_b <- bird$release_GPS.W
  
  if (!is.na(dt_b) & !is.na(lat_b) & !is.na(lon_b)) {
    clim_b <- extract_pt(dt_b, lat_b, lon_b, data_path)
    
    # put data into big list
    clim_data$prec_b[i, ] <- clim_b$prec
    clim_data$tmin_b[i, ] <- clim_b$tmin
    clim_data$tmax_b[i, ] <- clim_b$tmax
    land_coords$lat_b[i] <- clim_b$land_lat
    land_coords$lon_b[i] <- clim_b$land_lon
  }
  if (!is.na(dt_w) & !is.na(lat_w) & !is.na(lon_w)) {
    clim_w <- extract_pt(dt_w, lat_w, lon_w, data_path)
    
    # put data into big list
    clim_data$prec_w[i, ] <- clim_w$prec
    clim_data$tmin_w[i, ] <- clim_w$tmin
    clim_data$tmax_w[i, ] <- clim_w$tmax
    land_coords$lat_w[i] <- clim_w$land_lat
    land_coords$lon_w[i] <- clim_w$land_lon
  }
}

# rename cols
colnames(clim_data$prec_b) <- paste0('precb_offset_', seq(0, 6))
colnames(clim_data$tmin_b) <- paste0('tminb_offset_', seq(0, 6))
colnames(clim_data$tmax_b) <- paste0('tmaxb_offset_', seq(0, 6))
colnames(clim_data$prec_w) <- paste0('precw_offset_', seq(0, 6))
colnames(clim_data$tmin_w) <- paste0('tminw_offset_', seq(0, 6))
colnames(clim_data$tmax_w) <- paste0('tmaxw_offset_', seq(0, 6))

# add to main df
llg_env <- cbind(llg_sub, clim_data$prec_b, clim_data$tmin_b, clim_data$tmax_b, 
                 clim_data$prec_w, clim_data$tmin_w, clim_data$tmax_w)

# get daylengths
daylens_w <- clim_df
daylens_b <- clim_df
for (j in 1:nrow(llg_sub)) {
  bird <- llg_sub[j, ]
  
  # get lats and dates
  dt_b <- bird$fall_dep_new
  dt_w <- bird$spring_dep_new
  lat_b <- bird$release_GPS.N
  lat_w <- bird$wintering_lat
  
  # using photoperiod at breeding grounds for fall migration 
  # and at wintering grounds for spring migration
  if (!is.na(dt_b) & !is.na(lat_b)) {
    daylens_b[j, ] <- photoperiod(dt_b, lat_b)
  }
  
  if (!is.na(dt_w) & !is.na(lat_w)) {
    daylens_w[j, ] <- photoperiod(dt_w, lat_w)
  }
}

# rename columns
colnames(daylens_b) <- paste0('daylenb_offset_', seq(0, 6))
colnames(daylens_w) <- paste0('daylenw_offset_', seq(0, 6))

# add to main df
# can include land_coords vector here for land alternatives to ocean pixelss
llg_env <- cbind(llg_env, daylens_b, daylens_w) 

write.csv(llg_env, 'data/llg_env.csv')

# check that land coordinates are actually reasonable
# land_coords
# comp_df <- data.frame(lat = llg_sub$wintering_lat, lon = llg_sub$wintering_long,
#                       land_lat = land_coords$lat_w, land_lon = land_coords$lon_w)
# 
# comp_df <- comp_df |> filter(!is.na(land_lat) & (land_lat != 0))
# comp_df$id <- as.factor(1:nrow(comp_df))
# 
# # plot each point
# n_america <- ne_countries(continent = 'north america')
# nums <-34:36
# ggplot() +
#   geom_sf(data = n_america) + 
#   geom_point(data = comp_df[nums, ], aes(x = lon, y = lat), color = 'cornflowerblue') +
#   geom_point(data = comp_df[nums, ], aes(x = land_lon, y = land_lat), color = 'brown') +
#   facet_grid(id ~ .)

