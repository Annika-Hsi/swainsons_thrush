
# PURPOSE: run climwin sliding window on environmental data from Google Earth 
# Engine corresponding to Swainson's Thrush observations from the wintering and
# breeding seasons
# NOTES: have to remerge with original LLG dataset to get release site as covariate
# also there was a MODIS outage in April 2024 so have to remove some entries
# ISSUES: 'cohort' argument being ignored and still working on checking if best 
# models are overfitting since adj R^2 values are very close to 1

library(dplyr)
library(climwin)
library(tidyr)
library(lubridate)

# SET UP------------------------------------------------------------------------
# WINTERING
# load data
llg_env_land <- read.csv('data/llg_env_land.csv')
ndvi <- read.csv('env_data/gee/NDVI_wintering.csv')
clim <- read.csv('env_data/gee/clim_wintering.csv')
llg_orig <- read.csv('data/pheno_archival_raw_steph.csv')

# re-add IDs to old env data
llg_env_land$ID <- 1:nrow(llg_env_land)

# filter out the old climate variables that we no longer need
llg_land <- llg_env_land |> select('ID', 'release_date', 'spring_dep_new', 'release_GPS.N', 'release_GPS.W', 'fall_dep_new', 'wintering_lat', 'wintering_long', 'land_coords.lat_w', 'land_coords.lon_w')

# remove lines from original data with no fall and spring departure dates
llg_orig <- llg_orig |> 
  filter(!is.na(spring_dep_new) | !is.na(fall_dep_new)) |> 
  filter(spring_dep_new != '' | fall_dep_new != '')

# merge with data that has land coordinates
llg_merged <- merge(x = llg_land, y = llg_orig, by = c('release_GPS.N', 'release_GPS.W', 'wintering_lat', 'wintering_long'))

# convert all dates and check the duplicate columns all match
llg_merged$release_date.x <- as.Date(llg_merged$release_date.x, '%Y-%m-%d')
llg_merged$release_date.y <- as.Date(llg_merged$release_date.y, '%m/%d/%y')
sum(llg_merged$release_date.x != llg_merged$release_date.y)

llg_merged$spring_dep_new.x <- as.Date(llg_merged$spring_dep_new.x, '%Y-%m-%d')
llg_merged$spring_dep_new.y <- as.Date(llg_merged$spring_dep_new.y, '%m/%d/%Y')
sum(llg_merged$spring_dep_new.x != llg_merged$spring_dep_new.y, na.rm = TRUE)

llg_merged$fall_dep_new.x <- as.Date(llg_merged$fall_dep_new.x, '%Y-%m-%d')
llg_merged$fall_dep_new.y <- as.Date(llg_merged$fall_dep_new.y, '%m/%d/%Y')
sum(llg_merged$fall_dep_new.x != llg_merged$fall_dep_new.y, na.rm = TRUE)

# remove duplicate/excess columns
colnames(llg_merged)
llg_sub <- llg_merged[, c(1:10, 12)]

# TEST ON WINTERING NDVI--------------------------------------------------------

# get wintering data
llg_w <- llg_sub |> select('ID', 'wintering_lat', 'wintering_long', 'spring_dep_new.x', 'land_coords.lat_w', 'land_coords.lon_w', 'release_site')

# replace recorded coordinates with land coordinates where necessary
llg_w <- llg_w |> 
  mutate(lat_w = if_else(`land_coords.lat_w` == 0 | is.na(`land_coords.lat_w`), wintering_lat, `land_coords.lat_w`)) |> 
  mutate(lon_w = if_else(`land_coords.lon_w` == 0 | is.na(`land_coords.lon_w`), wintering_long, `land_coords.lon_w`))

# merge with NDVI data so we have release site 
ndvi_merged <- merge(x = llg_w, y = ndvi, by = c('ID'), all.x = FALSE)
ndvi_merged_clean <- ndvi_merged |> select('ID', 'spring_dep_new.x', 'lat_w', 'lon_w', 'date', 'ndvi', 'release_site')

# split into datasets for sliding window
dep_w <- llg_w |> select('spring_dep_new.x', 'ID', 'release_site') |> filter(!is.na(spring_dep_new.x))
ndvi_w <- ndvi_merged_clean |> select('date', 'ndvi', 'ID', 'release_site')

# convert dates
ndvi_w$date <- as.Date(ndvi_w$date)

# add cohort
dep_w$cohort <- year(dep_w$spring_dep_new.x)

# save IDs of birds with missing data
missing_birds <- (ndvi_merged_clean |>
                    group_by(ID, spring_dep_new.x) |>
                    summarise(earliest_ndvi = max(date)) |> 
                    filter(spring_dep_new.x != earliest_ndvi))$ID

# remove birds that had departure dates during the MODIS outage
dep_w_clean <- dep_w |> filter(!(ID %in% missing_birds))
ndvi_w_clean <- ndvi_w |> filter(!(ID %in% missing_birds))

# figure out reference dep date using earliest departure
(dep_w_clean |> 
    mutate(new_date = as.Date(paste0(month(spring_dep_new.x), '-', day(spring_dep_new.x), '-2018'), '%m-%d-%Y')))$new_date |> min() 
# N.B. I just picked a random year (2018) bc I needed to average the dates as if they were all the same year

# sliding window over 4 months
ndvi_w_window <- slidingwin(xvar = list(NDVI = ndvi_w_clean$ndvi),
                            cdate = ndvi_w_clean$date,
                            bdate = dep_w_clean$spring_dep_new.x,
                            baseline = lm(as.numeric(spring_dep_new.x) ~ release_site + cohort, data = dep_w_clean), # I'm confused on how to incorporate release site as a covariate...
                            cinterval = "day",
                            cmissing = 'method1',
                            #cohort = dep_w_clean$cohort,
                            range = c(120, 0),
                            type = "absolute", refday = c(03, 03),
                            stat = "mean",
                            func = "lin", 
                            spatial = list(dep_w_clean$ID, ndvi_w_clean$ID))

# examine model
summary(ndvi_w_window[[1]]$BestModel)

# TEST ON WINTERING CLIMATE DATA------------------------------------------------
# i.e., precipitation, wind speed, temperature

# merge with climate data so we have release site 
clim_merged <- merge(x = llg_w, y = clim, by = c('ID'), all.x = FALSE)
clim_merged_clean <- clim_merged |> select('ID', 'spring_dep_new.x', 'lat_w', 'lon_w', 'date', 'temp_2m', 'total_precip', 'u_wind', 'v_wind', 'release_site')

# calculate wind speed
clim_merged_clean <- clim_merged_clean |> mutate(wind_speed = sqrt(u_wind^2 + v_wind^2))

# split into datasets for sliding window
dep_w <- llg_w |> select('spring_dep_new.x', 'ID', 'release_site') |> filter(!is.na(spring_dep_new.x))
clim_w <- clim_merged_clean |> select('date', 'temp_2m', 'total_precip', 'wind_speed', 'ID', 'release_site')

# convert dates
clim_w$date <- as.Date(clim_w$date)

# add cohort
dep_w$cohort <- as.factor(year(dep_w$spring_dep_new.x))

# examine and remove NAs
clim_w |> filter(is.na(temp_2m)) |> group_by(ID) |> summarise(n_rows = n())
missing_clim <- (clim_w |> filter(is.na(temp_2m)))$ID |> unique()
clim_w_clean <- clim_w |> filter(!(ID %in% missing_clim))
dep_w_clean2 <- dep_w |> filter(!(ID %in% missing_clim))

# sliding window

# temperature
temp_w_window <- slidingwin(xvar = list(Temp = clim_w_clean$temp_2m),
                            cdate = clim_w_clean$date,
                            bdate = dep_w_clean2$spring_dep_new.x,
                            baseline = lm(as.numeric(spring_dep_new.x) ~ release_site + cohort, data = dep_w_clean2),
                            cinterval = "day",
                            cmissing = 'method2',
                            #cohort = dep_w_clean2$cohort,
                            range = c(184, 0),
                            type = "absolute", refday = c(03, 03),
                            stat = "mean",
                            func = "lin", spatial = list(dep_w_clean2$ID, clim_w_clean$ID))
# precipitation
precip_w_window <- slidingwin(xvar = list(Precip = clim_w_clean$total_precip),
                              cdate = clim_w_clean$date,
                              bdate = dep_w_clean2$spring_dep_new.x,
                              baseline = lm(as.numeric(spring_dep_new.x) ~ release_site + cohort, data = dep_w_clean2),
                              cinterval = "day",
                              cmissing = 'method1',
                              #cohort = dep_w_clean2$cohort,
                              range = c(184, 0),
                              type = "absolute", refday = c(03, 03),
                              stat = "mean",
                              func = "lin", spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# wind
wind_w_window <- slidingwin(xvar = list(Wind = clim_w_clean$wind_speed),
                            cdate = clim_w_clean$date,
                            bdate = dep_w_clean2$spring_dep_new.x,
                            baseline = lm(as.numeric(spring_dep_new.x) ~ release_site + cohort, data = dep_w_clean2),
                            cinterval = "day",
                            cmissing = 'method1',
                            #cohort = dep_w_clean2$cohort,
                            range = c(184, 0),
                            type = "absolute", refday = c(03, 03),
                            stat = "mean",
                            func = "lin", spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine models
summary(precip_w_window[[1]]$BestModel)
summary(wind_w_window[[1]]$BestModel)
summary(temp_w_window[[1]]$BestModel)

# BREEDING SEASON DATA----------------------------------------------------------
# update gee data to breeding season
ndvi <- read.csv('env_data/gee/NDVI_breeding.csv')
clim <- read.csv('env_data/gee/clim_breeding.csv')

# NDVI DATA---------------------------------------------------------------------
# examine data
ndvi$ID |> unique() |> length()
clim$ID |> unique() |> length()

(ndvi |> group_by(ID) |> summarise(n_rows = n()))$n_rows |> max()
(ndvi |> group_by(ID) |> summarise(n_rows = n()))$n_rows |> min()

(clim |> group_by(ID) |> summarise(n_rows = n()))$n_rows |> max()
(clim |> group_by(ID) |> summarise(n_rows = n()))$n_rows |> min()

# get breeding data
llg_b <- llg_sub |> select('ID', 'release_GPS.N', 'release_GPS.W', 'fall_dep_new.x', 'release_site')

# merge with NDVI data so we have release site 
ndvi_merged <- merge(x = llg_b, y = ndvi, by = c('ID'), all.x = FALSE)
ndvi_merged_clean <- ndvi_merged |> select('ID', 'fall_dep_new.x', 'lat', 'lon', 'date', 'ndvi', 'release_site')

# split into datasets for sliding window
dep_b <- llg_b |> select('fall_dep_new.x', 'ID', 'release_site') |> filter(!is.na(fall_dep_new.x))
ndvi_b <- ndvi_merged_clean |> select('date', 'ndvi', 'ID', 'release_site')

# convert dates
ndvi_b$date <- as.Date(ndvi_b$date)

# add cohort
dep_b$cohort <- year(dep_b$fall_dep_new.x)

# figure out reference dep date using earliest departure
(dep_b |> 
    mutate(new_date = as.Date(paste0(month(fall_dep_new.x), '-', day(fall_dep_new.x), '-2018'), '%m-%d-%Y')))$new_date |> min() 
# N.B. I just picked a random year (2018) bc I needed to average the dates as if they were all the same year

# sliding window over 4 months
ndvi_b_window <- slidingwin(xvar = list(NDVI = ndvi_b$ndvi),
                            cdate = ndvi_b$date,
                            bdate = dep_b$fall_dep_new.x,
                            baseline = lm(as.numeric(fall_dep_new.x) ~ release_site + cohort, data = dep_b), # I'm confused on how to incorporate release site as a covariate...
                            cinterval = "day",
                            cmissing = 'method1',
                            #cohort = dep_b$cohort,
                            range = c(184, 0),
                            type = "absolute", refday = c(18, 07),
                            stat = "mean",
                            func = "lin", spatial = list(dep_b$ID, ndvi_b$ID))

# examine model
summary(ndvi_b_window[[1]]$BestModel)

# CLIMATE DATA------------------------------------------------------------------
# merge with climate data so we have release site 
clim_merged <- merge(x = llg_b, y = clim, by = c('ID'), all.x = FALSE)
clim_merged_clean <- clim_merged |> select('ID', 'fall_dep_new.x', 'lat', 'lon', 'date', 'temp_2m', 'total_precip', 'u_wind', 'v_wind', 'release_site')

# calculate wind speed
clim_merged_clean <- clim_merged_clean |> mutate(wind_speed = sqrt(u_wind^2 + v_wind^2))

# split into datasets for sliding window
dep_b <- llg_b |> select('fall_dep_new.x', 'ID', 'release_site') |> filter(!is.na(fall_dep_new.x))
clim_b <- clim_merged_clean |> select('date', 'temp_2m', 'total_precip', 'wind_speed', 'ID', 'release_site')

# convert dates
clim_b$date <- as.Date(clim_b$date)

# add cohort
dep_b$cohort <- as.factor(year(dep_b$fall_dep_new.x))

# sliding window

# temperature
temp_b_window <- slidingwin(xvar = list(Temp = clim_b$temp_2m),
                            cdate = clim_b$date,
                            bdate = dep_b$fall_dep_new.x,
                            baseline = lm(as.numeric(fall_dep_new.x) ~ release_site + cohort, data = dep_b),
                            cinterval = "day",
                            cmissing = 'method2',
                            #cohort = dep_b$cohort,
                            range = c(184, 0),
                            type = "absolute", refday = c(18, 07),
                            stat = "mean",
                            func = "lin", spatial = list(dep_b$ID, clim_b$ID))

# precipitation
precip_b_window <- slidingwin(xvar = list(Temp = clim_b$total_precip),
                              cdate = clim_b$date,
                              bdate = dep_b$fall_dep_new.x,
                              baseline = lm(as.numeric(fall_dep_new.x) ~ release_site + cohort, data = dep_b),
                              cinterval = "day",
                              cmissing = 'method2',
                              #cohort = dep_b$cohort,
                              range = c(184, 0),
                              type = "absolute", refday = c(18, 07),
                              stat = "mean",
                              func = "lin", spatial = list(dep_b$ID, clim_b$ID))

# wind
wind_b_window <- slidingwin(xvar = list(Temp = clim_b$wind_speed),
                            cdate = clim_b$date,
                            bdate = dep_b$fall_dep_new.x,
                            baseline = lm(as.numeric(fall_dep_new.x) ~ release_site + cohort, data = dep_b),
                            cinterval = "day",
                            cmissing = 'method2',
                            #cohort = dep_b$cohort,
                            range = c(184, 0),
                            type = "absolute", refday = c(18, 07),
                            stat = "mean",
                            func = "lin", spatial = list(dep_b$ID, clim_b$ID))

# examine models
summary(precip_b_window[[1]]$BestModel)
summary(wind_b_window[[1]]$BestModel)
summary(temp_b_window[[1]]$BestModel)

