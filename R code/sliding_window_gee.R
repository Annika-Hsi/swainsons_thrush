
# PURPOSE: run climwin sliding window on environmental data from Google Earth 
# Engine corresponding to Swainson's Thrush observations from the wintering and
# breeding seasons
# NOTES: have to remerge with original LLG dataset to get release site as covariate
# also there was a MODIS outage in April 2024 so have to remove some entries
# ISSUES: cross validation does not always work because some release sites have 
# very few entries--just rerun until it does not give error

# load packages
library(dplyr)
library(climwin)
library(tidyr)
library(lubridate)
library(caret)
library(car)
library(geosphere)

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
sum(llg_merged$release_date.x != llg_merged$release_date.y) # should be 0

llg_merged$spring_dep_new.x <- as.Date(llg_merged$spring_dep_new.x, '%Y-%m-%d')
llg_merged$spring_dep_new.y <- as.Date(llg_merged$spring_dep_new.y, '%m/%d/%Y')
sum(llg_merged$spring_dep_new.x != llg_merged$spring_dep_new.y, na.rm = TRUE)

llg_merged$fall_dep_new.x <- as.Date(llg_merged$fall_dep_new.x, '%Y-%m-%d')
llg_merged$fall_dep_new.y <- as.Date(llg_merged$fall_dep_new.y, '%m/%d/%Y')
sum(llg_merged$fall_dep_new.x != llg_merged$fall_dep_new.y, na.rm = TRUE)

# remove duplicate/excess columns
colnames(llg_merged)
llg_sub <- llg_merged[, c(1:10, 12)]

llg_sub$release_site <- (llg_sub |> 
                           mutate(release_site_comb = case_when(release_site == 'Pacific Spirit' ~ 'Porpoise Bay',
                                                                release_site == 'Northern_BC' ~ 'Alaska',
                                                                ((release_site != 'Pacific Spirit') & (release_site != 'Northern_BC')) ~ release_site
                           )))$release_site_comb

# WINTERING NDVI----------------------------------------------------------------

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
set.seed(123)
ndvi_w_window <- slidingwin(xvar = list(NDVI = ndvi_w_clean$ndvi),
                            k = 3,
                            cdate = ndvi_w_clean$date,
                            bdate = dep_w_clean$spring_dep_new.x,
                            baseline = lm(as.numeric(spring_dep_new.x) ~ release_site, data = dep_w_clean), # I'm confused on how to incorporate release site as a covariate...
                            cinterval = "day",
                            cmissing = 'method1',
                            range = c(120, 0),
                            type = "absolute", refday = c(03, 03),
                            stat = "mean",
                            func = "lin", 
                            spatial = list(dep_w_clean$ID, ndvi_w_clean$ID))

# examine model
ndvi_w_mod <- ndvi_w_window[[1]]$BestModel
summary(ndvi_w_mod)
vif(ndvi_w_mod)

# function to extract and plot best data
window_plot <- function(window, varname, type) {
  # extract best data for plots
  bestdata <- window[[1]]$BestModelData
  bestdata$yvar <- as.Date(bestdata$yvar, origin = '1970-01-01')
  bestdata$yvar |> yday()
  
  # faceted plot
  if (type == 'facet year') {
    plot <- ggplot(bestdata, aes(x = yvar, y = climate)) +
      geom_point() +
      facet_wrap(as.factor(year(yvar))~., scales = 'free') +
      labs(x = 'Departure Date',
           y = varname,
           title = paste0(varname,' vs. Departure Date',' by Year')) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  }
  if (type == 'facet site') {
    plot <- ggplot(bestdata, aes(x = yday(yvar), y = climate)) +
      geom_point() +
      facet_wrap(release_site~., scales = 'free') +
      labs(x = 'Departure Day of Year',
           y = varname,
           title = paste0(varname,' vs. Departure Day of Year',' by Release Site'))
  }
  if (type == 'year') {
    plot <- ggplot(bestdata, aes(x = yday(yvar), y = climate, color = as.factor(year(yvar)))) +
      geom_point() +
      labs(x = 'Departure Day of Year',
           y = varname,
           title = paste0(varname,' vs. Departure Day of Year')) +
      scale_color_discrete(name = 'Year')
  }
  if (type == 'site') {
    plot <- ggplot(bestdata, aes(x = yday(yvar), y = climate, color = release_site)) +
      geom_point() +
      labs(x = 'Departure Day of Year',
           y = varname,
           title = paste0(varname,' vs. Departure Day of Year')) +
      scale_color_discrete(name = 'Release Site')
  }
  return(plot)
}

# plots
window_plot(ndvi_w_window, 'NDVI', 'facet year')
window_plot(ndvi_w_window, 'NDVI', 'facet site')
window_plot(ndvi_w_window, 'NDVI', 'year')
window_plot(ndvi_w_window, 'NDVI', 'site')

# WINTERING CLIMATE DATA--------------------------------------------------------
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

# examine and remove NAs
clim_w |> filter(is.na(temp_2m)) |> group_by(ID) |> summarise(n_rows = n())
missing_clim <- (clim_w |> filter(is.na(temp_2m)))$ID |> unique()
clim_w_clean <- clim_w |> filter(!(ID %in% missing_clim))
dep_w_clean2 <- dep_w |> filter(!(ID %in% missing_clim))

# sliding window

# temperature
set.seed(123)
temp_w_window <- slidingwin(xvar = list(Temp = clim_w_clean$temp_2m),
                            k = 3,
                            cdate = clim_w_clean$date,
                            bdate = dep_w_clean2$spring_dep_new.x,
                            baseline = lm(as.numeric(spring_dep_new.x) ~ release_site, data = dep_w_clean2),
                            cinterval = "day",
                            cmissing = 'method2',
                            range = c(184, 0),
                            type = "absolute", refday = c(03, 03),
                            stat = "mean",
                            func = "lin", spatial = list(dep_w_clean2$ID, clim_w_clean$ID))
# precipitation
set.seed(123)
precip_w_window <- slidingwin(xvar = list(Precip = clim_w_clean$total_precip),
                              k = 3,
                              cdate = clim_w_clean$date,
                              bdate = dep_w_clean2$spring_dep_new.x,
                              baseline = lm(as.numeric(spring_dep_new.x) ~ release_site, data = dep_w_clean2),
                              cinterval = "day",
                              cmissing = 'method1',
                              range = c(184, 0),
                              type = "absolute", refday = c(03, 03),
                              stat = "mean",
                              func = "lin", spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# wind
set.seed(123)
wind_w_window <- slidingwin(xvar = list(Wind = clim_w_clean$wind_speed),
                            k = 3,
                            cdate = clim_w_clean$date,
                            bdate = dep_w_clean2$spring_dep_new.x,
                            baseline = lm(as.numeric(spring_dep_new.x) ~ release_site, data = dep_w_clean2),
                            cinterval = "day",
                            cmissing = 'method1',
                            range = c(184, 0),
                            type = "absolute", refday = c(03, 03),
                            stat = "mean",
                            func = "lin", spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine models
summary(temp_w_window[[1]]$BestModel) 
summary(precip_w_window[[1]]$BestModel)
summary(wind_w_window[[1]]$BestModel)

vif(temp_w_window[[1]]$BestModel) 
vif(precip_w_window[[1]]$BestModel)
vif(wind_w_window[[1]]$BestModel)

# plots
window_plot(temp_w_window, 'Temperature', 'facet year')
window_plot(temp_w_window, 'Temperature', 'facet site')
window_plot(temp_w_window, 'Temperature', 'year')
window_plot(temp_w_window, 'Temperature', 'site')
window_plot(precip_w_window, 'Precipitation', 'facet year')
window_plot(precip_w_window, 'Precipitation', 'facet site')
window_plot(precip_w_window, 'Precipitation', 'year')
window_plot(precip_w_window, 'Precipitation', 'site')
window_plot(wind_w_window, 'Wind Speed', 'facet year')
window_plot(wind_w_window, 'Wind Speed', 'facet site')
window_plot(wind_w_window, 'Wind Speed', 'year')
window_plot(wind_w_window, 'Wind Speed', 'site')

# WINTERING DAYLENGTH-----------------------------------------------------------
# get rid of any rows missing lat or date
llg_w_clean <- llg_w |> 
  filter((!is.na(spring_dep_new.x)) & !is.na(lat_w)) |> 
  arrange(ID)

# figure out latest day of year departure
max_dt <- (llg_w_clean |> mutate(doy = yday(spring_dep_new.x)) |> arrange(-doy))[1, ]$spring_dep_new.x

# figure out 185 days back from earliest departure date (March 3rd)
as.Date('2018-03-03') - 185

# function to get sequence of dates from each departure date
get_all_dts_w <- function(dt) {
  yr <- year(dt)
  return(seq.Date(from = dt, to = (as.Date(paste0(yr, '-03-03'), '%Y-%m-%d')-185), by = '-1 days'))
}

all_dts <- lapply(llg_w_clean$spring_dep_new.x, get_all_dts_w)

# convert list of lists to one long data frame labeled by ID
dts_df <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dts_df) <- c('date', 'ID')
for (i in 1:nrow(llg_w_clean)) {
  curr_id <- llg_w_clean$ID[i]
  curr_list <- all_dts[[i]]
  df <- data.frame(date = curr_list)
  df$ID <- curr_id
  dts_df <- rbind(dts_df, df)
}

# check number of rows per bird look correct
(dts_df |> group_by(ID) |> 
    summarise(nrow = n()))$nrow |> range()

# append dates
llg_dts_w <- merge(x = llg_w_clean, y = dts_df, by = 'ID', all.y = TRUE)

# extract photoperiod
llg_dts_w <- llg_dts_w |> mutate(daylen = daylength(lat_w, date))

# split into 2 data frames
daylen_w <- llg_dts_w |> select('ID', 'release_site', 'date', 'daylen')
dep_daylen_w <- llg_dts_w |> select('ID', 'release_site', 'spring_dep_new.x')

# run sliding window
set.seed(123)
daylen_w_window <- slidingwin(xvar = list(Daylength = daylen_w$daylen),
                              k = 3,
                              cdate = daylen_w$date,
                              bdate = dep_daylen_w$spring_dep_new.x,
                              baseline = lm(as.numeric(spring_dep_new.x) ~ release_site, data = dep_daylen_w),
                              cinterval = "day",
                              cmissing = 'method2',
                              range = c(184, 0),
                              type = "absolute", refday = c(03, 03),
                              stat = "mean",
                              func = "lin", spatial = list(dep_daylen_w$ID, daylen_w$ID))

# examine models
summary(daylen_w_window[[1]]$BestModel)
vif(daylen_w_window[[1]]$BestModel)

# plots
window_plot(daylen_w_window, 'Day Length', 'facet year')
window_plot(daylen_w_window, 'Day Length', 'facet site')
window_plot(daylen_w_window, 'Day Length', 'year')
window_plot(daylen_w_window, 'Day Length', 'site')

# BREEDING SEASON DATA----------------------------------------------------------
# update gee data to breeding season
ndvi <- read.csv('env_data/gee/NDVI_breeding.csv')
clim <- read.csv('env_data/gee/clim_breeding.csv')

# BREEDING NDVI-----------------------------------------------------------------
# examine data

# number of unique birds
ndvi$ID |> unique() |> length()
clim$ID |> unique() |> length()

# minimum and maximum number of observations per bird
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

# figure out reference dep date using earliest departure
(dep_b |> 
    mutate(new_date = as.Date(paste0(month(fall_dep_new.x), '-', day(fall_dep_new.x), '-2018'), '%m-%d-%Y')))$new_date |> min() 
# N.B. I just picked a random year (2018) bc I needed to average the dates as if they were all the same year

# sliding window over 4 months
set.seed(123)
ndvi_b_window <- slidingwin(xvar = list(NDVI = ndvi_b$ndvi),
                            k = 3,
                            cdate = ndvi_b$date,
                            bdate = dep_b$fall_dep_new.x,
                            baseline = lm(as.numeric(fall_dep_new.x) ~ release_site, data = dep_b), # I'm confused on how to incorporate release site as a covariate...
                            cinterval = "day",
                            cmissing = 'method1',
                            range = c(184, 0),
                            type = "absolute", refday = c(18, 07),
                            stat = "mean",
                            func = "lin", spatial = list(dep_b$ID, ndvi_b$ID))

# examine model
summary(ndvi_b_window[[1]]$BestModel) 
vif(ndvi_b_window[[1]]$BestModel)

# plots
window_plot(ndvi_b_window, 'NDVI', 'facet year')
window_plot(ndvi_b_window, 'NDVI', 'facet site')
window_plot(ndvi_b_window, 'NDVI', 'year')
window_plot(ndvi_b_window, 'NDVI', 'site')

# BREEDING CLIMATE DATA---------------------------------------------------------
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

# sliding window

# temperature
set.seed(123)
temp_b_window <- slidingwin(xvar = list(Temp = clim_b$temp_2m),
                            k = 3,
                            cdate = clim_b$date,
                            bdate = dep_b$fall_dep_new.x,
                            baseline = lm(as.numeric(fall_dep_new.x) ~ release_site, data = dep_b),
                            cinterval = "day",
                            cmissing = 'method2',
                            range = c(184, 0),
                            type = "absolute", refday = c(18, 07),
                            stat = "mean",
                            func = "lin", spatial = list(dep_b$ID, clim_b$ID))

# precipitation
set.seed(123)
precip_b_window <- slidingwin(xvar = list(Precip = clim_b$total_precip),
                              k = 3,
                              cdate = clim_b$date,
                              bdate = dep_b$fall_dep_new.x,
                              baseline = lm(as.numeric(fall_dep_new.x) ~ release_site, data = dep_b),
                              cinterval = "day",
                              cmissing = 'method2',
                              range = c(184, 0),
                              type = "absolute", refday = c(18, 07),
                              stat = "mean",
                              func = "lin", spatial = list(dep_b$ID, clim_b$ID))

# wind
set.seed(123)
wind_b_window <- slidingwin(xvar = list(Wind = clim_b$wind_speed),
                            k = 3,
                            cdate = clim_b$date,
                            bdate = dep_b$fall_dep_new.x,
                            baseline = lm(as.numeric(fall_dep_new.x) ~ release_site, data = dep_b),
                            cinterval = "day",
                            cmissing = 'method2',
                            range = c(184, 0),
                            type = "absolute", refday = c(18, 07),
                            stat = "mean",
                            func = "lin", spatial = list(dep_b$ID, clim_b$ID))

# examine models
summary(temp_b_window[[1]]$BestModel)
summary(precip_b_window[[1]]$BestModel)
summary(wind_b_window[[1]]$BestModel)

vif(temp_b_window[[1]]$BestModel)
vif(precip_b_window[[1]]$BestModel)
vif(wind_b_window[[1]]$BestModel)

# plots
window_plot(temp_b_window, 'Temperature', 'facet year')
window_plot(temp_b_window, 'Temperature', 'facet site')
window_plot(temp_b_window, 'Temperature', 'year')
window_plot(temp_b_window, 'Temperature', 'site')
window_plot(precip_b_window, 'Precipitation', 'facet year')
window_plot(precip_b_window, 'Precipitation', 'facet site')
window_plot(precip_b_window, 'Precipitation', 'year')
window_plot(precip_b_window, 'Precipitation', 'site')
window_plot(wind_b_window, 'Wind Speed', 'facet year')
window_plot(wind_b_window, 'Wind Speed', 'facet site')
window_plot(wind_b_window, 'Wind Speed', 'year')
window_plot(wind_b_window, 'Wind Speed', 'site')

# BREEDING DAYLENGTH------------------------------------------------------------
# get rid of any rows missing lat or date
llg_b_clean <- llg_b |> 
  filter((!is.na(fall_dep_new.x)) & !is.na(release_GPS.N)) |> 
  arrange(ID)

# figure out latest day of year departure
max_dt <- (llg_b_clean |> mutate(doy = yday(fall_dep_new.x)) |> arrange(-doy))[1, ]$fall_dep_new.x

# figure out 185 days back from earliest departure date (March 3rd)
as.Date('2018-07-18') - 185

# function to get sequence of dates from each departure date
get_all_dts_b <- function(dt) {
  yr <- year(dt)
  return(seq.Date(from = dt, to = (as.Date(paste0(yr, '-07-18'), '%Y-%m-%d')-185), by = '-1 days'))
}

all_dts <- lapply(llg_b_clean$fall_dep_new.x, get_all_dts_b)

# convert list of lists to one long data frame labeled by ID
dts_df <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(dts_df) <- c('date', 'ID')
for (i in 1:nrow(llg_b_clean)) {
  curr_id <- llg_b_clean$ID[i]
  curr_list <- all_dts[[i]]
  df <- data.frame(date = curr_list)
  df$ID <- curr_id
  dts_df <- rbind(dts_df, df)
}

# check number of rows per bird look correct
(dts_df |> group_by(ID) |> 
    summarise(nrow = n()))$nrow |> range()

# append dates
llg_dts_b <- merge(x = llg_b_clean, y = dts_df, by = 'ID', all.y = TRUE)

# extract photoperiod
llg_dts_b <- llg_dts_b |> mutate(daylen = daylength(release_GPS.N, date))

# split into 2 data frames
daylen_b <- llg_dts_b |> select('ID', 'release_site', 'date', 'daylen')
dep_daylen_b <- llg_dts_b |> select('ID', 'release_site', 'fall_dep_new.x')

# run sliding window
set.seed(123)
daylen_b_window <- slidingwin(xvar = list(Daylength = daylen_b$daylen),
                              k = 3,
                              cdate = daylen_b$date,
                              bdate = dep_daylen_b$fall_dep_new.x,
                              baseline = lm(as.numeric(fall_dep_new.x) ~ release_site, data = dep_daylen_b),
                              cinterval = "day",
                              cmissing = 'method2',
                              range = c(184, 0),
                              type = "absolute", refday = c(18, 07),
                              stat = "mean",
                              func = "lin", spatial = list(dep_daylen_b$ID, daylen_b$ID))

# examine model
daylen_b_window
summary(daylen_b_window[[1]]$BestModel)
vif(daylen_b_window[[1]]$BestModel) # high multicollinearity

# plots
window_plot(daylen_b_window, 'Day Length', 'facet year')
window_plot(daylen_b_window, 'Day Length', 'facet site')
window_plot(daylen_b_window, 'Day Length', 'year')
window_plot(daylen_b_window, 'Day Length', 'site')
