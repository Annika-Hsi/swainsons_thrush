# ROUGH DRAFT VERSION: testing different baseline models
# PURPOSE: run climwin sliding window on environmental data from Google Earth 
# Engine corresponding to Swainson's Thrush observations from the wintering 
# season only
# NOTES: have to remerge with original LLG dataset to get release site as covariate;
# ENSO data not currently being used, but added in for future use

# load packages
library(dplyr)
library(climwin)
library(tidyr)
library(lubridate)
library(caret)
library(car)
library(geosphere)
library(lme4)
library(MuMIn)
library(survival)
library(lmtest)

# SET UP------------------------------------------------------------------------------------------------------------------------------------------------
# WINTERING
# load data
llg_env_land <- read.csv('data/llg_env_land.csv')
ndvi <- read.csv('env_data/gee/NDVI_wintering.csv')
clim <- read.csv('env_data/gee/clim_wintering.csv')
llg_orig <- read.csv('data/pheno_archival_raw_steph.csv')
llg_ancestry <- read.csv('data/AIMs_metadata.20250327.csv')

# re-add IDs to old env data
llg_env_land$ID <- 1:nrow(llg_env_land)

# filter out the old climate variables that we no longer need
llg_land <- llg_env_land |> select('ID', 'release_date', 'spring_dep_new', 'release_GPS.N', 'release_GPS.W', 'fall_dep_new', 'wintering_lat', 'wintering_long', 'land_coords.lat_w', 'land_coords.lon_w')

# remove lines from original data with no fall and spring departure dates
llg_orig <- llg_orig |> 
  filter(!is.na(spring_dep_new) | !is.na(fall_dep_new)) |> 
  filter(spring_dep_new != '' | fall_dep_new != '')

# add ancestry data in case we want to use it in the future
llg_ancestry <- llg_ancestry |> select('reference', 'aims_ancestry', 'aims_heterozygosity')
llg_orig <- merge(x = llg_orig, y = llg_ancestry, by = c('reference'))

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
llg_sub <- llg_merged[, c(1:5, 7:12, 46:47)]

llg_sub$release_site <- (llg_sub |> 
                           mutate(release_site_comb = case_when(release_site == 'Pacific Spirit' ~ 'Porpoise Bay',
                                                                release_site == 'Northern_BC' ~ 'Alaska',
                                                                ((release_site != 'Pacific Spirit') & (release_site != 'Northern_BC')) ~ release_site
                           )))$release_site_comb

# add transformed ancestry variable
llg_sub$transf_ancestry <- ifelse(llg_sub$aims_ancestry < .5, llg_sub$aims_ancestry, 1 - llg_sub$aims_ancestry)

# FUNCTIONS------------------------------------------------------------------------------------------------------------------------------------------
# helper function to check if window overlaps with other windows in list
overlaps <- function(win_opens, win_closes, o, c) {
  flag <- FALSE
  for (j in 1:length(win_opens)) {
    v1 <- c(win_opens[j], win_closes[j])
    v2 <- c(o, c)
    # swap so that v1 has greater start than v2
    if (v1[1] < v2[1]) {
      temp <- v1
      v1 <- v2
      v2 <- temp
    }
    if (v2[1] >= v1[2]) {
      flag <- TRUE
      break
    }
  }
  return(flag)
}

# function that takes climwin obj + finds all non-overlapping windows within 2 
# deltaAICs of best model and have climate as semi-significant predictor of 
# departure date (p <= .1)
# NOTE: climate variable must be renamed to 'climate' and departure date must be 
# renamed to 'dep_date'
# REFERENCE: Burnham KP, Anderson DR. 2002 Model selection and multimodel 
# inference: a practical information-theoretic approach. New York, NY: Springer.
all_win <- function(win_obj, clim_data, ref_day, varname, model, agg) {
  bestmod <- win_obj[[1]]$BestModel
  coefs <- summary(bestmod)$coefficients
  
  # if best model doesn't have climate as semi-significant predictor (p >= .1)
  if (coefs[nrow(coefs), ncol(coefs)] > .1) {
    return(paste0(varname, ' not significant predictor of departure date (p >= .1)'))
    
    # if linear model has poor R^2
    if (model == 'linear') {
      if (summary(bestmod)$adj.r.squared < 0.01) {
        return('poor model (adjusted R^2 < .01)')
      }
    }
    
    # if coxph model has poor R^2
    if (model == 'coxph') {
      if (summary(bestmod)$rsq[[1]]/summary(bestmod)$rsq[[2]] < 0.01) {
        return('poor model (Nagelkerke R^2 < .01)')
      }
    }
  }
  
  # get all windows within 2 delta AICs of best
  dataset <- win_obj[[1]]$Dataset
  maxdAIC <- dataset$deltaAICc[1]
  dataset <- dataset |> filter(deltaAICc < maxdAIC + 2)
  
  # list of window opens + closes
  win_opens <- c()
  win_closes <- c()
  for (i in 1:nrow(dataset)) {
    # current window info
    row <- dataset[i, ]
    open <- row$WindowOpen
    close <- row$WindowClose
    
    # check if overlaps with any recorded windows; go to next window if so
    if (length(win_opens) > 0) {
      if (overlaps(win_opens, win_closes, open, close) == TRUE) {
        next
      }
    }
    
    # get data
    # using mean as aggregate
    if (agg == 'mean') {
      temp_clim <- clim_data |>
        mutate(ref_date = as.Date(paste0(dep_year, ref_day)),
               open_date = ref_date - open,
               close_date = ref_date - close) |>
        group_by(ID, dep_date, release_site, dep_year, transf_ancestry, aims_heterozygosity) |>
        filter((date >= open_date) & (date <= close_date)) |>
        summarise(agg_clim = mean(climate, na.rm = TRUE))
    }
    
    # using slope as aggregate
    if (agg == 'slope') {
      temp_clim0 <- clim_data |>
        mutate(ref_date = as.Date(paste0(dep_year, ref_day)),
               open_date = ref_date - open,
               close_date = ref_date - close) 
      keep_clim <- (temp_clim1 |> 
                      group_by(ID) |> 
                      summarise(min = min(date), max = max(date)) |> 
                      filter(min != max))$ID
      temp_clim <- temp_clim0 |> 
        filter(ID %in% keep_clim) |> 
        group_by(ID, dep_date, release_site, dep_year, transf_ancestry, aims_heterozygosity) |>
        filter((date >= open_date) & (date <= close_date)) |>
        summarise(agg_clim = lm(climate ~ as.numeric(as.Date(date)))$coefficients[[2]])
    }
    
    # get model
    # if linear model
    if (model == 'linear') {
      # model
      mod <- lm(yday(dep_date) ~ transf_ancestry + agg_clim, temp_clim)
      mod_coefs <- summary(mod)$coefficients
    }
    
    # if coxph model
    if (model == 'coxph') {
      mod <- coxph(Surv(yday(dep_date), rep(1, nrow(temp_clim))) ~ transf_ancestry + agg_clim, data = temp_clim)
      mod_coefs <- summary(mod)$coefficients
    }
    
    # if climate IS semi-significant predictor add to list
    if (mod_coefs[nrow(mod_coefs), ncol(mod_coefs)] <= .1) {
      win_opens <- c(win_opens, open)
      win_closes <- c(win_closes, close)
      if (i == 1) {
        bestdat <- temp_clim
        colnames(bestdat)[which(colnames(bestdat) == 'agg_clim')] <- paste0(varname, '_', open, '_', close)
      }
      else {
        bestdat <- merge(bestdat, temp_clim, by = c('ID', 'dep_date', 'release_site', 'dep_year', 'transf_ancestry', 'aims_heterozygosity'), all.x = TRUE)
        colnames(bestdat)[which(colnames(bestdat) == 'agg_clim')] <- paste0(varname, '_', open, '_', close)
      }
    }
  }
  
  # final dataset of all kept windows
  return(list(data.frame(Opens = win_opens, Closes = win_closes), bestdat))
}

# function to extract and plot best data
window_plot <- function(window, varname) {
  # extract best data for plots
  bestdata <- window[[1]]$BestModelData
  bestdata$yvar <- as.Date(bestdata$yvar, origin = '1970-01-01')
  bestdata$yvar |> yday()
  
  plot <- ggplot(bestdata, aes(x = yvar, y = climate, color = transf_ancestry)) + 
    geom_point() + 
    #geom_smooth(method = 'lm', alpha = .25) +
    labs(x = 'Departure Day of Year',
         y = varname,
         title = paste0('Departure vs. ', varname))
  
  return(plot)
}
# WINTERING DATA--------------------------------------------------------------------------------------------------------------------------------

# get wintering data
llg_w <- llg_sub |> 
  select('ID', 'wintering_lat', 'wintering_long', 'spring_dep_new.x', 'land_coords.lat_w', 'land_coords.lon_w', 'release_site', 'aims_ancestry',  'transf_ancestry', 'aims_heterozygosity') |> 
  mutate(dep_year = as.factor(year(spring_dep_new.x)))

# replace recorded coordinates with land coordinates where necessary
llg_w <- llg_w |> 
  mutate(lat_w = if_else(`land_coords.lat_w` == 0 | is.na(`land_coords.lat_w`), wintering_lat, `land_coords.lat_w`)) |> 
  mutate(lon_w = if_else(`land_coords.lon_w` == 0 | is.na(`land_coords.lon_w`), wintering_long, `land_coords.lon_w`))

# add la nina occurrence (no el nino years in our range)
yrs <- sort(unique(llg_w$dep_year))
enso <- c(1, 1, 0, 0, 0, 1, 1)
enso_df <- data.frame(dep_year = yrs, enso = enso)
llg_w <- merge(x = llg_w, y = enso_df, by = c('dep_year'))

# merge with NDVI data so we have release site 
ndvi_merged <- merge(x = llg_w, y = ndvi, by = c('ID'), all.x = FALSE)
ndvi_merged_clean <- ndvi_merged |> 
  select('ID', 'spring_dep_new.x', 'lat_w', 'lon_w', 'date', 'ndvi', 'release_site', 'dep_year', 'transf_ancestry',  'aims_heterozygosity', 'enso')

# split into datasets for sliding window
dep_w <- llg_w |> 
  select('spring_dep_new.x', 'ID', 'release_site', 'dep_year', 'transf_ancestry', 'aims_heterozygosity', 'enso') |> 
  filter(!is.na(spring_dep_new.x))
ndvi_w <- ndvi_merged_clean |> select('date', 'ndvi', 'ID', 'release_site', 'dep_year', 'transf_ancestry', 'aims_heterozygosity', 'enso')

# convert dates
ndvi_w$date <- as.Date(ndvi_w$date)

# save IDs of birds with missing data
missing_table_ndvi <- ndvi_merged_clean |> 
  group_by(ID, spring_dep_new.x) |> 
  summarise(n_missing = sum(is.na(ndvi)))

# merge with climate data so we have release site 
clim_merged <- merge(x = llg_w, y = clim, by = c('ID'), all.x = FALSE)
clim_merged_clean <- clim_merged |> 
  select('ID', 'spring_dep_new.x', 'lat_w', 'lon_w', 'date', 'temp_2m', 'total_precip', 'surf_pressure', 'u_wind', 'v_wind', 'release_site', 'dep_year', 'transf_ancestry', 'aims_heterozygosity', 'enso')

# calculate wind speed
clim_merged_clean <- clim_merged_clean |> mutate(wind_speed = sqrt(u_wind^2 + v_wind^2))

# split into datasets for sliding window
dep_w <- llg_w |> 
  select('spring_dep_new.x', 'ID', 'release_site', 'dep_year', 'transf_ancestry', 'aims_heterozygosity', 'enso') |> 
  filter(!is.na(spring_dep_new.x))
clim_w <- clim_merged_clean |> select('date', 'temp_2m', 'total_precip', 'surf_pressure', 'wind_speed', 'ID', 'release_site', 'dep_year', 'transf_ancestry', 'aims_heterozygosity', 'enso')

# convert dates
clim_w$date <- as.Date(clim_w$date)

# examine and remove NAs
clim_w |> filter(is.na(temp_2m)) |> group_by(ID) |> summarise(n_rows = n())
missing_clim <- (clim_w |> filter(is.na(temp_2m)))$ID |> unique()
clim_w_clean <- clim_w |> filter(!(ID %in% missing_clim))
dep_w_clean2 <- dep_w |> filter(!(ID %in% missing_clim))


# figure out reference dep date using earliest departure
(dep_w |> 
    mutate(new_date = as.Date(paste0(month(spring_dep_new.x), '-', day(spring_dep_new.x), '-2018'), '%m-%d-%Y')))$new_date |> min() 

# figure out latest dep date using latest departure
(dep_w |> 
    mutate(new_date = as.Date(paste0(month(spring_dep_new.x), '-', day(spring_dep_new.x), '-2018'), '%m-%d-%Y')))$new_date |> max() 

# WINTERING PRECIPITATION--------------------------------------------------------------------------------------------------------------

# linear models
# precipitation sliding window using mean as aggregate and daily window
set.seed(123)
precip_w_win_mean <- slidingwin(xvar = list(Precip = clim_w_clean$total_precip),
                                k = 4,
                                cdate = clim_w_clean$date,
                                bdate = dep_w_clean2$spring_dep_new.x,
                                baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry, data = dep_w_clean2), 
                                cinterval = "day",
                                cmissing = 'method1',
                                range = c(184, 0),
                                type = "absolute", refday = c(03, 03),
                                stat = "mean",
                                func = "lin", 
                                spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
precip_w_win_mean
precip_w_mean_bestmod <- precip_w_win_mean[[1]]$BestModel
summary(precip_w_mean_bestmod)
vif(precip_w_mean_bestmod)
window_plot(precip_w_win_mean, 'Precipitation', NULL)

# check assumptions
plot(precip_w_mean_bestmod)
hist(resid(precip_w_mean_bestmod))
shapiro.test(resid(precip_w_mean_bestmod)) # non-normal

# precipitation sliding window using slope as aggregate
set.seed(123)
precip_w_win_sl <- slidingwin(xvar = list(Precip = clim_w_clean$total_precip),
                              k = 4,
                              cdate = clim_w_clean$date,
                              bdate = dep_w_clean2$spring_dep_new.x,
                              baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry, data = dep_w_clean2), 
                              cinterval = "day",
                              cmissing = 'method1',
                              range = c(184, 0),
                              type = "absolute", refday = c(03, 03),
                              stat = "slope",
                              func = "lin", 
                              spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
precip_w_win_sl
precip_w_sl_bestmod <- precip_w_win_sl[[1]]$BestModel
summary(precip_w_sl_bestmod)
vif(precip_w_sl_bestmod)
window_plot(precip_w_win_sl, 'Precipitation', NULL)

# check assumptions
plot(precip_w_sl_bestmod)
hist(resid(precip_w_sl_bestmod))
shapiro.test(resid(precip_w_sl_bestmod))
bptest(yvar ~ transf_ancestry + climate, data = precip_w_win_sl[[1]]$BestModelData)

# get all windows
clim_data <- clim_merged_clean |> filter(!ID %in% missing_clim)
colnames(clim_data)[c(2, 7)] <- c('dep_date', 'climate')
precip_w_sl_allwin <- all_win(precip_w_win_sl, clim_data, '-03-03', 'precip', 'linear', 'slope')

# check for overfitting using random window -> doesn't work bc wants all years of data for all birds
# set.seed(123)
# precip_w_rand_sl <- randwin(repeats = 5,
#                             window = 'sliding',
#                             xvar = list(Precip = clim_w_clean$total_precip),
#                             k = 4,
#                             cdate = clim_w_clean$date,
#                             bdate = dep_w_clean2$spring_dep_new.x,
#                             baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry, data = dep_w_clean2),
#                             cinterval = "day",
#                             cmissing = 'method1',
#                             range = c(184, 0),
#                             type = "absolute", refday = c(03, 03),
#                             stat = "slope",
#                             func = "lin" , 
#                             spatial = list(dep_w_clean2$ID, clim_w_clean$ID))
# 
# # probability of getting our deltaAICc if just random
# pvalue(dataset = precip_w_win_sl[[1]]$Dataset, datasetrand = precip_w_rand_sl[[1]], metric = "C", sample.size = length(yrs))

# Cox proportional hazards models
# using mean as aggregate
precip_w_cph_mean <- slidingwin(xvar = list(Precip = clim_w_clean$total_precip),
                                cdate = clim_w_clean$date,
                                bdate = dep_w_clean2$spring_dep_new.x,
                                baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean2))) ~ transf_ancestry, data = dep_w_clean2), 
                                cinterval = "day",
                                cmissing = 'method1',
                                range = c(184, 0),
                                type = "absolute", refday = c(03, 03),
                                stat = "mean",
                                func = "lin", 
                                spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
precip_w_cph_mean_bestmod <- precip_w_cph_mean[[1]]$BestModel
summary(precip_w_cph_mean_bestmod)
precip_w_cph_mean_survfit <- survfit(precip_w_cph_mean_bestmod)
AIC(precip_w_cph_mean_bestmod)

# Kaplan-Meier curve
plot(precip_w_cph_mean_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")


# get all windows
precip_w_cph_mean_allwin <- all_win(precip_w_cph_mean, clim_data, '-03-03', 'precip', 'cph', 'mean')

# using slope as aggregate
precip_w_cph_sl <- slidingwin(xvar = list(Precip = clim_w_clean$total_precip),
                              cdate = clim_w_clean$date,
                              bdate = dep_w_clean2$spring_dep_new.x,
                              baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean2))) ~ transf_ancestry, data = dep_w_clean2), 
                              cinterval = "day",
                              cmissing = 'method1',
                              range = c(184, 0),
                              type = "absolute", refday = c(03, 03),
                              stat = "slope",
                              func = "lin", 
                              spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
precip_w_cph_sl_bestmod <- precip_w_cph_sl[[1]]$BestModel
summary(precip_w_cph_sl_bestmod)
precip_w_cph_sl_survfit <- survfit(precip_w_cph_sl_bestmod)
AIC(precip_w_cph_sl_bestmod)

# Kaplan-Meier curve
plot(precip_w_cph_sl_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get all windows
precip_w_cph_sl_allwin <- all_win(precip_w_cph_sl, clim_data, '-03-03', 'precip', 'cph', 'slope')

# WINTERING TEMPERATURE------------------------------------------------------------------------------------------------------------------
# linear models
# temperature sliding window using mean as aggregate and daily window
set.seed(123)
temp_w_win_mean <- slidingwin(xvar = list(Temp = clim_w_clean$temp_2m),
                              k = 4,
                              cdate = clim_w_clean$date,
                              bdate = dep_w_clean2$spring_dep_new.x,
                              baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry, data = dep_w_clean2),
                              cinterval = "day",
                              cmissing = 'method1',
                              range = c(184, 0),
                              type = "absolute", refday = c(03, 03),
                              stat = "mean",
                              func = "lin", 
                              spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
temp_w_win_mean
temp_w_mean_bestmod <- temp_w_win_mean[[1]]$BestModel
summary(temp_w_mean_bestmod)
vif(temp_w_mean_bestmod)
window_plot(temp_w_win_mean, 'Temperature', NULL)

# check assumptions
plot(temp_w_mean_bestmod)
hist(resid(temp_w_mean_bestmod))
shapiro.test(resid(temp_w_mean_bestmod)) # non-normal

# temperature sliding window using slope as aggregate
set.seed(123)
temp_w_win_sl <- slidingwin(xvar = list(Temp = clim_w_clean$temp_2m),
                            k = 4,
                            cdate = clim_w_clean$date,
                            bdate = dep_w_clean2$spring_dep_new.x,
                            baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry, data = dep_w_clean2), 
                            cinterval = "day",
                            cmissing = 'method1',
                            range = c(184, 0),
                            type = "absolute", refday = c(03, 03),
                            stat = "slope",
                            func = "lin", 
                            spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
temp_w_win_sl
temp_w_sl_bestmod <- temp_w_win_sl[[1]]$BestModel
summary(temp_w_sl_bestmod)
vif(temp_w_sl_bestmod)
window_plot(temp_w_win_sl, 'Temperature', NULL)

# check assumptions
plot(temp_w_sl_bestmod)
hist(resid(temp_w_sl_bestmod))
shapiro.test(resid(temp_w_sl_bestmod)) # non-normal

# Cox proportional hazards models
# using mean as aggregate
temp_w_cph_mean <- slidingwin(xvar = list(Temp = clim_w_clean$temp_2m),
                              cdate = clim_w_clean$date,
                              bdate = dep_w_clean2$spring_dep_new.x,
                              baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean2))) ~ transf_ancestry, data = dep_w_clean2), 
                              cinterval = "day",
                              cmissing = 'method1',
                              range = c(184, 0),
                              type = "absolute", refday = c(03, 03),
                              stat = "mean",
                              func = "lin", 
                              spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

temp_w_cph_mean_bestmod <- temp_w_cph_mean[[1]]$BestModel
summary(temp_w_cph_mean_bestmod)
temp_w_cph_mean_survfit <- survfit(temp_w_cph_mean_bestmod)
AIC(temp_w_cph_mean_bestmod)

# Kaplan-Meier curve
plot(temp_w_cph_mean_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get all windows
clim_data <- clim_merged_clean |> filter(!ID %in% missing_clim)
colnames(clim_data)[c(2, 6)] <- c('dep_date', 'climate')
temp_w_cph_mean_allwin <- all_win(temp_w_cph_mean, clim_data, '-03-03', 'temp', 'cph', 'mean')

# using slope as aggregate
temp_w_cph_sl <- slidingwin(xvar = list(Temp = clim_w_clean$temp_2m),
                            cdate = clim_w_clean$date,
                            bdate = dep_w_clean2$spring_dep_new.x,
                            baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean2))) ~ transf_ancestry, data = dep_w_clean2), 
                            cinterval = "day",
                            cmissing = 'method1',
                            range = c(184, 0),
                            type = "absolute", refday = c(03, 03),
                            stat = "slope",
                            func = "lin", 
                            spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

temp_w_cph_sl_bestmod <- temp_w_cph_sl[[1]]$BestModel
summary(temp_w_cph_sl_bestmod)
temp_w_cph_sl_survfit <- survfit(temp_w_cph_sl_bestmod)
AIC(temp_w_cph_sl_bestmod)

# Kaplan-Meier curve
plot(temp_w_cph_sl_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get all windows
temp_w_cph_sl_allwin <- all_win(temp_w_cph_sl, clim_data, '-03-03', 'temp', 'cph', 'slope')

# WINTERING SURFACE/ATMOSPHERIC PRESSURE--------------------------------------------------------------------------------------------------------------

# linear models
# precipitation sliding window using mean as aggregate and daily window
set.seed(123)
press_w_win_mean <- slidingwin(xvar = list(Pressure = clim_w_clean$surf_pressure),
                               k = 4,
                               cdate = clim_w_clean$date,
                               bdate = dep_w_clean2$spring_dep_new.x,
                               baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry, data = dep_w_clean2), 
                               cinterval = "day",
                               cmissing = 'method1',
                               range = c(184, 0),
                               type = "absolute", refday = c(03, 03),
                               stat = "mean",
                               func = "lin", 
                               spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
press_w_win_mean
press_w_win_mean_bestmod <- press_w_win_mean[[1]]$BestModel
summary(press_w_win_mean_bestmod)
vif(press_w_win_mean_bestmod)
window_plot(press_w_win_mean, 'Surface Pressure', NULL)

# check assumptions
plot(press_w_win_mean_bestmod)
hist(resid(press_w_win_mean_bestmod))
shapiro.test(resid(press_w_win_mean_bestmod)) # non-normal

# precipitation sliding window using slope as aggregate
set.seed(123)
press_w_win_sl <- slidingwin(xvar = list(Pressure = clim_w_clean$surf_pressure),
                             k = 4,
                             cdate = clim_w_clean$date,
                             bdate = dep_w_clean2$spring_dep_new.x,
                             baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry, data = dep_w_clean2), 
                             cinterval = "day",
                             cmissing = 'method1',
                             range = c(184, 0),
                             type = "absolute", refday = c(03, 03),
                             stat = "slope",
                             func = "lin", 
                             spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
press_w_win_sl
press_w_sl_bestmod <- press_w_win_sl[[1]]$BestModel
summary(press_w_sl_bestmod)
vif(press_w_sl_bestmod)
window_plot(press_w_win_sl, 'Surface Pressure', NULL)

# check assumptions
plot(press_w_sl_bestmod)
hist(resid(press_w_sl_bestmod))
shapiro.test(resid(press_w_sl_bestmod)) # non-normal

# Cox proportional hazards models
# using mean as aggregate
press_w_cph_mean <- slidingwin(xvar = list(Pressure = clim_w_clean$surf_pressure),
                               cdate = clim_w_clean$date,
                               bdate = dep_w_clean2$spring_dep_new.x,
                               baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean2))) ~ transf_ancestry, data = dep_w_clean2), 
                               cinterval = "day",
                               cmissing = 'method1',
                               range = c(184, 0),
                               type = "absolute", refday = c(03, 03),
                               stat = "mean",
                               func = "lin", 
                               spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
press_w_cph_mean
press_w_cph_mean_bestmod <- press_w_cph_mean[[1]]$BestModel
summary(press_w_cph_mean_bestmod)
press_w_cph_mean_survfit <- survfit(press_w_cph_mean_bestmod)
AIC(press_w_cph_mean_bestmod)

# Kaplan-Meier curve
plot(press_w_cph_mean_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get all windows
press_w_cph_mean_allwin <- all_win(press_w_cph_mean, clim_data, '-03-03', 'press', 'cph', 'mean')

# using slope as aggregate
press_w_cph_sl <- slidingwin(xvar = list(Pressure = clim_w_clean$surf_pressure),
                             cdate = clim_w_clean$date,
                             bdate = dep_w_clean2$spring_dep_new.x,
                             baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean2))) ~ transf_ancestry, data = dep_w_clean2), 
                             cinterval = "day",
                             cmissing = 'method1',
                             range = c(184, 0),
                             type = "absolute", refday = c(03, 03),
                             stat = "slope",
                             func = "lin", 
                             spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
press_w_cph_sl
press_w_cph_sl_bestmod <- press_w_cph_sl[[1]]$BestModel
summary(press_w_cph_sl_bestmod)
press_w_cph_sl_survfit <- survfit(press_w_cph_sl_bestmod)
AIC(press_w_cph_sl_bestmod)

# Kaplan-Meier curve
plot(press_w_cph_sl_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get all windows
press_w_cph_sl_allwin <- all_win(press_w_cph_sl, clim_data, '-03-03', 'press', 'cph', 'slope')

# WINTERING NDVI------------------------------------------------------------------------------------------------------------------
# linear models
# NDVI sliding window using mean as aggregate and daily window
set.seed(123)
ndvi_w_win_mean <- slidingwin(xvar = list(NDVI = ndvi_w$ndvi),
                              k = 4,
                              cdate = ndvi_w$date,
                              bdate = dep_w$spring_dep_new.x,
                              baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry, data = dep_w), 
                              cinterval = "day",
                              cmissing = 'method1',
                              range = c(184, 0),
                              type = "absolute", refday = c(03, 03),
                              stat = "mean",
                              func = "lin",
                              spatial = list(dep_w$ID, ndvi_w$ID))

# examine model
ndvi_w_mean_bestmod <- ndvi_w_win_mean[[1]]$BestModel
summary(ndvi_w_mean_bestmod)
vif(ndvi_w_mean_bestmod)
window_plot(ndvi_w_win_mean, 'NDVI', NULL)

# check assumptions
plot(ndvi_w_mean_bestmod)
hist(resid(ndvi_w_mean_bestmod))
shapiro.test(resid(ndvi_w_mean_bestmod)) # non-normal

# temperature sliding window using slope as aggregate
set.seed(123)
ndvi_w_win_sl <- slidingwin(xvar = list(NDVI = ndvi_w$ndvi),
                            k = 4,
                            cdate = ndvi_w$date,
                            bdate = dep_w$spring_dep_new.x,
                            baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry, data = dep_w), 
                            cinterval = "day",
                            cmissing = 'method1',
                            range = c(184, 0),
                            type = "absolute", refday = c(03, 03),
                            stat = "slope",
                            func = "lin", 
                            spatial = list(dep_w$ID, ndvi_w$ID))

# examine model
ndvi_w_sl_bestmod <- ndvi_w_win_sl[[1]]$BestModel
summary(ndvi_w_sl_bestmod)
vif(ndvi_w_sl_bestmod)
window_plot(ndvi_w_win_sl, 'NDVI', NULL)

# check assumptions
plot(ndvi_w_sl_bestmod)
hist(resid(ndvi_w_sl_bestmod))
shapiro.test(resid(ndvi_w_sl_bestmod)) # non-normal

# Cox proportional hazards models
# using mean as aggregate
ndvi_w_cph_mean <- slidingwin(xvar = list(NDVI = ndvi_w$ndvi),
                              cdate = ndvi_w$date,
                              bdate = dep_w$spring_dep_new.x,
                              baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w))) ~ transf_ancestry, data = dep_w),
                              cinterval = "day",
                              cmissing = 'method1',
                              range = c(184, 0),
                              type = "absolute", refday = c(03, 03),
                              stat = "mean",
                              func = "lin", 
                              spatial = list(dep_w$ID, ndvi_w$ID))

# examine model
ndvi_w_cph_mean
ndvi_w_cph_mean_bestmod <- ndvi_w_cph_mean[[1]]$BestModel
summary(ndvi_w_cph_mean_bestmod)
ndvi_w_cph_mean_survfit <- survfit(ndvi_w_cph_mean_bestmod)
AIC(ndvi_w_cph_mean_bestmod)

# Kaplan-Meier curve
plot(ndvi_w_cph_mean_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get all windows
ndvi_w_cph_mean_allwin <- all_win(ndvi_w_cph_mean, clim_data, '-03-03', 'ndvi', 'cph', 'mean')

# using slope as aggregate
ndvi_w_cph_sl <- slidingwin(xvar = list(NDVI = ndvi_w$ndvi),
                            cdate = ndvi_w$date,
                            bdate = dep_w$spring_dep_new.x,
                            baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w))) ~ transf_ancestry, data = dep_w), # I'm confused on how to incorporate release site as a covariate...
                            cinterval = "day",
                            cmissing = 'method1',
                            range = c(184, 0),
                            type = "absolute", refday = c(03, 03),
                            stat = "slope",
                            func = "lin", 
                            spatial = list(dep_w$ID, ndvi_w$ID))

# examine model
ndvi_w_cph_sl
ndvi_w_cph_sl_bestmod <- ndvi_w_cph_sl[[1]]$BestModel
summary(ndvi_w_cph_sl_bestmod)
ndvi_w_cph_sl_survfit <- survfit(ndvi_w_cph_sl_bestmod)
AIC(ndvi_w_cph_sl_bestmod)

# Kaplan-Meier curve
plot(ndvi_w_cph_sl_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get all windows
ndvi_w_cph_sl_allwin <- all_win(ndvi_w_cph_sl, clim_data, '-03-03', 'ndvi', 'cph', 'slope')

# WINTERING DAYLENGTH-----------------------------------------------------------
# get rid of any rows missing lat or date
llg_w_clean <- llg_w |> 
  filter((!is.na(spring_dep_new.x)) & !is.na(lat_w)) |> 
  arrange(ID)

# figure out latest day of year departure
max_dt <- (llg_w_clean |> mutate(doy = yday(spring_dep_new.x)) |> arrange(-doy))[1, ]$spring_dep_new.x

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
daylen_w <- llg_dts_w |> select('ID', 'release_site', 'date', 'daylen', 'dep_year', 'transf_ancestry', 'aims_heterozygosity')
dep_daylen_w <- llg_w_clean |> select('ID', 'release_site', 'spring_dep_new.x', 'dep_year', 'transf_ancestry', 'aims_heterozygosity')

# linear models
# daylength sliding window using mean as aggregate
set.seed(123)
daylen_w_win_mean <- slidingwin(xvar = list(Daylength = daylen_w$daylen),
                                k = 4,
                                cdate = daylen_w$date,
                                bdate = dep_daylen_w$spring_dep_new.x,
                                baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry, data = dep_daylen_w),
                                cinterval = "day",
                                cmissing = 'method2',
                                range = c(184, 0),
                                type = "absolute", refday = c(03, 03),
                                stat = "mean",
                                func = "lin", 
                                spatial = list(dep_daylen_w$ID, daylen_w$ID))

# examine model
daylen_w_mean_bestmod <- daylen_w_win_mean[[1]]$BestModel
summary(daylen_w_mean_bestmod)
vif(daylen_w_mean_bestmod)
window_plot(daylen_w_win_mean, 'Daylength', NULL)

# check assumptions
plot(daylen_w_mean_bestmod)
hist(resid(daylen_w_mean_bestmod))
shapiro.test(resid(daylen_w_mean_bestmod)) # non-normal

# using slope as aggregate
set.seed(123)
daylen_w_win_sl <- slidingwin(xvar = list(Daylength = daylen_w$daylen),
                              k = 4,
                              cdate = daylen_w$date,
                              bdate = dep_daylen_w$spring_dep_new.x,
                              baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry, data = dep_daylen_w),
                              cinterval = "day",
                              cmissing = 'method2',
                              range = c(184, 0),
                              type = "absolute", refday = c(03, 03),
                              stat = "slope",
                              func = "lin", 
                              spatial = list(dep_daylen_w$ID, daylen_w$ID))

# examine model
daylen_w_sl_bestmod <- daylen_w_win_sl[[1]]$BestModel
summary(daylen_w_sl_bestmod)
vif(daylen_w_sl_bestmod)
window_plot(daylen_w_win_sl, 'Daylength', NULL)

# check assumptions
plot(daylen_w_sl_bestmod)
hist(resid(daylen_w_sl_bestmod))
shapiro.test(resid(daylen_w_sl_bestmod)) # non-normal

# Cox proportional hazards models
# using mean as aggregate
daylen_w_cph_mean <- slidingwin(xvar = list(Daylength = daylen_w$daylen),
                                cdate = daylen_w$date,
                                bdate = dep_daylen_w$spring_dep_new.x,
                                baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w))) ~ transf_ancestry, data = dep_daylen_w),
                                cinterval = "day",
                                cmissing = 'method1',
                                range = c(184, 0),
                                type = "absolute", refday = c(03, 03),
                                stat = "mean",
                                func = "lin", 
                                spatial = list(dep_daylen_w$ID, daylen_w$ID))

# examine model
daylen_w_cph_mean
daylen_w_cph_mean_bestmod <- daylen_w_cph_mean[[1]]$BestModel
summary(daylen_w_cph_mean_bestmod)
daylen_w_cph_mean_survfit <- survfit(daylen_w_cph_mean_bestmod)
AIC(daylen_w_cph_mean_bestmod)

# Kaplan-Meier curve
plot(daylen_w_cph_mean_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get all windows
daylen_w_cph_mean_allwin <- all_win(daylen_w_cph_mean, clim_data, '-03-03', 'daylen', 'cph', 'mean')

# using slope as aggregate
daylen_w_cph_sl <- slidingwin(xvar = list(Daylength = daylen_w$daylen),
                              cdate = daylen_w$date,
                              bdate = dep_daylen_w$spring_dep_new.x,
                              baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w))) ~ transf_ancestry, data = dep_daylen_w),
                              cinterval = "day",
                              cmissing = 'method1',
                              range = c(184, 0),
                              type = "absolute", refday = c(03, 03),
                              stat = "slope",
                              func = "lin", 
                              spatial = list(dep_daylen_w$ID, daylen_w$ID))

# examine model
daylen_w_cph_sl
daylen_w_cph_sl_bestmod <- daylen_w_cph_sl[[1]]$BestModel
summary(daylen_w_cph_sl_bestmod)
daylen_w_cph_sl_survfit <- survfit(daylen_w_cph_sl_bestmod)
AIC(daylen_w_cph_sl_bestmod)

# Kaplan-Meier curve
plot(daylen_w_cph_sl_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get all windows
daylen_w_cph_sl_allwin <- all_win(daylen_w_cph_sl, clim_data, '-03-03', 'daylen', 'cph', 'slope')