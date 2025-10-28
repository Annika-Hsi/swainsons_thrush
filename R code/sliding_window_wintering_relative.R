# ROUGH DRAFT VERSION: using relative window and testing linear vs. quadratic, 
# weekly vs. daily, slope vs. mean
# PURPOSE: run climwin sliding window on environmental data from Google Earth 
# Engine corresponding to Swainson's Thrush observations from the wintering 
# season only
# NOTES: have to remerge with original LLG dataset to get release site as covariate;
# ENSO data not currently being used, but added in for future use; can toggle
# 5-day window threshold

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
aqua <- read.csv('env_data/gee/NDVI_wintering.csv')
terra <- read.csv('env_data/gee/NDVI_wintering_terra.csv')
clim <- read.csv('env_data/gee/clim_wintering.csv')
llg_orig <- read.csv('data/pheno_archival_raw_steph.csv')
llg_ancestry <- read.csv('data/AIMs_metadata.20250327.csv')

# remove NAs from aqua data (terra already clean)
aqua <- aqua |> filter(!is.na(ndvi))

# merge and average MODIS data
ndvi <- merge(x = aqua, y = terra, by = c('ID', 'date', 'lat', 'lon'), all = TRUE)
ndvi <- ndvi |> 
  mutate(num_sat = (as.numeric(!is.na(ndvi.x) & !is.na(ndvi.y)) + 1))
ndvi$ndvi.x[which(is.na(ndvi$ndvi.x))] <- 0
ndvi$ndvi.y[which(is.na(ndvi$ndvi.y))] <- 0
ndvi <- ndvi |> mutate(ndvi = ((ndvi.x + ndvi.y)/num_sat))
ndvi <- ndvi |> select(all_of(c('ID', 'date', 'lat', 'lon', 'ndvi')))

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

# function that gets top window for each of the methods tested
top_wins <- function(win_obj, agg) {
  df <- as.data.frame(matrix(ncol = 7))
  combos <- win_obj$combos
  colnames(df) <- colnames(combos[,-1])
  if (grepl('Surv', combos$response[1]) == TRUE) {
    mod_type <- 'Coxph'
  }
  else {
    mod_type <- 'linear'
  }
  df <- combos[, -1]
  df$model <- mod_type
  df$season <- 'wintering/spring migration'
  df$aggregate <- agg
  
  r2 <- c()
  p_clim <- c()
  p_clim2 <- c()
  coef_clim <- c()
  coef_clim2 <- c()
  # get R^2 and p-values
  for (i in 1:4) {
    bestmod <- win_obj[[i]]$BestModel
    bestsum <- summary(bestmod)
    coefs <- as.data.frame(bestsum$coefficients)
    if (mod_type == 'linear') {
      r2[i] <- bestsum$adj.r.squared
      clim_coefs <- coefs[rownames(coefs) %in% c('climate', 'I(climate^2)'), ]
      p_clim[i] <- clim_coefs[1, ]$`Pr(>|t|)`
      coef_clim[i] <- clim_coefs[1, ]$Estimate
      if (combos[i, ]$func == 'quad') {
        p_clim2[i] <- clim_coefs[2, ]$`Pr(>|t|)`
        coef_clim2[i] <- clim_coefs[2, ]$Estimate
      }
      else {
        p_clim2[i] <- NULL
        coef_clim2[i] <- NULL
      }
    }
    if (mod_type == 'Coxph') {
      r2[i] <- bestsum$rsq[1]
      clim_coefs <- coefs[rownames(coefs) %in% c('climate', 'I(climate^2)'), ]
      p_clim[i] <- clim_coefs[1, ]$`Pr(>|z|)`
      coef_clim[i] <- clim_coefs[1, ]$coef
      if (combos[i, ]$func == 'quad') {
        p_clim2[i] <- clim_coefs[2, ]$`Pr(>|z|)`
        coef_clim2[i] <- clim_coefs[2, ]$coef
      }
      else {
        p_clim2[i] <- NULL
        coef_clim2[i] <- NULL
      }
    }
  } 
  df$r2 <- r2
  df$clim_pvalue <- p_clim
  df$clim2_pvalue <- p_clim2
  df$clim_coef <- coef_clim
  df$clim2_coefs <- coef_clim2
  return(df)
}

# function to extract and plot best data
window_plot <- function(window, m, varname) {
  # extract best data for plots
  bestdata <- window[[m]]$BestModelData
  bestmod <- window[[m]]$BestModel
  preds <- predict(bestmod)
  df <- data.frame(yvar = preds, climate = bestdata$climate, type = 'pred')
  bestdata$type <- 'true'
  df <- rbind(df, bestdata |> select(all_of(c('yvar', 'climate', 'type'))))
  plot <- ggplot(df, aes(x = climate, y = yvar, color = type)) +
    geom_point() +
    labs(y = 'Departure Day of Year',
         x = varname,
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

# figure out average dep date
(dep_w |> 
    mutate(new_date = as.Date(paste0(month(spring_dep_new.x), '-', day(spring_dep_new.x), '-2018'), '%m-%d-%Y')))$new_date |> mean() 

# WINTERING PRECIPITATION--------------------------------------------------------------------------------------------------------------

# linear models
# precipitation sliding window using mean as aggregate and daily window

# temporary fix to 2019 being "new" during cross validation
dep_w_clean2 <- dep_w_clean2 |> filter(dep_year != 2019)
clim_w_clean <- clim_w_clean |> filter(dep_year != 2019)

# linear models
# precipitation sliding window using linear model and daily window
set.seed(123)
precip_w_lm_d <- slidingwin(xvar = list(Precip = clim_w_clean$total_precip),
                            k = 4,
                            cdate = clim_w_clean$date,
                            bdate = dep_w_clean2$spring_dep_new.x,
                            baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry + dep_year, data = dep_w_clean2), 
                            cinterval = "day",
                            cmissing = 'method1',
                            range = c(92, 0),
                            type = "relative",
                            stat = c("mean", "slope"),
                            func = c("lin", "quad"), 
                            spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
precip_w_lm_d
m <- which(precip_w_lm_d$combos$DeltaAICc == min(precip_w_lm_d$combos$DeltaAICc))
precip_w_lm_d_bestmod <- precip_w_lm_d[[m]]$BestModel
summary(precip_w_lm_d_bestmod)
vif(precip_w_lm_d_bestmod)
AIC(precip_w_lm_d_bestmod)
window_plot(precip_w_lm_d, m, 'Precipitation')

# check assumptions
plot(precip_w_lm_d_bestmod)
hist(resid(precip_w_lm_d_bestmod))
shapiro.test(resid(precip_w_lm_d_bestmod)) # not normal

# weekly linear model
set.seed(123)
precip_w_lm_w <- slidingwin(xvar = list(Precip = clim_w_clean$total_precip),
                            k = 4,
                            cdate = clim_w_clean$date,
                            bdate = dep_w_clean2$spring_dep_new.x,
                            baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry + dep_year, data = dep_w_clean2), 
                            cinterval = "week",
                            cmissing = 'method1',
                            range = c(13, 0),
                            type = "relative",
                            stat = c("mean", "slope"),
                            func = c("lin", "quad"), 
                            spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
precip_w_lm_w
m <- which(precip_w_lm_w$combos$DeltaAICc == min(precip_w_lm_w$combos$DeltaAICc))
precip_w_lm_w_bestmod <- precip_w_lm_w[[m]]$BestModel
summary(precip_w_lm_w_bestmod)
vif(precip_w_lm_w_bestmod)

# best model overfit, try simpler model
m <- m - 2
precip_w_lm_w_bestmod <- precip_w_lm_w[[m]]$BestModel
summary(precip_w_lm_w_bestmod)
vif(precip_w_lm_w_bestmod)
window_plot(precip_w_lm_w, 2, 'Precipitation')

# check assumptions
plot(precip_w_lm_w_bestmod)
hist(resid(precip_w_lm_w_bestmod))
shapiro.test(resid(precip_w_lm_w_bestmod)) # not normal

# daily Coxph
precip_w_cph_d <- slidingwin(xvar = list(Precip = clim_w_clean$total_precip),
                             cdate = clim_w_clean$date,
                             bdate = dep_w_clean2$spring_dep_new.x,
                             baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean2))) ~ transf_ancestry + dep_year, data = dep_w_clean2), 
                             cinterval = "day",
                             cmissing = 'method1',
                             range = c(92, 0),
                             type = "relative",
                             stat = c("mean", "slope"),
                             func = c("lin", "quad"), 
                             spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
precip_w_cph_d
m <- which(precip_w_cph_d$combos$DeltaAICc == min(precip_w_cph_d$combos$DeltaAICc))
precip_w_cph_d_bestmod <- precip_w_cph_d[[m]]$BestModelData
summary(precip_w_cph_d_bestmod)
summary(precip_w_cph_d_bestmod)$rsq
precip_w_cph_d_survfit <- survfit(precip_w_cph_d_bestmod)
AIC(precip_w_cph_d_bestmod)

# Kaplan-Meier curve
plot(precip_w_cph_d_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get top windows
top_windows <- top_wins(precip_w_cph_d, 'daily')

# weekly Coxph
precip_w_cph_w <- slidingwin(xvar = list(Precip = clim_w_clean$total_precip),
                             cdate = clim_w_clean$date,
                             bdate = dep_w_clean2$spring_dep_new.x,
                             baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean2))) ~ transf_ancestry + dep_year, data = dep_w_clean2), 
                             cinterval = "week",
                             cmissing = 'method1',
                             range = c(13, 0),
                             type = "relative",
                             stat = c("mean", "slope"),
                             func = c("lin", "quad"), 
                             spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
precip_w_cph_w
m <- which(precip_w_cph_w$combos$DeltaAICc == min(precip_w_cph_w$combos$DeltaAICc))
m <- m[1] # multiple models with same AIC -> pick simpler one
precip_w_cph_w_bestmod <- precip_w_cph_w[[m]]$BestModel
summary(precip_w_cph_w_bestmod)
summary(precip_w_cph_w_bestmod)$rsq
precip_w_cph_w_survfit <- survfit(precip_w_cph_w_bestmod)
AIC(precip_w_cph_w_bestmod)

# Kaplan-Meier curve
plot(precip_w_cph_w_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get top windows
top_windows <- rbind(top_windows, top_wins(precip_w_cph_w, 'weekly'))

# WINTERING TEMPERATURE------------------------------------------------------------------------------------------------------------------
# linear models
# temperature sliding window using linear model and daily window
set.seed(123)
temp_w_lm_d <- slidingwin(xvar = list(Temp = clim_w_clean$temp_2m),
                          k = 4,
                          cdate = clim_w_clean$date,
                          bdate = dep_w_clean2$spring_dep_new.x,
                          baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry + dep_year, data = dep_w_clean2),
                          cinterval = "day",
                          cmissing = 'method1',
                          range = c(92, 0),
                          type = "relative",
                          stat = c("mean", "slope"),
                          func = c("lin", "quad"), 
                          spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
temp_w_lm_d
m <- which(temp_w_lm_d$combos$DeltaAICc == min(temp_w_lm_d$combos$DeltaAICc))
temp_w_lm_d_bestmod <- temp_w_lm_d[[m]]$BestModel
summary(temp_w_lm_d_bestmod)
vif(temp_w_lm_d_bestmod) 

# best model overfit, try simpler model
m <- m - 2
temp_w_lm_d_bestmod <- temp_w_lm_d[[m]]$BestModel
summary(temp_w_lm_d_bestmod)
vif(temp_w_lm_d_bestmod) 
window_plot(temp_w_lm_d, m, 'Temperature')

# check assumptions
plot(temp_w_lm_d_bestmod)
hist(resid(temp_w_lm_d_bestmod))
shapiro.test(resid(temp_w_lm_d_bestmod))
bptest(yvar ~ ., data = temp_w_lm_d[[m]]$BestModelData)

# get top windows
top_windows <- rbind(top_windows, top_wins(temp_w_lm_d, 'daily'))

# weekly linear model
set.seed(123)
temp_w_lm_w <- slidingwin(xvar = list(Temp = clim_w_clean$temp_2m),
                          k = 4,
                          cdate = clim_w_clean$date,
                          bdate = dep_w_clean2$spring_dep_new.x,
                          baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry + dep_year, data = dep_w_clean2), 
                          cinterval = "week",
                          cmissing = 'method1',
                          range = c(13, 0),
                          type = "relative",
                          stat = c("mean", "slope"),
                          func = c("lin", "quad"), 
                          spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
temp_w_lm_w
m <- which(temp_w_lm_w$combos$DeltaAICc == min(temp_w_lm_w$combos$DeltaAICc))
temp_w_lm_w_bestmod <- temp_w_lm_w[[m]]$BestModel
summary(temp_w_lm_w_bestmod)
vif(temp_w_lm_w_bestmod) 

# best model overfit, try simpler model
m <- m - 2
temp_w_lm_w_bestmod <- temp_w_lm_w[[m]]$BestModel
summary(temp_w_lm_w_bestmod)
vif(temp_w_lm_w_bestmod) 
window_plot(temp_w_lm_w, m, 'Temperature')

# check assumptions
plot(temp_w_lm_w_bestmod)
hist(resid(temp_w_lm_w_bestmod))
shapiro.test(resid(temp_w_lm_w_bestmod))
bptest(yvar ~ ., data = temp_w_lm_w[[m]]$BestModelData)

# get top windows
top_windows <- rbind(top_windows, top_wins(temp_w_lm_w, 'weekly'))

# daily Coxph
temp_w_cph_d <- slidingwin(xvar = list(Temp = clim_w_clean$temp_2m),
                           cdate = clim_w_clean$date,
                           bdate = dep_w_clean2$spring_dep_new.x,
                           baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean2))) ~ transf_ancestry + dep_year, data = dep_w_clean2), 
                           cinterval = "day",
                           cmissing = 'method1',
                           range = c(92, 0),
                           type = "relative",
                           stat = c("mean", "slope"),
                           func = c("lin", "quad"), 
                           spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
temp_w_cph_d
m <- which(temp_w_cph_d$combos$DeltaAICc == min(temp_w_cph_d$combos$DeltaAICc))
temp_w_cph_d_bestmod <- temp_w_cph_d[[m]]$BestModel
summary(temp_w_cph_d_bestmod)
summary(temp_w_cph_d_bestmod)$rsq
temp_w_cph_d_survfit <- survfit(temp_w_cph_d_bestmod)
AIC(temp_w_cph_d_bestmod)

# Kaplan-Meier curve
plot(temp_w_cph_d_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

top_windows <- rbind(top_windows, top_wins(temp_w_cph_d, 'daily'))

# weekly Coxph
temp_w_cph_w <- slidingwin(xvar = list(Temp = clim_w_clean$temp_2m),
                           cdate = clim_w_clean$date,
                           bdate = dep_w_clean2$spring_dep_new.x,
                           baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean2))) ~ transf_ancestry + dep_year, data = dep_w_clean2), 
                           cinterval = "week",
                           cmissing = 'method1',
                           range = c(13, 0),
                           type = "relative",
                           stat = c("mean", "slope"),
                           func = c("lin", "quad"), 
                           spatial = list(dep_w_clean2$ID, clim_w_clean$ID))


# examine model
temp_w_cph_w
m <- which(temp_w_cph_w$combos$DeltaAICc == min(temp_w_cph_w$combos$DeltaAICc))
temp_w_cph_w_bestmod <- temp_w_cph_w[[m]]$BestModel
summary(temp_w_cph_w_bestmod)
summary(temp_w_cph_w_bestmod)$rsq
temp_w_cph_w_survfit <- survfit(temp_w_cph_w_bestmod)
AIC(temp_w_cph_w_bestmod)

# Kaplan-Meier curve
plot(temp_w_cph_w_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# WINTERING WIND SPEED------------------------------------------------------------------------------------------------------------------
# linear models
# wind speed sliding window using linear model and daily window
set.seed(123)
wind_w_lm_d <- slidingwin(xvar = list(Wind = clim_w_clean$wind_speed),
                          k = 4,
                          cdate = clim_w_clean$date,
                          bdate = dep_w_clean2$spring_dep_new.x,
                          baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry + dep_year, data = dep_w_clean2),
                          cinterval = "day",
                          cmissing = 'method1',
                          range = c(92, 0),
                          type = "relative",
                          stat = c("mean", "slope"),
                          func = c("lin", "quad"), 
                          spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
wind_w_lm_d
m <- which(wind_w_lm_d$combos$DeltaAICc == min(wind_w_lm_d$combos$DeltaAICc))
wind_w_lm_d_bestmod <- wind_w_lm_d[[m]]$BestModel
summary(wind_w_lm_d_bestmod) # climate not significant predictor
vif(wind_w_lm_d_bestmod)
window_plot(wind_w_lm_d, m, 'Wind Speed')

# check assumptions
plot(wind_w_lm_d_bestmod)
hist(resid(wind_w_lm_d_bestmod))
shapiro.test(resid(wind_w_lm_d_bestmod)) # not normal

# weekly linear model
set.seed(123)
wind_w_lm_w <- slidingwin(xvar = list(Wind = clim_w_clean$wind_speed),
                          k = 4,
                          cdate = clim_w_clean$date,
                          bdate = dep_w_clean2$spring_dep_new.x,
                          baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry + dep_year, data = dep_w_clean2), 
                          cinterval = "week",
                          cmissing = 'method1',
                          range = c(13, 0),
                          type = "relative",
                          stat = c("mean", "slope"),
                          func = c("lin", "quad"),
                          spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
wind_w_lm_w
m <- which(wind_w_lm_w$combos$DeltaAICc == min(wind_w_lm_w$combos$DeltaAICc))
wind_w_lm_w_bestmod <- wind_w_lm_w[[m]]$BestModel
summary(wind_w_lm_w_bestmod) # climate not significant predictor
vif(wind_w_lm_w_bestmod)
window_plot(wind_w_lm_w, m, 'Wind Speed')

# check assumptions
plot(wind_w_lm_w_bestmod)
hist(resid(wind_w_lm_w_bestmod))
shapiro.test(resid(wind_w_lm_w_bestmod))
bptest(yvar ~ transf_ancestry + dep_year + climate, data = wind_w_lm_w[[m]]$BestModelData)

# get top windows
top_windows <- rbind(top_windows, top_wins(wind_w_lm_w, 'weekly'))

# daily Coxph
wind_w_cph_d <- slidingwin(xvar = list(Wind = clim_w_clean$wind_speed),
                           cdate = clim_w_clean$date,
                           bdate = dep_w_clean2$spring_dep_new.x,
                           baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean2))) ~ transf_ancestry + dep_year, data = dep_w_clean2), 
                           cinterval = "day",
                           cmissing = 'method1',
                           range = c(92, 0),
                           type = "relative",
                           stat = c("mean", "slope"),
                           func = c("lin", "quad"), 
                           spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

wind_w_cph_d
m <- which(wind_w_cph_d$combos$DeltaAICc == min(wind_w_cph_d$combos$DeltaAICc))
wind_w_cph_d_bestmod <- wind_w_cph_d[[m]]$BestModel
summary(wind_w_cph_d_bestmod)
summary(wind_w_cph_d_bestmod)$rsq
wind_w_cph_d_survfit <- survfit(wind_w_cph_d_bestmod)
AIC(wind_w_cph_d_bestmod)

# Kaplan-Meier curve
plot(wind_w_cph_d_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get top windows
top_windows <- rbind(top_windows, top_wins(wind_w_cph_d, 'daily'))

# weekly Coxph
wind_w_cph_w <- slidingwin(xvar = list(Wind = clim_w_clean$wind_speed),
                           cdate = clim_w_clean$date,
                           bdate = dep_w_clean2$spring_dep_new.x,
                           baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean2))) ~ transf_ancestry + dep_year, data = dep_w_clean2), 
                           cinterval = "week",
                           cmissing = 'method1',
                           range = c(13, 0),
                           type = "relative",
                           stat = c("mean", "slope"),
                           func = c("lin", "quad"), 
                           spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

wind_w_cph_w
m <- which(wind_w_cph_w$combos$DeltaAICc == min(wind_w_cph_w$combos$DeltaAICc))
wind_w_cph_w_bestmod <- wind_w_cph_w[[m]]$BestModel
summary(wind_w_cph_w_bestmod)
summary(wind_w_cph_w_bestmod)$rsq
wind_w_cph_w_survfit <- survfit(wind_w_cph_w_bestmod)
AIC(wind_w_cph_w_bestmod)

# Kaplan-Meier curve
plot(wind_w_cph_w_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get top windows
top_windows <- rbind(top_windows, top_wins(wind_w_cph_w, 'weekly'))

# WINTERING SURFACE/ATMOSPHERIC PRESSURE--------------------------------------------------------------------------------------------------------------

# linear models
# surface pressure linear model using daily window
set.seed(123)
press_w_lm_d <- slidingwin(xvar = list(Pressure = clim_w_clean$surf_pressure),
                           k = 4,
                           cdate = clim_w_clean$date,
                           bdate = dep_w_clean2$spring_dep_new.x,
                           baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry + dep_year, data = dep_w_clean2), 
                           cinterval = "day",
                           cmissing = 'method1',
                           range = c(92, 0),
                           type = "relative",
                           stat = c("mean", "slope"),
                           func = c("lin", "quad"),
                           spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
press_w_lm_d
m <- which(press_w_lm_d$combos$DeltaAICc == min(press_w_lm_d$combos$DeltaAICc))
press_w_lm_d_bestmod <- press_w_lm_d[[m]]$BestModel
summary(press_w_lm_d_bestmod)
vif(press_w_lm_d_bestmod)
window_plot(press_w_lm_d, m, 'Surface Pressure')

# check assumptions
plot(press_w_lm_d_bestmod)
hist(resid(press_w_lm_d_bestmod))
shapiro.test(resid(press_w_lm_d_bestmod)) # not normal

# weekly linear model
set.seed(123)
press_w_lm_w <- slidingwin(xvar = list(Pressure = clim_w_clean$surf_pressure),
                           k = 4,
                           cdate = clim_w_clean$date,
                           bdate = dep_w_clean2$spring_dep_new.x,
                           baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry + dep_year, data = dep_w_clean2), 
                           cinterval = "week",
                           cmissing = 'method1',
                           range = c(13, 0),
                           type = "relative", 
                           stat = c("mean", "slope"),
                           func = c("lin", "quad"),
                           spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
press_w_lm_w
m <- which(press_w_lm_w$combos$DeltaAICc == min(press_w_lm_w$combos$DeltaAICc))
press_w_lm_w_bestmod <- press_w_lm_w[[m]]$BestModel
summary(press_w_lm_w_bestmod)
vif(press_w_lm_w_bestmod)
window_plot(press_w_lm_w, m, 'Surface Pressure')

# check assumptions
plot(press_w_lm_w_bestmod)
hist(resid(press_w_lm_w_bestmod))
shapiro.test(resid(press_w_lm_w_bestmod)) # not normal

# Cox proportional hazards models
press_w_cph_d <- slidingwin(xvar = list(Pressure = clim_w_clean$surf_pressure),
                            cdate = clim_w_clean$date,
                            bdate = dep_w_clean2$spring_dep_new.x,
                            baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean2))) ~ transf_ancestry + dep_year, data = dep_w_clean2), 
                            cinterval = "day",
                            cmissing = 'method1',
                            range = c(92, 0),
                            type = "relative",
                            stat = c("mean", "slope"),
                            func = c("lin", "quad"),
                            spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
press_w_cph_d
m <- which(press_w_cph_d$combos$DeltaAICc == min(press_w_cph_d$combos$DeltaAICc))
press_w_cph_d_bestmod <- press_w_cph_d[[m]]$BestModel
summary(press_w_cph_d_bestmod)
summary(press_w_cph_d_bestmod)$rsq
press_w_cph_d_survfit <- survfit(press_w_cph_d_bestmod)
AIC(press_w_cph_d_bestmod)

# Kaplan-Meier curve
plot(press_w_cph_d_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get top window
top_windows <- rbind(top_windows, top_wins(press_w_cph_d, 'daily'))

# using slope as aggregate
press_w_cph_w <- slidingwin(xvar = list(Pressure = clim_w_clean$surf_pressure),
                            cdate = clim_w_clean$date,
                            bdate = dep_w_clean2$spring_dep_new.x,
                            baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean2))) ~ transf_ancestry + dep_year, data = dep_w_clean2), 
                            cinterval = "week",
                            cmissing = 'method1',
                            range = c(13, 0),
                            type = "relative",
                            stat = c("mean", "slope"),
                            func = c("lin", "quad"),
                            spatial = list(dep_w_clean2$ID, clim_w_clean$ID))

# examine model
press_w_cph_w
m <- which(press_w_cph_w$combos$DeltaAICc == min(press_w_cph_w$combos$DeltaAICc))
press_w_cph_w_bestmod <- press_w_cph_w[[m]]$BestModel
summary(press_w_cph_w_bestmod)
summary(press_w_cph_w_bestmod)$rsq
press_w_cph_w_survfit <- survfit(press_w_cph_w_bestmod)
AIC(press_w_cph_w_bestmod)

# Kaplan-Meier curve
plot(press_w_cph_w_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get top window
top_windows <- rbind(top_windows, top_wins(press_w_cph_w, 'weekly'))

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

# split into 2 data frames and remove 2019
daylen_w <- llg_dts_w |> 
  select('ID', 'release_site', 'date', 'daylen', 'dep_year', 'transf_ancestry', 'aims_heterozygosity', 'enso') |> 
  filter(dep_year != 2019)
dep_daylen_w <- llg_w_clean |> 
  select('ID', 'release_site', 'spring_dep_new.x', 'dep_year', 'transf_ancestry', 'aims_heterozygosity', 'enso') |> 
  filter(dep_year != 2019)

# linear models
# daylength sliding window using linear model and daily window
set.seed(123)
daylen_w_lm_d <- slidingwin(xvar = list(Daylength = daylen_w$daylen),
                            k = 4,
                            cdate = daylen_w$date,
                            bdate = dep_daylen_w$spring_dep_new.x,
                            baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry + dep_year, data = dep_daylen_w),
                            cinterval = "day",
                            cmissing = 'method2',
                            range = c(92, 0),
                            type = "relative",
                            stat = c("mean", "slope"),
                            func = c("lin", "quad"),
                            spatial = list(dep_daylen_w$ID, daylen_w$ID))

# examine model
daylen_w_lm_d
m <- which(daylen_w_lm_d$combos$DeltaAICc == min(daylen_w_lm_d$combos$DeltaAICc))
daylen_w_lm_d_bestmod <- daylen_w_lm_d[[m]]$BestModel
summary(daylen_w_lm_d_bestmod)
vif(daylen_w_lm_d_bestmod) 

# best model severely overfit, try simpler model
m <- m - 2
daylen_w_lm_d_bestmod <- daylen_w_lm_d[[m]]$BestModel
summary(daylen_w_lm_d_bestmod)
vif(daylen_w_lm_d_bestmod)
window_plot(daylen_w_lm_d, m, 'Daylength')

# check assumptions
plot(daylen_w_lm_d_bestmod)
hist(resid(daylen_w_lm_d_bestmod))
shapiro.test(resid(daylen_w_lm_d_bestmod)) # not normal

# weekly linear model
set.seed(123)
daylen_w_lm_w <- slidingwin(xvar = list(Daylength = daylen_w$daylen),
                            k = 4,
                            cdate = daylen_w$date,
                            bdate = dep_daylen_w$spring_dep_new.x,
                            baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry + dep_year, data = dep_daylen_w),
                            cinterval = "week",
                            cmissing = 'method2',
                            range = c(13, 0),
                            type = "relative",
                            stat = c("mean", "slope"),
                            func = c("lin", "quad"),
                            spatial = list(dep_daylen_w$ID, daylen_w$ID))

# examine model
daylen_w_lm_w
m <- which(daylen_w_lm_w$combos$DeltaAICc == min(daylen_w_lm_w$combos$DeltaAICc))
daylen_w_lm_w_bestmod <- daylen_w_lm_w[[m]]$BestModel
summary(daylen_w_lm_w_bestmod)
vif(daylen_w_lm_w_bestmod)

# best model severely overfit, try simpler model
m <- m - 2
daylen_w_lm_w_bestmod <- daylen_w_lm_w[[m]]$BestModel
summary(daylen_w_lm_w_bestmod)
vif(daylen_w_lm_w_bestmod)
window_plot(daylen_w_lm_w, m, 'Daylength')

# check assumptions
plot(daylen_w_lm_w_bestmod)
hist(resid(daylen_w_lm_w_bestmod))
shapiro.test(resid(daylen_w_lm_w_bestmod)) # not normal

# daily Coxph
daylen_w_cph_d <- slidingwin(xvar = list(Daylength = daylen_w$daylen),
                             cdate = daylen_w$date,
                             bdate = dep_daylen_w$spring_dep_new.x,
                             baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_daylen_w))) ~ transf_ancestry + dep_year, data = dep_daylen_w),
                             cinterval = "day",
                             cmissing = 'method1',
                             range = c(92, 0),
                             type = "relative", 
                             stat = c("mean", "slope"),
                             func = c("lin", "quad"),
                             spatial = list(dep_daylen_w$ID, daylen_w$ID))

# examine model
daylen_w_cph_d
m <- which(daylen_w_cph_d$combos$DeltaAICc == min(daylen_w_cph_d$combos$DeltaAICc))
daylen_w_cph_d_bestmod <- daylen_w_cph_d[[m]]$BestModel
summary(daylen_w_cph_d_bestmod)
summary(daylen_w_cph_d_bestmod)$rsq
daylen_w_cph_d_survfit <- survfit(daylen_w_cph_d_bestmod)
AIC(daylen_w_cph_d_bestmod)

# Kaplan-Meier curve
plot(daylen_w_cph_d_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get top window
top_windows <- rbind(top_windows, top_wins(daylen_w_cph_d, 'daily'))

# using slope as aggregate
daylen_w_cph_w <- slidingwin(xvar = list(Daylength = daylen_w$daylen),
                             cdate = daylen_w$date,
                             bdate = dep_daylen_w$spring_dep_new.x,
                             baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_daylen_w))) ~ transf_ancestry + dep_year, data = dep_daylen_w),
                             cinterval = "week",
                             cmissing = 'method1',
                             range = c(13, 0),
                             type = "relative",
                             stat = c("mean", "slope"),
                             func = c("lin", "quad"),
                             spatial = list(dep_daylen_w$ID, daylen_w$ID))

# examine model
daylen_w_cph_w
m <- which(daylen_w_cph_w$combos$DeltaAICc == min(daylen_w_cph_w$combos$DeltaAICc))
daylen_w_cph_w_bestmod <- daylen_w_cph_w[[m]]$BestModel
summary(daylen_w_cph_w_bestmod)
summary(daylen_w_cph_w_bestmod)$rsq
daylen_w_cph_w_survfit <- survfit(daylen_w_cph_w_bestmod)
AIC(daylen_w_cph_w_bestmod)

# Kaplan-Meier curve
plot(daylen_w_cph_w_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get top window
top_windows <- rbind(top_windows, top_wins(daylen_w_cph_w, 'weekly'))

# WINTERING NDVI------------------------------------------------------------------------------------------------------------------
# remove 2019
ndvi_w_clean <- ndvi_w |> filter(dep_year != 2019)
dep_w_clean <- dep_w |> filter(dep_year != 2019)

# linear models
# NDVI sliding window using mean as aggregate and daily window
set.seed(123)
ndvi_w_lm_d <- slidingwin(xvar = list(NDVI = ndvi_w_clean$ndvi),
                          k = 4,
                          cdate = ndvi_w_clean$date,
                          bdate = dep_w_clean$spring_dep_new.x,
                          baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry + dep_year, data = dep_w_clean), 
                          cinterval = "day",
                          cmissing = 'method1',
                          range = c(92, 0),
                          type = "relative", 
                          stat = c("mean", "slope"),
                          func = c("lin", "quad"),
                          spatial = list(dep_w_clean$ID, ndvi_w_clean$ID))


# examine model
ndvi_w_lm_d
m <- which(ndvi_w_lm_d$combos$DeltaAICc == min(ndvi_w_lm_d$combos$DeltaAICc))
ndvi_w_lm_d_bestmod <- ndvi_w_lm_d[[m]]$BestModel
summary(ndvi_w_lm_d_bestmod)
vif(ndvi_w_lm_d_bestmod)
window_plot(ndvi_w_lm_d, m, 'NDVI')

# check assumptions
plot(ndvi_w_lm_d_bestmod)
hist(resid(ndvi_w_lm_d_bestmod))
shapiro.test(resid(ndvi_w_lm_d_bestmod))
bptest(yvar ~., data = ndvi_w_lm_d[[m]]$BestModelData)

# get top windows
top_windows <- rbind(top_windows, top_wins(ndvi_w_lm_d, 'daily'))

# temperature sliding window using slope as aggregate
set.seed(123)
ndvi_w_lm_w <- slidingwin(xvar = list(NDVI = ndvi_w_clean$ndvi),
                          k = 4,
                          cdate = ndvi_w_clean$date,
                          bdate = dep_w_clean$spring_dep_new.x,
                          baseline = lm(yday(spring_dep_new.x) ~ transf_ancestry + dep_year, data = dep_w_clean), 
                          cinterval = "week",
                          cmissing = 'method1',
                          range = c(13, 0),
                          type = "relative",
                          stat = c("mean", "slope"),
                          func = c("lin", "quad"),
                          spatial = list(dep_w_clean$ID, ndvi_w_clean$ID))

# examine model
ndvi_w_lm_w
m <- which(ndvi_w_lm_w$combos$DeltaAICc == min(ndvi_w_lm_w$combos$DeltaAICc))
ndvi_w_lm_w_bestmod <- ndvi_w_lm_w[[m]]$BestModel
summary(ndvi_w_lm_w_bestmod)
vif(ndvi_w_lm_w_bestmod)
window_plot(ndvi_w_lm_w, m, 'NDVI')

# check assumptions
plot(ndvi_w_lm_w_bestmod)
hist(resid(ndvi_w_lm_w_bestmod))
shapiro.test(resid(ndvi_w_lm_w_bestmod))
bptest(yvar ~., data = ndvi_w_lm_w[[m]]$BestModelData)

# get top windows
top_windows <- rbind(top_windows, top_wins(ndvi_w_lm_w, 'weekly'))

# Cox proportional hazards models
# using mean as aggregate
ndvi_w_cph_d <- slidingwin(xvar = list(NDVI = ndvi_w_clean$ndvi),
                           cdate = ndvi_w_clean$date,
                           bdate = dep_w_clean$spring_dep_new.x,
                           baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean))) ~ transf_ancestry + dep_year, data = dep_w_clean),
                           cinterval = "day",
                           cmissing = 'method1',
                           range = c(92, 0),
                           type = "relative",
                           stat = c("mean", "slope"),
                           func = c("lin", "quad"),
                           spatial = list(dep_w_clean$ID, ndvi_w_clean$ID))

# examine model
ndvi_w_cph_d
m <- which(ndvi_w_cph_d$combos$DeltaAICc == min(ndvi_w_cph_d$combos$DeltaAICc))
ndvi_w_cph_d_bestmod <- ndvi_w_cph_d[[m]]$BestModel
summary(ndvi_w_cph_d_bestmod)
summary(ndvi_w_cph_d_bestmod)$rsq
ndvi_w_cph_d_survfit <- survfit(ndvi_w_cph_d_bestmod)
AIC(ndvi_w_cph_d_bestmod)

# Kaplan-Meier curve
plot(ndvi_w_cph_d_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get top windows
top_windows <- rbind(top_windows, top_wins(ndvi_w_cph_d, 'daily'))

# using slope as aggregate
ndvi_w_cph_w <- slidingwin(xvar = list(NDVI = ndvi_w_clean$ndvi),
                           cdate = ndvi_w_clean$date,
                           bdate = dep_w_clean$spring_dep_new.x,
                           baseline = coxph(Surv(yday(spring_dep_new.x), rep(1, nrow(dep_w_clean))) ~ transf_ancestry + dep_year, data = dep_w_clean), # I'm confused on how to incorporate release site as a covariate...
                           cinterval = "week",
                           cmissing = 'method1',
                           range = c(13, 0),
                           type = "relative",
                           stat = c("mean", "slope"),
                           func = c("lin", "quad"), 
                           spatial = list(dep_w_clean$ID, ndvi_w_clean$ID))

# examine model
ndvi_w_cph_w
m <- which(ndvi_w_cph_w$combos$DeltaAICc == min(ndvi_w_cph_w$combos$DeltaAICc))
ndvi_w_cph_w_bestmod <- ndvi_w_cph_w[[m]]$BestModel
summary(ndvi_w_cph_w_bestmod)
summary(ndvi_w_cph_w_bestmod)$rsq
ndvi_w_cph_w_survfit <- survfit(ndvi_w_cph_w_bestmod)
AIC(ndvi_w_cph_w_bestmod)

# Kaplan-Meier curve
plot(ndvi_w_cph_w_survfit, xlab = "Day",
     ylab = "Estimated Probability of Not Departing")

# get top windows
top_windows <- rbind(top_windows, top_wins(ndvi_w_cph_w, 'weekly'))

write.csv(top_windows, 'env_data/wintering_relative_all_top_windows.csv')

# COMBINE WINTERING WINDOWS-----------------------------------------------------------

# for use later
extract_win <- function(win_obj, win_type, n) {
  data <- win_obj[[n]]$BestModelData
  data$yvar <- as.numeric(data$yvar)[1:nrow(data)]
  data$climate <- as.numeric(data$climate)
  close <- win_obj$combos[n, ]$WindowClose
  open <- win_obj$combos[n, ]$WindowOpen
  clim <- win_obj$combos[n, ]$climate
  agg <- win_obj$combos[n, ]$stat
  if (win_obj$combos[n, ]$func == 'quad') {
    data$`I(climate^2)` <- as.numeric(data$`I(climate^2)`)
    colnames(data)[5] <- paste0(clim, '2_', win_type, '_', agg, '_', open, '_', close)
  }
  colnames(data)[4] <- paste0(clim, '_', win_type, '_', agg, '_', open, '_', close)
  return(data)
}
