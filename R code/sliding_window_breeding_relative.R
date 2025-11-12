# ROUGH DRAFT VERSION: using relative window and testing linear vs. quadratic, 
# weekly vs. daily, slope vs. mean
# CHANGES: using Coxph only, and attempting randomization to validate results; 
# only using weekly windows for speed and to prevent super short windows
# PURPOSE: run climwin sliding window on environmental data from Google Earth 
# Engine corresponding to Swainson's Thrush observations from the wintering 
# season only
# NOTES: have to remerge with original LLG dataset to get release site as covariate;
# STILL NEED TO ADD IN RANDOMIZATION FOR VARIABLES BESIDES DAYLENGTH

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
library(coxme)

# SET UP------------------------------------------------------------------------------------------------------------------------------------------------
# WINTERING
# load data
llg_env_land <- read.csv('data/llg_env_land.csv')
ndvi <- read.csv('env_data/gee/NDVI_breeding.csv')
clim <- read.csv('env_data/gee/clim_breeding.csv')
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
llg_ancestry <- llg_ancestry |> mutate(cat_ancestry = case_when(aims_heterozygosity >= 0.75~'F1',
                                                                aims_ancestry <= 0.1~'pure_inland',
                                                                aims_ancestry >= 0.9~'pure_coastal',
                                                                aims_ancestry >= 0.1 & aims_ancestry <= 0.3 & aims_heterozygosity <= 0.5~'inland_BC',
                                                                aims_ancestry >= 0.7 & aims_ancestry <= 0.9 & aims_heterozygosity <= 0.5~'coastal_BC',
                                                                aims_ancestry >= 0.25 & aims_ancestry <= 0.75 & aims_heterozygosity < 0.75~'later'))
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
llg_sub <- llg_merged[, c(1:5, 7:12, 15, 48)]

llg_sub$release_site <- (llg_sub |> 
                           mutate(release_site_comb = case_when(release_site == 'Pacific Spirit' ~ 'Porpoise Bay',
                                                                release_site == 'Northern_BC' ~ 'Alaska',
                                                                ((release_site != 'Pacific Spirit') & (release_site != 'Northern_BC')) ~ release_site
                           )))$release_site_comb

llg_b <- llg_sub |> 
  select('ID', 'release_GPS.N', 'release_GPS.W', 'fall_dep_new.x', 'release_site', 'age', 'cat_ancestry') |> 
  mutate(dep_year = as.factor(year(fall_dep_new.x)))

# merge with NDVI data so we have release site 
ndvi_merged <- merge(x = llg_b, y = ndvi, by = c('ID'), all.x = FALSE)
ndvi_merged_clean <- ndvi_merged |> 
  select('ID', 'fall_dep_new.x', 'release_GPS.N', 'release_GPS.W', 'date', 'ndvi', 'release_site', 'dep_year', 'cat_ancestry',  'age')

# split into datasets for sliding window
dep_b <- llg_b |> 
  select('fall_dep_new.x', 'ID', 'release_site', 'dep_year', 'cat_ancestry', 'age') |> 
  filter(!is.na(fall_dep_new.x))
ndvi_b <- ndvi_merged_clean |> select('date', 'ndvi', 'ID', 'release_site', 'dep_year', 'cat_ancestry', 'age')

# convert dates
ndvi_b$date <- as.Date(ndvi_b$date)

# save IDs of birds with missing data
missing_table_ndvi <- ndvi_merged_clean |> 
  group_by(ID, fall_dep_new.x) |> 
  summarise(n_missing = sum(is.na(ndvi)))

# merge with climate data so we have release site 
clim_merged <- merge(x = llg_b, y = clim, by = c('ID'), all.x = FALSE)
clim_merged_clean <- clim_merged |> 
  select('ID', 'fall_dep_new.x', 'release_GPS.N', 'release_GPS.W', 'date', 'temp_2m', 'total_precip', 'surf_pressure', 'u_wind', 'v_wind', 'release_site', 'dep_year', 'cat_ancestry', 'age')

# calculate wind speed
clim_merged_clean <- clim_merged_clean |> mutate(wind_speed = sqrt(u_wind^2 + v_wind^2))

# split into datasets for sliding window
dep_b <- llg_b |> 
  select('fall_dep_new.x', 'ID', 'release_site', 'dep_year', 'cat_ancestry', 'age') |> 
  filter(!is.na(fall_dep_new.x))
clim_b <- clim_merged_clean |> select('date', 'temp_2m', 'total_precip', 'surf_pressure', 'wind_speed', 'ID', 'release_site', 'dep_year', 'cat_ancestry', 'age')

# convert dates
clim_b$date <- as.Date(clim_b$date)

# examine and remove NAs
clim_b |> filter(is.na(temp_2m)) |> group_by(ID) |> summarise(n_rows = n())
missing_clim <- (clim_b |> filter(is.na(temp_2m)))$ID |> unique() # no missing data, but keeping this for consistency with wintering code
clim_b_clean <- clim_b |> filter(!(ID %in% missing_clim))
dep_b_clean2 <- dep_b |> filter(!(ID %in% missing_clim))


# figure out reference dep date using earliest departure
(dep_b |> 
    mutate(new_date = as.Date(paste0(month(fall_dep_new.x), '-', day(fall_dep_new.x), '-2018'), '%m-%d-%Y')))$new_date |> min() 

# figure out latest dep date using latest departure
(dep_b |> 
    mutate(new_date = as.Date(paste0(month(fall_dep_new.x), '-', day(fall_dep_new.x), '-2018'), '%m-%d-%Y')))$new_date |> max() 

# figure out average dep date
(dep_b |> 
    mutate(new_date = as.Date(paste0(month(fall_dep_new.x), '-', day(fall_dep_new.x), '-2018'), '%m-%d-%Y')))$new_date |> mean() 

# BREEDING DAYLENGTH----------------------------------------------------------------------------------------------------------------------
# get rid of any rows missing lat or date
llg_b_clean <- llg_b |> 
  filter((!is.na(fall_dep_new.x)) & !is.na(release_GPS.N)) |> 
  arrange(ID)

# figure out latest day of year departure
max_dt <- (llg_b_clean |> mutate(doy = yday(fall_dep_new.x)) |> arrange(-doy))[1, ]$fall_dep_new.x
max_dt_str <- substr(as.character(max_dt), 5, 11)
min_dt <- (llg_b_clean |> mutate(doy = yday(fall_dep_new.x)) |> arrange(doy))[1, ]$fall_dep_new.x
min_dt_str <- substr(as.character(min_dt), 5, 11)

# get years
yrs <- llg_b_clean$dep_year |> unique()

# get sequence of dates for full range of departure days across all years
all_dts <- c()
for (y in yrs) {
  yr <- as.character(y)
  dt_seq <- seq.Date(from = (as.Date(paste0(yr, max_dt_str), '%Y-%m-%d')+7), to = (as.Date(paste0(yr, min_dt_str), '%Y-%m-%d')-100), by = '-1 days')
  all_dts <- c(all_dts, dt_seq)
}
# convert to dates
all_dts <- as.Date(all_dts)

# convert list of lists to one long data frame labeled by ID
ids <- llg_b_clean$ID |> unique()
dts_df <- as.data.frame(cbind(rep(ids, each = length(all_dts)), rep(all_dts,  times = length(ids))))
colnames(dts_df) <- c('ID', 'date')
# check number of rows per bird look correct
(dts_df |> group_by(ID) |> 
    summarise(nrow = n()))$nrow |> range()

# append dates
llg_dts_b <- merge(x = llg_b_clean, y = dts_df, by = 'ID', all.y = TRUE)

# extract photoperiod
llg_dts_b <- llg_dts_b |> mutate(daylen = daylength(release_GPS.N, date))
llg_dts_b$date <- as.Date(llg_dts_b$date)

# split into 2 data frames
daylen_b <- llg_dts_b |> 
  select('ID', 'release_site', 'date', 'daylen', 'dep_year', 'age', 'cat_ancestry')
dep_daylen_b <- llg_b_clean |> 
  select('ID', 'release_site', 'fall_dep_new.x', 'dep_year', 'age', 'cat_ancestry')

# split into 2 data frames
daylen_b <- llg_dts_b |> 
  select('ID', 'release_site', 'date', 'daylen', 'dep_year', 'cat_ancestry', 'age')
dep_daylen_b <- llg_b_clean |> 
  select('ID', 'release_site', 'fall_dep_new.x', 'dep_year', 'cat_ancestry', 'age')

# weekly Coxph
daylen_b_cph_w <- slidingwin(xvar = list(Daylength = daylen_b$daylen),
                             cdate = daylen_b$date,
                             bdate = dep_daylen_b$fall_dep_new.x,
                             baseline = coxph(Surv(yday(fall_dep_new.x), rep(1, nrow(dep_daylen_b))) ~ cat_ancestry + age, data = dep_daylen_b),
                             cinterval = "week",
                             cmissing = 'method1',
                             range = c(13, 0),
                             type = "relative", 
                             stat = c("mean", "slope"),
                             func = c("lin", "quad"),
                             spatial = list(dep_daylen_b$ID, daylen_b$ID))

daylen_b_rand_w <- randwin(repeats = 100, xvar = list(Daylength = daylen_b$daylen),
                           cdate = daylen_b$date,
                           bdate = dep_daylen_b$fall_dep_new.x,
                           baseline = coxph(Surv(yday(fall_dep_new.x), rep(1, nrow(dep_daylen_b))) ~ cat_ancestry + age, data = dep_daylen_b),
                           cinterval = "week",
                           cmissing = 'method1',
                           range = c(13, 0),
                           type = "relative", 
                           stat = "mean",
                           func = "lin",
                           spatial = list(dep_daylen_b$ID, daylen_b$ID))

# examine best model
daylen_b_cph_w$combos
m <- which(daylen_b_cph_w$combos$DeltaAICc == min(daylen_b_cph_w$combos$DeltaAICc))
daylen_b_cph_w_bestmod <- daylen_b_cph_w[[m]]$BestModel
AIC(daylen_b_cph_w_bestmod)
summary(daylen_b_cph_w_bestmod)$rsq

# check p value from comparison with random windows
pvalue(daylen_b_cph_w[[m]]$Dataset, daylen_b_rand_w[[1]], metric = 'AIC')

# best day length dataset
best_daylen_data <- daylen_b_cph_w[[m]]$BestModelData
best_daylen_info <- daylen_b_cph_w$combos[m, ]

# refit model
best_daylen_data$dep_year <- dep_daylen_b$dep_year
best_daylen_data$event <- 1
best_daylen_data$yvar <- as.numeric(best_daylen_data$yvar)[1:nrow(best_daylen_data)]
best_daylen_data$ID <- dep_daylen_b$ID
daylen_refit_mod <- coxme(Surv(yvar, event) ~ cat_ancestry + age + (1|dep_year) + climate, data = best_daylen_data)

# examine refit model
summary(daylen_refit_mod)
AIC(daylen_refit_mod)

# plot data
ggplot(best_daylen_data, aes(x = yvar, y = climate)) + geom_point()

# BREEDING PRECIPITATION--------------------------------------------------------------------------------------------------------------

# weekly Coxph
precip_b_cph_w <- slidingwin(xvar = list(Precip = clim_b_clean$total_precip),
                             cdate = clim_b_clean$date,
                             bdate = dep_b_clean2$fall_dep_new.x,
                             baseline = coxph(Surv(yday(fall_dep_new.x), rep(1, nrow(dep_b_clean2))) ~ cat_ancestry + age, data = dep_b_clean2), 
                             cinterval = "week",
                             cmissing = 'method1',
                             range = c(13, 0),
                             type = "relative",
                             stat = c("mean", "slope"),
                             func = c("lin", "quad"), 
                             spatial = list(dep_b_clean2$ID, clim_b_clean$ID))

# examine best model
precip_b_cph_w$combos
m <- which(precip_b_cph_w$combos$DeltaAICc == min(precip_b_cph_w$combos$DeltaAICc))
m <- m -  2
precip_b_cph_w_bestmod <- precip_b_cph_w[[m]]$BestModel
AIC(precip_b_cph_w_bestmod)
summary(precip_b_cph_w_bestmod)$rsq

# best day length dataset
best_precip_data <- precip_b_cph_w[[m]]$BestModelData
best_precip_info <- precip_b_cph_w$combos[m, ]

# refit model
best_precip_data$dep_year <- dep_b_clean2$dep_year
best_precip_data$event <- 1
best_precip_data$yvar <- as.numeric(best_precip_data$yvar)[1:nrow(best_precip_data)]
best_precip_data$ID <- dep_b_clean2$ID
precip_refit_mod <- coxme(Surv(yvar, event) ~ cat_ancestry + age + (1|dep_year) + climate, data = best_precip_data)

# examine refit model
summary(precip_refit_mod)
AIC(precip_refit_mod)

# plot data to get general sense
ggplot(best_precip_data, aes(x = yvar, y = climate)) + geom_point()

# WINTERING TEMPERATURE--------------------------------------------------------------------------------------------------------------

# weekly Coxph
temp_b_cph_w <- slidingwin(xvar = list(Temp = clim_b_clean$temp_2m),
                           cdate = clim_b_clean$date,
                           bdate = dep_b_clean2$fall_dep_new.x,
                           baseline = coxph(Surv(yday(fall_dep_new.x), rep(1, nrow(dep_b_clean2))) ~ cat_ancestry + age, data = dep_b_clean2), 
                           cinterval = "week",
                           cmissing = 'method1',
                           range = c(13, 0),
                           type = "relative",
                           stat = c("mean", "slope"),
                           func = c("lin", "quad"), 
                           spatial = list(dep_b_clean2$ID, clim_b_clean$ID))

# examine best model
temp_b_cph_w$combos
m <- which(temp_b_cph_w$combos$DeltaAICc == min(temp_b_cph_w$combos$DeltaAICc))
m <- m - 2
temp_b_cph_w_bestmod <- temp_b_cph_w[[m]]$BestModel
AIC(temp_b_cph_w_bestmod)
summary(temp_b_cph_w_bestmod)$rsq

# best day length dataset
best_temp_data <- temp_b_cph_w[[m]]$BestModelData
best_temp_info <- temp_b_cph_w$combos[m, ]

# refit model
best_temp_data$dep_year <- dep_b_clean2$dep_year
best_temp_data$event <- 1
best_temp_data$yvar <- as.numeric(best_temp_data$yvar)[1:nrow(best_temp_data)]
best_temp_data$ID <- dep_b_clean2$ID
temp_refit_mod <- coxme(Surv(yvar, event) ~ cat_ancestry + age + (1|dep_year) + climate, data = best_temp_data)

# examine refit model
summary(temp_refit_mod)
AIC(temp_refit_mod)

# plot data to get general sense
ggplot(best_temp_data, aes(x = yvar, y = climate)) + geom_point()

# WINTERING SURFACE PRESSURE--------------------------------------------------------------------------------------------------------------

# weekly Coxph
press_b_cph_w <- slidingwin(xvar = list(Press = clim_b_clean$surf_pressure),
                            cdate = clim_b_clean$date,
                            bdate = dep_b_clean2$fall_dep_new.x,
                            baseline = coxph(Surv(yday(fall_dep_new.x), rep(1, nrow(dep_b_clean2))) ~ cat_ancestry + age, data = dep_b_clean2), 
                            cinterval = "week",
                            cmissing = 'method1',
                            range = c(13, 0),
                            type = "relative",
                            stat = c("mean", "slope"),
                            func = c("lin", "quad"), 
                            spatial = list(dep_b_clean2$ID, clim_b_clean$ID))

# examine best model
press_b_cph_w$combos
m <- which(press_b_cph_w$combos$DeltaAICc == min(press_b_cph_w$combos$DeltaAICc))
press_b_cph_w_bestmod <- press_b_cph_w[[m]]$BestModel
AIC(press_b_cph_w_bestmod)
summary(press_b_cph_w_bestmod)$rsq

# best day length dataset
best_press_data <- press_b_cph_w[[m]]$BestModelData
best_press_info <- press_b_cph_w$combos[m, ]

# refit model
best_press_data$dep_year <- dep_b_clean2$dep_year
best_press_data$event <- 1
best_press_data$yvar <- as.numeric(best_press_data$yvar)[1:nrow(best_press_data)]
best_press_data$ID <- dep_b_clean2$ID
press_refit_mod <- coxme(Surv(yvar, event) ~ cat_ancestry + age + (1|dep_year) + climate, data = best_press_data)

# examine refit model
summary(press_refit_mod)
AIC(press_refit_mod)

# plot data to get general sense
ggplot(best_press_data, aes(x = yvar, y = climate)) + geom_point()

# WINTERING WIND SPEED--------------------------------------------------------------------------------------------------------------

# weekly Coxph
wind_b_cph_w <- slidingwin(xvar = list(Wind = clim_b_clean$wind_speed),
                           cdate = clim_b_clean$date,
                           bdate = dep_b_clean2$fall_dep_new.x,
                           baseline = coxph(Surv(yday(fall_dep_new.x), rep(1, nrow(dep_b_clean2))) ~ cat_ancestry + age, data = dep_b_clean2), 
                           cinterval = "week",
                           cmissing = 'method1',
                           range = c(13, 0),
                           type = "relative",
                           stat = c("mean", "slope"),
                           func = c("lin", "quad"), 
                           spatial = list(dep_b_clean2$ID, clim_b_clean$ID))

# examine best model
wind_b_cph_w$combos
m <- which(wind_b_cph_w$combos$DeltaAICc == min(wind_b_cph_w$combos$DeltaAICc))
wind_b_cph_w_bestmod <- wind_b_cph_w[[m]]$BestModel
AIC(wind_b_cph_w_bestmod)
summary(wind_b_cph_w_bestmod)$rsq

# best day length dataset
best_wind_data <- wind_b_cph_w[[m]]$BestModelData
best_wind_info <- wind_b_cph_w$combos[m, ]

# refit model
best_wind_data$dep_year <- dep_b_clean2$dep_year
best_wind_data$event <- 1
best_wind_data$yvar <- as.numeric(best_wind_data$yvar)[1:nrow(best_wind_data)]
best_wind_data$ID <- dep_b_clean2$ID
wind_refit_mod <- coxme(Surv(yvar, event) ~ cat_ancestry + age + (1|dep_year) + climate + `I(climate^2)`, data = best_wind_data)

# examine refit model
summary(wind_refit_mod)
AIC(wind_refit_mod)

# plot data to get general sense
ggplot(best_wind_data, aes(x = yvar, y = climate)) + geom_point()


# WINTERING NDVI-------------------------------------------------------------------------------------------------------------

# weekly Coxph
ndvi_b_cph_w <- slidingwin(xvar = list(NDVI = ndvi_b$ndvi),
                           cdate = ndvi_b$date,
                           bdate = dep_b$fall_dep_new.x,
                           baseline = coxph(Surv(yday(fall_dep_new.x), rep(1, nrow(dep_b))) ~ cat_ancestry + age, data = dep_b), # I'm confused on how to incorporate release site as a covariate...
                           cinterval = "week",
                           cmissing = 'method1',
                           range = c(13, 0),
                           type = "relative",
                           stat = c("mean", "slope"),
                           func = c("lin", "quad"), 
                           spatial = list(dep_b$ID, ndvi_b$ID))

# examine best model
ndvi_b_cph_w$combos
m <- which(ndvi_b_cph_w$combos$DeltaAICc == min(ndvi_b_cph_w$combos$DeltaAICc))
ndvi_b_cph_w_bestmod <- ndvi_b_cph_w[[m]]$BestModel
AIC(ndvi_b_cph_w_bestmod)
summary(ndvi_b_cph_w_bestmod)$rsq

# best day length dataset
best_ndvi_data <- ndvi_b_cph_w[[m]]$BestModelData
best_ndvi_info <- ndvi_b_cph_w$combos[m, ]

# refit model
best_ndvi_data$dep_year <- dep_b$dep_year
best_ndvi_data$event <- 1
best_ndvi_data$yvar <- as.numeric(best_ndvi_data$yvar)[1:nrow(best_ndvi_data)]
best_ndvi_data$ID <- dep_b$ID
ndvi_refit_mod <- coxme(Surv(yvar, event) ~ cat_ancestry + age + (1|dep_year) + climate, data = best_ndvi_data)

# examine refit model
summary(ndvi_refit_mod)
AIC(ndvi_refit_mod)

# plot data to get general sense
ggplot(best_ndvi_data, aes(x = yvar, y = climate)) + geom_point()

# SAVING DATA------------------------------------------------------------------------------------------------------------
rename_cols <- function(best_data, best_info) {
  colnames(best_data)[which(colnames(best_data) == 'climate')] <- paste0(best_info$climate, '_', best_info$WindowOpen, '_', best_info$WindowClose)
  if (best_info$func == 'quad') {
    colnames(best_data)[which(colnames(best_data) == 'I(climate^2)')] <- paste0(best_info$climate, '2_', best_info$WindowOpen, '_', best_info$WindowClose)
  }
  return(best_data)
}

best_data_b <- rename_cols(best_daylen_data, best_daylen_info)
best_data_b <- merge(best_data_b, rename_cols(best_precip_data, best_precip_info), by = c('ID', 'yvar', 'cat_ancestry', 'age', 'dep_year', 'event'))
best_data_b <- merge(best_data_b, rename_cols(best_temp_data, best_temp_info), by = c('ID', 'yvar', 'cat_ancestry', 'age', 'dep_year', 'event'))
best_data_b <- merge(best_data_b, rename_cols(best_press_data, best_press_info), by = c('ID', 'yvar', 'cat_ancestry', 'age', 'dep_year', 'event'))
best_data_b <- merge(best_data_b, rename_cols(best_wind_data, best_wind_info), by = c('ID', 'yvar', 'cat_ancestry', 'age', 'dep_year', 'event'))
best_data_b <- merge(best_data_b, rename_cols(best_ndvi_data, best_ndvi_info), by = c('ID', 'yvar', 'cat_ancestry', 'age', 'dep_year', 'event'))

write.csv(best_data_b, 'env_data/breeding_rel_slidingwin_v2_windows.csv')

write.csv(rbind(best_daylen_info, best_precip_info, best_temp_info, best_press_info, best_wind_info, best_ndvi_info), 'env_data/breeding_rel_slidingwin_v2_info.csv')
