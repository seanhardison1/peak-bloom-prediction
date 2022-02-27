library(tidyverse)
library(sf)
library(rEDM)
library(lubridate)
library(tsibble)
library(forecast)
library(GGally)
library(mgcv)
library(INLA)
library(gratia)
library(rnoaa)
library(raster)
library(sdmTMB)
library(fpp3)
library(rnaturalearth)
library(readxl)

# do environmental forecasts or pre-load
fc <- F

#load data
load(here::here("data/all_cb_data.rdata"))

# process environmental data----

# query raw
if (fc){
  dc_temps <- ghcnd(stationid = "USC00186350", refresh = TRUE) 
  save(dc_temps, file = here::here("data/dc_temps_raw.rdata"))
} else {
  load(here::here("data/dc_temps_raw.rdata"))
}

dc_temps2 <-
  dc_temps  %>% 
  dplyr::select_at(vars(year, month, element, contains("VALUE"))) %>% 
  rowwise() %>% 
  mutate(val = mean(c_across(VALUE1:VALUE31), na.rm = T)) %>% 
  ungroup() %>% 
  dplyr::select(year, month, element, val) %>% 
  spread(.,element, val) %>% 
  dplyr::filter(month %in% 1:3) %>% 
  dplyr::select(tmax = TMAX,
                tmin = TMIN,
                temp = TOBS,
                precip = PRCP,
                year, month)

jan <- dc_temps2 %>% filter(month == 1) %>% 
  rename_at(vars(1:4), function(x)paste0("j_",x)) %>% 
  dplyr::select(-month)
feb <- dc_temps2 %>% filter(month == 2) %>% 
  rename_at(vars(1:4), function(x)paste0("f_",x)) %>% 
  dplyr::select(-month)
mar <- dc_temps2 %>% filter(month == 3) %>% 
  rename_at(vars(1:4), function(x)paste0("m_",x)) %>% 
  dplyr::select(-month)
dc_temps3 <- jan %>% 
  left_join(.,feb, by = c("year")) %>% 
  left_join(.,mar, by = c("year"))

ggplot(dc) +
  geom_point(aes(y = bloom_doy, x = year)) +
  geom_line(aes(y = bloom_doy, x = year))

# project environmental data----
base <- dc_temps3 %>%
  tsibble(index = "year") %>% 
  fill_gaps()

if (fc){
  output_mtmax <- 
    base %>% 
    model(
      m_tmax = NNETAR(m_tmax)
    ) %>% 
    forecast(h = 10) %>% 
    tibble()
  
  output_ftmax <- 
    base %>% 
    model(
      f_tmax = NNETAR(f_tmax)
    ) %>% 
    forecast(h = 10) %>% 
    tibble()
  save(output_ftmax, output_mtmax, file = here::here("data/dc_env_fc.rdata"))
} else {
  load(here::here("data/dc_env_fc.rdata"))
}

proj_df <- 
  output_ftmax %>% 
  dplyr::select(year, f_tmax = .mean) %>% 
  left_join(.,output_mtmax %>% 
              dplyr::select(year, m_tmax = .mean))

# fit the bloom DOY model----
dc_sample <- dc %>% 
  left_join(.,dc_temps3) %>% 
  tsibble(index = "year")

m <- gam(bloom_doy ~
           s(year) +
           s(m_tmax) +
           s(f_tmax),
    data = dc_sample)

summary(m)
acf(m$residuals)
gratia::appraise(m)
gratia::draw(m)

# new data for prediction/projection
ndf <- 
  tibble(f_tmax = dc_sample$f_tmax,
         m_tmax = dc_sample$m_tmax,
         bloom_doy = dc_sample$bloom_doy,
         year = dc_sample$year)  %>% 
  dplyr::select(year, m_tmax, f_tmax) %>% 
  bind_rows(.,proj_df) %>% 
  # enter f_tmax for most of feb 2022
  mutate(f_tmax = ifelse(year == 2022, 111.5384616, f_tmax))

pred <- 
      predict(m, newdata = ndf, se.fit = T)

pred_df <- tibble(fit = pred$fit,
                  se = pred$se.fit,
                  year = ndf$year) 

ggplot() +
  geom_point(data = dc_sample, aes(y = bloom_doy, x = year)) +
  geom_line(data = pred_df, aes(y = fit, x = year))

# write out
dc_proj <- pred_df %>% 
  filter(year > 2021) %>% 
  mutate(location = "washington, dc",
         bloom_doy = round(fit)) %>% 
  dplyr::select(bloom_doy, year, location)
write.csv(dc_proj, file = here::here("data/dc_projection.csv"), row.names = F)
