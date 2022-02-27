library(tidyverse)
library(sf)
library(rEDM)
library(lubridate)
library(tsibble)
library(forecast)
library(mgcv)
library(INLA)
library(gratia)
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
  bc_temps <- ghcnd(stationid = "CA001108395", refresh = TRUE) 
  save(bc_temps, file = here::here("data/bc_temps_raw.rdata"))
} else {
  load(here::here("data/bc_temps_raw.rdata"))
}

bc_temps2 <-
  bc_temps  %>% 
  dplyr::select_at(vars(year, month, element, contains("VALUE"))) %>% 
  rowwise() %>% 
  mutate(val = mean(c_across(VALUE1:VALUE31), na.rm = T)) %>% 
  ungroup() %>% 
  dplyr::select(year, month, element, val) %>% 
  spread(.,element, val) %>% 
  dplyr::filter(month %in% 1:3) %>% 
  dplyr::select(tmax = TMAX,
                tmin = TMIN,
                temp = TAVG,
                precip = PRCP,
                year, month)

jan <- bc_temps2 %>% filter(month == 1) %>% 
  rename_at(vars(1:4), function(x)paste0("j_",x)) %>% 
  dplyr::select(-month)
feb <- bc_temps2 %>% filter(month == 2) %>% 
  rename_at(vars(1:4), function(x)paste0("f_",x)) %>% 
  dplyr::select(-month)
mar <- bc_temps2 %>% filter(month == 3) %>% 
  rename_at(vars(1:4), function(x)paste0("m_",x)) %>% 
  dplyr::select(-month)
bc_temps3 <- jan %>% 
  left_join(.,feb, by = c("year")) %>% 
  left_join(.,mar, by = c("year"))


#load data
load(here::here("data/all_cb_data.rdata"))

bc_blooms <- read_excel(here::here("data/bc_blooms.xlsx")) %>% 
  mutate(bloom_doy = yday(date)) %>% 
  left_join(.,bc_temps3)

ggplot(bc_blooms) +
  geom_point(aes(y = bloom_doy, x = year))

m <- gam(bloom_doy ~ s(f_tmax, k = 3) + s(m_tmax), 
              data = bc_blooms)
summary(m)

base <- bc_blooms %>% 
  tsibble(index = "year")

if (fc){
  output_m_tmax <- 
    base %>% 
    model(
      j_tmax = NNETAR(m_tmax)
    ) %>% 
    forecast(h = 10) %>% 
    tibble()
  
  output_f_tmax <- 
    base %>% 
    model(
      f_tmax = NNETAR(f_tmax)
    ) %>% 
    forecast(h = 10) %>% 
    tibble()
  
  bc_proj <- 
    output_m_tmax %>% 
    dplyr::select(year, m_tmax = .mean) %>% 
    left_join(.,output_f_tmax %>% 
                dplyr::select(year, f_tmax = .mean))
  save(bc_proj, file = here::here("data/bc_env_fc.rdata"))
} else {
  load(here::here("data/bc_env_fc.rdata"))
}

ndf <- bc_blooms %>% 
  dplyr::select(m_tmax, f_tmax, year) %>% 
  bind_rows(.,bc_proj)

pred <- predict(m, newdata = ndf)
pred_df <- tibble(bloom_doy = pred,
                  year = ndf$year) 

ggplot() +
  geom_point(data = bc_blooms, aes(y = bloom_doy, x = year)) +
  geom_line(data = pred_df, aes(y = bloom_doy, x = year))

# write out
bc_proj <- pred_df %>% 
  filter(year > 2021) %>% 
  mutate(location = "vancouver, bc",
         bloom_doy = round(bloom_doy)) %>% 
  dplyr::select(bloom_doy, year, location)
write.csv(bc_proj, file = here::here("data/bc_projection.csv"), row.names = F)
