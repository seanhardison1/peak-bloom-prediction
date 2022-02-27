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

fc <- F

#load data
load(here::here("data/all_cb_data.rdata"))

bas_temps <- 
  read_excel(here::here("data/basel_temps2.xlsx")) %>% 
  dplyr::rename(temp = 2,
                gdd = 3,
                precip = 4,
                soil_moist = 5,
                wind_speed = 5,
                wind_dir = 6)

y <- str_sub(bas_temps$timestamp,1,4)
m <- str_sub(bas_temps$timestamp,5,6)
d <- str_sub(bas_temps$timestamp,7,8)
hr <- str_sub(bas_temps$timestamp,9,13)

bas_temps2 <- bas_temps %>% 
  mutate(date = as.Date(paste(y, m, d, sep = "-")))

mmtemps <- bas_temps2 %>% 
  group_by(date) %>% 
  dplyr::summarise(tmax = max(temp),
                   tmin = min(temp)) %>% 
  group_by(ymon = yearmonth(date)) %>% 
  dplyr::summarise(m_tmax = (mean(tmax) - 32)*(5/9),
                   m_tmin = (mean(tmin) - 32)*(5/9))

bas_temps3 <- bas_temps2 %>% 
  group_by(ymon = yearmonth(date)) %>% 
  dplyr::summarise(m_temp = mean(temp),
                   m_precip = mean(precip)) %>% 
  left_join(.,mmtemps) %>%  
  filter(month(ymon) %in% 1:3) %>% 
  mutate(year = year(ymon),
         month = month(ymon)) %>%
  dplyr::select(-ymon)

jan <- bas_temps3 %>% 
  filter(month == 1) %>% 
  rename_at(vars(1:4), function(x)paste0("j_",x)) %>% 
  dplyr::select(-month)

feb <- bas_temps3 %>% 
  filter(month == 2) %>% 
  rename_at(vars(1:4), function(x)paste0("f_",x)) %>% 
  dplyr::select(-month)

mar <- bas_temps3 %>% 
  filter(month == 3) %>% 
  rename_at(vars(1:4), function(x)paste0("m_",x)) %>% 
  dplyr::select(-month)

bas_temps4 <- left_join(feb,
                        jan) %>% 
  left_join(.,mar)

# switzerland
sw_sf <- ne_countries(country  = "switzerland",
                      returnclass = "sf",
                      scale = "medium") %>% 
  st_transform(st_crs("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")) %>% 
  st_union() 

sw_bf <- meteoswiss %>% 
  filter(str_detect(location, "Switzerland/Liestal")) %>% 
  dplyr::select(long, lat) %>% 
  distinct() %>% 
  st_as_sf(.,coords = c("long","lat"), 
           crs = 4326) %>% 
  st_transform(st_crs("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")) %>% 
  st_buffer(.,dist = 200) 

# sw first
sw <- 
  meteoswiss %>% 
  mutate(bloom_date = as.Date(bloom_date)) %>%
  st_as_sf(.,coords = c("long","lat"), 
           crs = 4326) %>% 
  st_transform(st_crs("+proj=utm +zone=32 +datum=WGS84 +units=km +no_defs")) %>% 
  dream::sfc_as_cols(names = c("longitude","latitude")) %>% 
  st_set_geometry(NULL)

y <- 1984
sw_lats <- sw %>% 
  filter(year >= y) %>%
  group_by(location) %>% 
  summarise(longitude = mean(longitude),
            latitude = mean(latitude))

sw_sample <- sw %>% 
  filter(year >= y) %>%
  group_by(location, bloom_date) %>% 
  dplyr::summarise(bloom_doy = mean(bloom_doy, na.rm = T)) %>% 
  left_join(.,sw_lats) %>% 
  mutate(year = year(bloom_date)) 

sw_sample2 <- sw_sample %>% 
  dplyr::select(location, bloom_doy, year,
                latitude, longitude) %>%
  tsibble(index = "year", key = "location") %>% 
  fill_gaps() %>% 
  mutate(bloom_doy_l1 = lag(bloom_doy),
         bloom_doy_l2 = lag(bloom_doy, 2)) %>% 
  left_join(.,bas_temps4) %>%
  na.exclude()

m <- gam(bloom_doy ~
           # s(year, bs = "fs") +
           s(bloom_doy_l1) +
           s(bloom_doy_l2) +
           s(f_m_precip) +
           s(f_m_tmax) +
           s(j_m_tmax) +
           s(f_m_temp) +
           te(longitude, latitude, year),
         data = sw_sample2)
# plot(m)
# draw(m)
acf(m$residuals)
summary(m)
appraise(m)
y <- 1988
pred_df1 <- sw_sample2 %>% 
  filter(location == "Switzerland/Liestal",
         year >= y) 

prec_2022 <- bas_temps4 %>% filter(year == 2022) %>% pull(f_m_precip)
tmax_2022 <- bas_temps4 %>% filter(year == 2022) %>% pull(f_m_tmax)
temp_2022 <- bas_temps4 %>% filter(year == 2022) %>% pull(f_m_temp)
jmax_2022 <- bas_temps4 %>% filter(year == 2022) %>% pull(j_m_tmax)


ndf <- tibble(longitude = 404.3579,
              latitude = 5259.444,
              year = y:2022,
              bloom_doy_l1 = c(pred_df1 %>% pull(bloom_doy_l1), 87),
              bloom_doy_l2 = c(pred_df1 %>% pull(bloom_doy_l2), 79),
              f_m_precip = c(pred_df1 %>% pull(f_m_precip), prec_2022),
              f_m_tmax = c(pred_df1 %>% pull(f_m_tmax), tmax_2022),
              f_m_temp = c(pred_df1 %>% pull(f_m_temp), temp_2022),
              j_m_tmax = c(pred_df1 %>% pull(j_m_tmax), jmax_2022))

if (fc){
  base <- ndf %>% tsibble(index = "year")
  
  output_bdl1 <- 
    base %>% 
    model(
      bloom_doy_l1 = NNETAR(bloom_doy_l1)
    ) %>% 
    forecast(h = 10) %>% 
    tibble()
  
  output_fmpr <- 
    base %>% 
    model(
      f_m_precip = NNETAR(f_m_precip)
    ) %>% 
    forecast(h = 10) %>% 
    tibble()
  
  output_jmtmax <- 
    base %>% 
    model(
      j_m_tmax = NNETAR(j_m_tmax)
    ) %>% 
    forecast(h = 10) %>% 
    tibble()
  
  output_fmtmax <- 
    base %>% 
    model(
      f_m_tmax = NNETAR(f_m_tmax)
    ) %>% 
    forecast(h = 10) %>% 
    tibble()
  
  output_fmtem <- 
    base %>% 
    model(
      f_m_temp = NNETAR(f_m_temp)
    ) %>% 
    forecast(h = 10) %>% 
    tibble()
  
  sw_proj <- 
    output_bdl1 %>% 
    dplyr::select(year, bloom_doy_l1 = .mean) %>% 
    left_join(.,output_fmpr %>% 
                dplyr::select(year, f_m_precip = .mean)) %>% 
    left_join(.,output_jmtmax %>% 
                dplyr::select(year, j_m_tmax = .mean)) %>% 
    left_join(.,output_fmtmax %>% 
                dplyr::select(year, f_m_tmax = .mean)) %>% 
    left_join(.,output_fmtem %>% 
                dplyr::select(year, f_m_temp = .mean)) 
    
  save(sw_proj, file = here::here("data/sw_env_fc.rdata"))
} else {
  load(here::here("data/sw_env_fc.rdata"))
}


liestal <-
  sw_sample2 %>% 
  filter(str_detect(location, "Switzerland/Liestal"),
         year >= y)

ndf2 <- ndf %>% 
  bind_rows(
  sw_proj %>% 
    mutate(bloom_doy_l2 = lag(bloom_doy_l1),
         latitude = unique(ndf$latitude),
         longitude = unique(ndf$longitude)) 
  ) %>% 
  mutate(bloom_doy_l2 = ifelse(year == 2023, 87, bloom_doy_l2))

pred <- 
  predict(m, newdata = ndf2, se.fit = T)

pred_df <- tibble(bloom_doy = pred$fit,
                  year = ndf2$year) 

ggplot(pred_df) +
  geom_point(aes(x = year, y = bloom_doy))

# write out
sw_proj <- pred_df %>% 
  filter(year > 2021) %>% 
  mutate(location = "liestal, ch",
         bloom_doy = round(bloom_doy)) %>% 
  dplyr::select(bloom_doy, year, location)
write.csv(sw_proj, file = here::here("data/sw_projection.csv"), row.names = F)
