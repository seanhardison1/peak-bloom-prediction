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

ky_temps <- 
  read_excel(here::here('data/kyoto_temps.xlsx')) %>% 
  dplyr::select(year = Year,
                mar = Mar,
                feb = Feb,
                jmax = Jan_max,
                fmax = Feb_max,
                mmax = Mar_max,
                jpre = Jan_precip,
                fpre = Feb_precip,
                mpre = Mar_precip)

#load data
load(here::here("data/all_cb_data.rdata"))

# japan 
jp_sf <- ne_countries(country  = "japan",
                      returnclass = "sf",
                      scale = "medium") %>% 
  st_transform(st_crs("+proj=utm +zone=54 +datum=WGS84 +units=km +no_defs")) %>% 
  st_union() 

ky_bf <- japan %>% 
  filter(str_detect(location, "Kyoto")) %>% 
  dplyr::select(long, lat) %>% 
  distinct() %>% 
  st_as_sf(.,coords = c("long","lat"), 
           crs = 4326) %>% 
  st_transform(st_crs("+proj=utm +zone=54 +datum=WGS84 +units=km +no_defs")) %>% 
  st_buffer(.,dist = 200)

# japan first
jp <- 
  japan %>% 
  mutate(bloom_date = as.Date(bloom_date)) %>%
  st_as_sf(.,coords = c("long","lat"), 
           crs = 4326) %>% 
  st_transform(st_crs("+proj=utm +zone=54 +datum=WGS84 +units=km +no_defs")) %>% 
  st_intersection(.,ky_bf) %>% 
  dream::sfc_as_cols(names = c("longitude","latitude")) %>% 
  st_set_geometry(NULL)

y <- 1950
jp_lats <- jp %>% 
  filter(year >= y) %>% 
  group_by(location) %>% 
  summarise(longitude = mean(longitude),
            latitude = mean(latitude))


jp_sample <- jp %>% 
  filter(year >= y) %>% 
  group_by(location, bloom_date) %>% 
  dplyr::summarise(bloom_doy = mean(bloom_doy, na.rm = T)) %>% 
  left_join(.,jp_lats) %>% 
  mutate(year = year(bloom_date)) %>% 
  tsibble(key = "location", index = "year") %>% 
  fill_gaps()

sites <- 
  jp_sample %>% 
  as_tibble() %>% 
  group_by(location) %>% 
  dplyr::summarise(n = n()) %>% 
  filter(n == 69) %>% 
  pull(location)

jp_sample2 <- jp_sample %>% 
  filter(location %in% sites) %>% 
  dplyr::select(location, bloom_doy, year,
                latitude, longitude) %>% 
  mutate(bloom_doy_l1 = lag(bloom_doy),
         bloom_doy_l2 = lag(bloom_doy, 2)) %>% 
  # filter(location == "Japan/Kyoto") %>% 
  left_join(.,ky_temps) %>% 
  mutate(lag_mar = lag(mar),
         lag_mmax = lag(mmax))

m <- gam(bloom_doy ~ 
           s(fpre) + 
           s(fmax) +
           s(jmax) +
           s(feb) +
           te(longitude, latitude, year), 
         data = jp_sample2)
draw(m)
acf(m$residuals)
summary(m)
appraise(m)
ndf <- tibble(longitude = 19.2106,
              latitude = 3887.377,
              year = 1953:2022,
              fpre = c(jp_sample2 %>% filter(location == "Japan/Kyoto") %>% pull(fpre),17.0),
              fmax = c(jp_sample2 %>% filter(location == "Japan/Kyoto") %>% pull(fmax),9.0),
              feb = c(jp_sample2 %>% filter(location == "Japan/Kyoto") %>% pull(feb),4.2),
              jmax = c(jp_sample2 %>% filter(location == "Japan/Kyoto") %>% pull(jmax),8.2))
pred_df <- bind_cols(jp_sample2 %>% 
                       filter(str_detect(location, "Kyoto")) %>% 
                       add_row(year = 2022),
                     predict(m, newdata = ndf,
                             se.fit = T))
ggplot(pred_df) +
  geom_point(aes(x = year, y = bloom_doy)) +
  geom_line(aes(x = year, y = fit))


pred_df %>% filter(year == 2022) %>% pull(fit)
