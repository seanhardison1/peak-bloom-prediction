library(tidyverse)
library(sf)

japan <-read.csv(here::here("data/japan.csv"))
kyoto <-read.csv(here::here("data/kyoto.csv"))
meteoswiss <- read.csv(here::here("data/meteoswiss.csv"))
sk <- read.csv(here::here("data/south_korea.csv"))
usa_npn_ind <- read.csv(here::here("data/USA-NPN_individual_phenometrics_data.csv"))
usa_npn_ind_desc <- read.csv(here::here("data/USA-NPN_individual_phenometrics_datafield_descriptions.csv"))
usa_npn_sint <- read.csv(here::here("data/USA-NPN_status_intensity_datafield_descriptions.csv"))
usa_npb_sint_obs <- read.csv(here::here("data/USA-NPN_status_intensity_observations_data.csv"))
dc <- read.csv(here::here("data/washingtondc.csv"))

save(japan, kyoto, meteoswiss, sk, usa_npn_ind, 
     usa_npn_ind_desc, usa_npn_sint, usa_npb_sint_obs,
     dc, file = here::here("data/all_cb_data.rdata"))
