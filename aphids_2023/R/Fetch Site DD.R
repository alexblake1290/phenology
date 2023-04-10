#### Setup ####

# Load libraries
library(tidyverse)
library(lubridate)
library(mgcv)
library(leaflet)
source("./aphids_2023/src/cdd_density_functions.R") #used for download_daymet_v()

# Get data (site lat/long)
# Original (pan trap) sites used for training the GAM, here: https://resjournals.onlinelibrary.wiley.com/doi/full/10.1111/afe.12564

# Sites * location isn't actually distinct, lat/long keep jittering. Use averages
sites = read.csv("./aphids_2023/osf/Aggregated Aphid Data/Idaho pan trap with degree days.csv") %>%
  mutate(
    Site = str_remove(SiteName, "[:digit:]{4}-")
  ) %>%
  select(Site,Latitude,Longitude) %>%
  distinct() %>%
  group_by(Site) %>%
  summarize(Lat = mean(Latitude), Long = mean(Longitude))

# New sites where we want predictions (Vetch/sweep dataset)
# duplicate locations with different names
dups = c("Hells Gate State Park", "Nisqually Canyon", "Blyton Landing", "Chris", "Wawawai", "Golf")
newsites = read.csv("./aphids_2023/osf/Aggregated Aphid Data/Vetch 2019 and 2020 with degree days.csv") %>%
  rename(
    Site = SiteID
  ) %>%
  filter(!grepl("^[0-9][A-z]$|^[A-z][0-9]$",Site)) %>%
  select(Site,Lat=Latitude,Long=Longitude) %>%
  filter(!(Site %in% dups)) %>%
  distinct()

# comment out to get predictions on the training sites instead of the new sites
sites = newsites

# Set aphid DD min and max in F
lb = 41.9; ub = 82.4

# Set date range
end = 2022; start = 2000
ref = as.Date("2023-01-01")

# Set upper, mean, and lower DD bounds for aphid peak
# Extracted from Pan Trap m1 (data aggregated across place and time)
#peak = data_frame(cdd = c(1201,1338,1452))
# Extracted from sweeps data, upper bound is crap due to sparse data at upper extreme
 peak = data_frame(cdd = c(1578, 1654, 1771))

# Confidence interval for timing predictions. 1.645 for 90%, 1.96 for 95%, 2.58 for 99%
cilvl = 2.58

#### Fetch historical temperature data ####
# dat should only contain a label/name column (Site), and Lat and Long columns, one row per site

dat = sites %>%
  select(-Site) %>%
  mutate(cl = pmap(., .f = function(Lat, Long){
    # Fetch historical temp data by site and year from daymet
    download_daymet_v(lat = Lat, lon = Long, start=start, end=end, vars="tmax,tmin", internal=T, force=F) %>%
    .$data %>%
    rename(
      tmax = `tmax..deg.c.`,
      tmin = `tmin..deg.c.`
    )
  }
  # Reattach site names, calculate CDD and julian days
  )) %>%
  right_join(sites,by=c("Lat","Long")) %>%
  unnest(cl) %>%
  mutate(tavg = ((tmax + tmin)/2) * (9/5) + 32) %>%
  rename(julian_day = yday) %>%
  # select(-tmax, -tmin) %>%
  mutate(tavg = if_else(tavg < lb | tavg > ub, lb, tavg),
         dd = tavg - lb) %>%
  group_by(Site, Lat, Long, year) %>%
  mutate(cdd = cumsum(dd)) %>%
  select(-tavg, -dd) %>%
  ungroup()


#### Intermediary visualizations ####

# Visualise the relationship by site
dat %>%
  ggplot(aes(y=julian_day,x=cdd)) +
    geom_smooth(size=.8) +
    geom_vline(xintercept=1339,color="red") +
    facet_wrap(~Site)

# Visualise tmax and tmin all sites
dat %>%
  ggplot() +
    geom_smooth(aes(x=julian_day,y=tmin,group=Site),color="lightblue", size=.1, se=F) +
    geom_smooth(aes(x=julian_day,y=tmax,group=Site),color="red", size=.1, se=F)

# Visualise tmax and tmin all years (one site)
dat %>%
  filter(Site == "Alpowa Creek") %>%
  ggplot() +
    geom_smooth(aes(x=julian_day,y=tmin,group=year,color=year),size=.2, se=F) +
    geom_smooth(aes(x=julian_day,y=tmax,group=year,color=year),size=.2, se=F) +
    scale_color_gradient(low="yellow",high="red", na.value=NA)

# Visualise CDD all sites
dat %>%
  ggplot(aes(x=julian_day,y=cdd,group=Site)) +
    geom_smooth(size=.1,se=F) +
    geom_hline(yintercept=1339,color="red")

# Visualise CDD all years (one site)
dat %>%
  filter(Site == "Alpowa Creek") %>%
    ggplot() +
    geom_point(aes(x=julian_day,y=cdd,group=year,color=year),size=0.2,alpha=0.5) +
    scale_color_gradient(low="yellow",high="red", na.value=NA) +
    geom_hline(yintercept=1339,color="black")

# Map sites to validate coordinates
leaflet(sites) %>%
  addTiles() %>%
  addMarkers(lng= ~Long, lat= ~Lat, popup= ~Site)


#### Fit GAM for each site ####

# purrr the gam fits using mgcv
fits = dat %>%
  group_by(Site,Lat,Long) %>%
  nest() %>%
  mutate(
    mod = map(data, ~ mgcv::gam(julian_day ~ s(cdd), data= .x))
  )

# unnest all the outputs, then purrr predicted dates based on known CDDs for aphid peaks, get CIs from se
preds = fits %>%
  mutate(
    pred = map(mod, ~ predict(.x, se.fit=T, newdata=peak))
  ) %>%
  select(-mod,-data) %>%
  unnest_wider(pred) %>%
  unnest_wider(fit) %>%
  rename(
    julian_peak_lower = "1",
    julian_peak_median = "2",
    julian_peak_upper = "3"
  ) %>%
  unnest_wider(se.fit) %>%
  rename(
    ci_lower = "1",
    ci_median = "2",
    ci_upper = "3"
  ) %>%
  mutate(
    ci_lower = ci_lower*cilvl,
    ci_median = ci_median*cilvl,
    ci_upper = ci_upper*cilvl,
    date_lower = as.character(as_date(julian_peak_lower, origin=ref)),
    date_median = as.character(as_date(julian_peak_median, origin=ref)),
    date_upper = as.character(as_date(julian_peak_upper, origin=ref))
  ) %>%
  mutate(
    across(where(is.array), as.numeric) # avoid a column of 1D arrays from the gam outputs, which causes write_csv issues
  ) %>%
  ungroup()

# Rename according to training data
write_csv(preds,"./aphids_2023/outputs/peak_predictions_2023_sweeps.csv")
