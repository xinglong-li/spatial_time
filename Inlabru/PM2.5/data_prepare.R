library(jsonlite)
library(readr)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(sf)
library(sp)
library(raster)

# rely on the R package "raqs" to automatically download data from the EPA website
# this is equivalent to the previous method using the API in the EPA website:
# https://aqs.epa.gov/aqsweb/documents/data_api.html#param

library(raqs)
# setup the account for access to the EPA data
# store your AQS credentials once per session
set_aqs_user(
  email = "xinglong.li@stat.ubc.ca",
  key   = "khakicrane79"
)

# ---- Whole-history pull: 1999 (first PM2.5 year) to present ----
ca_pm25_annual <- aqs_annualdata(                     # main “annual” service
  aqs_filter = "byState",                             # grab *all* CA sites
  aqs_variables = list(
    param = "88101",         # PM2.5 – Local Conditions
    bdate = "19990101",      # begin date – only the *year* matters
    edate = "20251231",      # end date  – adjust as needed
    state = "06"             # CA FIPS
  )
)

# peek
names(ca_pm25_annual)
head(ca_pm25_annual[c("year", "site_number", "arithmetic_mean")])

# Some information of the data
print(sprintf("Total number of records: %s", dim(ca_pm25_annual)[1]))
print(sprintf("Total number of site_numbers in California: %s", length(unique(ca_pm25_annual$site_number))))
print(sprintf("Total number of local_site names: %s", length(unique(ca_pm25_annual$local_site_name))))


# Prepare the data for model fitting ===============================================================

# Descriptions of data each column is in 
# https://aqs.epa.gov/aqsweb/airdata/FileFormats.html#_annual_summary_files

# Remove Extreme events, there are sites with only one record "Concurred Events Excluded",
# and the total number of sites will reduce if we exclude such records. 
# So we keep such records at the moments.

PM25 <- dplyr::select(ca_pm25_annual, 
                      c("site_number", "poc", "latitude", "longitude", "datum", "year", "arithmetic_mean"))


# Convert coordinates of the data ==================================================================

# It is strange to see that some 'site_number's has multiple 'local_site_name's, and each local name
# has different latitude and longitude.
# But some site has no local_site_name.
# For example, site 0007 in year 2021 has 2 local_site_names, even stranger, the poc #1 of this site
# has 2 local_site_name and they have different locations.

# The approach:
# Take the mean value of the PM10 arithmetic means over all pocs of one site;
# Take the mean value of all locations of one site, but since the locations are measured under different
# coordinate reference system (CRS), we firstly convert all positions into same projected CRS.

# The CRS used in the dataset:
# Datum: WGS84 (EPSG:4326)
# Datum: NAD83 (EPSG:4296)
# By some referece: Transforming between NAD83 and WGS84 typically isn’t recommended for most applications
# because standard transformations can introduce error that is large relative to their difference. 
# Complicating the matter, the difference between NAD83 and WGS84 varies with time and location.
# Both systems have frequent new realizations due to more data and improved techniques.

# We want to convert all locations into Projected CRS (Easting/Northing) coordinates:
# Datum: UTM, Zone 10 (EPSG: 32610) // Zone 10 is used in the Pacific Northwest

# Locations with CRS NAD83
PM25_nad83 <- filter(PM25, datum == "NAD83")
loc_nad83 <- dplyr::select(PM25_nad83, c("longitude", "latitude"))

# Locations with CRS WGS84
PM25_wgs84 <- filter(PM25, datum == "WGS84")
loc_wgs84 <- dplyr::select(PM25_wgs84, c("longitude", "latitude"))

stopifnot(
  dim(loc_nad83)[1] + dim(loc_wgs84)[1] == dim(PM25)[1]
)

# Project the map into UTM zone 10 with unit kilometer
km_proj <- CRS("+proj=utm +zone=10 + ellps=WGS84 +units=km")

loc_nad83_spt <- SpatialPoints(loc_nad83, proj4string=CRS("+init=epsg:4296"))
loc_nad83_to_utm <- spTransform(loc_nad83_spt, km_proj) %>% coordinates()
loc_wgs84_spt <- SpatialPoints(loc_wgs84, proj4string=CRS("+init=epsg:4326"))
loc_wgs84_to_utm <- spTransform(loc_wgs84_spt, km_proj) %>% coordinates()

PM25_nad83_utm <- dplyr::mutate(PM25_nad83, "N" = loc_nad83_to_utm[, 2],
                                "E" = loc_nad83_to_utm[, 1])
PM25_wgs84_utm <- dplyr::mutate(PM25_wgs84, "N" = loc_wgs84_to_utm[, 2],
                                "E" = loc_wgs84_to_utm[, 1])

PM25_utm <- bind_rows(PM25_nad83_utm, PM25_wgs84_utm) %>% 
  dplyr::select(-c("latitude", "longitude", "datum"))


# SOCAB border

SOCAB_border <- readRDS("./PM10/SOCAB/SOCAB_border.rds")

# Filter the sites in SOCAB ------------------------------------------------------------------------
site_locs <- data.frame("N" = PM25_utm$N, "E" = PM25_utm$E) %>%
  st_as_sf(coords = c("E", "N"), crs=km_proj)

xy_in <- st_contains(SOCAB_border, site_locs) %>% as.matrix() %>% c()

PM25_SOCAB_utm <- PM25_utm[xy_in, ]

ggplot(PM25_SOCAB_utm) + geom_sf(data = SOCAB_border) +
  geom_point(aes(x = E, y = N), color = "blue") + 
  xlab('East / km') +
  ylab('North / km') +
  coord_sf()

# Aggregate by site --------------------------------------------------------------------------------
# (Should we aggregate monitors belong to the same sites together?)
PM25_SOCAB_summary <- group_by(PM25_SOCAB_utm, site_number, year) %>%
  summarise(annual_mean = mean(arithmetic_mean),
            north = mean(N),
            east = mean(E)) %>%
  mutate(north = mean(north[site_number == site_number]),
         east = mean(east[site_number == site_number]))

ggplot(PM25_SOCAB_summary) + geom_sf(data = SOCAB_border) +
  geom_point(aes(x = east, y = north), color = "blue") + 
  xlab('East / km') +
  ylab('North / km') +
  coord_sf()


# Import the border map of California ==============================================================

# # The border file of California downloaded have associated projection data,such as ESRI shapefiles. 
# # ATTENTION: function readOGR from package rgdal is no longer available, should use read_sf from package sf instead 
# ca_boundary <- readOGR(dsn = "~/Downloads/ca-state-boundary", layer = "CA_State_TIGER2016")
# ca_boundary_km <- spTransform(ca_boundary, km_proj)
# proj4string(ca_boundary)
# plot(ca_boundary)
# saveRDS(ca_boundary_km, sprintf("%sCA_border.rds", result_path))
# 
# ca_boundary_km <- readRDS("PM2.5/CA_border.rds")

# Rescale the coordinates of sites and border, the new unit is 10km ===============================

sd_north <- sd(PM25_SOCAB_summary$north) 
sd_east <- sd(PM25_SOCAB_summary$east) 

PM25_SOCAB_scaled <- PM25_SOCAB_summary
PM25_SOCAB_scaled$north <- PM25_SOCAB_summary$north / 10
PM25_SOCAB_scaled$east <- PM25_SOCAB_summary$east / 10

SOCAB_border_scaled <- (SOCAB_border / matrix(data = c(10, 10), ncol = 2)) %>%
  as("Spatial")

ggplot(PM25_SOCAB_scaled) + gg(SOCAB_border_scaled) + 
  geom_point(aes(x = east, y = north), color = "blue") + 
  coord_sf()


saveRDS(PM25_SOCAB_scaled, "PM2.5/PM25_SOCAB_scaled.rds")
saveRDS(SOCAB_border_scaled, "PM2.5/SOCAB_border_scaled.rds")


