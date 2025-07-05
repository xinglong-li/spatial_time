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
  key   = "sandbird73"
)

o3_ca_2024 <- aqs_sampledata(
  aqs_filter   = "byState",
  aqs_variables = list(
    state     = "06",       # California FIPS
    param     = "44201",    # O₃
    bdate     = "20240101",
    edate     = "20241231"
  )
) 

o3_ca_2025 <- aqs_sampledata(
  aqs_filter   = "byState",
  aqs_variables = list(
    state     = "06",       # California FIPS
    param     = "44201",    # O₃
    bdate     = "20250101",
    edate     = "20250531"
  )
) 

Ozone <- o3_ca_2025 %>% dplyr::select(site_number, latitude, longitude, 
                                      datum, date_local, time_local, 
                                      sample_measurement)

head(Ozone)

# Some information of the data
print(sprintf("Total number of records: %s", dim(Ozone)[1]))
print(sprintf("Total number of site_numbers in California: %s", length(unique(Ozone$site_number))))

# Convert coordinates of the data ==================================================================
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
Ozone_nad83 <- filter(Ozone, datum == "NAD83")
loc_nad83 <- dplyr::select(Ozone_nad83, c("longitude", "latitude"))

# Locations with CRS WGS84
Ozone_wgs84 <- filter(Ozone, datum == "WGS84")
loc_wgs84 <- dplyr::select(Ozone_wgs84, c("longitude", "latitude"))

stopifnot(
  dim(loc_nad83)[1] + dim(loc_wgs84)[1] == dim(Ozone)[1]
)

# Project the map into UTM zone 10 with unit kilometer
km_proj <- CRS("+proj=utm +zone=10 + ellps=WGS84 +units=km")

loc_nad83_spt <- SpatialPoints(loc_nad83, proj4string=CRS("+init=epsg:4296"))
loc_nad83_to_utm <- spTransform(loc_nad83_spt, km_proj) %>% coordinates()
loc_wgs84_spt <- SpatialPoints(loc_wgs84, proj4string=CRS("+init=epsg:4326"))
loc_wgs84_to_utm <- spTransform(loc_wgs84_spt, km_proj) %>% coordinates()

Ozone_nad83_utm <- dplyr::mutate(Ozone_nad83, "N" = loc_nad83_to_utm[, 2],
                                "E" = loc_nad83_to_utm[, 1])
Ozone_wgs84_utm <- dplyr::mutate(Ozone_wgs84, "N" = loc_wgs84_to_utm[, 2],
                                "E" = loc_wgs84_to_utm[, 1])

Ozone_utm <- bind_rows(Ozone_nad83_utm, Ozone_wgs84_utm) %>% 
  dplyr::select(-c("latitude", "longitude", "datum"))


# SOCAB border

SOCAB_border <- readRDS("./PM10/SOCAB/SOCAB_border.rds")

# Filter the sites in SOCAB ------------------------------------------------------------------------
site_locs <- data.frame("N" = Ozone_utm$N, "E" = Ozone_utm$E) %>%
  st_as_sf(coords = c("E", "N"), crs=km_proj)

xy_in <- st_contains(SOCAB_border, site_locs) %>% as.matrix() %>% c()

Ozone_SOCAB_utm <- Ozone_utm[xy_in, ]

ggplot(Ozone_SOCAB_utm) + geom_sf(data = SOCAB_border) +
  geom_point(aes(x = E, y = N), color = "blue") + 
  xlab('East / km') +
  ylab('North / km') +
  coord_sf()

# # Aggregate by site --------------------------------------------------------------------------------
# # (Should we aggregate monitors belong to the same sites together?)
# Ozone_SOCAB_summary <- group_by(Ozone_SOCAB_utm, 
#                                 site_number, date_local, time_local) %>%
#   summarise(hourly_mean = mean(sample_measurement),
#             north = mean(N),
#             east = mean(E)) %>%
#   mutate(north = mean(north[site_number == site_number]),
#          east = mean(east[site_number == site_number]))

Ozone_SOCAB_summary <- Ozone_SOCAB_utm %>% rename(north = N, east = E)

ggplot(Ozone_SOCAB_summary) + geom_sf(data = SOCAB_border) +
  geom_point(aes(x = east, y = north), color = "blue") + 
  xlab('East / km') +
  ylab('North / km') +
  coord_sf()

# Rescale the coordinates of sites and border, the new unit is 10km ===============================

sd_north <- sd(Ozone_SOCAB_summary$north) 
sd_east <- sd(Ozone_SOCAB_summary$east) 

Ozone_SOCAB_scaled <- Ozone_SOCAB_summary
Ozone_SOCAB_scaled$north <- Ozone_SOCAB_summary$north / 10
Ozone_SOCAB_scaled$east <- Ozone_SOCAB_summary$east / 10

SOCAB_border_scaled <- (SOCAB_border / matrix(data = c(10, 10), ncol = 2)) %>%
  as("Spatial")

saveRDS(Ozone_SOCAB_scaled, "Ozone/Ozone_SOCAB_scaled.rds")
saveRDS(SOCAB_border_scaled, "Ozone/SOCAB_border_scaled.rds")


