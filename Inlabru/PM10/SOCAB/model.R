library(dplyr)
library(ggplot2)
library(sf)
library(sp)
library(inlabru)


# The raw data downloaded from EPA
PM10s_CA_raw <- readRDS("./PM10/PM10s_CA_raw.rds")

# The data after filtering out the extreme events
PM10s_CA_selected <- dplyr::select(PM10s_CA_raw, c("site_number",
                                                   "poc",
                                                   "latitude",
                                                   "longitude",
                                                   "datum",
                                                   "year",
                                                   "event_type",
                                                   "arithmetic_mean")) %>% 
  filter(event_type %in% c("No Events",
                           "Events Excluded", 
                           "Concurred Events Excluded"))

# Preprocess the data and border of SOCAB ==========================================================

# Project the sites to UTM CRS ---------------------------------------------------------------------
# Locations with CRS NAD83
PM10s_nad83 <- filter(PM10s_CA_selected, datum == "NAD83")
loc_nad83 <- dplyr::select(PM10s_nad83, c("longitude", "latitude"))

# Locations with CRS WGS84
PM10s_wgs84 <- filter(PM10s_CA_selected, datum == "WGS84")
loc_wgs84 <- dplyr::select(PM10s_wgs84, c("longitude", "latitude"))

stopifnot(
  dim(loc_nad83)[1] + dim(loc_wgs84)[1] == dim(PM10s_CA_selected)[1]
)

# Project the map into UTM zone 10 with unit kilometer
crs_utm_km <- CRS("+proj=utm +zone=10 + ellps=WGS84 +units=km")

loc_nad83_spt <- SpatialPoints(loc_nad83, proj4string=CRS("+init=epsg:4296"))
loc_nad83_to_utm <- spTransform(loc_nad83_spt, crs_utm_km) %>% coordinates()
loc_wgs84_spt <- SpatialPoints(loc_wgs84, proj4string=CRS("+init=epsg:4326"))
loc_wgs84_to_utm <- spTransform(loc_wgs84_spt, crs_utm_km) %>% coordinates()

PM10s_nad83_utm <- dplyr::mutate(PM10s_nad83, "N" = loc_nad83_to_utm[, 2],
                                 "E" = loc_nad83_to_utm[, 1])
PM10s_wgs84_utm <- dplyr::mutate(PM10s_wgs84, "N" = loc_wgs84_to_utm[, 2],
                                 "E" = loc_wgs84_to_utm[, 1])

PM10s_CA_utm <- bind_rows(PM10s_nad83_utm, PM10s_wgs84_utm) %>% 
  dplyr::select(-c("latitude", "longitude", "datum"))


# Load the border of SOCAB (South Coast Air Basin) -------------------------------------------------
CA_AirBasins <- st_read(dsn = "./Data/Air_Basins_SCAG_Region", layer = "Air_Basins_SCAG_Region")
SOCAB_border <- CA_AirBasins[CA_AirBasins$AirBasins == "South Coast Air Basin",] %>%
  st_transform(crs = crs_utm_km) %>%
  # st_polygonize() %>%
  st_union() %>%
  st_geometry()

# Filter the sites in SOCAB ------------------------------------------------------------------------

site_locs <- data.frame("N" = PM10s_CA_utm$N, "E" = PM10s_CA_utm$E) %>%
  st_as_sf(coords = c("E", "N"), crs = crs_utm_km)

xy_in <- st_contains(SOCAB_border, site_locs) %>% as.matrix() %>% c()

PM10s_SOCAB_utm <- PM10s_CA_utm[xy_in, ]

# Convert the border to 'SP' object, which is accepted object for gg function in inlabru
SOCAB_border <- as(SOCAB_border, "Spatial")
ggplot(PM10s_SOCAB_utm) + gg(SOCAB_border) + geom_point(aes(x = E, y = N), color = "blue") + coord_equal()


