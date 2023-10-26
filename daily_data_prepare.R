library(dplyr)
library(readr)
library(sp)
library(ggplot2)


# The daily records of PM10 are downloaded from 
# https://aqs.epa.gov/aqsweb/airdata/download_files.html
# ONLY keep records from 1990 to 2022, since records in 1980 - 1989 are too small, for example,
# the file of 1980 and 1981 contains 0 row?

data_path <- "/home/xinglong/Downloads/daily_summary/"
result_path <- "/home/xinglong/git_local/spatial_time/PM10/"


years <- seq(1990, 2022)
summary_daily <- NULL

for (year in years) {
  record <- read_csv(sprintf("%sdaily_81102_%s.csv", data_path, year))
  summary_daily <- bind_rows(summary_daily, record)
}

names(summary_daily) <- gsub(" ", "_", tolower(names(summary_daily)))

print(sprintf("Total number of records: %s", dim(summary_daily)[1]))
print(sprintf("Total number of site_numbers in California: %s", length(unique(summary_daily$site_num))))
print(sprintf("Total number of local_site names: %s", length(unique(summary_daily$local_site_name))))

saveRDS(summary_daily, sprintf("%sPM10s_CA_daily_raw.rds", result_path))


# Prepare the data for model fitting ===============================================================

# Descriptions of data each column is in 
# https://aqs.epa.gov/aqsweb/airdata/FileFormats.html#_annual_summary_files

# Remove Extreme events.
PM10s_daily <- dplyr::select(summary_daily, c("site_num",
                                              "poc",
                                              "latitude",
                                              "longitude",
                                              "datum",
                                              "date_local",
                                              "event_type",
                                              "arithmetic_mean")) %>% 
  filter(event_type %in% c("None")) %>%
  select(-c("event_type"))

# saveRDS(PM10s_daily, sprintf("%sPM10s_CA__daily_selected.rds", result_path))


# Convert coordinates of the data ==================================================================

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
PM10s_nad83 <- filter(PM10s_daily, datum == "NAD83")
loc_nad83 <- dplyr::select(PM10s_nad83, c("longitude", "latitude"))

# Locations with CRS WGS84
PM10s_wgs84 <- filter(PM10s_daily, datum == "WGS84")
loc_wgs84 <- dplyr::select(PM10s_wgs84, c("longitude", "latitude"))

stopifnot(
  dim(loc_nad83)[1] + dim(loc_wgs84)[1] == dim(PM10s_daily)[1]
)

# Project the map into UTM zone 10 with unit kilometer
km_proj <- CRS("+proj=utm + zone=10 + ellps=WGS84 +units=km")

loc_nad83_spt <- SpatialPoints(loc_nad83, proj4string=CRS("+init=epsg:4296"))
loc_nad83_to_utm <- spTransform(loc_nad83_spt, km_proj) %>% coordinates()
loc_wgs84_spt <- SpatialPoints(loc_wgs84, proj4string=CRS("+init=epsg:4326"))
loc_wgs84_to_utm <- spTransform(loc_wgs84_spt, km_proj) %>% coordinates()

PM10s_nad83_utm <- dplyr::mutate(PM10s_nad83, "N" = loc_nad83_to_utm[, 2],
                                       "E" = loc_nad83_to_utm[, 1])
PM10s_wgs84_utm <- dplyr::mutate(PM10s_wgs84, "N" = loc_wgs84_to_utm[, 2],
                                       "E" = loc_wgs84_to_utm[, 1])

PM10s_daily_utm <- bind_rows(PM10s_nad83_utm, PM10s_wgs84_utm) %>% 
  dplyr::select(-c("latitude", "longitude", "datum"))


# Aggregate by site --------------------------------------------------------------------------------

# Before aggregating by sites, it can be seen that some of the records are negative ?


PM10s_daily_summary <- group_by(PM10s_daily_utm, site_num, date_local) %>%
  summarise(daily_mean = mean(arithmetic_mean),
            north = mean(N),
            east = mean(E)) %>%
  mutate(north = mean(north[site_num == site_num]),
         east = mean(east[site_num == site_num]))

print(PM10s_daily_summary, n = 100)

# saveRDS(PM10s_daily_summary, sprintf("%sPM10s_CA_daily_summary.rds", result_path))


# Rescale the coordinates of sites and border, the new unit is 100km ===============================

sd_north <- sd(PM10s_daily_summary$north) # = 181.9436
sd_east <- sd(PM10s_daily_summary$east) # = 149.8781

PM10s_daily_summary_scaled <- PM10s_daily_summary
PM10s_daily_summary_scaled$north <- PM10s_daily_summary$north / 100
PM10s_daily_summary_scaled$east <- PM10s_daily_summary$east / 100

saveRDS(PM10s_daily_summary_scaled, sprintf("%sPM10s_CA_daily_summary_scaled.rds", result_path))


# Preprocess the data to remove outliers ===========================================================

# It can be seen from this histogram that most values are smaller than 200 or even 150.
# And there are obvious outliers.
hist(PM10s_daily_summary_scaled$daily_mean) 

# Remove the outliers?


# Prepare monthly measurements =====================================================================



# Prepare weekly measurements ======================================================================






