library(jsonlite)
library(readr)
library(rgdal)
library(dplyr)
library(ggplot2)
library(ggspatial)
library(sf)
library(sp)
library(raster)


# Read data from each year (from 1980 to 2022) and combine the PM10 records ========================

data_path <- "/home/xinglong/Downloads/annual_by_monitor/"
result_path <- "/home/xinglong/git_local/spacial_time/PM10/"

annul_records <- list.files(data_path, pattern="*.json")

PM10s_raw <- NULL

for (annual_record in annul_records) {
  record <- read_json(sprintf("%s%s", data_path, annual_record))$Data %>%
    bind_rows()
  PM10s_raw <- bind_rows(PM10s_raw, record)
}


# Check data source ================================================================================

# Make sure that the time is from 1985 --- 2022
stopifnot(
  length(
    setdiff(unique(PM10s_raw$year), seq(1985, 2022))
    ) == 0
  )

# Make sure that all measurements is from California
stopifnot(
  all(PM10s_raw$state == 'California')
)

# Make sure that all measurements is PM10
stopifnot(
  all(PM10s_raw$parameter == "PM10 Total 0-10um STP")
)

# Make sure that all records are measured in the same unit
stopifnot(
  length(unique(PM10s_raw$units_of_measure)) == 1
)

# Some information of the data
print(sprintf("Total number of records: %s", dim(PM10s_raw)[1]))
print(sprintf("Total number of site_numbers in California: %s", length(unique(PM10s_raw$site_number))))
print(sprintf("Total number of local_site names: %s", length(unique(PM10s_raw$local_site_name))))

saveRDS(PM10s_raw, sprintf("%sPM10s_CA_raw.rds", result_path))


# Prepare the data for model fitting ===============================================================

# Descriptions of data each column is in 
# https://aqs.epa.gov/aqsweb/airdata/FileFormats.html#_annual_summary_files

# Remove Extreme events, there are sites with only one record "Concurred Events Excluded",
# and the total number of sites will reduce if we exclude such records. 
# So we keep such records at the moments.

PM10s <- dplyr::select(PM10s_raw, c("site_number",
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

saveRDS(PM10s, sprintf("%sPM10s_CA_selected.rds", result_path))


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
PM10s_nad83 <- filter(PM10s, datum == "NAD83")
loc_nad83 <- dplyr::select(PM10s_nad83, c("longitude", "latitude"))

# Locations with CRS WGS84
PM10s_wgs84 <- filter(PM10s, datum == "WGS84")
loc_wgs84 <- dplyr::select(PM10s_wgs84, c("longitude", "latitude"))

stopifnot(
  dim(loc_nad83)[1] + dim(loc_wgs84)[1] == dim(PM10s)[1]
)

# Project the map into UTM zone 10 with unit kilometer
km_proj <- CRS("+proj=utm +zone=10 + ellps=WGS84 +units=km")

loc_nad83_spt <- SpatialPoints(loc_nad83, proj4string=CRS("+init=epsg:4296"))
loc_nad83_to_utm <- spTransform(loc_nad83_spt, km_proj) %>% coordinates()
loc_wgs84_spt <- SpatialPoints(loc_wgs84, proj4string=CRS("+init=epsg:4326"))
loc_wgs84_to_utm <- spTransform(loc_wgs84_spt, km_proj) %>% coordinates()

PM10s_nad83_utm <- dplyr::mutate(PM10s_nad83, "N" = loc_nad83_to_utm[, 2],
                                              "E" = loc_nad83_to_utm[, 1])
PM10s_wgs84_utm <- dplyr::mutate(PM10s_wgs84, "N" = loc_wgs84_to_utm[, 2],
                                              "E" = loc_wgs84_to_utm[, 1])

PM10s_utm <- bind_rows(PM10s_nad83_utm, PM10s_wgs84_utm) %>% 
  dplyr::select(-c("latitude", "longitude", "datum"))


# Aggregate by site --------------------------------------------------------------------------------

PM10s_summary <- group_by(PM10s_utm, site_number, year) %>%
  summarise(annual_mean = mean(arithmetic_mean),
            north = mean(N),
            east = mean(E))
print(PM10s_summary, n = 100)

saveRDS(PM10s_summary, sprintf("%sPM10s_CA_summary.rds", result_path))


# Import the border map of California ==============================================================

# The border file of California downloaded have associated projection data,such as ESRI shapefiles. 

ca_boundary <- readOGR(dsn = "~/Downloads/ca-state-boundary", layer = "CA_State_TIGER2016")
ca_boundary_km <- spTransform(ca_boundary, km_proj)
proj4string(ca_boundary)
plot(ca_boundary)
saveRDS(ca_boundary_km, sprintf("%sCA_border.rds", result_path))


# Rescale the coordinates of sites and border, the new unit is 100km ===============================

sd_north <- sd(PM10s_summary$north) # = 181.9436
sd_east <- sd(PM10s_summary$east) # = 149.8781

PM10s_summary_scaled <- PM10s_summary
PM10s_summary_scaled$north <- PM10s_summary$north / 100
PM10s_summary_scaled$east <- PM10s_summary$east / 100

saveRDS(PM10s_summary_scaled, sprintf("%sPM10s_CA_summary_scaled.rds", result_path))

# Re-scale coordinates of CA boundary
#Extract coordinates and put them in a list
extractCoords <- function(sp.df)
{
  results <- list()
  for(i in 1:length(sp.df@polygons[[1]]@Polygons))
  {
    results[[i]] <- sp.df@polygons[[1]]@Polygons[[i]]@coords
  }
  results
}

vertices<-extractCoords(ca_boundary_km)

scaled.vertices <- lapply(vertices, function(x) x / 100)

#Create Polygons and Spatial Polygons Data Frame from scaled list of coordinates
Polys <- list()
for (i in 1:length(scaled.vertices)) {
  Polys[i] <- sp::Polygon(scaled.vertices[[i]])
}
Polys.plural <- sp::Polygons(Polys, ID = "0")
Polys.sp <- sp::SpatialPolygons(list(Polys.plural), proj4string = km_proj)
ca_boundary_km_scaled <- sp::SpatialPolygonsDataFrame(Polys.sp, data = ca_boundary_km@data)

saveRDS(ca_boundary_km_scaled, sprintf("%sCA_border_scaled.rds", result_path))

plot(ca_boundary_km_scaled)
points(PM10s_summary_scaled$east, PM10s_summary_scaled$north, pch=24, col='blue', cex=0.6)


# Preprocess the data to remove outliers ===========================================================

# It can be seen from this histogram that most values are smaller than 200 or even 150.
# And there are obvious outliers.
hist(PM10s_summary_scaled$annual_mean) 

# Remove the outliers?
dim(PM10s)
# Only 5 out of 8566 records are larger than 150, which means taht we can safely discard these points.
print(PM10s[PM10s$arithmetic_mean > 150, ], n=100)
# These 5 extreme points come from only 2 sites:
# The site 0030 has 22 records in total, starting from 2011 to 2022
print(arrange(filter(PM10s, site_number == "0030"), year), n=100)
# The site 0022 has 50 records in total, starts from 1999 to 2022
print(arrange(filter(PM10s, site_number == "0022"), year), n=100)

# So we remove these 2 sites at all, since they only have 50 + 22 = 72 records out of 8566 records

PM10s <- filter(PM10s, ! site_number %in% c("0030", "0022"))
saveRDS(PM10s, sprintf("%sPM10s_CA_selected.rds", result_path))

PM10s_summary <- filter(PM10s_summary, ! site_number %in% c("0030", "0022"))
saveRDS(PM10s_summary, sprintf("%sPM10s_CA_summary.rds", result_path))

PM10s_summary_scaled <- filter(PM10s_summary_scaled, ! site_number %in% c("0030", "0022"))
saveRDS(PM10s_summary_scaled, sprintf("%sPM10s_CA_summary_scaled.rds", result_path))

# The histgram now looks fine!
hist(PM10s_summary_scaled$annual_mean)





















# ##################################################
# # Use Annual Summary Data Files From The Website #
# ##################################################
# 
# # The annual summary data file is similar to that from downloaded json files, but are not exactly same.
# # After filtering out the extreme events, there are 8451 raws in total from 1985 to 2022 in the annual summary data.
# # But there are 8566 raws in total from 1985 to 2022 in the PM10s data. 
# # The set of sites is the same in the two data. 
# # So for the moment we use the PM10s dataset.But we can check later the annual summary data.
# 
# data_path <- "/home/xinglong/Downloads/annual_summary/"
# years <- seq(1980, 2022)
# summary_raw <- NULL
# 
# for (year in years) {
#   record <- read_csv(sprintf("%sannual_conc_by_monitor_%s.csv", data_path, year))
#   summary_raw <- bind_rows(summary_raw, record)
# }
# 
# names(summary_raw) <- gsub(" ", "_", tolower(names(summary_raw)))
# 
# annual_summary <- filter(summary_raw,
#                          state_code == '06',
#                          parameter_code == 81102,
#                          event_type %in% c("No Events",
#                                            "Events Excluded",
#                                            "Concurred Events Excluded")) %>%
#   rename(site_number = site_num) %>%
#   select(c("site_number",
#            "poc",
#            "latitude",
#            "longitude",
#            "datum",
#            "year",
#            "event_type",
#            "arithmetic_mean"))


# # Load site/monitor information of California / PM10 ===============================================
# 
# # For the purpose of this project, the site file is not very useful ,since one site can measure more
# # pollutants than PM10 and in this file each site also have multiple locations, even more than these
# # from the PM10s dataset, since here site has more locations for other pollutants.
# 
# sites <- read_csv("PM10/aqs_sites.csv")
# names(sites) <- gsub(" ", "_", tolower(names(sites)))
# ca_sites <- filter(sites, state_code == '06')
# 
# monitors <- read_csv("PM10/aqs_monitors.csv")
# names(monitors) <- gsub(" ", "_", tolower(names(monitors)))
# ca_monitors <- filter(monitors, 
#                       state_code == '06',
#                       parameter_code == 81102)
