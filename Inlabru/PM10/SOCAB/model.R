library(dplyr)
library(ggplot2)
library(sf)
library(sp)
library(inlabru)
library(reshape2)


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

ggplot(PM10s_SOCAB_utm) + gg(as(SOCAB_border, "Spatial")) + 
  geom_point(aes(x = E, y = N), color = "blue") + coord_equal()

# Aggregate by site --------------------------------------------------------------------------------
# (Should we aggregate monitors belong to the same sites together?)
PM10s_SOCAB_summary <- group_by(PM10s_SOCAB_utm, site_number, year) %>%
  summarise(annual_mean = mean(arithmetic_mean),
            north = mean(N),
            east = mean(E)) %>%
  mutate(north = mean(north[site_number == site_number]),
         east = mean(east[site_number == site_number]))

ggplot(PM10s_SOCAB_summary) + gg(as(SOCAB_border, "Spatial")) + 
  geom_point(aes(x = east, y = north), color = "blue") + coord_equal()

# Rescale the coordinates of sites and border, the new unit is 10km --------------------------------

sd_north <- sd(PM10s_SOCAB_summary$north) # = 20.45
sd_east <- sd(PM10s_SOCAB_summary$east) # = 41.68

PM10s_SOCAB_scaled <- PM10s_SOCAB_summary
PM10s_SOCAB_scaled$north <- PM10s_SOCAB_summary$north / 10
PM10s_SOCAB_scaled$east <- PM10s_SOCAB_summary$east / 10

SOCAB_border_scaled <- (SOCAB_border / matrix(data = c(10, 10), ncol = 2)) %>%
  as("Spatial")

ggplot(PM10s_SOCAB_scaled) + gg(SOCAB_border_scaled) + 
  geom_point(aes(x = east, y = north), color = "blue") + coord_equal()


# Fit the Preferential sampling model ==============================================================

# Data preprocessing -------------------------------------------------------------------------------

PM10s0 <- PM10s_SOCAB_scaled
print(sprintf("Total number of annual measurements: %s", dim(PM10s0)[1]))
print(sprintf("Total number of sites: %s", length(unique(PM10s0$site_number))))
print(sprintf("The sd of East coordinates is: %s", sd(PM10s0$east)))
print(sprintf("The sd of North coordinates is: %s", sd(PM10s0$north)))

mean_annually <- group_by(PM10s0, year) %>%
  summarise(mean_pm = mean(annual_mean)) 
mean_pm_annually <- mean(mean_annually$mean_pm)
PM10s0$annual_mean <- log(PM10s0$annual_mean) - log(mean_pm_annually)

no_sites = length(unique(PM10s0$site_number))
no_T = length(unique(PM10s0$year))

# Reform the data so that each site has num_of_year rows, the empty records are filled with NAs
PM10s_flat <- dcast(PM10s0, site_number + north + east ~ year, value.var = "annual_mean")
PM10s = melt(PM10s_flat, id.vars = c(1,2,3), variable.name = 'year', value.name = 'annual_mean')
PM10s$year <- as.numeric(as.character(factor(PM10s$year, labels = 1985:2022)))

stopifnot(dim(PM10s)[1] == no_sites * no_T)

means_plot = ggplot(data = PM10s, aes(x = year, y = annual_mean)) +
  geom_line(aes(group = site_number, colour = site_number)) + xlab("Year") + ylab("log(PM10)")
means_plot

var_annual <- group_by(PM10s, year) %>% 
  summarise(var_pm = var(annual_mean, na.rm = T))
variance_plot =  ggplot(data = var_annual, aes(x = year, y = var_pm)) +
  geom_line() + geom_smooth() + xlab("Year") + ylab("Variance log(PM10)") 
ggtitle('A plot of the variance of the log annual means with fitted smoother')
variance_plot

# Prepare variables for the model ==================================================================

PM10s$time <- (PM10s$year - min(PM10s$year)) / (max(PM10s$year) - min(PM10s$year))
PM10s$locs <- coordinates(PM10s[, c("east", "north")])
PM10s$site_number <- as.numeric(as.factor(PM10s$site_number))

# Site-selection indicator
PM10s$R <- as.numeric(!is.na(PM10s$annual_mean))
PM10s$R_lag <- c(rep(NA, no_sites), PM10s$R[1:(dim(PM10s)[1]-no_sites)])

# Response variable for the auxiliary model
PM10s$zero <- 0

# Compute Euclidean distances between all the sites
dists <- spDists(cbind(PM10s$east[1:no_sites], PM10s$north[1:no_sites]))

r <- 0.1 # The maximum radius of interest to be 10 km
PM10s$repulsion_ind <- 0

counter <- no_sites + 1
for (i in sort(unique(PM10s$year))[-1]) {
  # First extract the data at time i
  data_i <- PM10s[PM10s$year == i, ]
  # Compute the repulsion indicator. Was there a site at year i-1 within radius r of it?
  PM10s$repulsion_ind[counter:(counter+(no_sites - 1))] = rowSums(dists[, which(data_i$R_lag==1)] < r) > 0
  counter <- counter + no_sites
}

# Create mesh
edge_in = 0.15 # 20km
edge_out = 2 * edge_in
mesh = fm_mesh_2d_inla(loc = cbind(PM10s$east, PM10s$north),
                       boundary = CA_border,
                       offset = c(2*edge_in, edge_out), 
                       max.edge = 2*c(edge_in, edge_out),
                       cutoff = edge_in,
                       min.angle = 21)
mesh$n
ggplot(PM10s) + gg(mesh) + geom_point(aes(x = east, y = north), color = "blue") + coord_fixed()

# Maybe we should consider set the PC prior using data info
spde_obj <- inla.spde2.pcmatern(mesh = mesh, 
                                alpha = 2, 
                                prior.range = c(edge_in, 0.01),
                                prior.sigma = c(sqrt(mean(var_annual$var_pm)), 0.1),
                                constr = T)


