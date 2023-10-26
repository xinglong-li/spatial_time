library(dplyr)
library(sf)
library(sp)
library(reshape2)
library(INLA)
library(ggplot2)
library(inlabru)



# Prepare the data =============================================================

# Load data 
load("./Reproducibility/Data2Joe.RData")
# Load the indicator telling us if the mesh vertices lie in GB
xy_in = readRDS("./Data/Reproducibility/xy_in.rds") 

# Standardize the location
BS <- BlackSmokePrefData
BS[, c(2,3)] <- BS[, c(2,3)]/sd(BS[, 2])
ggplot(BS) + geom_point(aes(x = east, y = north)) + coord_equal()

# Reshape the data with one observation per row (required by INLA)
BS2 <- melt(BS, id.vars = c(1, 2, 3), variable.name = 'year', value.name = 'bsmoke')
BS2$year = as.numeric(as.character(factor(BS2$year, labels =66:96 )))

ggplot(BS2) + geom_histogram(aes(x = bsmoke), bins = 40) 
# Right skew - take natural log and make it unitless by minus log of mean value
BS2$bsmoke = log(BS2$bsmoke / mean(colMeans(BS[,4:34], na.rm = T)))
ggplot(BS2) + geom_histogram(aes(x = bsmoke), bins = 40) 


subsmaple_plotting_data = BS2[which(BS2$site %in% levels(BS2$site)[
  sample.int(1466, size = 30)]), ]

ggplot(data = subsmaple_plotting_data, aes(x = year, y = bsmoke)) +
  geom_line(aes(group = site, colour = site))

var_annual <- group_by(BS2, year) %>% 
  summarise(var_bs = var(bsmoke, na.rm = T))
ggplot(data = var_annual, aes(x = year, y = var_bs)) +
  geom_line() + geom_smooth() + 
  ggtitle('A plot of the variance of the log annual means with fitted smoother')


# Prepare the GB border shape file =============================================

# gb_boundary <- read_sf(dsn = "/Users/Xinglong/Downloads/gb_2011", layer = "infuse_gb_2011_clipped")
# gb_boundary <- gb_boundary$geometry
# km_proj <- CRS("+proj=utm +zone=30 + ellps=WGS84 +units=km")
# gb_boundary_km <- st_transform(gb_boundary, km_proj)
# plot(gb_boundary_km)
# # saveRDS(ca_boundary_km, sprintf("%sCA_border.rds", result_path))
# # Re-scale coordinates of CA boundary
# # Extract coordinates and put them in a list
# extractCoords <- function(sp.df)
# {
#   results <- list()
#   for(i in 1:length(sp.df@polygons[[1]]@Polygons))
#   {
#     results[[i]] <- sp.df@polygons[[1]]@Polygons[[i]]@coords
#   }
#   results
# }
# 
# vertices <- st_coordinates(gb_boundary_km)
# 
# scaled.vertices <- lapply(vertices, function(x) x / 100)
# 
# #Create Polygons and Spatial Polygons Data Frame from scaled list of coordinates
# Polys <- list()
# for (i in 1:length(scaled.vertices)) {
#   Polys[i] <- sp::Polygon(scaled.vertices[[i]])
# }
# Polys.plural <- sp::Polygons(Polys, ID = "0")
# Polys.sp <- sp::SpatialPolygons(list(Polys.plural), proj4string = km_proj)
# gb_boundary_km_scaled <- sp::SpatialPolygonsDataFrame(Polys.sp, data = gb_boundary_km@data)
# 
# saveRDS(ca_boundary_km_scaled, sprintf("%sCA_border_scaled.rds", result_path))



# Prepare variables for the model ==============================================

# 

# The mesh grid used by Joe 
mesh_jw <- readRDS('./Data/Reproducibility/mesh_5_7.rds')

ggplot(BS) + gg(mesh_jw) + geom_point(aes(x = east, y = north)) + coord_fixed()
mesh_jw$n


