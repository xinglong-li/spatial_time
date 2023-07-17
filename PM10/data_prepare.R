library(jsonlite)
library(dplyr)


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

saveRDS(PM10s_raw, sprintf("%sPM10s_raw.rds", result_path))


# Prepare the data for model fitting ===============================================================

# Descriptions of data each column is in 
# https://aqs.epa.gov/aqsweb/airdata/FileFormats.html#_annual_summary_files

PM10s <- select(PM10s_raw, c("site_number",
                             "local_site_name",
                             "poc",
                             "latitude",
                             "longitude",
                             "year",
                             "event_type",
                             "arithmetic_mean"))

saveRDS(PM10s, sprintf("%sPM10s.rds", result_path))

# Remove Extreme events, there are sites with only one record "Concurred Events Excluded",
# and the total number of sites will reduce if we exclude such records. 
# So we keep such records at the moments.
PM10s <- filter(PM10s, event_type %in% c("No Events", 
                                         "Events Excluded", 
                                         "Concurred Events Excluded"))

# It is strange to see that some 'site_number's has multiple 'local_site_name's, and each local name
# has different latitude and longitude.
# But some site has no local_site_name.
# For example, site 0007 in year 2021 has 2 local_site_names, even stranger, the poc #1 of this site
# has 2 local_site_name and they have different locations.


# Aggregate by site --------------------------------------------------------------------------------

# But how to determine the the location of each site in this case? 
PM10s_by_site <- group_by(PM10s, year, site_number) %>%
  summarise(arithmetic_mean = mean(arithmetic_mean))
print(PM10s_by_site)

# Aggregate by "local_site_name".
PM10s_by_name <- group_by(PM10s, year, site_number, local_site_name) %>%
  summarise(arithmetic_mean = mean(arithmetic_mean))
print(PM10s_by_name)



