library(jsonlite)
library(dplyr)


# Read data from each year (from 1980 to 2022) and combine the PM10 records ========================

data_path <- "/home/xinglong/Downloads/annual_by_monitor/"

annul_records <- list.files(data_path, pattern="*.json")

PM10s <- NULL

for (annual_record in annul_records) {
  record <- read_json(sprintf("%s%s", data_path, annual_record))$Data %>%
    bind_rows()
  PM10s <- bind_rows(PM10_all, record)
}


# Check data source ================================================================================

# Make sure that the time is from 1985 --- 2022
stopifnot(
  length(
    setdiff(unique(PM10s$year), seq(1985, 2022))
    ) == 0
  )

# Make sure that all measurements is from California
stopifnot(
  all(PM10s$state == 'California')
)

# Make sure that all measurements is PM10
stopifnot(
  all(PM10s$parameter == "PM10 Total 0-10um STP")
)

# Make sure that all records are measured in the same unit
stopifnot(
  length(unique(PM10s$units_of_measure)) == 1
)

# Some information of the data
print(sprintf("Total number of records: %s", dim(PM10s)[1]))
print(sprintf("Total number of sites in California: %s", length(unique(PM10s$site_number))))


# Prepare the data for model fitting ===============================================================

PM10s <- select(PM10s, c("site_number",
                         "poc",
                         "latitude",
                         "longitude",
                         "year",
                         "event_type",
                         "observation_count",
                         "observation_percent",
                         "validity_indicator",
                         "valid_day_count",
                         "required_day_count",
                         "exceptional_data_count",
                         "null_observation_count",
                         "certification_indicator",
                         "arithmetic_mean",
                         "standard_deviation",
                         "local_site_name",
                         "county",
                         "city",
                         "cbsa_code",
                         "cbsa"))

result_path <- "/home/xinglong/git_local/spacial_time/PM10/PM10s.rds"
saveRDS(PM10s, result_path)


