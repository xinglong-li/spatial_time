library(dplyr)
library(ggplot2)

path <- "/home/xinglong/git_local/spacial_time/PM10/"

PM10s <- readRDS(sprintf("%sPM10s_CA_summary.rds", path))
print(PM10s, n=30)

