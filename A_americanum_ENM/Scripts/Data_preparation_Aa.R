################################################################################
# Project: RII Track-2 FEC: Marshalling Diverse Big Data Streams to Understand 
#          Risk of Tick-Borne Diseases in the Great Plains
# Title: Data preparation for Ecological niche modeling
# Authors: Marlon E. Cobos, Ismari Martinez, Taylor Winters, A. Townsend Peterson
# Date updated: 21/03/20234(dd/mm/yyyy)
# Funding: NSF EPSCoR (IIA-1920946)
################################################################################



# Description ------------------------------------------------------------------
# This script serves to prepare tick presence-absence data for ecological niche
# modeling. 
#
# 1. The whole data set is processed to obtain presence-absence records from the
#    complete set of records gathered in this project.
# 2. Environmental data, previously prepared and downloaded using GEE is 
#    processed (PCA) and associated with the records.
#    - Make sure the Envrionmental layers are organized as follows:
#      - Data/Daymet/
#        - Variable name (e.g., dayl)/
#          - Year (e.g., 2000/)
#            - All files for weakly environmental averages for that year.
# 3. Files with data ready for analyses are written to working directory.
# ------------------------------------------------------------------------------



# Packages ---------------------------------------------------------------------
# install packages
#install.packages("terra")

# load packages
library(terra)
# ------------------------------------------------------------------------------



# Data to be used in all ENM related projects ----------------------------------
# working directory
setwd("YOUR/DIRECTORY")  # change as needed (contains the sub folder Data)

# tick data
tick_data <- read.csv("Data/KS_OK_tickdata_Aa.csv")
# ------------------------------------------------------------------------------



# General filtering and fixing of data -----------------------------------------
# column names
colnames(tick_data)


# fixing data
## fix life stage
tick_data$Life_stage <- gsub("Male", "Adult", tick_data$Life_stage)
tick_data$Life_stage <- gsub("Female", "Adult", tick_data$Life_stage)


# filter
## columns to keep
cols <- c("Species", "Longitude", "Latitude", "Day", "Month", "Year", 
          "Life_stage", "site.code")

## species to keep
sps <- c("Amblyomma americanum", "Other species")

## filter
row_filter <- tick_data$Species %in% sps & !is.na(tick_data$Species)  
adv_ticks <- tick_data[row_filter, cols]

## erase rows that contain NA in any of the columns
adv_ticks <- na.omit(adv_ticks)

unique(adv_ticks[, c(2:3, 8)])  # check unique site name and xy coordinate
# ------------------------------------------------------------------------------



# Presences and absences per species -------------------------------------------
# erase duplicated records
dup <- paste(adv_ticks$Species, adv_ticks$Day, adv_ticks$Month, adv_ticks$Year,
             adv_ticks$site.code)

tick_unique <- adv_ticks[!duplicated(dup), ]


# presences per species (columns 2:6 contain relevant information)
aa_ticks <- tick_unique[tick_unique$Species == sps[1], 2:6]


# absences
## absences for all species
abs_filter <- unique(c(grep("no_ticks", tick_data$Your.ID.Number),
                       grep("No_ticks", tick_data$Def.ID.number)))

all_abs <- tick_data[abs_filter, cols[2:6]]  # keep relevant rows and columns
all_abs <- unique(all_abs)  # erase duplicates just in case

## absences per species
### function to make it easier to find absences per species
get_abs <- function(unique_data, species_col, species_name, vtf_duplicates,
                    all_absences) {
  aaa_rule1 <- unique_data[, species_col] != species_name
  aa_nodup1 <- vtf_duplicates[aaa_rule1]
  aa_dup1 <- vtf_duplicates[unique_data[, species_col] == species_name]
  
  rule_aa <- aa_nodup1[!aa_nodup1 %in% aa_dup1]
  
  no_ticks <- unique_data[aaa_rule1 & vtf_duplicates %in% rule_aa, 2:6]
  no_ticks <- rbind(no_ticks, all_absences)
  no_ticks <- unique(no_ticks)
}

### vector to find absences
dup1 <- paste(tick_unique$Day, tick_unique$Month, tick_unique$Year, 
              tick_unique$site.code)

### all absences species vy species
aa_noticks <- get_abs(unique_data = tick_unique, species_col = "Species", 
                      species_name = sps[1], vtf_duplicates = dup1, 
                      all_absences = all_abs)


# combine presence-absence data per species
aa_all <- rbind(data.frame(Presence_absence = 1, aa_ticks),
                data.frame(Presence_absence = 0, aa_noticks))


# write files
write.csv(aa_all, "Data/aa_presence_absence.csv", row.names = FALSE)
# ------------------------------------------------------------------------------



# Presences and absences per species and life stage ----------------------------
# erase duplicated records
dup <- paste(adv_ticks$Species, adv_ticks$Day, adv_ticks$Month, adv_ticks$Year,
             adv_ticks$Life_stage, adv_ticks$site.code)

tick_uniquels <- adv_ticks[!duplicated(dup), ]


# presences per species (columns 2:7 also include life stage)
aa_ticksls <- tick_uniquels[tick_uniquels$Species == sps[1], 2:7]


# absences per species and life stage and only per species are the same
# only adding extra column for compativility
aa_noticks$Life_stage <- NA_character_

# combining data
aa_allls <- rbind(data.frame(Presence_absence = 1, aa_ticksls),
                  data.frame(Presence_absence = 0, aa_noticks))


# write files
write.csv(aa_allls, "Data/aa_presence_absence_ls.csv", row.names = FALSE)
# ------------------------------------------------------------------------------



# Environmental data preparation -----------------------------------------------
# variables
varnames1 <- c("dayl", "prcp", "srad", "tmax", "tmin", "vp")

# get values (2020-2022) relevant considering tick sampling
varval <- lapply(varnames1, function(x) {
  folders <- paste0("Data/Daymet/", x, "/", c(2020:2022))
  
  layers <- c(list.files(folders[1], pattern = ".tif$", full.names = TRUE),
              list.files(folders[2], pattern = ".tif$", full.names = TRUE),
              list.files(folders[3], pattern = ".tif$", full.names = TRUE))
  message(x)
  c(as.matrix(rast(layers)))
})

## put all information together
varval1 <- do.call(cbind, varval)
varval1 <- na.omit(varval1)
dim(varval1)

colnames(varval1) <- paste0(varnames1, "_mean")  # to match layer names


# PCA
vpca1 <- prcomp(varval1, retx = FALSE, center = TRUE, scale. = TRUE)

summary(vpca1)

## save PCA results
### folder
pcadir <- "Data/PCA_results"
dir.create(pcadir)

### pca as RData
save(vpca1, file = "Data/PCA_results/PCA_all.RData")

### summary and rotations as txt
sink("Data/PCA_results/PCA_summary.txt")
summary(vpca1)
cat("\n\n\n")
print(vpca1$rotation)
sink()


# predictions to all years (1980-2022)
## all weeks (according to Julian days)
ejdays <- seq(0, 360, 8)

## year week combination for names
yjday <- paste0(rep(1980:2022, each = length(ejdays)), "_", ejdays)

pcnames1 <- paste0("Data/PCA_results/pcs_", yjday, ".tif")

## all names of variables to be used
all_layers <- list.files("Data/Daymet", pattern = ".tif$", full.names = TRUE, 
                         recursive = TRUE)

### exclude swe layers. They were downloaded with GEE but not used in analysis
ex_swe <- grep("/swe/", all_layers)
all_layers <- all_layers[-ex_swe]

## prediction in loop
pcpred <- lapply(1:length(pcnames1), function(x) {
  varpc <- grep(paste0(yjday[x], ".tif"), all_layers, value = TRUE)
  pcs <- terra::predict(rast(varpc), vpca1, filename = pcnames1[x])
  x
})
# ------------------------------------------------------------------------------



# Environmental data extraction (tick data per species) ------------------------
# extraction of Edata according to dates in data and environmental layers
## function to make it easier to extract variables per species
extract_matching <- function(pa_data, day_col, month_col, year_col, date_format,
                             lonlat_cols, jdays_tomatch, layer_names) {
  ## dates in data
  dates <- as.Date(paste(pa_data[, day_col], pa_data[, month_col], 
                         pa_data[, year_col], sep = "-"), format = date_format)
  
  jdays <- as.POSIXlt(dates)$yday
  years <- pa_data[, year_col]
  
  ## matching date from data and environmental layers 
  matched_days <- sapply(jdays, function(y) {
    aday <- jdays_tomatch[jdays_tomatch <= y]
    aday[length(aday)]
  })
  
  ## year and Julian day combination
  year_jday <- paste(years, matched_days, sep = "_")
  uyear_jday <- unique(year_jday)
  
  
  # data extraction
  ## extracting raw variable values according to dates in a loop
  edata_extracted <- lapply(uyear_jday, function(y) {
    wweek <- which(year_jday == y)
    
    vars <-  grep(paste0(y, ".tif$"), layer_names, value = TRUE)
    varr <- rast(vars)
    
    setrec <- pa_data[wweek, ]
    varset <- extract(varr, setrec[, lonlat_cols])[, -1]
    cbind(setrec, varset)
  })
  
  edata_extracted <- do.call(rbind, edata_extracted)
  
  return(edata_extracted)
}

## defining repeated arguments
head(aa_all)

dy <- "Day"
mt <- "Month"
yr <- "Year"
dt_format <- "%d-%B-%Y"
lonlat <- c("Longitude", "Latitude")

## extraction of raw environmental data per species
aa_all_env <- extract_matching(pa_data = aa_all, day_col = dy, month_col = mt, 
                               year_col = yr, date_format = dt_format, 
                               lonlat_cols = lonlat, jdays_tomatch = ejdays, 
                               layer_names = all_layers)

## extraction of PCs per species
### names of all pcs
all_pcs <- list.files("Data/PCA_results", pattern = ".tif$", full.names = TRUE)

### extraction
aa_all_pcs <- extract_matching(pa_data = aa_all, day_col = dy, month_col = mt, 
                               year_col = yr, date_format = dt_format, 
                               lonlat_cols = lonlat, jdays_tomatch = ejdays, 
                               layer_names = all_pcs)


# write files
write.csv(aa_all_env, "Data/aa_presence_absence_env.csv", row.names = FALSE)
write.csv(aa_all_pcs, "Data/aa_presence_absence_pcs.csv", row.names = FALSE)
# ------------------------------------------------------------------------------




# Environmental data extraction (tick data per species and life stage) ---------
# extraction of Edata according to dates in data and environmental layers
## extraction of raw environmental data per species and life stage
aa_allls_env <- extract_matching(pa_data = aa_allls, day_col = dy, month_col = mt, 
                                 year_col = yr, date_format = dt_format, 
                                 lonlat_cols = lonlat, jdays_tomatch = ejdays, 
                                 layer_names = all_layers)

## extraction of PCs per species and life stage
### extraction
aa_allls_pcs <- extract_matching(pa_data = aa_allls, day_col = dy, month_col = mt, 
                                 year_col = yr, date_format = dt_format, 
                                 lonlat_cols = lonlat, jdays_tomatch = ejdays, 
                                 layer_names = all_pcs)


# write files
write.csv(aa_allls_env, "Data/aa_presence_absence_env_ls.csv", row.names = FALSE)
write.csv(aa_allls_pcs, "Data/aa_presence_absence_pcs_ls.csv", row.names = FALSE)
# ------------------------------------------------------------------------------
