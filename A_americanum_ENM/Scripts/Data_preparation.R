################################################################################
# Project: RII Track-2 FEC: Marshalling Diverse Big Data Streams to Understand 
#          Risk of Tick-Borne Diseases in the Great Plains
# Title: Data preparation for Ecological niche modeling
# Authors: Marlon E. Cobos, Ismari Martinez, Taylor Winters, A. Townsend Peterson
# Date: 10/08/2023 (dd/mm/yyyy)
# Funding: NSF EPSCoR (IIA-1920946)
################################################################################



# Description ------------------------------------------------------------------

# ------------------------------------------------------------------------------



# Packages ---------------------------------------------------------------------
# install packages
#install.packages("terra")

# load packages
#library(terra)
library(geodata)
# ------------------------------------------------------------------------------



# Data to be used in all ENM related projects ----------------------------------
# state borders (Download, store in folder Data)
states <- gadm("USA", level = 1, path = "Data")


# tick data
tick_data <- read.csv("Data/KS_OK_tickdata.csv")
# ------------------------------------------------------------------------------



# General filtering and fixing of data -----------------------------------------
# column names
colnames(tick_data)


# fixing data
## fix life stage
tick_data$Life.Stage <- gsub("Male", "Adult", tick_data$Life.Stage)
tick_data$Life.Stage <- gsub("Female", "Adult", tick_data$Life.Stage)

colnames(tick_data)[35] <- "Life_stage"


# filter
## columns to keep
cols <- c("Species", "Longitude", "Latitude", "Day", "Month", "Year", 
          "Life_stage", "site.code")

## species to keep
sps <- c("Amblyomma americanum", "Amblyomma maculatum", 
         "Dermacentor variabilis", "Dermacentor albipictus", 
         "Ixodes scapularis")

## filter
row_filter <- tick_data$Species %in% sps & !is.na(tick_data$Species)  

adv_ticks <- tick_data[row_filter, cols]

sum(is.na(adv_ticks$Year))

adv_ticks[is.na(adv_ticks$Year), ]

## erase NA rows
adv_ticks <- na.omit(adv_ticks)

unique(adv_ticks[, c(2:3, 8)])
# ------------------------------------------------------------------------------



# Presences and absences per species -------------------------------------------
# erase duplicated records
dup <- paste(adv_ticks$Species, adv_ticks$Day, adv_ticks$Month, adv_ticks$Year,
             adv_ticks$site.code)

tick_unique <- adv_ticks[!duplicated(dup), ]


# presences per species (columns 2:6 contain relevant information)
aa_ticks <- tick_unique[tick_unique$Species == sps[1], 2:6]
am_ticks <- tick_unique[tick_unique$Species == sps[2], 2:6]
dv_ticks <- tick_unique[tick_unique$Species == sps[3], 2:6]
da_ticks <- tick_unique[tick_unique$Species == sps[4], 2:6]
is_ticks <- tick_unique[tick_unique$Species == sps[5], 2:6]


# absences
## absences for all species
grep("no_ticks", tick_data$Your.ID.Number)

grep("No_ticks", tick_data$Def.ID.number)

grep("no_ticks", tick_data$Your.ID.Number) %in% 
  grep("No_ticks", tick_data$Def.ID.number)


abs_filter <- c(grep("no_ticks", tick_data$Your.ID.Number),
                grep("No_ticks", tick_data$Def.ID.number))

all_abs <- tick_data[abs_filter, cols]
all_abs <- all_abs[, 2:6]

all_abs <- unique(all_abs)

## absences per species
### vector to find absences
dup1 <- paste(tick_unique$Day, tick_unique$Month, tick_unique$Year, 
              tick_unique$site.code)

### A. americanum absences
aaa_rule1 <- tick_unique$Species != sps[1]
aa_nodup1 <- dup1[aaa_rule1]
aa_dup1 <- dup1[tick_unique$Species == sps[1]]

rule_aa <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

aa_noticks <- tick_unique[aaa_rule1 & dup1 %in% rule_aa, 2:6]
aa_noticks <- rbind(aa_noticks, all_abs)
aa_noticks <- unique(aa_noticks)

### A. maculatum absences
aaa_rule1 <- tick_unique$Species != sps[2]
aa_nodup1 <- dup1[aaa_rule1]
aa_dup1 <- dup1[tick_unique$Species == sps[2]]

rule_am <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

am_noticks <- tick_unique[aaa_rule1 & dup1 %in% rule_am, 2:6]
am_noticks <- rbind(am_noticks, all_abs)
am_noticks <- unique(am_noticks)

### D. variabilis absences
aaa_rule1 <- tick_unique$Species != sps[3]
aa_nodup1 <- dup1[aaa_rule1]
aa_dup1 <- dup1[tick_unique$Species == sps[3]]

rule_dv <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

dv_noticks <- tick_unique[aaa_rule1 & dup1 %in% rule_dv, 2:6]
dv_noticks <- rbind(dv_noticks, all_abs)
dv_noticks <- unique(dv_noticks)

### D. albipictus absences
aaa_rule1 <- tick_unique$Species != sps[4]
aa_nodup1 <- dup1[aaa_rule1]
aa_dup1 <- dup1[tick_unique$Species == sps[4]]

rule_da <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

da_noticks <- tick_unique[aaa_rule1 & dup1 %in% rule_da, 2:6]
da_noticks <- rbind(da_noticks, all_abs)
da_noticks <- unique(da_noticks)

### I. scapularis absences
aaa_rule1 <- tick_unique$Species != sps[5]
aa_nodup1 <- dup1[aaa_rule1]
aa_dup1 <- dup1[tick_unique$Species == sps[5]]

rule_is <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

is_noticks <- tick_unique[aaa_rule1 & dup1 %in% rule_is, 2:6]
is_noticks <- rbind(is_noticks, all_abs)
is_noticks <- unique(is_noticks)


# combine presence-absence data per species
aa_all <- rbind(data.frame(Presence_absence = 1, aa_ticks),
                data.frame(Presence_absence = 0, aa_noticks))

am_all <- rbind(data.frame(Presence_absence = 1, am_ticks),
                data.frame(Presence_absence = 0, am_noticks))

dv_all <- rbind(data.frame(Presence_absence = 1, dv_ticks),
                data.frame(Presence_absence = 0, dv_noticks))

da_all <- rbind(data.frame(Presence_absence = 1, da_ticks),
                data.frame(Presence_absence = 0, da_noticks))

is_all <- rbind(data.frame(Presence_absence = 1, is_ticks),
                data.frame(Presence_absence = 0, is_noticks))


# write files
write.csv(aa_all, "Data/aa_presence_absence.csv", row.names = FALSE)
write.csv(am_all, "Data/am_presence_absence.csv", row.names = FALSE)
write.csv(dv_all, "Data/dv_presence_absence.csv", row.names = FALSE)
write.csv(da_all, "Data/da_presence_absence.csv", row.names = FALSE)
write.csv(is_all, "Data/is_presence_absence.csv", row.names = FALSE)
# ------------------------------------------------------------------------------



# Presences and absences per species and life stage ----------------------------
# erase duplicated records
dup <- paste(adv_ticks$Species, adv_ticks$Day, adv_ticks$Month, adv_ticks$Year,
             adv_ticks$Life_stage, adv_ticks$site.code)

tick_uniquels <- adv_ticks[!duplicated(dup), ]


# presences per species (columns 2:7 also include life stage)
aa_ticksls <- tick_uniquels[tick_uniquels$Species == sps[1], 2:7]
am_ticksls <- tick_uniquels[tick_uniquels$Species == sps[2], 2:7]
dv_ticksls <- tick_uniquels[tick_uniquels$Species == sps[3], 2:7]
da_ticksls <- tick_uniquels[tick_uniquels$Species == sps[4], 2:7]
is_ticksls <- tick_uniquels[tick_uniquels$Species == sps[5], 2:7]


# absences per species (adding life stage column)
aa_noticks$Life_stage <- NA_character_
am_noticks$Life_stage <- NA_character_
dv_noticks$Life_stage <- NA_character_
da_noticks$Life_stage <- NA_character_
is_noticks$Life_stage <- NA_character_

# combining data
aa_allls <- rbind(data.frame(Presence_absence = 1, aa_ticksls),
                  data.frame(Presence_absence = 0, aa_noticks))

am_allls <- rbind(data.frame(Presence_absence = 1, am_ticksls),
                  data.frame(Presence_absence = 0, am_noticks))

dv_allls <- rbind(data.frame(Presence_absence = 1, dv_ticksls),
                  data.frame(Presence_absence = 0, dv_noticks))

da_allls <- rbind(data.frame(Presence_absence = 1, da_ticksls),
                  data.frame(Presence_absence = 0, da_noticks))

is_allls <- rbind(data.frame(Presence_absence = 1, is_ticksls),
                  data.frame(Presence_absence = 0, is_noticks))


# write files
write.csv(aa_allls, "Data/aa_presence_absence_ls.csv", row.names = FALSE)
write.csv(am_allls, "Data/am_presence_absence_ls.csv", row.names = FALSE)
write.csv(dv_allls, "Data/dv_presence_absence_ls.csv", row.names = FALSE)
write.csv(da_allls, "Data/da_presence_absence_ls.csv", row.names = FALSE)
write.csv(is_allls, "Data/is_presence_absence_ls.csv", row.names = FALSE)
# ------------------------------------------------------------------------------



# Presences and absences for A. americanum pathogens ---------------------------
# column names
colnames(tick_data)


# filter data
## columns to keep
cols <- c("Species", "Longitude", "Latitude", "Day", "Month", "Year", 
          "tested", "pathogen1_ricketsia", "pathogen2_e_chaf",        
          "pathogen3_e_ewingii", "pathogen4_borrelia", "site.code")

## species to keep
unique(tick_data$Species)

sp <- "Amblyomma americanum"

row_filter <- tick_data$Species %in% sp & !is.na(tick_data$Species) 

## filter
aa_ticks <- tick_data[row_filter, cols]


# counts 
## keep only tested ticks
table(aa_ticks$tested)
aa_tested <- aa_ticks[aa_ticks$tested == "yes", ]

## records per test and pathogen
table(aa_tested$pathogen1_ricketsia)
table(aa_tested$pathogen2_e_chaf)
table(aa_tested$pathogen3_e_ewingii)
table(aa_tested$pathogen4_borrelia)


# fix name inconsistencies
aa_tested$pathogen1_ricketsia <- gsub("Ricketsia spp.", "Rickettsia sp.", 
                                      aa_tested$pathogen1_ricketsia)

aa_tested$pathogen1_ricketsia <- gsub("R.amblyommatis", "R. amblyommatis", 
                                      aa_tested$pathogen1_ricketsia)

aa_tested$pathogen2_e_chaf <- gsub("neg", "NEG", aa_tested$pathogen2_e_chaf)

aa_tested$pathogen3_e_ewingii <- gsub("neg", "NEG", 
                                      aa_tested$pathogen3_e_ewingii)

aa_tested$pathogen4_borrelia <- gsub("neg", "NEG", aa_tested$pathogen4_borrelia)


# create new columns with test results per pathogen (0 = Negative; 1 = Positive)
aa_tested$R_amblyommatis <- (aa_tested$pathogen1_ricketsia == 
                               "R. amblyommatis") * 1
table(aa_tested$R_amblyommatis)

aa_tested$B_lonestari <- (aa_tested$pathogen4_borrelia == "B. lonestari") * 1
table(aa_tested$B_lonestari)

aa_tested$E_ewingii <- (aa_tested$pathogen3_e_ewingii == "E. ewingii") * 1
table(aa_tested$E_ewingii)

aa_tested$E_chaffeensis <- (aa_tested$pathogen2_e_chaf == "E. chaff") * 1
table(aa_tested$E_chaffeensis)


# keep only useful data
colnames(aa_tested)
col_keep <- c(2:6, 13:16)

aa_pathogens <- aa_tested[, col_keep]

# save prepared data
write.csv(aa_pathogens, "Data/aa_pathogens.csv", row.names = FALSE)
# ------------------------------------------------------------------------------



# Environmental data preparation -----------------------------------------------
# variables
all_layers <- list.files("Data/raster", pattern = ".tif$", full.names = TRUE)

# get values
varnames1 <- c("dayl", "prcp", "srad", "tmax", "tmin", "vp")

varval <- lapply(varnames1, function(x) {
  layvar <- grep(x, all_layers, value = TRUE)
  message(x)
  c(as.matrix(rast(layvar)))
})

varval1 <- do.call(cbind, varval)
varval1 <- na.omit(varval1)
dim(varval1)

colnames(varval1) <- paste0(varnames1, "_mean")

# PCA
vpca1 <- prcomp(varval1, retx = TRUE, center = TRUE, scale. = TRUE)

summary(vpca1)

# save initial results
pcadir <- "Data/PCA_results_noswe"
dir.create(pcadir)

save(vpca1, file = "Data/PCA_results_noswe/PCA_all.RData")

sink("Data/PCA_results_noswe/PCA_summary.txt")
summary(vpca1)
cat("\n\n\n")
print(vpca1$rotation)
sink()

# predictions
## all Julian e days
ejdays <- seq(0, 360, 8)

yjday <- paste0(rep(c(2020, 2021, 2022), each = length(ejdays)), "_", ejdays)

pcnames1 <- paste0("Data/PCA_results_noswe/pcs_", yjday, ".tif")

pcpred <- lapply(1:length(pcnames1), function(x) {
  varpc <- grep(paste0(yjday[x], ".tif"), all_layers, value = TRUE)
  pcs <- terra::predict(rast(varpc), vpca1, filename = pcnames1[x])
})
# ------------------------------------------------------------------------------
