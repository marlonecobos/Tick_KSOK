################################################################################
# Project: ENM tick species (KSOK)
# Title: Complete process
# Authors: Marlon E. Cobos, Ismari Martinez, Taylor Winters, A. Townsend Peterson
# Date: 10/08/2023 (dd/mm/yyyy)
################################################################################


# Description ------------------------------------------------------------------

# ------------------------------------------------------------------------------


# Packages ---------------------------------------------------------------------
# install packages
#install.packages("terra")
#remotes::install_github("marlonecobos/enmpa")

# load packages
library(terra)
library(geodata)
library(enmpa)
# ------------------------------------------------------------------------------


# General filtering and fixing of data -----------------------------------------
# state borders
states <- gadm("USA", level = 1, path = "Data")

# reading the data
tick_data <- read.csv("Data/KS_OK_uniontickdata_fix_20230811.csv")

## column names
colnames(tick_data)


# fixing data
## fix life stage
tick_data$Life.Stage <- gsub("Male", "Adult", tick_data$Life.Stage)
tick_data$Life.Stage <- gsub("Female", "Adult", tick_data$Life.Stage)

colnames(tick_data)[35] <- "Life_stage"


# filtering
## columns to keep
cols <- c("Species", "Longitude", "Latitude", "Day", "Month", "Year", 
          "Life_stage", "site.code")

## species to keep
sps <- c("Amblyomma americanum", "Ixodes scapularis", "Dermacentor variabilis",
         "Amblyomma maculatum", "Dermacentor albipictus")

## filtering
row_filter <- tick_data$Species %in% sps & !is.na(tick_data$Species)  

adv_ticks <- tick_data[row_filter, cols]

sum(is.na(adv_ticks$Year))

adv_ticks[is.na(adv_ticks$Year), ]

## erasing NA rows
adv_ticks <- na.omit(adv_ticks)

unique(adv_ticks[, c(2:3, 8)])
# ------------------------------------------------------------------------------


# Presences and absences per species -------------------------------------------
# erase duplicated records
dup <- paste(adv_ticks$Species, adv_ticks$Day, adv_ticks$Month, adv_ticks$Year,
             adv_ticks$site.code)

tick_unique <- adv_ticks[!duplicated(dup), ]

# separate files
aa_ticks <- tick_unique[tick_unique$Species == "Amblyomma americanum", 2:6]
am_ticks <- tick_unique[tick_unique$Species == "Amblyomma maculatum", 2:6]
dv_ticks <- tick_unique[tick_unique$Species == "Dermacentor variabilis", 2:6]
da_ticks <- tick_unique[tick_unique$Species == "Dermacentor albipictus", 2:6]
is_ticks <- tick_unique[tick_unique$Species == "Ixodes scapularis", 2:6]


# absences
## absences for all
grep("no_ticks", tick_data$Your.ID.Number) %in% 
  grep("No_ticks", tick_data$Def.ID.number)

abs_filter <- c(grep("no_ticks", tick_data$Your.ID.Number),
                  grep("No_ticks", tick_data$Def.ID.number))

all_abs <- tick_data[abs_filter, cols]
all_abs <- all_abs[, 2:6]

all_abs <- unique(all_abs)

## to find absences
dup1 <- paste(tick_unique$Day, tick_unique$Month, tick_unique$Year, 
              tick_unique$site.code)

# aa
aaa_rule1 <- tick_unique$Species != "Amblyomma americanum"
aa_dup1 <- dup1[tick_unique$Species == "Amblyomma americanum"]
aa_nodup1 <- dup1[aaa_rule1]

rule_aa <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

aa_noticks <- tick_unique[aaa_rule1 & dup1 %in% rule_aa, 2:6]
aa_noticks <- rbind(aa_noticks, all_abs)
aa_noticks <- unique(aa_noticks)

# am
aaa_rule1 <- tick_unique$Species != "Amblyomma maculatum"
aa_dup1 <- dup1[tick_unique$Species == "Amblyomma maculatum"]
aa_nodup1 <- dup1[aaa_rule1]

rule_am <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

am_noticks <- tick_unique[aaa_rule1 & dup1 %in% rule_am, 2:6]
am_noticks <- rbind(am_noticks, all_abs)
am_noticks <- unique(am_noticks)

# dv
aaa_rule1 <- tick_unique$Species != "Dermacentor variabilis"
aa_dup1 <- dup1[tick_unique$Species == "Dermacentor variabilis"]
aa_nodup1 <- dup1[aaa_rule1]

rule_dv <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

dv_noticks <- tick_unique[aaa_rule1 & dup1 %in% rule_dv, 2:6]
dv_noticks <- rbind(dv_noticks, all_abs)
dv_noticks <- unique(dv_noticks)

# da
aaa_rule1 <- tick_unique$Species != "Dermacentor albipictus"
aa_dup1 <- dup1[tick_unique$Species == "Dermacentor albipictus"]
aa_nodup1 <- dup1[aaa_rule1]

rule_da <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

da_noticks <- tick_unique[aaa_rule1 & dup1 %in% rule_da, 2:6]
da_noticks <- rbind(da_noticks, all_abs)
da_noticks <- unique(da_noticks)

# is
aaa_rule1 <- tick_unique$Species != "Ixodes scapularis"
aa_dup1 <- dup1[tick_unique$Species == "Ixodes scapularis"]
aa_nodup1 <- dup1[aaa_rule1]

rule_is <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

is_noticks <- tick_unique[aaa_rule1 & dup1 %in% rule_is, 2:6]
is_noticks <- rbind(is_noticks, all_abs)
is_noticks <- unique(is_noticks)

# combining data
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

# separate files
aa_ticksls <- tick_uniquels[tick_uniquels$Species == "Amblyomma americanum", 2:7]
am_ticksls <- tick_uniquels[tick_uniquels$Species == "Amblyomma maculatum", 2:7]
dv_ticksls <- tick_uniquels[tick_uniquels$Species == "Dermacentor variabilis", 2:7]
da_ticksls <- tick_uniquels[tick_uniquels$Species == "Dermacentor albipictus", 2:7]
is_ticksls <- tick_uniquels[tick_uniquels$Species == "Ixodes scapularis", 2:7]


# absences
## absences for all
grep("no_ticks", tick_data$Your.ID.Number) %in% 
  grep("No_ticks", tick_data$Def.ID.number)

abs_filter <- c(grep("no_ticks", tick_data$Your.ID.Number),
                grep("No_ticks", tick_data$Def.ID.number))

all_abs <- tick_data[abs_filter, cols]
all_abs <- all_abs[, 2:7]

all_abs <- unique(all_abs)

## to find absences
dup1ls <- paste(tick_uniquels$Day, tick_uniquels$Month, tick_uniquels$Year, 
                tick_uniquels$Life_stage, tick_uniquels$site.code)

# aa
aaa_rule1 <- tick_uniquels$Species != "Amblyomma americanum"
aa_dup1 <- dup1ls[tick_uniquels$Species == "Amblyomma americanum"]
aa_nodup1 <- dup1ls[aaa_rule1]

rule_aa <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

aa_noticksls <- tick_uniquels[aaa_rule1 & dup1ls %in% rule_aa, 2:7]
aa_noticksls <- rbind(aa_noticksls, all_abs)
aa_noticksls$Life_stage <- NA_character_
aa_noticksls <- unique(aa_noticksls)

# am
aaa_rule1 <- tick_uniquels$Species != "Amblyomma maculatum"
aa_dup1 <- dup1ls[tick_uniquels$Species == "Amblyomma maculatum"]
aa_nodup1 <- dup1ls[aaa_rule1]

rule_am <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

am_noticksls <- tick_uniquels[aaa_rule1 & dup1ls %in% rule_am, 2:7]
am_noticksls <- rbind(am_noticksls, all_abs)
am_noticksls$Life_stage <- NA_character_
am_noticksls <- unique(am_noticksls)

# dv
aaa_rule1 <- tick_uniquels$Species != "Dermacentor variabilis"
aa_dup1 <- dup1ls[tick_uniquels$Species == "Dermacentor variabilis"]
aa_nodup1 <- dup1ls[aaa_rule1]

rule_dv <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

dv_noticksls <- tick_uniquels[aaa_rule1 & dup1ls %in% rule_dv, 2:7]
dv_noticksls <- rbind(dv_noticksls, all_abs)
dv_noticksls$Life_stage <- NA_character_
dv_noticksls <- unique(dv_noticksls)

# da
aaa_rule1 <- tick_uniquels$Species != "Dermacentor albipictus"
aa_dup1 <- dup1ls[tick_uniquels$Species == "Dermacentor albipictus"]
aa_nodup1 <- dup1ls[aaa_rule1]

rule_da <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

da_noticksls <- tick_uniquels[aaa_rule1 & dup1ls %in% rule_da, 2:7]
da_noticksls <- rbind(da_noticksls, all_abs)
da_noticksls$Life_stage <- NA_character_
da_noticksls <- unique(da_noticksls)

# is
aaa_rule1 <- tick_uniquels$Species != "Ixodes scapularis"
aa_dup1 <- dup1ls[tick_uniquels$Species == "Ixodes scapularis"]
aa_nodup1 <- dup1ls[aaa_rule1]

rule_is <- aa_nodup1[!aa_nodup1 %in% aa_dup1]

is_noticksls <- tick_uniquels[aaa_rule1 & dup1ls %in% rule_is, 2:7]
is_noticksls <- rbind(is_noticksls, all_abs)
is_noticksls$Life_stage <- NA_character_
is_noticksls <- unique(is_noticksls)

# combining data
aa_allls <- rbind(data.frame(Presence_absence = 1, aa_ticksls),
                  data.frame(Presence_absence = 0, aa_noticksls))

am_allls <- rbind(data.frame(Presence_absence = 1, am_ticksls),
                  data.frame(Presence_absence = 0, am_noticksls))

dv_allls <- rbind(data.frame(Presence_absence = 1, dv_ticksls),
                  data.frame(Presence_absence = 0, dv_noticksls))

da_allls <- rbind(data.frame(Presence_absence = 1, da_ticksls),
                  data.frame(Presence_absence = 0, da_noticksls))

is_allls <- rbind(data.frame(Presence_absence = 1, is_ticksls),
                  data.frame(Presence_absence = 0, is_noticksls))


# write files
write.csv(aa_allls, "Data/aa_presence_absence_ls.csv", row.names = FALSE)
write.csv(am_allls, "Data/am_presence_absence_ls.csv", row.names = FALSE)
write.csv(dv_allls, "Data/dv_presence_absence_ls.csv", row.names = FALSE)
write.csv(da_allls, "Data/da_presence_absence_ls.csv", row.names = FALSE)
write.csv(is_allls, "Data/is_presence_absence_ls.csv", row.names = FALSE)
# ------------------------------------------------------------------------------


# Environmental data preparation -----------------------------------------------
# variables
all_layers <- list.files("Data/raster", pattern = ".tif$", full.names = TRUE)

# get values
varnames <- c("dayl", "prcp", "srad", "swe", "tmax", "tmin", "vp")

varval <- lapply(varnames, function(x) {
  layvar <- grep(x, all_layers, value = TRUE)
  message(x)
  c(as.matrix(rast(layvar)))
})

varval1 <- do.call(cbind, varval)
varval1 <- na.omit(varval1)
dim(varval1)

colnames(varval1) <- paste0(varnames, "_mean")

# PCA
vpca <- prcomp(varval1, retx = TRUE, center = TRUE, scale. = TRUE)

summary(vpca)

# save initial results
pcadir <- "Data/PCA_results"
dir.create(pcadir)
  
save(vpca, file = "Data/PCA_results/PCA_all.RData")

sink("Data/PCA_results/PCA_summary.txt")
summary(vpca)
cat("\n\n\n")
print(vpca$rotation)
sink()

# predictions
## all Julian e days
ejdays <- seq(0, 360, 8)

yjday <- paste0(rep(c(2020, 2021, 2022), each = length(ejdays)), "_", ejdays)

pcnames <- paste0("Data/PCA_results/pcs_", yjday, ".tif")

pcpred <- lapply(1:length(pcnames), function(x) {
  varpc <- grep(paste0(yjday[x], ".tif"), all_layers, value = TRUE)
  pcs <- terra::predict(rast(varpc), vpca, filename = pcnames[x])
})
# ------------------------------------------------------------------------------


# Environmental data preparation no swe ----------------------------------------
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



# Environmental data extraction ------------------------------------------------
# dates in data
dates <- as.Date(paste(aa_all$Day, aa_all$Month, aa_all$Year, sep = "-"), 
                 format = "%d-%B-%Y")

jdays <- as.POSIXlt(dates)$yday
years <- aa_all$Year

# environmental date matching
matched_days <- sapply(jdays, function(x) {
  aday <- ejdays[ejdays <= x]
  aday[length(aday)]
})


# data extraction
## year and Julian day combination
year_jday <- paste(years, matched_days, sep = "_")

t(t(table(year_jday)))

## extracting raw variable values according to dates in a loop
uyear_jday <- unique(year_jday)

edata_extracted <- lapply(uyear_jday, function(x) {
  wweek <- which(year_jday == x)
  
  vars <- paste0("Data/raster/", varnames, "_", x, ".tif")
  varr <- rast(vars)
  
  setrec <- aa_all[wweek, ]
  varset <- extract(varr, setrec[, 2:3])[, -1]
  cbind(setrec, varset)
})

aa_all_env <- do.call(rbind, edata_extracted)

pairs(aa_all_env[, 7:13])

## extracting PCs values according to dates in a loop
edata_extracted <- lapply(uyear_jday, function(x) {
  wweek <- which(year_jday == x)
  
  vars <- paste0("Data/PCA_results/", "pcs_", x, ".tif")
  varr <- rast(vars)
  
  setrec <- aa_all[wweek, ]
  varset <- extract(varr, setrec[, 2:3])[, -1]
  cbind(setrec, varset)
})

aa_all_pcs <- do.call(rbind, edata_extracted)

pairs(aa_all_pcs[, 7:13])



## extracting PCs values according to dates in a loop (NO SWE)
edata_extracted <- lapply(uyear_jday, function(x) {
  wweek <- which(year_jday == x)
  
  vars <- paste0("Data/PCA_results_noswe/", "pcs_", x, ".tif")
  varr <- rast(vars)
  
  setrec <- aa_all[wweek, ]
  varset <- extract(varr, setrec[, 2:3])[, -1]
  cbind(setrec, varset)
})

aa_all_pcs1 <- do.call(rbind, edata_extracted)

pairs(aa_all_pcs1[, 7:12])
# ------------------------------------------------------------------------------



# Niche signal explorations ----------------------------------------------------
# multivariate test comparing all data vs positive records
colnames(aa_all_env)

vnames <- paste0(varnames, "_mean")

perma4 <- niche_signal(data = aa_all_env[, -(2:6)], 
                       condition = "Presence_absence", 
                       variables = vnames, 
                       method = "permanova")

perma4$permanova_results

# univariate test comparing all data vs positive records
uni4 <- lapply(vnames, function(x) {
  niche_signal(data = aa_all_env[, -(2:6)], condition = "Presence_absence", 
               variables = x, method = "univariate")
})

names(uni4) <- vnames


# save results
dir.create("Results")

## save objects 
save(perma4, uni4, file = "Results/aa_niche_signal_results.RData")

### univariate results as tables
u_res4 <- lapply(vnames, function(x) {
  c(variable = x, uni4[[x]]$univariate_results$hypothesis_test)
})

u_res4 <- do.call(rbind, u_res4)


### multivariate results as text'
sink("Results/aa_summary_permanova.txt", append = FALSE)
print(perma4$permanova_results)
sink()


## writing results
write.csv(u_res4, "Results/aa_summary_univariate.csv", row.names = FALSE)
# ------------------------------------------------------------------------------


# GLMs to predict risks --------------------------------------------------------
# preparing data
## check variable correlation
vcor4 <- cor(aa_all_env[, -(1:6)])

# model calibration
head(aa_all_pcs)

pcnam <- paste0("PC", 1:7)

cal_res4 <- calibration_glm(data = aa_all_pcs[, -(2:6)], 
                            dependent = "Presence_absence",
                            independent = pcnam, response_type = "lq", 
                            selection_criterion = "TSS", cv_kfolds = 10, 
                            parallel = TRUE, n_cores = 32)

cal_res4$selected

save(cal_res4, aa_all_pcs, file = "Results/glm_calibration_lq.RData")

cal_resb <- enmpa:::calibration_glmb(data = aa_all_pcs[, -(2:6)], 
                                     dependent = "Presence_absence",
                                     independent = pcnam, response_type = "lq", 
                                     selection_criterion = "TSS", cv_kfolds = 10, 
                                     exclude_bimodal = TRUE,
                                     parallel = T, n_cores = 32)

cal_resb$selected

save(cal_resb, file = "Results/to_bisec.RData")


pcnam1 <- paste0("PC", 1:6)

cal_resb1 <- enmpa:::calibration_glmb(data = aa_all_pcs1[, -(2:6)], 
                                     dependent = "Presence_absence",
                                     independent = pcnam1, response_type = "lq", 
                                     selection_criterion = "TSS", cv_kfolds = 10, 
                                     exclude_bimodal = TRUE,
                                     parallel = TRUE, n_cores = 32)

cal_resb1$selected

save(cal_resb1, aa_all_pcs1, file = "Results/glm_calibration_lq_noswe.RData")



# model predictions
# function for weighted averages
waver <- function(var_names, cal_results) {
  pcsw <- rast(var_names)
  pred4 <- predict_selected(x = cal_results, newdata = pcsw)
  
  #### Weighted average based on Akaike weights (wAIC)
  wAIC <- cal_results$selected$AIC_weight
  app(pred4$predictions*wAIC, sum)
}

## new data
pcsss <- paste0("Data/PCA_results/", "pcs_", yjday, ".tif")

## predictions per week
### figure features
mains <- paste0("A. americanum activity (Y/JD): ", gsub("_", "/",  yjday))

aa_res <- "Results/Aa_predictions"
dir.create(aa_res)

figname <-  paste0(aa_res, "/aa_wpred_", yjday, ".png")

### figures in loop
plot_gif <- lapply(1:length(mains), function(x) {
  #### Weighted average based on Akaike weights (wAIC)
  c_wmean <- waver(pcsss[(x)], cal_res4)
  
  #### plot
  png(figname[x], width = 120, height = 120, units = "mm", res = 300)
  plot(c_wmean, main = mains[x])
  plot(states, add = TRUE)
  dev.off()
})


## predictions per week (moving windoew, three week averages)
### figure features
mainsmw <- paste0("A. americanum activity MW (Y/JD): ", gsub("_", "/",  yjday))

aa_resmw <- "Results/Aa_predictions_mw"
dir.create(aa_resmw)

fignamemw <-  paste0(aa_resmw, "/aa_wpredmw_", yjday, ".png")

last <- length(mainsmw)

### figures in loop
plot_gifmw <- lapply(1:last, function(x) {
  if (x == 1) {
    wm_wmean <- app(c(waver(pcsss[(x)], cal_res4),
                      waver(pcsss[(x + 1)], cal_res4)), mean)
  } else {
    if (x == last) {
      wm_wmean <- app(c(waver(pcsss[(x - 1)], cal_res4),
                        waver(pcsss[(x)], cal_res4)), mean)
    } else {
      wm_wmean <- app(c(waver(pcsss[(x - 1)], cal_res4),
                        waver(pcsss[(x)], cal_res4),
                        waver(pcsss[(x + 1)], cal_res4)), mean)
    }
  }
  
  #### plot
  png(fignamemw[x], width = 120, height = 120, units = "mm", res = 300)
  plot(wm_wmean, main = mainsmw[x])
  plot(states, add = TRUE)
  dev.off()
})






## predictions per week (moving windoew, three week averages) no bimodal
### figure features
mainsmw <- paste0("A. americanum activity MW (Y/JD): ", gsub("_", "/",  yjday))

aa_resmw <- "Results/Aa_predictions_nobim_mw"
dir.create(aa_resmw)

fignamemw <-  paste0(aa_resmw, "/aa_wpredmw_", yjday, ".png")

last <- length(mainsmw)

### figures in loop
plot_gifmw <- lapply(1:last, function(x) {
  if (x == 1) {
    wm_wmean <- app(c(waver(pcsss[(x)], cal_resb),
                      waver(pcsss[(x + 1)], cal_resb)), mean)
  } else {
    if (x == last) {
      wm_wmean <- app(c(waver(pcsss[(x - 1)], cal_resb),
                        waver(pcsss[(x)], cal_resb)), mean)
    } else {
      wm_wmean <- app(c(waver(pcsss[(x - 1)], cal_resb),
                        waver(pcsss[(x)], cal_resb),
                        waver(pcsss[(x + 1)], cal_resb)), mean)
    }
  }
  
  #### plot
  png(fignamemw[x], width = 120, height = 120, units = "mm", res = 300)
  plot(wm_wmean, main = mainsmw[x])
  plot(states, add = TRUE)
  dev.off()
})







## predictions per week (moving window, three week averages) no bimodal no swe
### figure features
mainsmw <- paste0("A. americanum activity MW (Y/JD): ", gsub("_", "/",  yjday))

aa_resmw <- "Results/Aa_predictions_nobim_noswe_mw"
dir.create(aa_resmw)

fignamemw <-  paste0(aa_resmw, "/aa_wpredmw_", yjday, ".png")

last <- length(mainsmw)

pcsss1 <- paste0("Data/PCA_results_noswe/", "pcs_", yjday, ".tif")


### figures in loop
plot_gifmw <- lapply(1:last, function(x) {
  if (x == 1) {
    wm_wmean <- app(c(waver(pcsss1[(x)], cal_resb1),
                      waver(pcsss1[(x + 1)], cal_resb1)), mean)
  } else {
    if (x == last) {
      wm_wmean <- app(c(waver(pcsss1[(x - 1)], cal_resb1),
                        waver(pcsss1[(x)], cal_resb1)), mean)
    } else {
      wm_wmean <- app(c(waver(pcsss1[(x - 1)], cal_resb1),
                        waver(pcsss1[(x)], cal_resb1),
                        waver(pcsss1[(x + 1)], cal_resb1)), mean)
    }
  }
  
  #### plot
  png(fignamemw[x], width = 120, height = 120, units = "mm", res = 300)
  plot(wm_wmean, range = c(0, 1), main = mainsmw[x])
  plot(states, add = TRUE)
  dev.off()
})






# response curves
## variable names 
pcnam <- paste0("PC", 1:7)
pcranges <- apply(vpca$x, 2, range)
preds <- predict_selected(x = cal_res4, newdata = as.data.frame(pcranges))

x11()
par(mfrow = c(9, 14), mar = c(4, 4, .2, .2), cex = 0.3)

plots <- lapply(preds$fitted_models, function(x) {
  lapply(pcnam, function(y) {
    response_curve(model = x, variable = y, new_data = pcranges)
  })
})

enmpa:::model_selection()


preds <- predict_selected(x = cal_resb, newdata = as.data.frame(pcranges))

x11()
par(mfrow = c(5, 14), mar = c(4, 4, .2, .2), cex = 0.3)

plots <- lapply(preds$fitted_models, function(x) {
  lapply(pcnam, function(y) {
    response_curve(model = x, variable = y, new_data = pcranges)
  })
})



pcranges <- apply(vpca1$x, 2, range)
preds <- predict_selected(x = cal_resb1, newdata = as.data.frame(pcranges))

x11()
par(mfrow = c(5, 12), mar = c(4, 4, .2, .2), cex = 0.3)

plots <- lapply(preds$fitted_models, function(x) {
  lapply(pcnam1, function(y) {
    response_curve(model = x, variable = y, new_data = pcranges)
  })
})



# variable contribution
# Relative contribution of the deviance explained
varimport <- lapply(preds$fitted_models, var_importance)
# ------------------------------------------------------------------------------

