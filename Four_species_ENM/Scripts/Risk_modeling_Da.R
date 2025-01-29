################################################################################
# Project: RII Track-2 FEC: Marshalling Diverse Big Data Streams to Understand 
#          Risk of Tick-Borne Diseases in the Great Plains
# Title: Risk modeling for D. albipictus
# Authors: Marlon E. Cobos
# Date: 29/01/2025 (dd/mm/yyyy)
# Funding: NSF EPSCoR (IIA-1920946)
################################################################################



# Description ------------------------------------------------------------------
# This script serves to perform ecological niche modeling analyses and model
# projections, as well as MOP analyses. 
#
# 1. Model calibration and selection is performed.
# 2. Selected models are fitted and projected to all weeks 2020-2022.
# 3. Month averages of projections are obtained.
# 4. MOP analyses are performed for all weeks 2020-2022.
# 5. MOP results are averaged and monthly representations are obtained. 
# ------------------------------------------------------------------------------



# Packages ---------------------------------------------------------------------
# install packages
install.packages("geodata")
install.packages("enmpa")
install.packages("biosurvey") # for now only in github # remotes::install_github("claununez/biosurvey")

# load packages
library(geodata)  # it loads terra which we use in many lines
library(enmpa)
library(mop)
library(biosurvey)  # only for color palettes, alpha colors, and legend bars 
# ------------------------------------------------------------------------------



# Data to be used in ENM and plots ---------------------------------------------
# world countries
wrd <- world(resolution = 2, path = ".")

# state borders (Download, store in folder Data)
states <- gadm("USA", level = 1, path = "Data")

## select local states for outlines shown on final plot
selected <- c("Kansas", "Oklahoma", "Missouri")
pstates <- states[states$NAME_1 %in% selected, ]

sastates <- states[states$NAME_1 %in% selected[-3], ]

# tick data
da_all_env <- read.csv("Data/da_presence_absence_env.csv")
da_all_pcs <- read.csv("Data/da_presence_absence_pcs.csv")
# ------------------------------------------------------------------------------


# Niche signal explorations ----------------------------------------------------
# multivariate test comparing all data vs positive records
colnames(da_all_env)

vnames <- colnames(da_all_env)[7:12]

da_perma <- niche_signal(data = da_all_env, condition = "Presence_absence", 
                         variables = vnames, method = "permanova")

da_perma$permanova_results


# univariate test comparing all data vs positive records
da_univ <- lapply(vnames, function(x) {
  niche_signal(data = da_all_env, condition = "Presence_absence", 
               variables = x, method = "univariate")
})

names(da_univ) <- vnames


# save results
#dir.create("Results")

## save objects 
save(da_perma, da_univ, file = "Results/da_niche_signal_results.RData")

### univariate results as tables
da_univa <- lapply(vnames, function(x) {
  c(variable = x, da_univ[[x]]$univariate_results$hypothesis_test)
})

da_univa <- do.call(rbind, da_univa)

### multivariate results as text'
sink("Results/da_summary_permanova.txt", append = FALSE)
print(da_perma$permanova_results)
sink()

## writing results
write.csv(da_univa, "Results/da_summary_univariate.csv", row.names = FALSE)
# ------------------------------------------------------------------------------



# GLMs to predict risks --------------------------------------------------------
# model calibration
head(da_all_pcs)

pcnad <- colnames(da_all_pcs)[7:12]

da_glm_cal <- calibration_glm(data = da_all_pcs, dependent = "Presence_absence",
                              independent = pcnad, response_type = "lq", 
                              selection_criterion = "TSS", cv_kfolds = 9, 
                              exclude_bimodal = TRUE, parallel = TRUE, 
                              n_cores = 32)

da_glm_cal$selected


# glms of selected models
da_sel_glms <- fit_selected(da_glm_cal)


# variable contribution
# Relative contribution of the deviance explained
da_var_imp <- var_importance(da_sel_glms)

plot_importance(da_var_imp)

# save glm results
save(da_glm_cal, da_sel_glms, da_var_imp, file = "Results/da_glm_results.RData")


# variable response curves
x11()
par(mfrow = c(2, 3), mar = c(4, 4, .5, .2), cex = 0.8)

plots <- lapply(pcnad, function(y) {
  response_curve(fitted = da_sel_glms, variable = y)
})
# ------------------------------------------------------------------------------


# GLMs prediction and averaging -----------------------------------------------
# model predictions
## new data (layers where models will be projected, only names now)
pcsss <- list.files("Data/PCA_results", pattern = ".tif$", full.names = TRUE)

## sorting files by year and day
yrs <- as.numeric(gsub(".*_(\\d{4})_\\d*.tif$", "\\1", pcsss))
day <- as.numeric(gsub(".*_\\d{4}_(\\d*).tif$", "\\1", pcsss))

pcsss <- pcsss[order(yrs, day)]

## new folder for predictions
da_res <- "Results/Da_predictions"
dir.create(da_res)

## prediction file names
yweeks <- gsub("Data/PCA_results/pcs_", "", pcsss)

pred_name <- paste0(da_res, "/da_pred_", yweeks)

## predictions per week
pred_yweeks <- lapply(1:length(pcsss), function(x) {
  ### Weighted average based on AIC weights (wAIC)
  c_wmean <- predict_selected(fitted = da_sel_glms, newdata = rast(pcsss[x]), 
                              consensus = TRUE)
  
  ### write
  writeRaster(c_wmean$consensus$Weighted_average, pred_name[x])
  
  message(x, " of ", length(pcsss))
})


# three week moving window averages
## new folder for predictions
da_resmw <- "Results/Da_predictions_mw"
dir.create(da_resmw)

## prediction file names
pred_namemw <- paste0(da_resmw, "/da_pred_", yweeks)

## weekly averages using moving window in loop
last <- length(pcsss)

pred_yweeksmw <- lapply(1:last, function(x) {
  if (x == 1) {
    ### average when it is the first week
    app(rast(c(pred_name[x], pred_name[x + 1])), mean,
        filename = pred_namemw[x])
  } else {
    if (x == last) {
      ### average when it is the last week
      app(rast(c(pred_name[x - 1], pred_name[x])), mean,
          filename = pred_namemw[x])
    } else {
      ### average when it is every other week
      app(rast(c(pred_name[x - 1], pred_name[x], pred_name[x + 1])), mean,
          filename = pred_namemw[x])
    }
  }
  
  message(x, " of ", length(pcsss))
})


# monthly averages
## years and days to be considered in names
years <- 1980:2022
wda <- seq(0, 364, 8)

## days in months
dmon <- c(31, NA, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

## pattern of species name
sppatt <- "da"

## directory for results
mondir <- "Results/Da_predictions_monthly"
dir.create(mondir)

## monthly averages in loop
res <- lapply(years, function(x) {
  ### predictions for a complete year
  predx <- rast(grep(x, pred_name, value = TRUE, fixed = TRUE))
  names(predx) <- wda
  
  ### running averages for each month in each year
  mres <- lapply(1:12, function(y) {
    ### number of day in the month
    if (y == 2) {
      mdays <- as.numeric(difftime(as.Date(paste0(x, "-03-01")), 
                                   as.Date(paste0(x, "-02-01"))))
    } else {
      mdays <- dmon[y]
    }
    
    ### number indicating Julian day at the initial and last day of the month
    iday <- as.POSIXlt(paste0(x, "-", y, "-", 1), format = "%Y-%m-%d")$yday
    fday <- as.POSIXlt(paste0(x, "-", y, "-", mdays), format = "%Y-%m-%d")$yday
    
    ### detecting layers at the initial, middly, and last parts of the month
    ones <- wda[wda >= iday & wda <= (fday - 8)]
    ini <- wda[wda >= (iday - 8) & wda < iday]
    fin <- wda[wda > (fday - 8) & wda <= fday]
    
    ### grouping middle and initial predictions (weighting initial if needed)
    if (length(ini) > 0) {
      iweight <- (ini + 8 - iday) / 8
      if (iweight == 1) {
        morast <- c(predx[[as.character(ini)]], predx[[as.character(ones)]])
      } else {
        morast <- c(predx[[as.character(ini)]] * iweight, 
                    predx[[as.character(ones)]])
      }
    } else {
      morast <- predx[[as.character(ones)]]
    }
    
    ### grouping previous ones and final predictions (weighting final if needed)
    if (length(fin) > 0) {
      fweight <- (fday - fin + 1) / 8
      if (fweight == 1) {
        morast <- c(morast, predx[[as.character(fin)]])
      } else {
        morast <- c(morast, predx[[as.character(fin)]] * fweight)
      }
    } 
    
    ### preparing monthly averages
    moname <- paste0(mondir, "/", sppatt, "_pred_", x, '_', y, '.tif')
    app(morast, mean, filename = moname)
  })
  
  message(x)
})


# normal monthly averages for every 10 years
## new directory
mndir <- "Results/Da_predictions_monthly_normals"
dir.create(mndir)

## decades to consider
decs <- seq(1980, 2020, 10)

## average calculation in loo
norms <- lapply(1:(length(decs) - 1), function(x) {
  mres <- lapply(1:12, function(y) {
    ### file names for averaging
    filesm <- paste0(mondir, "/", sppatt, "_pred_", decs[x]:(decs[x + 1] - 1), 
                     "_", y, ".tif")
    
    ### new file name
    avname <- paste0(mndir, "/", sppatt, "_pred_", decs[x], "s_", y, 
                     "_normal.tif")
    
    ### average and writing
    app(rast(filesm), mean, filename = avname)
  })
  
  message(x, " of ", (length(decs) - 1))
})


# normal monthly standard deviations for every 10 years
## new directory
mndirsd <- "Results/Da_predictions_monthly_normals_SD"
dir.create(mndirsd)

## average calculation in loo
normssd <- lapply(1:(length(decs) - 1), function(x) {
  mres <- lapply(1:12, function(y) {
    ### file names for averaging
    filesm <- paste0(mondir, "/", sppatt, "_pred_", decs[x]:(decs[x + 1] - 1), 
                     "_", y, ".tif")
    
    ### new file name
    sdname <- paste0(mndirsd, "/", sppatt, "_pred_", decs[x], "s_", y, 
                     "_normal_SD.tif")
    
    ### average and writing
    app(rast(filesm), sd, filename = sdname)
  })
  
  message(x, " of ", (length(decs) - 1))
})
# ------------------------------------------------------------------------------



# MOP: extrapolation areas -----------------------------------------------------
## new folder for MOP results
da_mop <- "Results/Da_MOP"
dir.create(da_mop)

## mop file names
mop_name <- paste0(da_mop, "/da_mop_", yweeks)

## mop per week
mop_yweeks <- lapply(1:length(pcsss), function(x) {
  ### Weighted average based on AIC weights (wAIC)
  mopb <- mop(m = da_all_pcs[, -(1:6)], g = rast(pcsss[x]), type = "basic", 
              calculate_distance = FALSE)
  
  ### write
  writeRaster(mopb$mop_basic, mop_name[x])
  
  message(x, " of ", length(pcsss))
})


# monthly averages
## years and days to be considered in names
years <- 1980:2022
wda <- seq(0, 364, 8)

## days in months
dmon <- c(31, NA, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)

## pattern of species name
sppatt <- "da"

## directory for results
mopmdir <- "Results/Da_MOP_monthly"
dir.create(mopmdir)

## monthly averages in loop
res <- lapply(years, function(x) {
  ### predictions for a complete year
  predx <- rast(grep(x, mop_name, value = TRUE, fixed = TRUE))
  names(predx) <- wda
  
  ### replace NA by 0 in mops
  predx <- subst(predx, NA, 0)
  
  ### running averages for each month in each year
  mres <- lapply(1:12, function(y) {
    ### number of day in the month
    if (y == 2) {
      mdays <- as.numeric(difftime(as.Date(paste0(x, "-03-01")), 
                                   as.Date(paste0(x, "-02-01"))))
    } else {
      mdays <- dmon[y]
    }
    
    ### number indicating Julian day at the initial and last day of the month
    iday <- as.POSIXlt(paste0(x, "-", y, "-", 1), format = "%Y-%m-%d")$yday
    fday <- as.POSIXlt(paste0(x, "-", y, "-", mdays), format = "%Y-%m-%d")$yday
    
    ### detecting layers at the initial, middly, and last parts of the month
    ones <- wda[wda >= iday & wda <= (fday - 8)]
    ini <- wda[wda >= (iday - 8) & wda < iday]
    fin <- wda[wda > (fday - 8) & wda <= fday]
    
    ### grouping middle and initial predictions (weighting initial if needed)
    if (length(ini) > 0) {
      iweight <- (ini + 8 - iday) / 8
      if (iweight == 1) {
        morast <- c(predx[[as.character(ini)]], predx[[as.character(ones)]])
      } else {
        morast <- c(predx[[as.character(ini)]] * iweight, 
                    predx[[as.character(ones)]])
      }
    } else {
      morast <- predx[[as.character(ones)]]
    }
    
    ### grouping previous ones and final predictions (weighting final if needed)
    if (length(fin) > 0) {
      fweight <- (fday - fin + 1) / 8
      if (fweight == 1) {
        morast <- c(morast, predx[[as.character(fin)]])
      } else {
        morast <- c(morast, predx[[as.character(fin)]] * fweight)
      }
    } 
    
    ### preparing monthly averages
    moname <- paste0(mopmdir, "/", sppatt, "_mop_", x, '_', y, '.tif')
    
    morast <- app(morast, mean)
    morast <- subst(morast, 0, NA)
    
    writeRaster(morast, filename = moname)
  })
  
  message(x)
})
# ------------------------------------------------------------------------------