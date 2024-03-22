################################################################################
# Project: RII Track-2 FEC: Marshalling Diverse Big Data Streams to Understand 
#          Risk of Tick-Borne Diseases in the Great Plains
# Title: Risk modeling for A. americanum life stages
# Authors: Marlon E. Cobos
# Date: 21/03/2024 (dd/mm/yyyy)
# Funding: NSF EPSCoR (IIA-1920946)
################################################################################



# Description ------------------------------------------------------------------
# This script serves to perform ecological niche modeling analyses and model
# projections, as well as MOP analyses. This is done for all live stages of 
# A. americanum: Larva, Nymph, and Adult.
#
# 1. Model calibration and selection is performed.
# 2. Selected models are fitted and projected to all weeks 2020-2022.
# 3. Month averages of projections are obtained.
# 4. MOP analyses are performed for all weeks 2020-2022.
# 5. MOP results are averaged and monthly representations are obtained. 
# ------------------------------------------------------------------------------



# Packages ---------------------------------------------------------------------
# install packages
#install.packages("geodata")
#install.packages("enmpa")

# load packages
library(geodata)  # it loads terra which we use in many lines
library(enmpa)
# ------------------------------------------------------------------------------



# Data to be used in ENM and plots ---------------------------------------------
# working directory
setwd("YOUR/DIRECTORY")  # change as needed (contains the sub folder Data)

# state borders (Download, store in folder Data)
states <- gadm("USA", level = 1, path = "Data")

## select local states for outlines shown on final plot
selected <- c("Kansas", "Oklahoma", "Missouri")
pstates <- states[states$NAME_1 %in% selected, ]

# tick data
aals_all_env <- read.csv("Data/aa_presence_absence_env_ls.csv")
aals_all_pcs <- read.csv("Data/aa_presence_absence_pcs_ls.csv")

## separate data for analyses according to life stage
l_stage <- na.omit(unique(aals_all_env$Life_stage))

### data with raw variables
aa_stage_env <- lapply(l_stage, function(x) {
  cond <- aals_all_env$Life_stage == x | is.na(aals_all_env$Life_stage)
  newdata <- aals_all_env[cond, ]
  
  exabs <- aals_all_env$Life_stage != x & !is.na(aals_all_env$Life_stage)
  exabs <- aals_all_env[exabs, ]
  
  exabsf <- !do.call(paste, exabs[, 2:6]) %in% do.call(paste, newdata[, 2:6])
  exabs <- exabs[exabsf, ]
  exabs$Presence_absence <- 0
  
  newdata <- rbind(newdata, exabs)
})
names(aa_stage_env) <- l_stage

### data with PCs
aa_stage_pcs <- lapply(l_stage, function(x) {
  cond <- aals_all_pcs$Life_stage == x | is.na(aals_all_pcs$Life_stage)
  newdata <- aals_all_pcs[cond, ]
  
  exabs <- aals_all_pcs$Life_stage != x & !is.na(aals_all_pcs$Life_stage)
  exabs <- aals_all_pcs[exabs, ]
  
  exabsf <- !do.call(paste, exabs[, 2:6]) %in% do.call(paste, newdata[, 2:6])
  exabs <- exabs[exabsf, ]
  exabs$Presence_absence <- 0
  
  newdata <- rbind(newdata, exabs)
})

names(aa_stage_pcs) <- l_stage
# ------------------------------------------------------------------------------


# Niche signal explorations ----------------------------------------------------
# multivariate test comparing all data vs positive records
colnames(aals_all_env)

vnames <- colnames(aals_all_env)[8:13]

aals_perma <- lapply(aa_stage_env, function(x) {
  niche_signal(data = x, condition = "Presence_absence", 
               variables = vnames, method = "permanova")
})


aals_perma$Nymph$permanova_results


# univariate test comparing all data vs positive records
aals_univ <- lapply(aa_stage_env, function(w) {
  uni <- lapply(vnames, function(x) {
    niche_signal(data = w, condition = "Presence_absence", 
                 variables = x, method = "univariate")
  })
  names(uni) <- vnames
  uni
})

aals_univ$Nymph$dayl_mean$univariate_results$hypothesis_test


# save results
#dir.create("Results")  # folder created in "Risk_modeling_Aa.R"

## save objects 
save(aals_perma, aals_univ, file = "Results/aals_niche_signal_results.RData")

### univariate results as tables
aals_univa <- lapply(1:length(aals_univ), function(w) {
  sumuni <- lapply(vnames, function(x) {
    c(variable = x, aals_univ[[w]][[x]]$univariate_results$hypothesis_test)
  })
  
  cbind(life_stage = l_stage[w], do.call(rbind, sumuni))
})

aals_univa <- do.call(rbind, aals_univa)

### multivariate results as text'
sink("Results/aals_summary_permanova.txt", append = FALSE)
lapply(aals_perma, function(x) {print(x$permanova_results)})
sink()

## writing results
write.csv(aals_univa, "Results/aals_summary_univariate.csv", row.names = FALSE)
# ------------------------------------------------------------------------------



# GLMs to predict risks --------------------------------------------------------
# model calibration
head(aals_all_pcs)

pcnam <- colnames(aals_all_pcs)[8:13]

aals_glm_cal <- lapply(aa_stage_pcs, function(x) {
  calibration_glm(data = x, dependent = "Presence_absence",
                  independent = pcnam, response_type = "lq", 
                  formula_mode = "intensive",
                  selection_criterion = "TSS", cv_kfolds = 10, 
                  exclude_bimodal = TRUE, parallel = TRUE, 
                  n_cores = 32)  # change number of cores as needed
})


# glms of selected models
aals_sel_glms <- lapply(aals_glm_cal, function(x) {
  fit_selected(x)
})

# variable contribution
# Relative contribution of the deviance explained
aals_var_imp <- lapply(aals_sel_glms, var_importance)

x11()
plot_importance(aals_var_imp$Larva)
plot_importance(aals_var_imp$Nymph)
plot_importance(aals_var_imp$Adult)

# save glm results
save(aals_glm_cal, aals_sel_glms, aals_var_imp, 
     file = "Results/aals_glm_results.RData")


# variable response curves
plt <- lapply(aals_sel_glms, function(x) {
  x11()
  par(mfrow = c(2, 3), mar = c(4, 4, .5, .2), cex = 0.8)
  
  plots <- lapply(pcnam, function(y) {
    response_curve(x, variable = y, extrapolate = TRUE, ylim = c(0, 1))
  })
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

## looooop for all life stages

## pattern of species and life stage name and directories
sppwrt <- paste0("/aa_", names(aals_glm_cal), "_pred_")
sppdir <- paste0("Results/Aa_", names(aals_glm_cal), "_predictions")


predss <- lapply(1:length(aals_glm_cal), function(w) {
  ## new folder for predictions
  aals_res <- sppdir[w]
  dir.create(aals_res)
  
  ## prediction file names
  yweeks <- gsub("Data/PCA_results/pcs_", "", pcsss)
  
  pred_name <- paste0(aals_res, sppwrt[w], yweeks)
  
  ## predictions per week
  pred_yweeks <- lapply(1:length(pcsss), function(x) {
    ### Weighted average based on AIC weights (wAIC)
    c_wmean <- waver(pcsss[x], aals_glm_cal[[w]])
    
    ### write
    writeRaster(c_wmean, pred_name[x])
    
    message(x, " of ", length(pcsss))
  })
  
  
  # three week moving window averages
  ## new folder for predictions
  aals_resmw <- paste0(sppdir[w], "_mw")
  dir.create(aals_resmw)
  
  ## prediction file names
  pred_namemw <- paste0(aals_resmw, sppwrt[w], yweeks)
  
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

  ## directory for results
  mondir <- paste0(sppdir[w], "_monthly")
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
      moname <- paste0(mondir, sppwrt[w], x, '_', y, '.tif')
      app(morast, mean, filename = moname)
    })
    
    message(x)
  })
  
  
  # normal monthly averages for every 10 years
  ## new directory
  mndir <- paste0(sppdir[w], "_monthly_normals")
  dir.create(mndir)
  
  ## decades to consider
  decs <- seq(1980, 2020, 10)
  
  ## average calculation in loo
  norms <- lapply(1:(length(decs) - 1), function(x) {
    mres <- lapply(1:12, function(y) {
      ### file names for averaging
      filesm <- paste0(mondir, sppwrt[w], decs[x]:(decs[x + 1] - 1), 
                       "_", y, ".tif")
      
      ### new file name
      avname <- paste0(mndir, sppwrt[w], decs[x], "s_", y, 
                       "_normal.tif")
      
      ### average and writing
      app(rast(filesm), mean, filename = avname)
    })
    
    message(x, " of ", (length(decs) - 1))
  })
  
  
  # normal monthly standard deviations for every 10 years
  ## new directory
  mndirsd <- paste0(sppdir[w], "_monthly_normals_SD")
  dir.create(mndirsd)
  
  ## average calculation in loo
  normssd <- lapply(1:(length(decs) - 1), function(x) {
    mres <- lapply(1:12, function(y) {
      ### file names for averaging
      filesm <- paste0(mondir, sppwrt[w], decs[x]:(decs[x + 1] - 1), 
                       "_", y, ".tif")
      
      ### new file name
      sdname <- paste0(mndirsd, sppwrt[w], decs[x], "s_", y, 
                       "_normal_SD.tif")
      
      ### average and writing
      app(rast(filesm), sd, filename = sdname)
    })
    
    message(x, " of ", (length(decs) - 1))
  })
})
# ------------------------------------------------------------------------------
