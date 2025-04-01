###########################################################################################
##################### Scripts for single-trial analysis####################################
## This work © 2025 by Nick Fradgley, CSIRO is licensed under CC BY 4.0 ###################

#   'plot.data' is a data frame with all trial plot data. Column names include:
#   'Env' Factor or environment names
#   'Grain.Yield..t.ha.' Grain yield of each plot numeric
#   'ROW and COL row and column coordinates factors
#   'GID' is the ID factor for genotypes

library(ASRtriala)
library(asreml)
library(data.table)
library(stringr)

# Make empty lists to fill with fitted model outputs for each trial
All.BLUPs.and.BLUEs <- list()
All.best.mods <- list()

envs <- unique(plot.data$Env)

for (i in 1:length(envs)) { # For each trial environment
  trial_sub <- plot.data[plot.data$Env == envs[i], ] # subset data for one trial environment
  trial_sub <- fill.grid(data = trial_sub, row = "ROW", col = "COL", order = TRUE, message = TRUE) # Fill in any missing plots in row column grid
  # Make list of all possible row and column random effects combinations to try
  Rcovs <- list(
    "GRC" = as.formula(~ GID + ROW + COL),
    "GR" = as.formula(~ GID + ROW),
    "Gc" = as.formula(~ GID + COL),
    "G" = as.formula(~GID)
  )
  # Make list of all possible row and column residual spatial effects combinations to try with id, ar1 and ar2
  Rescovs <- list(
    "idar1" = as.formula(~ id(ROW):ar1(COL)),
    "ar1id" = as.formula(~ ar1(ROW):id(COL)),
    "ar1ar1" = as.formula(~ ar1(ROW):ar1(COL)),
    "ar2ar1" = as.formula(~ ar2(ROW):ar1(COL)),
    "ar1ar2" = as.formula(~ ar1(ROW):ar2(COL)),
    "ar2ar2" = as.formula(~ ar2(ROW):ar2(COL)),
    "units" = as.formula(~units)
  )

  asreml.options(workspace = "1028mb", trace = TRUE) # Set up ASREML workspace
  # Fit random genotype effects model with every combination of random and residual covariance effects combinations
  all.lmms <- list()
  for (R in 1:length(Rcovs)) {
    for (r in 1:length(Rescovs)) {
      asmod <- asreml(
        fixed = Grain.Yield..t.ha. ~ 1,
        random = Rcovs[[R]],
        residual = Rescovs[[r]],
        data = trial_sub,
        na.action = list(x = "include", y = "include")
      )
      max.updates <- 10
      while (asmod$converge == F & max.updates > 1) { # make sure they converge ok
        asmod <- update(asmod)
        max.updates <- max.updates - 1
      }
      all.lmms[[length(all.lmms) + 1]] <- asmod
      names(all.lmms)[length(all.lmms)] <- paste(names(Rcovs)[R], names(Rescovs)[r], sep = "_")
    }
  }
  # Pick the best model with the lowest AIC value
  best.mod <- all.lmms[[which.min(unlist(lapply(all.lmms, function(x) summary(x)$aic[1])))]]
  # Get genotype BLUPs from model coefficients
  preds <- summary(best.mod, coef = TRUE)$coef.random
  preds <- as.data.frame(preds[grepl(pattern = "GID", rownames(preds)), ])
  preds$solution <- preds$solution + as.numeric(best.mod$coefficients$fixed)
  Gpred <- preds
  Gpred$GID <- factor(gsub(pattern = "GID_", replacement = "", x = rownames(Gpred)))
  # Get genotype BLUPs from model predictions
  predictions <- predict.asreml(best.mod, classify = "GID", only = "GID", vcov = T)
  # Calculate accuracy and reliability of BLUPs
  Gpred$weights <- diag(solve(predictions$vcov))
  Gpred$SEMwts <- 1 / (Gpred$std.error^2)
  Gpred$PEV <- Gpred$std.error^2
  Va <- summary(best.mod)$varcomp["GID", "component"]
  Gpred$Accuracy <- sqrt(1 - Gpred$PEV / Va)
  Gpred$Reliability <- 1 - Gpred$PEV / Va
  colnames(Gpred)[colnames(Gpred) == "solution"] <- "BLUPs"


  # Fit fixed genotype effects model with the same random and residual covariance effects as for the best random model
  which.res.cov <- c(sapply(1:length(Rescovs), function(x) rep(x, length(Rcovs))))[which.min(unlist(lapply(all.lmms, function(x) summary(x)$aic[1])))]
  which.R.cov <- c(sapply(1:length(Rcovs), function(x) rep(x, length(Rescovs))))[which.min(unlist(lapply(all.lmms, function(x) summary(x)$aic[1])))]
  Rcovs.noGID <- lapply(Rcovs, function(x) as.formula(paste("~", gsub("GID \\+ ", "", as.character(x)[2]))))

  if (Rcovs.noGID[[which.R.cov]] == "~GID") {
    fixasmod <- asreml(
      fixed = Grain.Yield..t.ha. ~ 1 + GID,
      data = trial_sub,
      na.action = list(x = "include", y = "include")
    )
    max.updates <- 10
    while (fixasmod$converge == F & max.updates > 1) {
      fixasmod <- update(fixasmod)
      max.updates <- max.updates - 1
    }
  } else {
    fixasmod <- asreml(
      fixed = Grain.Yield..t.ha. ~ 1 + GID,
      random = Rcovs.noGID[[which.R.cov]],
      residual = Rescovs[[which.res.cov]],
      data = trial_sub,
      na.action = list(x = "include", y = "include")
    )
    max.updates <- 10
    while (fixasmod$converge == F & max.updates > 1) {
      fixasmod <- update(fixasmod)
      max.updates <- max.updates - 1
    }
  }
  # Get the BLUEs for model predictions
  predictions <- predict.asreml(fixasmod, classify = "GID", vcov = T)
  blues <- predictions$pvals$predicted.value
  Gpred$BLUEsweights <- diag(solve(predictions$vcov))
  Gpred$BLUEsSEMwts <- 1 / (Gpred$std.error^2)
  Gpred$BLUEsPEV <- Gpred$std.error^2
  Gpred$BLUEs <- blues
  Gpred <- cbind(
    "Year" = substr(envs[i], 1, 4),
    "Site" = substr(envs[i], 6, nchar(as.character(envs[i]))),
    "Env" = envs[i], Gpred
  )

  All.best.mods[[i]] <- best.mod # Add best model to full list
  names(All.best.mods)[i] <- as.character(envs[i])
  All.BLUPs.and.BLUEs[[i]] <- Gpred # Add all BLUPs and BLUEs to full list
}
All.BLUPs.and.BLUEs <- as.data.frame(rbindlist(All.BLUPs.and.BLUEs))

# Get mean reliability of BLUPs at each environment
mean.reliability <- sapply(unique(All.BLUPs.and.BLUEs$Env), function(x) mean(All.BLUPs.and.BLUEs$Reliability[All.BLUPs.and.BLUEs$Env == x]))
names(mean.reliability) <- unique(All.BLUPs.and.BLUEs$Env)

# Get variance components at each environment
all.env.varcomps <- lapply(All.best.mods, function(x) summary(x)$varcomp)
names(all.env.varcomps) <- names(All.best.mods)

# Get a subset of bad trials with low reliability to exclude from further analysis
bad.trials <- names(which(mean.reliability < .3))
All.BLUPs.and.BLUEs.skimmed <- All.BLUPs.and.BLUEs[!All.BLUPs.and.BLUEs$Env %in% bad.trials, ]
