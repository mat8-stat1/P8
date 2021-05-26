
############################
# Training Weibull models  #
############################

library(tidyverse)
library(splines)
library(survival)
library(MASS)
library(glmnet)

source("data_cleaning.R")

covars <- colnames(data.imp)[c(3:21, 50)]
full.formula <- formula(paste0("Surv(OS, status) ~ ", paste(covars, collapse = " + ")))
null.model <- survreg(Surv(OS, status) ~ 1, data = data.imp)
full.model <- survreg(full.formula, data = data.imp)


# AIC forward/backward/both selection

mod.cand <- list()
mod.cand[[1]] <- stepAIC(null.model, full.formula, direction = "forward")
mod.cand[[2]] <- stepAIC(full.model, full.formula, direction = "backward")
mod.cand[[3]] <- stepAIC(null.model, full.formula, direction = "both")
mod.cand[[4]] <- stepAIC(full.model, full.formula, direction = "both")

min.aic <- which.min(lapply(mod.cand, AIC))

weibull.aic <- mod.cand[[min.aic]]


# BIC forward/backward/both selection

mod.cand <- list()
mod.cand[[1]] <- stepAIC(null.model, full.formula, direction = "forward", k = log(nrow(data.imp)))
mod.cand[[2]] <- stepAIC(full.model, full.formula, direction = "backward", k = log(nrow(data.imp)))
mod.cand[[3]] <- stepAIC(null.model, full.formula, direction = "both", k = log(nrow(data.imp)))
mod.cand[[4]] <- stepAIC(full.model, full.formula, direction = "both", k = log(nrow(data.imp)))

min.bic <- which.min(lapply(mod.cand, BIC))

weibull.bic <- mod.cand[[min.bic]]


# Cross validation and LASSO
library(penaltyLearning)
library(iregnet)
library(Rcpp)

## Data is scaled (for LASSO) and transformed to matrix-format (required by iregnet)
covars.std <- data.matrix(data.imp[ ,covars])
covars.std[, c("sex", "bsym", "LDH")] <- covars.std[, c("sex", "bsym", "LDH")] - 1
covars.std[, !(colnames(covars.std) %in% c("sex", "bsym", "LDH"))] <- scale(covars.std[, !(colnames(covars.std) %in% c("sex", "bsym", "LDH"))])

## LASSO for Weibull
Weibull.lasso <- iregnet(covars.std, Surv(data.imp$OS, data.imp$status), family = "weibull", standardize = F)

## Cross validation
set.seed(5)
lambda      <- seq(0.01,0.0001, by = -0.0001)
nfold       <- 10
sample.data <- sample(rep(1:nfold, l=nrow(covars.std)))
loglik.list <- as.data.frame(matrix(NA, nrow=10, ncol=length(lambda)))
for(i in 1:10){
  testIndex <- which(sample.data==i,arr.ind=TRUE)
  testData  <- covars.std[testIndex, ]
  trainData <- covars.std[-testIndex, ]
  
  data.imp.train <- data.imp[-testIndex,]
  data.imp.test  <- data.imp[testIndex,]
  lasso          <- iregnet(trainData, Surv(data.imp.train$OS, data.imp.train$status), family = "weibull", standardize = F, lambda=lambda)
  for(j in 1:length(lambda)){
    loglik <- rep(NA,nrow(testData))
    for (k in 1:nrow(testData)) {
      loglik[k] <- data.imp.test$status[k] * (lasso$beta[1,j] + sum(lasso$beta[-1,j] * testData[k,]) - log(lasso$scale[j]) + (1 / lasso$scale[j] - 1) * log(data.imp.test$OS[k])) - 
        exp(lasso$beta[1,j] + sum(lasso$beta[-1,j] * testData[k,])) * data.imp.test$OS[k]^(1 / lasso$scale[j])
    }
  loglik.list[i,j] <- sum(loglik)
  }
}

best.lambda <- lambda[which.max(base::colMeans(loglik.list))]

weibull.lasso <- iregnet(covars.std, Surv(data.imp$OS, data.imp$status), family = "weibull", standardize = F, lambda=best.lambda)


# Function for estimating survival based on fitted model
S_weibull <- function(t){ 
  exp(-exp(-weibull.lasso$beta[1,]/weibull.lasso$scale-sum(weibull.lasso$beta[-1,]/weibull.lasso$scale*covars.std[5,]))*t^(1/weibull.lasso$scale))
}












########################################
# SPLINE BASED MODELS (AIC and BIC only)
########################################


# ----------------
# AIC forward selection
# ----------------

# Defining covariates for spline and non-spline modelling. All spline-vars are checked for as non-spline covars.
covars         <- colnames(data.imp)[c(3:21, 50)]
spline.vars    <- c("Alder", "nodal", "Hemoglobin", "Albumin", "Lymfocytter.mia.L.", "Leukocytter.mia.L.", "Beta2Microglobulin",
                    "IgG", "IgM", "IgA", "CreatineMikroMol", "Thrombocytter.mia.L.", "organ.inv", "tumordiameter")
nonspline.vars <- c("sex", "Stadium..Ann.Arbor.", "ECOG", "bsym", "ENODAL", "LDH")
max.spline     <- 4 # largest spline df's allowed

# spline.df (for keeping track of each of splinevars df):
  # 0: covariate is not part of current.formula
  # 1: covariate is part of current.formula, but not modelled by spline
  # 2-4: covariate is part of current.formula, with spline of degree spline.df
spline.df           <- as.data.frame(matrix(rep(0, length(spline.vars), ncol=1)))
rownames(spline.df) <- spline.vars

# Initialisation: null model
  # cand.formula: The current model at each forward step
  # cand.formula: The candidate formula, one for each additional variable at each step
  # best.formula: The formula of the model with lowest AIC of all cand.formula at each step
current.formula <- cand.formula <- best.formula <- "Surv(OS, status) ~ 1"
current.model   <- coxph(formula(current.formula), data = data.imp)

best.models <- list()
best.aics   <- list()

#for (i in (1:(length(nonspline.vars) + length(spline.vars)*max.spline))){
for (i in 1:20){  
  
  current.formula <- best.formula
  best.aic        <- Inf # Initialising: Best AIC at each step i
  
  # Run through all covariates to check for improved fit:
  for (j in 1:length(covars)) {
    
    # Non-spline vars: check whether adding other covariates improves fit.
    if ((covars[j] %in% nonspline.vars) & !(str_detect(current.formula, fixed(covars[j])))){ # ensures that the covariate is not already part of the current model
      cand.formula <- paste0(current.formula, " + ", covars[j])
    }
    
    # Spline-vars: check whether adding one df improves fit.
    if (covars[j] %in% spline.vars){
      if (spline.df[covars[j], ] == 0){ # If not part of current model, model as non-spline covar.
        cand.formula <- paste0(current.formula, " + ", covars[j])
      }
      if (spline.df[covars[j], ] == 1){ # If part of current model as non-spline covar, model as spline with degree 2.
        cand.formula <- str_remove(current.formula, fixed(paste0(" + ", covars[j])))
        cand.formula <- paste0(cand.formula, " + ns(", covars[j], ", df = 2)")
      }
      if (spline.df[covars[j], ] > 1 & spline.df[covars[j], ] < max.spline){ # If part of current model, model as spline with degree + 1 (if not already of degree 4).
        cand.formula <- str_remove(current.formula, fixed(paste0(" + ns(", covars[j], ", df = ", spline.df[covars[j], ], ")")))
        cand.formula <- paste0(cand.formula, " + ns(", covars[j], ", df = ", spline.df[covars[j], ] + 1, ")")
      }
    }
    
    # Defining the candidate model (where j'th covariate is added)
    cand.model <- survreg(formula(cand.formula), data = data.imp)
    cand.aic   <- AIC(cand.model)
    # cand.aic   <- AIC(cand.model, k = log(nrow(data.imp))) # for BIC
    
    # Check whether adding the covariate improves fit; if true, define this as the new current model.
    if (cand.aic < best.aic){
      
      best.aic     <- cand.aic
      best.model   <- cand.model
      best.formula <- cand.formula
      
      spline.df.temp           <- as.data.frame(matrix(rep(0, length(spline.vars), ncol=1)))
      rownames(spline.df.temp) <- spline.vars
      
      # Adding 1 df to the list of df's of variables in the model. Only temp since better model might show up.
      if (covars[j] %in% spline.vars){
        spline.df.temp[covars[j], ] <- 1
      }
    }
    
  }
  
  # Updating the list of df's of spline variables in the model
  spline.df <- spline.df + spline.df.temp
  
  # Save the current best model and AIC's in separate list
  best.models[[i]] <- best.model
  best.aics[[i]]   <- best.aic
}


weibull.aic.spl <- best.models[[which.min(unlist(best.aics))]]



# ----------------
# BIC forward selection (same procedure as above)
# ----------------

# Defining covariates for spline and non-spline modelling. All spline-vars are checked for as non-spline covars.
covars         <- colnames(data.imp)[c(3:21, 50)]
spline.vars    <- c("Alder", "nodal", "Hemoglobin", "Albumin", "Lymfocytter.mia.L.", "Leukocytter.mia.L.", "Beta2Microglobulin",
                    "IgG", "IgM", "IgA", "CreatineMikroMol", "Thrombocytter.mia.L.", "organ.inv", "tumordiameter")
nonspline.vars <- c("sex", "Stadium..Ann.Arbor.", "ECOG", "bsym", "ENODAL", "LDH")
max.spline     <- 4 # largest spline df's allowed

# spline.df (for keeping track of each of splinevars df):
  # 0: covariate is not part of current.formula
  # 1: covariate is part of current.formula, but not modelled by spline
  # 2-4: covariate is part of current.formula, with spline of degree spline.df
spline.df           <- as.data.frame(matrix(rep(0, length(spline.vars), ncol=1)))
rownames(spline.df) <- spline.vars

# Initialisation: null model
  # cand.formula: The current model at each forward step
  # cand.formula: The candidate formula, one for each additional variable at each step
  # best.formula: The formula of the model with lowest AIC of all cand.formula at each step
current.formula <- cand.formula <- best.formula <- "Surv(OS, status) ~ 1"
current.model   <- coxph(formula(current.formula), data = data.imp)

best.models <- list()
best.aics   <- list()

#for (i in (1:(length(nonspline.vars) + length(spline.vars)*max.spline))){
for (i in 1:20){  
  
  current.formula <- best.formula
  best.aic        <- Inf # Initialising: Best AIC at each step i
  
  # Run through all covariates to check for improved fit:
  for (j in 1:length(covars)) {
    
    # Non-spline vars: check whether adding other covariates improves fit.
    if ((covars[j] %in% nonspline.vars) & !(str_detect(current.formula, fixed(covars[j])))){ # ensures that the covariate is not already part of the current model
      cand.formula <- paste0(current.formula, " + ", covars[j])
    }
    
    # Spline-vars: check whether adding one df improves fit.
    if (covars[j] %in% spline.vars){
      if (spline.df[covars[j], ] == 0){ # If not part of current model, model as non-spline covar.
        cand.formula <- paste0(current.formula, " + ", covars[j])
      }
      if (spline.df[covars[j], ] == 1){ # If part of current model as non-spline covar, model as spline with degree 2.
        cand.formula <- str_remove(current.formula, fixed(paste0(" + ", covars[j])))
        cand.formula <- paste0(cand.formula, " + ns(", covars[j], ", df = 2)")
      }
      if (spline.df[covars[j], ] > 1 & spline.df[covars[j], ] < max.spline){ # If part of current model, model as spline with degree + 1 (if not already of degree 4).
        cand.formula <- str_remove(current.formula, fixed(paste0(" + ns(", covars[j], ", df = ", spline.df[covars[j], ], ")")))
        cand.formula <- paste0(cand.formula, " + ns(", covars[j], ", df = ", spline.df[covars[j], ] + 1, ")")
      }
    }
    
    # Defining the candidate model (where j'th covariate is added)
    cand.model <- survreg(formula(cand.formula), data = data.imp)
    # cand.aic   <- AIC(cand.model)
    cand.aic   <- AIC(cand.model, k = log(nrow(data.imp)))
    
    # Check whether adding the covariate improves fit; if true, define this as the new current model.
    if (cand.aic < best.aic){
      
      best.aic     <- cand.aic
      best.model   <- cand.model
      best.formula <- cand.formula
      
      spline.df.temp           <- as.data.frame(matrix(rep(0, length(spline.vars), ncol=1)))
      rownames(spline.df.temp) <- spline.vars
      
      # Adding 1 df to the list of df's of variables in the model. Only temp since better model might show up.
      if (covars[j] %in% spline.vars){
        spline.df.temp[covars[j], ] <- 1
      }
    }
    
  }
  
  # Updating the list of df's of spline variables in the model
  spline.df <- spline.df + spline.df.temp
  
  # Save the current best model and AIC's in separate list
  best.models[[i]] <- best.model
  best.aics[[i]]   <- best.aic
}


weibull.bic.spl <- best.models[[which.min(unlist(best.aics))]]





# Predictions for Weibull models (here 5-y survival estimation for patient 5)

pct <- 1:99000/100000
ptime <- predict(weibull.bic.spl, newdata = data.imp[5,], type='quantile', p=pct, se=F)
plot(x=ptime, y=1-pct, type="l")

wh <- which(ptime == ptime[ptime<5][length(ptime[ptime<5])])
1-pct[wh] # 5y survival estimate based on model weibull.bic.spl


