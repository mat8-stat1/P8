
##########################################################
# Calculating logarithmic and Brier score for all models #
##########################################################

library(survival)
library(splines)
library(DescTools)
library(bayestestR)
library(tidyverse)
library(MASS)
library(glmnet)

source("data_cleaning.R")
source("IPS.R")
source("RSF.R")
source("coxph.R")
source("weibull.R")

#### Function for calculating Briern scores using IPCW (Censoring distribution calculated using Kaplan-Meier)

Brier_IPCW <- function(pred.surv, data, t){
  # data$cens <- 1 - data$status
  
  n <- nrow(data)
  cens.mod <- survfit(Surv(OS, cens) ~ 1, data = data)
  
  S_t   <- summary(cens.mod, t)$surv
  S_y   <- sapply(data[,"OS"], function(x) summary(cens.mod, x)$surv)
  S_y[data[, "status"] == 0] <- 1 # to avoid dividing by 0. these S_y values are not used.
  ## fordi vi for de censurerede efter t skal bruge s_c(t)
  
  w_inv <- (data[, "OS"] > t) * 1/S_t + (data[, "OS"] <= t & data[, "status"] == 1) * 1/S_y
  delta <- (data[, "OS"] > t) + (data[, "OS"] <= t & data[, "status"] == 1)
  
  resp <- data %>% dplyr::select(status, OS) %>% 
    mutate(resp = case_when(OS <= t & status == 1 ~ 0, OS > t ~ 1, OS <= t & status == 0 ~ 100)) %>%
    pull(resp)
  ## Vi har at prædikteret ssh. er ssh. for at overleve, så vi koder observeret værdi som 1, hvis vi overlever længere
  ## end 5 år.
  
  Brier_score <- 1-2*(resp-pred.surv)^2
  ## er pred.surv ssh for at overleve, altså respons er 1 ellers skal det ombyttes ovenfor
  
  IPCW <- (1/n * delta * w_inv * Brier_score) %>% sum
}


log_IPCW <- function(pred.surv, data, t){
  # data$cens <- 1 - data$status
  
  n <- nrow(data)
  cens.mod <- survfit(Surv(OS, cens) ~ 1, data = data)
  
  S_t   <- summary(cens.mod, t)$surv
  S_y   <- sapply(data[,"OS"], function(x) summary(cens.mod, x)$surv)
  S_y[data[, "status"] == 0] <- 1 # to avoid dividing by 0. these S_y values are not used.
  ## fordi vi for de censurerede efter t skal bruge s_c(t)
  
  w_inv <- (data[, "OS"] > t) * 1/S_t + (data[, "OS"] <= t & data[, "status"] == 1) * 1/S_y
  delta <- (data[, "OS"] > t) + (data[, "OS"] <= t & data[, "status"] == 1)
  
  resp <- data %>% dplyr::select(status, OS) %>% 
    mutate(resp = case_when(OS <= t & status == 1 ~ 0, OS > t ~ 1, OS <= t & status == 0 ~ 100)) %>%
    pull(resp)
  ## Vi har at prædikteret ssh. er ssh. for at overleve, så vi koder observeret værdi som 1, hvis vi overlever længere
  ## end 5 år.
  
  log_score <- data.frame(resp, pred.surv) %>% 
    mutate(log_score = case_when(resp == 1 ~ log(pred.surv), 
                                 resp == 0 ~ log(1-pred.surv),
                                 resp == 100 ~ 1000)) %>% 
    pull(log_score)
  ## er pred.surv ssh for at overleve, altså respons er 1 ellers skal det ombyttes ovenfor
  
  IPCW <- (1/n * delta * w_inv * log_score) %>% sum
}


####################### RSF model:
wh       <- length(pred_RSF$time.interest[pred_RSF$time.interest <= 5]) # this time point corresponds to estimated 5-year survival
pred.RSF <- pred_RSF$survival[, wh] # estimated 5-year surivival

Brier.RSF  <- Brier_IPCW(pred.RSF, test.imp, 5)
log.RSF <- log_IPCW(pred.RSF, test.imp, 5)

######################## IPS model
test.IPS <- case_when(test.imp$IPS %in% 5:7 ~ 5L, TRUE ~ test.imp$IPS)

sum_IPS  <- summary(IPS, times = 5)
pred.IPS <- sum_IPS$surv[1]*(test.IPS==0) + sum_IPS$surv[2]*(test.IPS==1) + sum_IPS$surv[3]*(test.IPS==2) +
  sum_IPS$surv[4]*(test.IPS==3) + sum_IPS$surv[5]*(test.IPS==4) + sum_IPS$surv[6]*(test.IPS==5) 

Brier.IPS <- Brier_IPCW(pred.IPS, test.imp, 5)
log.IPS <- log_IPCW(pred.IPS, test.imp, 5)

######################## Cox
## AIC + BIC models
cox.models <- list(coxph.aic, coxph.aic.spl, coxph.bic, coxph.bic.spl)
Brier.coxph  <- list()
log.coxph <- list()

k <- 1
for (j in cox.models){
  pred.coxph <- rep(NA, nrow(test.imp))
  for (i in 1:nrow(test.imp)){
    pred          <- survfit(j, newdata = test.imp[i,])
    pred.coxph[i] <- summary(pred, 5)$surv
  }
  
  Brier.coxph[[k]] <- Brier_IPCW(pred.coxph, test.imp, 5)
  log.coxph[[k]] <- log_IPCW(pred.coxph, test.imp, 5)
  k <- k + 1
}

## CV + LASSO model
# Standardising training data (and saving transformations for test data)
covars          <- colnames(test.imp)[c(3:21, 50)]

covars.std <- data.matrix(data.imp[ ,covars])
covars.std[, c("sex", "bsym", "LDH")] <- covars.std[, c("sex", "bsym", "LDH")] - 1

centers <- attr(scale(covars.std[, !(colnames(covars.std) %in% c("sex", "bsym", "LDH"))]), "scaled:center")
scales  <- attr(scale(covars.std[, !(colnames(covars.std) %in% c("sex", "bsym", "LDH"))]), "scaled:scale")

covars.std[, !(colnames(covars.std) %in% c("sex", "bsym", "LDH"))] <- scale(covars.std[, !(colnames(covars.std) %in% c("sex", "bsym", "LDH"))])

# Using same transformation on test data
covars.std.test <- data.matrix(test.imp[ ,covars])
covars.std.test[, c("sex", "bsym", "LDH")] <- covars.std.test[, c("sex", "bsym", "LDH")] - 1
for (i in names(scales)){
  covars.std.test[, i] <- sapply(covars.std.test[, i], function(x) (x - centers[i]) / scales[i])
}

# Making predictions from model
pred.coxph <- rep(NA, nrow(test.imp))
for (i in 1:nrow(test.imp)){
  pred          <- survfit(coxph.lasso, s = "lambda.min", x = covars.std, y = Surv(data.imp$OS, data.imp$status), newx = covars.std.test[i, ])
  pred.coxph[i] <- summary(pred, 5)$surv
}
Brier.coxph[[5]] <- Brier_IPCW(pred.coxph, test.imp, 5)
log.coxph[[5]] <- log_IPCW(pred.coxph, test.imp, 5)

######################## Weibull models:
## Predictions with the lasso model (estimated t-year survival):
S_wei.lasso <- function(newdata, t){ 
  wei_sum <- sapply(1:nrow(newdata), function(x){sum(weibull.lasso$beta[-1,]/weibull.lasso$scale*newdata[x,])})
  exp(-exp(-weibull.lasso$beta[1,]/weibull.lasso$scale - wei_sum)*t^(1/weibull.lasso$scale))
}

## Function for estimating t-year survival fct. for survreg Weibull objects (remaining models)

pred.weibull <- function(mod, data, t){
  
  out <- rep(NA, nrow(data))
  pct <- 1:99000/100000
  
  for (i in 1:nrow(data)){
    ptime  <- predict(mod, newdata = data[i,], type='quantile', p=pct, se=F)
    wh     <- which(ptime == ptime[ptime<5][length(ptime[ptime<5])])
    out[i] <- 1 - pct[wh]
  }
  
  out
}

## Using this fct. to save predictions for AIC/BIC Weibull models
weibull.models <- list(weibull.aic, weibull.aic.spl, weibull.bic, weibull.bic.spl)
Brier.wei <- list()
log.wei <- list()

k <- 1
for (j in weibull.models){
  pred.wei <- pred.weibull(j, test.imp, 5)
  Brier.wei[[k]] <- Brier_IPCW(pred.wei, test.imp, 5)
  log.wei[[k]] <- log_IPCW(pred.wei, test.imp, 5)
  k <- k + 1
}

# Making predictions based on CV+LASSO model (using fct. from "weibull.R")
pred.wei     <- S_wei.lasso(covars.std.test, 5)
Brier.wei[[5]] <- Brier_IPCW(pred.wei, test.imp, 5)
log.wei[[5]] <- log_IPCW(pred.wei, test.imp, 5)
