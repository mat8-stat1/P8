
####################################################################################
# Functions for fitting proposed models and predicting based on these models, and  #
# loss function (misclassification error) and IPCW estimated AUC                   #
####################################################################################


library(survival)
library(randomForestSRC)
library(bayestestR)

#### Sourcing necessary files (loading all model objects)
source("data_cleaning.R")
source("IPS.R")
source("weibull.R")
source("coxph.R")
source("RSF.R")


#### Defining IPC weights based on full data.imp / test.imp (for t=5)

t <- 5

## Training data weights/delta
cens.mod <- survfit(Surv(OS, cens) ~ 1, data = data.imp)

S_t   <- summary(cens.mod, t)$surv
S_y   <- sapply(data.imp[,"OS"], function(x) summary(cens.mod, x)$surv)
S_y[data.imp[, "status"] == 0] <- 1 # to avoid dividing by 0. these S_y values are not used.

w_inv.data <- (data.imp[, "OS"] > t) * 1/S_t + (data.imp[, "OS"] <= t & data.imp[, "status"] == 1) * 1/S_y
delta.data <- (data.imp[, "OS"] > t) + (data.imp[, "OS"] <= t & data.imp[, "status"] == 1)

## Test data weights/delta
cens.mod <- survfit(Surv(OS, cens) ~ 1, data = test.imp)

S_t   <- summary(cens.mod, t)$surv
S_y   <- sapply(test.imp[,"OS"], function(x) summary(cens.mod, x)$surv)
S_y[test.imp[, "status"] == 0] <- 1 # to avoid dividing by 0. these S_y values are not used.

w_inv.test <- (test.imp[, "OS"] > t) * 1/S_t + (test.imp[, "OS"] <= t & test.imp[, "status"] == 1) * 1/S_y
delta.test <- (test.imp[, "OS"] > t) + (test.imp[, "OS"] <= t & test.imp[, "status"] == 1)




#### Functions for efficiently fitting and predicting 5-year survival

## IPS

fit_pred.IPS <- function(data.fit, data.test, formula = NULL, model = NULL){
  if( !(is.null(formula)) ){
    model <- survfit(formula(paste0("Surv(OS, status) ~ ", formula)), data = data.fit)
  }
  sum_IPS  <- summary(model, times = 5)
  test.IPS <- data.test[, "IPS2"]
  pred.IPS <- sum_IPS$surv[1]*(test.IPS==0) + sum_IPS$surv[2]*(test.IPS==1) + sum_IPS$surv[3]*(test.IPS==2) +
              sum_IPS$surv[4]*(test.IPS==3) + sum_IPS$surv[5]*(test.IPS==4) + sum_IPS$surv[6]*(test.IPS==5)
  pred.IPS
} 



## Weibull
fit_pred.weibull <- function(data.fit, data.test, formula=NULL,  model=NULL, t){
  if(!(is.null(formula))){
    model <- survreg(formula(paste0("Surv(OS, status) ~ ", formula)), data = data.fit)
  }
  pct <- 1:990/1000
  
  ptime <- sapply(1:nrow(data.test), function(i){predict(model, newdata = data.test[i,], type='quantile', p=pct, se=F)} )
  wh <-  sapply(1:ncol(ptime), function(i){which(ptime[,i] == ptime[,i][(ptime[,i])<t][length(ptime[,i][(ptime[,i])<t])])}) 
  
  out <- 1 - pct[wh]
  out
}

# Example:
# fit_pred.weibull(model=weibull.aic.spl, data.test=test.imp, t=5)



### Weibull LASSO
# NB: Needs standardised data. Model should always be specified for extracting lambda.

fit_pred.weibull.lasso <- function(data.fit, covars.std.fit = covars.std, covars.std.test = covars.std.test, fit=F,  model=weibull.lasso, t){
  if(fit){
    best.lambda <- model$lambda
    model <- iregnet(covars.std.fit, Surv(data.fit$OS, data.fit$status), family = "weibull", standardize = F, lambda=best.lambda)
  }
  wei_sum <- sapply(1:nrow(covars.std.test), function(x){sum(model$beta[-1,]/model$scale*covars.std.test[x,])})
  
  exp(-exp(-model$beta[1,]/model$scale - wei_sum)*t^(1/model$scale))
}

# Example: Using standardised test data
# ------------
# 
# covars <- colnames(test.imp)[c(3:21, 50)]
# 
# covars.std <- data.matrix(data.imp[ ,covars])
# covars.std[, c("sex", "bsym", "LDH")] <- covars.std[, c("sex", "bsym", "LDH")] - 1
# 
# centers <- attr(scale(covars.std[, !(colnames(covars.std) %in% c("sex", "bsym", "LDH"))]), "scaled:center")
# scales  <- attr(scale(covars.std[, !(colnames(covars.std) %in% c("sex", "bsym", "LDH"))]), "scaled:scale")
# 
# covars.std[, !(colnames(covars.std) %in% c("sex", "bsym", "LDH"))] <- scale(covars.std[, !(colnames(covars.std) %in% c("sex", "bsym", "LDH"))])
# 
# # Using same transformation on test data
# covars.std.test <- data.matrix(test.imp[ ,covars])
# covars.std.test[, c("sex", "bsym", "LDH")] <- covars.std.test[, c("sex", "bsym", "LDH")] - 1
# for (i in names(scales)){
#   covars.std.test[, i] <- sapply(covars.std.test[, i], function(x) (x - centers[i]) / scales[i])
# }
# 
# fit_pred.weibull.lasso(data.imp,covars.std.fit = covars.std, covars.std.test = covars.std.test, fit=T,  model=weibull.lasso, 5)
# ------------



## Coxph

fit_pred.coxph <- function(data.fit, data.test, formula = NULL, model = NULL){
  if( !(is.null(formula))){
    model <- coxph(formula(paste0("Surv(OS, status) ~ ", formula)), data = data.fit)
  }
  pred.coxph <- as.vector(summary(survfit(model, newdata = data.test), 5)$surv)
  pred.coxph
}



## Coxph LASSO
# NB: Needs standardised data. Model should always be specified for extracting lambda.

fit_pred.coxph.lasso <- function(data.fit, covars.std.fit = covars.std, covars.std.test = covars.std.test, fit = T, model = coxph.lasso, t = 5){
  
  if( fit ){
    best.lambda <- model$lambda
    model       <- glmnet(covars.std.fit, Surv(data.fit$OS, data.fit$status), family = "cox", standardize = F, lambda = best.lambda)
  }
  pred.coxph.lasso <- as.vector(summary(survfit(model, s = "lambda.min", x = covars.std.fit, y = Surv(data.fit$OS, data.fit$status), newx = covars.std.test), t)$surv)
  pred.coxph.lasso
}

# Example: Using standardised test data
# ------------
# 
# covars <- colnames(test.imp)[c(3:21, 50)]
# 
# covars.std <- data.matrix(data.imp[ ,covars])
# covars.std[, c("sex", "bsym", "LDH")] <- covars.std[, c("sex", "bsym", "LDH")] - 1
# 
# centers <- attr(scale(covars.std[, !(colnames(covars.std) %in% c("sex", "bsym", "LDH"))]), "scaled:center")
# scales  <- attr(scale(covars.std[, !(colnames(covars.std) %in% c("sex", "bsym", "LDH"))]), "scaled:scale")
# 
# covars.std[, !(colnames(covars.std) %in% c("sex", "bsym", "LDH"))] <- scale(covars.std[, !(colnames(covars.std) %in% c("sex", "bsym", "LDH"))])
# 
# # Using same transformation on test data
# covars.std.test <- data.matrix(test.imp[ ,covars])
# covars.std.test[, c("sex", "bsym", "LDH")] <- covars.std.test[, c("sex", "bsym", "LDH")] - 1
# for (i in names(scales)){
#   covars.std.test[, i] <- sapply(covars.std.test[, i], function(x) (x - centers[i]) / scales[i])
# }
# 
# fit_pred.coxph(data.fit = data.imp, covars.std.fit = covars.std, covars.std.test = covars.std.test, fit = T, model = coxph.lasso)
#
# ------------




## RSF

# Data preparation: Hardcoding NA to be a factor
# ------------
full_data <- full_data %>% 
  mutate(Stadium..Ann.Arbor. = as.character(Stadium..Ann.Arbor.)) %>% 
  mutate(Stadium..Ann.Arbor. = case_when(is.na(Stadium..Ann.Arbor.) ~ "Not Available", T ~ Stadium..Ann.Arbor.)) %>% 
  mutate(Stadium..Ann.Arbor. = as.factor(Stadium..Ann.Arbor.))
full_data <- full_data %>% 
  mutate(ECOG = as.character(ECOG)) %>% 
  mutate(ECOG = case_when(is.na(ECOG) ~ "Not Available", T ~ ECOG)) %>% 
  mutate(ECOG = as.factor(ECOG))
full_data <- full_data %>% 
  mutate(bsym = as.character(bsym)) %>% 
  mutate(bsym = case_when(is.na(bsym) ~ "Not Available", T ~ bsym)) %>% 
  mutate(bsym = as.factor(bsym))

data <- full_data[(1:n_train),]
test <- full_data[((n_train+1):(n_train+n_test)),]
# ------------

# NB: data.fit and data.test needs to be non-imputed data.
fit_pred.RSF <- function(data.fit, data.test, fit = T, model = RSF){
  if( fit ){
    model <- rfsrc(Surv(OS, status) ~ . - OS - status, data = data.fit, na.action = "na.impute")
  }
  pred.RSF <- predict(model, newdata = data.test, na.action = "na.impute")
  wh       <- length(pred.RSF$time.interest[pred.RSF$time.interest <= 5]) # this time point corresponds to estimated 5-year survival
  pred.RSF <- pred.RSF$survival[, wh] # estimated 5-year surivival
  pred.RSF
}

# Example:
#fit_pred.RSF(data, test)









#### Loss functions (AUC and misclassification error)

## AUC function

AUC <- function(pred.surv, w_inv, delta, data, t){
  n <- nrow(data)
  alpha <- c(-1,sort(pred.surv) + c(diff(sort(pred.surv))/2, 1))
  
  TPR_numerator <- sapply(alpha, function(x) sum(1/n * delta * w_inv * (pred.surv > x & data[, "OS"] > t)))
  FPR_numerator <- sapply(alpha, function(x) sum(1/n * delta * w_inv * (pred.surv > x & data[, "OS"] <= t)))
  
  TPR_denominator <- sum(1/n * delta * w_inv * (data[, "OS"] > t))
  FPR_denominator <- sum(1/n * delta * w_inv * (data[, "OS"] <= t))
  
  TPR <- TPR_numerator / TPR_denominator
  FPR <- FPR_numerator / FPR_denominator
  
  AUC <- area_under_curve(FPR, TPR, method = "step")
  AUC
}

## Misclassification error
mis_error <- function(pred.surv, test.data, w_inv, delta, t, alpha=0.5){
  true_val <- test.data[,"OS"]< t & test.data[,"status"] == 1
  pred_t <- pred.surv < alpha
  
  1/nrow(test.data)*sum(w_inv*delta*(true_val!=pred_t))
}

