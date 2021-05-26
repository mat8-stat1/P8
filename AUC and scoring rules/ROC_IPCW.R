
#################################################
# IPCW estimated ROC-curves for all models, and #
# IPCW estimated scores for all models          #
#################################################


library(survival)
library(bayestestR)
library(splines)
library(glmnet)
library(ggplot2)
library(reshape2)
library(viridis)
library(scales)


#### Sourcing necessary files and loading all model objects
source("data_cleaning.R")
source("Scoring rule.R")
source("IPS.R")
source("coxph.R")
source("weibull.R")
source("RSF.R")

#### Function for calculating ROC curves using IPCW (Censoring distribution calculated using Kaplan-Meier)

ROC_IPCW <- function(pred.surv, data, t){
  n <- nrow(data)
  cens.mod <- survfit(Surv(OS, cens) ~ 1, data = data)
  
  S_t   <- summary(cens.mod, t)$surv
  S_y   <- sapply(data[,"OS"], function(x) summary(cens.mod, x)$surv)
  S_y[data[, "status"] == 0] <- 1 # to avoid dividing by 0. these S_y values are not used.
  
  w_inv <- (data[, "OS"] > t) * 1/S_t + (data[, "OS"] <= t & data[, "status"] == 1) * 1/S_y
  delta <- (data[, "OS"] > t) + (data[, "OS"] <= t & data[, "status"] == 1)
  alpha <- c(-1,sort(pred.surv) + c(diff(sort(pred.surv))/2, 1))
  
  TPR_numerator <- sapply(alpha, function(x) sum(1/n * delta * w_inv * (pred.surv > x & data[, "OS"] > t)))
  FPR_numerator <- sapply(alpha, function(x) sum(1/n * delta * w_inv * (pred.surv > x & data[, "OS"] <= t)))
  
  TPR_denominator <- sum(1/n * delta * w_inv * (data[, "OS"] > t))
  FPR_denominator <- sum(1/n * delta * w_inv * (data[, "OS"] <= t))
  
  TPR <- TPR_numerator / TPR_denominator
  FPR <- FPR_numerator / FPR_denominator
  
  AUC <- area_under_curve(FPR, TPR, method = "step")
  
  out <- list(TPR = TPR, FPR = FPR, AUC = AUC, alpha = alpha)
  out
}



#### IPS

test.IPS <- case_when(test.imp$IPS %in% 5:7 ~ 5L, TRUE ~ test.imp$IPS)

sum_IPS  <- summary(IPS, times = 5)
pred.IPS <- sum_IPS$surv[1]*(test.IPS==0) + sum_IPS$surv[2]*(test.IPS==1) + sum_IPS$surv[3]*(test.IPS==2) +
            sum_IPS$surv[4]*(test.IPS==3) + sum_IPS$surv[5]*(test.IPS==4) + sum_IPS$surv[6]*(test.IPS==5) 

roc.IPS <- ROC_IPCW(pred.IPS, test.imp, 5)




#### Cox

## AIC + BIC models
cox.models <- list(coxph.aic, coxph.aic.spl, coxph.bic, coxph.bic.spl)
roc.coxph  <- list()

k <- 1
for (j in cox.models){
  pred.coxph <- rep(NA, nrow(test.imp))
  for (i in 1:nrow(test.imp)){
    pred          <- survfit(j, newdata = test.imp[i,])
    pred.coxph[i] <- summary(pred, 5)$surv
  }
  
  roc.coxph[[k]] <- ROC_IPCW(pred.coxph, test.imp, 5)
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
roc.coxph[[5]] <- ROC_IPCW(pred.coxph, test.imp, 5)


#Plotting ROC curves
cols <- c("dodgerblue3", "lightblue", "darkorange", "goldenrod1", "forestgreen")

plot(roc.coxph[[1]]$FPR, roc.coxph[[1]]$TPR, type = "l", lwd = 2, col = cols[1],
     xlab = "False positive rate", ylab = "True positive rate")
for (i in 2:5){
  lines(roc.coxph[[i]]$FPR, roc.coxph[[i]]$TPR, lwd = 2, col = cols[i])
}
abline(0,1)
legend(.6, .4,
       legend = c("coxph.aic", "coxph.aic.spl", "coxph.bic", "coxph.bic.spl", "coxph.lasso"),
       col = cols,
       lty = 1,
       lwd = 2)



#### Weibull

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
    wh     <- which(ptime == ptime[ptime<t][length(ptime[ptime<t])])
    out[i] <- 1 - pct[wh]
  }
  
  out
}

## Using this fct. to save predictions for AIC/BIC Weibull models
weibull.models <- list(weibull.aic, weibull.aic.spl, weibull.bic, weibull.bic.spl)
roc.wei        <- list()

k <- 1
for (j in weibull.models){
  pred.wei <- pred.weibull(j, test.imp, 5)
  roc.wei[[k]] <- ROC_IPCW(pred.wei, test.imp, 5)
  k <- k + 1
}

# Making predictions based on CV+LASSO model (using fct. from "weibull.R")
pred.wei     <- S_wei.lasso(covars.std.test, 5)
roc.wei[[5]] <- ROC_IPCW(pred.wei, test.imp, 5)


#Plotting ROC curves
cols <- c("dodgerblue3", "lightblue", "darkorange", "goldenrod1", "forestgreen")

plot(roc.wei[[1]]$FPR, roc.wei[[1]]$TPR, type = "l", lwd = 2, col = cols[1],
     xlab = "False positive rate", ylab = "True positive rate")
for (i in 2:5){
  lines(roc.wei[[i]]$FPR, roc.wei[[i]]$TPR, lwd = 2, col = cols[i])
}
abline(0,1)
legend(.6, .4,
       legend = c("weibull.aic", "weibull.aic.spl", "weibull.bic", "weibull.bic.spl", "weibull.lasso"),
       col = cols,
       lty = 1,
       lwd = 2)


### Weibull + Cox ROC plot

par(mfrow = c(1,2))

plot(roc.wei[[1]]$FPR, roc.wei[[1]]$TPR, type = "l", lwd = 2, col = cols[1],
     xlab = "False positive rate", ylab = "True positive rate", main = "a)")
for (i in 2:5){
  lines(roc.wei[[i]]$FPR, roc.wei[[i]]$TPR, lwd = 2, col = cols[i])
}
abline(0,1)
legend(.5, .4,
       legend = c("weibull.aic", "weibull.aic.spl", "weibull.bic", "weibull.bic.spl", "weibull.lasso"),
       col = cols,
       lty = 1,
       lwd = 2)

plot(roc.coxph[[1]]$FPR, roc.coxph[[1]]$TPR, type = "l", lwd = 2, col = cols[1],
     xlab = "False positive rate", ylab = "True positive rate", main = "b)")
for (i in 2:5){
  lines(roc.coxph[[i]]$FPR, roc.coxph[[i]]$TPR, lwd = 2, col = cols[i])
}
abline(0,1)
legend(.5, .4,
       legend = c("coxph.aic", "coxph.aic.spl", "coxph.bic", "coxph.bic.spl", "coxph.lasso"),
       col = cols,
       lty = 1,
       lwd = 2)

par(mfrow = c(1,1))



#### RSF

wh       <- length(pred_RSF$time.interest[pred_RSF$time.interest <= 5]) # this time point corresponds to estimated 5-year survival
pred.RSF <- pred_RSF$survival[, wh] # estimated 5-year surivival
roc.RSF  <- ROC_IPCW(pred.RSF, test.imp, 5)



#### Plot: Best AIC/BIC models + RSF + IPS
rocs <- list(roc.IPS, roc.wei[[2]], roc.coxph[[1]], roc.RSF)
cols <- c("dodgerblue3", "darkorange", "forestgreen", "blueviolet")

plot(rocs[[1]]$FPR, rocs[[1]]$TPR, type = "l", lwd = 2, col = cols[1],
     xlab = "False positive rate", ylab = "True positive rate")
for (i in 2:4){
  lines(rocs[[i]]$FPR, rocs[[i]]$TPR, lwd = 2, col = cols[i])
}
abline(0,1)
legend(.6, .4,
       legend = c("IPS model", "weibull.aic.spl", "coxph.aic", "RSF model"),
       col = cols,
       lty = 1,
       lwd = 2)






##### Table: AUC + scores for all models

mod.names <- c("IPS model",
               "coxph.aic", "coxph.aic.spl", "coxph.bic", "coxph.bic.spl", "coxph.lasso",
               "weibull.aic", "weibull.aic.spl", "weibull.bic", "weibull.bic.spl", "weibull.lasso",
               "RSF model")

AUC.table           <- matrix(nrow = length(mod.names), ncol = 4)
colnames(AUC.table) <- c("Model", "AUC", "Logarithmic score", "Brier score")
AUCs                <- c(roc.IPS$AUC,
                         roc.coxph[[1]]$AUC, roc.coxph[[2]]$AUC, roc.coxph[[3]]$AUC, roc.coxph[[4]]$AUC, roc.coxph[[5]]$AUC,
                         roc.wei[[1]]$AUC, roc.wei[[2]]$AUC, roc.wei[[3]]$AUC, roc.wei[[4]]$AUC, roc.wei[[5]]$AUC,
                         roc.RSF$AUC)
logs                <- c(log.IPS,
                         log.coxph[[1]], log.coxph[[2]], log.coxph[[3]], log.coxph[[4]], log.coxph[[5]],
                         log.wei[[1]], log.wei[[2]], log.wei[[3]], log.wei[[4]], log.wei[[5]],
                         log.RSF)
Briers              <- c(Brier.IPS,
                          Brier.coxph[[1]], Brier.coxph[[2]], Brier.coxph[[3]], Brier.coxph[[4]], Brier.coxph[[5]],
                          Brier.wei[[1]], Brier.wei[[2]], Brier.wei[[3]], Brier.wei[[4]], Brier.wei[[5]],
                          Brier.RSF)

AUC.table[, 1] <- mod.names
AUC.table[, 2] <- sapply(AUCs, function(x) round(x, 4))
AUC.table[, 3] <- sapply(logs, function(x) round(x, 4))
AUC.table[, 4] <- sapply(Briers, function(x) round(x, 4))


##### Plot: AUC + scores for all models
AUC.table.plot <- matrix(nrow = length(mod.names), ncol = 4)
colnames(AUC.table.plot) <- c("Model", "AUC", "Logarithmic score + 1", "Brier score")

AUC.table.plot <- AUC.table.plot %>%
  as.data.frame %>% 
  mutate("Model" = mod.names) %>% 
  mutate("AUC" = AUCs) %>% 
  mutate("Logarithmic score + 1" = logs + 1) %>% 
  mutate("Brier score" = Briers)

AUC.table2 <- melt(AUC.table.plot, id = "Model", variable.name = "Score")
AUC.table2$Model <- factor(AUC.table2$Model, levels = mod.names[c(12, 6:2, 11:7, 1)])
AUC.table2$Score <- factor(AUC.table2$Score, levels = c("AUC", "Logarithmic score + 1", "Brier score"))

cols <- viridis(8)

ggplot(data=AUC.table2, aes(x = Model, y = value, fill = Score)) + 
  geom_bar(stat="identity", position = "dodge", width = .8) + 
  #scale_fill_manual(" ", values = c("True error" = "black", "Training error" = "forestgreen", "CV_10" = "darkorange", "Err^(1)" = "dodgerblue3", "Err^(.632)" = "deepskyblue", "Err^(.632+)" = "lightblue"), limits = c("True error", "Training error", "CV_10", "Err^(1)", "Err^(.632)", "Err^(.632+)")) +
  scale_fill_manual(" ", values = c("AUC" = cols[3], "Logarithmic score + 1" = cols[5], "Brier score" = cols[7]), limits = c("Brier score", "Logarithmic score + 1", "AUC")) +
  xlab("") + ylab("Score") +
  scale_y_continuous(limits=c(0.63,0.85),oob = rescale_none) +
  coord_flip() +
  theme_minimal()


