
#############################################################
# Error estimation tables and figures (using training data) #
#############################################################

##### Training error, cross validation error, Err^(1), Err^(.632) and Err^(.632+)
##### Loss functions used for error estimates: Misclassfification error and AUC (only training error and CV_10)

library(ggplot2)
library(reshape2)
library(viridis)
library(grid)
library(gridExtra)
library(xtable)

##### Sourcing necessary functions for calculating error estimates
source("fcts-fit-pred-loss.R")
source("error-fcts.R")

##### Defining model objects in lists
models <- list(IPS,
               weibull.aic, weibull.aic.spl, weibull.bic, weibull.bic.spl, weibull.lasso,
               coxph.aic, coxph.aic.spl, coxph.bic, coxph.bic.spl, coxph.lasso,
               RSF)
mod.names <- c("IPS model",
               "weibull.aic", "weibull.aic.spl", "weibull.bic", "weibull.bic.spl", "weibull.lasso",
               "coxph.aic", "coxph.aic.spl", "coxph.bic", "coxph.bic.spl", "coxph.lasso",
               "RSF model")

##### Table: Prediction error estimates of each model (AUC)

### Error estimates are calculated
## Defining standardised covariates
# ----------
covars <- colnames(test.imp)[c(3:21, 50)]

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
# ----------

true.AUCs     <- training_error(data = test.imp, data.RSF = test, w_inv = w_inv.test, delta = delta.test, models = models, type = "AUC", 
                                t = 5, alpha = 0.5, covars.std = covars.std.test, covars.std.test = covars.std.test)

train.err.AUC <- training_error(data = data.imp, data.RSF = data, w_inv = w_inv.data, delta = delta.data, models = models, type = "AUC", 
                                t = 5, alpha = 0.5, covars.std = covars.std, covars.std.test = covars.std)

CV_10.AUC     <- CV_K(data = data.imp, data.RSF = data, nfold = 10, delta = delta.data, models = models, type = "AUC", t = 5, alpha = 0.5)

predtable.AUC           <- matrix(nrow = length(mod.names), ncol = 4)
colnames(predtable.AUC) <- c("Model", "True AUC", "$\\overline{\\text{err}}$", "$\\text{CV}_{10}$")

predtable.AUC[, 1] <- mod.names
predtable.AUC[, 2] <- sapply(true.AUCs, function(x) round(x,4))
predtable.AUC[, 3] <- sapply(train.err.AUC, function(x) round(x, 4))
predtable.AUC[, 4] <- sapply(CV_10.AUC, function(x) round(x, 4))




##### Table: Prediction error estimates of each model (Misclassification error)

### Estimates are calculated and placed in table

true.misc      <- training_error(data = test.imp, data.RSF = test, w_inv = w_inv.test, delta = delta.test, models = models, type = "mis_error",
                                 t = 5, alpha = 0.5, covars.std = covars.std.test, covars.std.test = covars.std.test)

train.err.misc <- training_error(data = data.imp, data.RSF = data, w_inv = w_inv.data, delta = delta.data, models = models, type = "mis_error", 
                                 t = 5, alpha = 0.5, covars.std = covars.std, covars.std.test = covars.std)

CV_10.misc     <- CV_K(data = data.imp, data.RSF = data, nfold = 10, delta = delta.data, models = models, type = "mis_err", t = 5, alpha = 0.5)

Err_1.misc     <- Err_1(data = data.imp, data.RSF = data, w_inv = w_inv.data, delta = delta.data, models = models, type = "mis_err", t = 5, alpha = 0.5, B = 1000)

Err_632.misc   <- .368 * train.err.misc + .632 * Err_1.misc

gamma <- gamma(data = data.imp, data.RSF = data, w_inv = w_inv.data, delta = delta.data, models = models, type = "mis_error", t = 5, alpha = 0.5)
Rhat  <- (Err_1.misc - train.err.misc) / (gamma - train.err.misc)
omega <- .632 / (1 - .368 * Rhat)

Err_632p.misc  <- (1 - omega) * train.err.misc + omega * Err_1.misc


predtable.misc           <- matrix(nrow = length(mod.names), ncol = 7)
colnames(predtable.misc) <- c("Model", "True classification error", "$\\overline{\\text{err}}$", "$\\text{CV}_{10}$", "$\\widehat{\\text{Err}}^{(1)}$", "$\\widehat{\\text{Err}}^{(0.632)}$", "$\\widehat{\\text{Err}}^{(0.632+)}$")

predtable.misc[, 1] <- mod.names
predtable.misc[, 2] <- sapply(true.misc, function(x) round(x, 4))
predtable.misc[, 3] <- sapply(train.err.misc, function(x) round(x, 4))
predtable.misc[, 4] <- sapply(CV_10.misc, function(x) round(x, 4))
predtable.misc[, 5] <- sapply(Err_1.misc, function(x) round(x, 4))
predtable.misc[, 6] <- sapply(Err_632.misc, function(x) round(x, 4))
predtable.misc[, 7] <- sapply(Err_632p.misc, function(x) round(x, 4))




##### Figure for tables

cols  <- viridis(5)
cols2 <- viridis(10)

## AUC
error.est <- c("True error", "Training data estimate", "CV_10")

df <- data.frame("Model" = predtable.AUC[,1],
                 "True error" = as.numeric(predtable.AUC[,2]),
                 "Training data estimate" = as.numeric(predtable.AUC[,3]),
                 "CV_10" = as.numeric(predtable.AUC[,4])
)
df <- melt(df, id = "Model", variable.name = "Error")
df$Model                 <- factor(df$Model, levels = mod.names[21:1])
levels(df$Error)         <- error.est[1:3]
df$Error                 <- factor(df$Error, levels = error.est[3:1])

p1 <- ggplot(data = df, aes(x = Model, y = value, fill = Error)) + 
  geom_bar(stat="identity", position = "dodge", width = .7) + 
  scale_fill_manual(" ", values = c("True error" = "black", "Training data estimate" = cols2[3], "CV_10" = cols2[6]), limits = c("True error", "Training data estimate", "CV_10")) +
  xlab("") + ylab("AUC estimate") + 
  coord_flip(ylim = c(.6,.95)) +
  theme_bw()

## Misclassification error
error.est <- c("True error", "Training error", "CV_10", "Err^(1)", "Err^(.632)", "Err^(.632+)")

df <- data.frame("Model" = predtable.misc[,1],
                 "True error" = as.numeric(predtable.misc[,2]),
                 "Training error" = as.numeric(predtable.misc[,3]),
                 "CV_10" = as.numeric(predtable.misc[,4]),
                 "Err^(1)" = as.numeric(predtable.misc[,5]),
                 "Err^(.632)" = as.numeric(predtable.misc[,6]), 
                 "Err^(.632+)" = as.numeric(predtable.misc[,7])
)
df <- melt(df, id = "Model", variable.name = "Error")
df$Model                 <- factor(df$Model, levels = mod.names[21:1])
levels(df$Error)         <- error.est
df$Error                 <- factor(df$Error, levels = error.est[6:1])


p2 <- ggplot(data = df, aes(x = Model, y = value, fill = Error)) + 
  geom_bar(stat="identity", position = "dodge", width = .85) + 
  scale_fill_manual(" ", values = c("True error" = "black", "Training error" = cols2[3], "CV_10" = cols2[6], "Err^(1)" = cols2[8], "Err^(.632)" = cols2[10], "Err^(.632+)" = "darkorange"), limits = c("True error", "Training error", "CV_10", "Err^(1)", "Err^(.632)", "Err^(.632+)")) +
  xlab("") + ylab("Misclassification error estimate") + 
  coord_flip(ylim = c(.05,.2)) +
  theme_bw()



