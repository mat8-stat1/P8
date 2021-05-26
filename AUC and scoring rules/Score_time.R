
##################################################################
# Plot AUC, log and Brier score for t=1:15 and integrated scores #
##################################################################

library(pracma)
library(data.table)
library(reshape2)
library(viridis)
library(scales)

source("data_cleaning.R")
source("ROC_IPCW.R")


IPCW_loop_RSF <- function(t_values){
  dat <- data.frame(Model = rep("RSF", length(t_values)), t_value = t_values, AUC = rep(NA, length(t_values)),
                             log = rep(NA, length(t_values)), Brier = rep(NA, length(t_values)))
  k <- 1
  for (t in t_values){
    wh       <- length(pred_RSF$time.interest[pred_RSF$time.interest <= t]) # this time point corresponds to estimated 5-year survival
    pred.RSF <- pred_RSF$survival[, wh] # estimated 5-year surivival
    
    dat$Brier[k]  <- Brier_IPCW(pred.RSF, test.imp, t)
    dat$log[k] <- log_IPCW(pred.RSF, test.imp, t)
    dat$AUC[k]  <- ROC_IPCW(pred.RSF, test.imp, t)$AUC
    
    k <- k+1
  }
  dat
}

IPCW_loop_IPS <- function(t_values){
  dat <- data.frame(Model = rep("IPS", length(t_values)), t_value = t_values, AUC = rep(NA, length(t_values)),
                    log = rep(NA, length(t_values)), Brier = rep(NA, length(t_values)))
  k <- 1
  for (t in t_values){
    sum_IPS  <- summary(IPS, times = t)
    pred.IPS <- sum_IPS$surv[1]*(test.IPS==0) + sum_IPS$surv[2]*(test.IPS==1) + sum_IPS$surv[3]*(test.IPS==2) +
      sum_IPS$surv[4]*(test.IPS==3) + sum_IPS$surv[5]*(test.IPS==4) + sum_IPS$surv[6]*(test.IPS==5) 
    
    dat$Brier[k]  <- Brier_IPCW(pred.IPS, test.imp, t)
    dat$log[k] <- log_IPCW(pred.IPS, test.imp, t)
    dat$AUC[k]  <- ROC_IPCW(pred.IPS, test.imp, t)$AUC
    
    k <- k+1
  }
  dat
}

IPCW_loop_Cox <- function(t_values){
  dat <- data.frame(Model = rep(NA, length(t_values)*5), t_value = rep(t_values, rep(5, length(t_values))), AUC = rep(NA, length(t_values)*5),
                    log = rep(NA, length(t_values)*5), Brier = rep(NA, length(t_values)*5))
  l <- 1
  cox.names <- rep(c("coxph.aic", "coxph.aic.spl", "coxph.bic", "coxph.bic.spl", "coxph.lasso"), length(t_values))
  
  for (t in t_values){
    for (j in cox.models){
      pred.coxph <- rep(NA, nrow(test.imp))
      for (i in 1:nrow(test.imp)){
        pred          <- survfit(j, newdata = test.imp[i,])
        pred.coxph[i] <- summary(pred, t)$surv
      }
      
      dat$AUC[l] <- ROC_IPCW(pred.coxph, test.imp, t)$AUC
      dat$log[l] <- log_IPCW(pred.coxph, test.imp, t)
      dat$Brier[l] <- Brier_IPCW(pred.coxph, test.imp, t)
      dat$Model[l] <- cox.names[l]
      
      l <- l + 1
    }
    
    pred.coxph <- rep(NA, nrow(test.imp))
    for (i in 1:nrow(test.imp)){
      pred          <- survfit(coxph.lasso, s = "lambda.min", x = covars.std, y = Surv(data.imp$OS, data.imp$status), newx = covars.std.test[i, ])
      pred.coxph[i] <- summary(pred, t)$surv
    }
    
    dat$Brier[l]  <- Brier_IPCW(pred.coxph, test.imp, t)
    dat$log[l] <- log_IPCW(pred.coxph, test.imp, t)
    dat$AUC[l]  <- ROC_IPCW(pred.coxph, test.imp, t)$AUC
    dat$Model[l] <- cox.names[l]
    
    l <- l + 1
    
  }
  dat
}

IPCW_loop_Weibull <- function(t_values){
  dat <- data.frame(Model = rep(NA, length(t_values)*5), t_value = rep(t_values, rep(5, length(t_values))), AUC = rep(NA, length(t_values)*5),
                    log = rep(NA, length(t_values)*5), Brier = rep(NA, length(t_values)*5))
  l <- 1
  weibull.names <- rep(c("weibull.aic", "weibull.aic.spl", "weibull.bic", "weibull.bic.spl", "weibull.lasso"), length(t_values))
  
  for (t in t_values){
    for (j in weibull.models){
      pred.wei <- pred.weibull(j, test.imp, t)
      
      dat$AUC[l] <- ROC_IPCW(pred.wei, test.imp, t)$AUC
      dat$log[l] <- log_IPCW(pred.wei, test.imp, t)
      dat$Brier[l] <- Brier_IPCW(pred.wei, test.imp, t)
      dat$Model[l] <- weibull.names[l]
      
      l <- l + 1
    }
    
    pred.wei <- S_wei.lasso(covars.std.test, t)
    
    dat$Brier[l]  <- Brier_IPCW(pred.wei, test.imp, t)
    dat$log[l] <- log_IPCW(pred.wei, test.imp, t)
    dat$AUC[l]  <- ROC_IPCW(pred.wei, test.imp, t)$AUC
    dat$Model[l] <- weibull.names[l]
    
    l <- l + 1
  }
  dat
}

RSF_loop <- IPCW_loop_RSF(1:15)
IPS_loop <- IPCW_loop_IPS(1:15)
Cox_loop <- IPCW_loop_Cox(1:15)
Weibull_loop <- IPCW_loop_Weibull(1:15)

######## Plot AUC across t=1:15
loop_dat_AUC <- rbind(IPS_loop, RSF_loop, Cox_loop, Weibull_loop) %>% 
  dplyr::select(Model, t_value, AUC)

Orange <- scales::seq_gradient_pal("yellow", "firebrick2", "Lab")(seq(0,1,length.out=5))
Blue <- scales::seq_gradient_pal("lightskyblue1", "mediumblue", "Lab")(seq(0,1,length.out=5))
Mod_col <- c("darkorange4", "chartreuse3", Orange, Blue)

ggplot(loop_dat_AUC, aes(x = t_value, y = AUC, colour = forcats::fct_inorder(Model))) +
  geom_line(size = 1.1) +
  xlab("Time (year)") +
  scale_colour_manual(name = "Model", values = Mod_col) +
  guides(color = FALSE, size = FALSE) +
  theme_minimal()

######## Plot log score across t=1:15
loop_dat_log <- rbind(IPS_loop, RSF_loop, Cox_loop, Weibull_loop) %>% 
  dplyr::select(Model, t_value, log)

ggplot(loop_dat_log, aes(x = t_value, y = log, color = forcats::fct_inorder(Model))) +
  geom_line(size = 1.1) +
  xlab("Time (year)") +
  ylab("Logarithmic score") +
  scale_colour_manual(name = "Model", values = Mod_col) +
  guides(color = FALSE, size = FALSE) +
  theme_minimal()

######## Plot Brier score across t=1:15
loop_dat_Brier <- rbind(IPS_loop, RSF_loop, Cox_loop, Weibull_loop) %>% 
  dplyr::select(Model, t_value, Brier)

ggplot(loop_dat_Brier, aes(x = t_value, y = Brier, color = forcats::fct_inorder(Model))) +
  geom_line(size = 1.1) +
  xlab("Time (year)") +
  ylab("Brier score") +
  scale_colour_manual(name = "Model", values = Mod_col) +
  guides(color = FALSE, size = FALSE) +
  theme_minimal()

############### Area under curves
IPS_AUC_int <- trapz(IPS_loop$t_value, IPS_loop$AUC)
IPS_log_int <- trapz(IPS_loop$t_value, IPS_loop$log)
IPS_Brier_int <- trapz(IPS_loop$t_value, IPS_loop$Brier)

RSF_AUC_int <- trapz(RSF_loop$t_value, RSF_loop$AUC)
RSF_log_int <- trapz(RSF_loop$t_value, RSF_loop$log)
RSF_Brier_int <- trapz(RSF_loop$t_value, RSF_loop$Brier)

coxph.aic_loop <- Cox_loop %>% filter(Model == "coxph.aic")
coxph.aic.spl_loop <- Cox_loop %>% filter(Model == "coxph.aic.spl")
coxph.bic_loop <- Cox_loop %>% filter(Model == "coxph.bic")
coxph.bic.spl_loop <- Cox_loop %>% filter(Model == "coxph.bic.spl")
coxph.lasso_loop <- Cox_loop %>% filter(Model == "coxph.lasso")

weibull.aic_loop <- Weibull_loop %>% filter(Model == "weibull.aic")
weibull.aic.spl_loop <- Weibull_loop %>% filter(Model == "weibull.aic.spl")
weibull.bic_loop <- Weibull_loop %>% filter(Model == "weibull.bic")
weibull.bic.spl_loop <- Weibull_loop %>% filter(Model == "weibull.bic.spl")
weibull.lasso_loop <- Weibull_loop %>% filter(Model == "weibull.lasso")

coxph.aic_AUC_int <- trapz(coxph.aic_loop$t_value, coxph.aic_loop$AUC)
coxph.aic_log_int <- trapz(coxph.aic_loop$t_value, coxph.aic_loop$log)
coxph.aic_Brier_int <- trapz(coxph.aic_loop$t_value, coxph.aic_loop$Brier)

coxph.aic.spl_AUC_int <- trapz(coxph.aic.spl_loop$t_value, coxph.aic.spl_loop$AUC)
coxph.aic.spl_log_int <- trapz(coxph.aic.spl_loop$t_value, coxph.aic.spl_loop$log)
coxph.aic.spl_Brier_int <- trapz(coxph.aic.spl_loop$t_value, coxph.aic.spl_loop$Brier)

coxph.bic_AUC_int <- trapz(coxph.bic_loop$t_value, coxph.bic_loop$AUC)
coxph.bic_log_int <- trapz(coxph.bic_loop$t_value, coxph.bic_loop$log)
coxph.bic_Brier_int <- trapz(coxph.bic_loop$t_value, coxph.bic_loop$Brier)

coxph.bic.spl_AUC_int <- trapz(coxph.bic.spl_loop$t_value, coxph.bic.spl_loop$AUC)
coxph.bic.spl_log_int <- trapz(coxph.bic.spl_loop$t_value, coxph.bic.spl_loop$log)
coxph.bic.spl_Brier_int <- trapz(coxph.bic.spl_loop$t_value, coxph.bic.spl_loop$Brier)

coxph.lasso_AUC_int <- trapz(coxph.lasso_loop$t_value, coxph.lasso_loop$AUC)
coxph.lasso_log_int <- trapz(coxph.lasso_loop$t_value, coxph.lasso_loop$log)
coxph.lasso_Brier_int <- trapz(coxph.lasso_loop$t_value, coxph.lasso_loop$Brier)

weibull.aic_AUC_int <- trapz(weibull.aic_loop$t_value, weibull.aic_loop$AUC)
weibull.aic_log_int <- trapz(weibull.aic_loop$t_value, weibull.aic_loop$log)
weibull.aic_Brier_int <- trapz(weibull.aic_loop$t_value, weibull.aic_loop$Brier)

weibull.aic.spl_AUC_int <- trapz(weibull.aic.spl_loop$t_value, weibull.aic.spl_loop$AUC)
weibull.aic.spl_log_int <- trapz(weibull.aic.spl_loop$t_value, weibull.aic.spl_loop$log)
weibull.aic.spl_Brier_int <- trapz(weibull.aic.spl_loop$t_value, weibull.aic.spl_loop$Brier)

weibull.bic_AUC_int <- trapz(weibull.bic_loop$t_value, weibull.bic_loop$AUC)
weibull.bic_log_int <- trapz(weibull.bic_loop$t_value, weibull.bic_loop$log)
weibull.bic_Brier_int <- trapz(weibull.bic_loop$t_value, weibull.bic_loop$Brier)

weibull.bic.spl_AUC_int <- trapz(weibull.bic.spl_loop$t_value, weibull.bic.spl_loop$AUC)
weibull.bic.spl_log_int <- trapz(weibull.bic.spl_loop$t_value, weibull.bic.spl_loop$log)
weibull.bic.spl_Brier_int <- trapz(weibull.bic.spl_loop$t_value, weibull.bic.spl_loop$Brier)

weibull.lasso_AUC_int <- trapz(weibull.lasso_loop$t_value, weibull.lasso_loop$AUC)
weibull.lasso_log_int <- trapz(weibull.lasso_loop$t_value, weibull.lasso_loop$log)
weibull.lasso_Brier_int <- trapz(weibull.lasso_loop$t_value, weibull.lasso_loop$Brier)

#-------------# Table
##### Table: AUC + scores for all models

mod.names <- c("IPS model",
               "coxph.aic", "coxph.aic.spl", "coxph.bic", "coxph.bic.spl", "coxph.lasso",
               "weibull.aic", "weibull.aic.spl", "weibull.bic", "weibull.bic.spl", "weibull.lasso",
               "RSF model")

AUC.table           <- matrix(nrow = length(mod.names), ncol = 4)
colnames(AUC.table) <- c("Model", "Integrated AUC", "Integrated logarithmic score", "Integrated Brier score")
AUCs                <- c(IPS_AUC_int, coxph.aic_AUC_int, coxph.aic.spl_AUC_int, coxph.bic_AUC_int, coxph.bic.spl_AUC_int,
                         coxph.lasso_AUC_int, weibull.aic_AUC_int, weibull.aic.spl_AUC_int, weibull.bic_AUC_int,
                         weibull.bic.spl_AUC_int, weibull.lasso_AUC_int, RSF_AUC_int)

logs                <- c(IPS_log_int, coxph.aic_log_int, coxph.aic.spl_log_int, coxph.bic_log_int, coxph.bic.spl_log_int,
                         coxph.lasso_log_int, weibull.aic_log_int, weibull.aic.spl_log_int, weibull.bic_log_int,
                         weibull.bic.spl_log_int, weibull.lasso_log_int, RSF_log_int)

Briers              <- c(IPS_Brier_int, coxph.aic_Brier_int, coxph.aic.spl_Brier_int, coxph.bic_Brier_int, coxph.bic.spl_Brier_int,
                         coxph.lasso_Brier_int, weibull.aic_Brier_int, weibull.aic.spl_Brier_int, weibull.bic_Brier_int,
                         weibull.bic.spl_Brier_int, weibull.lasso_Brier_int, RSF_Brier_int)

AUC.table[, 1] <- mod.names
AUC.table[, 2] <- sapply(AUCs, function(x) round(x, 4))
AUC.table[, 3] <- sapply(logs, function(x) round(x, 4))
AUC.table[, 4] <- sapply(Briers, function(x) round(x, 4))

AUC.table

########## plot
AUC.table.plot <- matrix(nrow = length(mod.names), ncol = 4)
colnames(AUC.table.plot) <- c("Model", "Integrated AUC", "Integrated logarithmic score + 15", "Integrated Brier score")

AUC.table.plot <- AUC.table.plot %>%
  as.data.frame %>% 
  mutate("Model" = mod.names) %>% 
  mutate("Integrated AUC" = AUCs) %>% 
  mutate("Integrated logarithmic score + 15" = logs + 15) %>% 
  mutate("Integrated Brier score" = Briers)

AUC.table2 <- melt(AUC.table.plot, id = "Model", variable.name = "Score")
AUC.table2$Model <- factor(AUC.table2$Model, levels = mod.names[c(12, 6:2, 11:7, 1)])
AUC.table2$Score <- factor(AUC.table2$Score, levels = c("Integrated AUC",
                                                        "Integrated logarithmic score + 15", "Integrated Brier score"))

cols <- viridis(8)

ggplot(data=AUC.table2, aes(x = Model, y = value, fill = Score)) + 
  geom_bar(stat="identity", position = "dodge", width = .8) + 
  #scale_fill_manual(" ", values = c("True error" = "black", "Training error" = "forestgreen", "CV_10" = "darkorange", "Err^(1)" = "dodgerblue3", "Err^(.632)" = "deepskyblue", "Err^(.632+)" = "lightblue"), limits = c("True error", "Training error", "CV_10", "Err^(1)", "Err^(.632)", "Err^(.632+)")) +
  scale_fill_manual(" ", values = c("Integrated AUC" = cols[3], "Integrated logarithmic score + 15" = cols[5],
                                    "Integrated Brier score" = cols[7]), limits = c("Integrated Brier score", "Integrated logarithmic score + 15", "Integrated AUC")) +
  xlab("") + ylab("Integrated score") +
  coord_flip() +
  scale_y_continuous(limits=c(8.7,12),oob = rescale_none) +
  theme_minimal()
