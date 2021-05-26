
##############################################################
# Fitting RSF, getting 1 tree from the forest and prediction #
##############################################################

library(tidyverse)
library(xtable)
library(randomForestSRC)
library(survival)
library(viridis)
library(DescTools)

source("data_cleaning.R")
source("fcts-fit-pred-loss.R")
source("ROC_IPCW.R")

############### Checking NA values for covariates on training data
NA_var <- data.frame(names(data), colSums(is.na(data)))
colnames(NA_var) <- c("Variable name", "Number of NA")

################# Checking NA values for covariates on test data
NA_var_test <- data.frame(names(test), colSums(is.na(test)))
colnames(NA_var_test) <- c("Variable name", "Number of NA")

################# Hardcoding NA to be a factor
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

############### Random survival forest
set.seed(123)
RSF <- rfsrc(Surv(OS, status) ~ . - OS - status, data = data, na.action = "na.impute", importance = "TRUE", block.size = 1000)

############### Impotance measure (VIMP)
VIMP <- RSF$importance %>% as.numeric()

VIMP_dat <- data.frame(Variable = colnames(data)[-(3:4)], VIMP = VIMP, stringsAsFactors = FALSE) %>%
  arrange(VIMP) %>% 
  mutate(Variable = fct_reorder(Variable, VIMP))

ggplot(VIMP_dat, aes(VIMP, Variable)) +
  geom_point(aes(colour = "dodgerblue2", size = VIMP), data = VIMP_dat) +
  theme_minimal() +
  theme(legend.position = "none", axis.text.y = element_text(size = 9, angle = 20)) +
  ylab("") +
  xlim(c(-0.001, 0.065)) +
  scale_y_discrete(limits = c(levels(VIMP_dat$Variable), " "))

##################### Error for different size forests using AUC
pred_ntree_plot_cv_AUC <- function(start, slut, data, nfold){
  list <- list(err_rate_dat = data.frame(ntree = seq(start, slut, start), err_rate = rep(NA, (slut/start))), plot = NA)
  data$cens <- 1 - data$status
  
  n <- nrow(data)
  sample.data <- sample(rep(1:nfold, l = n))
  
  temp <- rep(NA, nfold)
  
  for (j in 1:(slut/start)) {
    for(i in 1:nfold){
      testIndex <- which(sample.data==i,arr.ind=TRUE)
      testData  <- data[testIndex, ]
      trainData <- data[-testIndex, ]
      
      mod_j <- rfsrc(Surv(OS, status) ~ . - OS - status, data = trainData, ntree = start*j,
                     na.action = "na.impute")
      pred_mod_j <- predict(mod_j, newdata = testData, na.action = "na.impute")
      wh       <- length(pred_mod_j$time.interest[pred_mod_j$time.interest <= 5]) # this time point corresponds to estimated 5-year survival
      pred.RSF <- pred_mod_j$survival[, wh]
      
      temp[i] <- ROC_IPCW(pred.RSF, testData, 5)$AUC
    }
    list$err_rate_dat$err_rate[j] <- mean(temp)
  }
  list$plot <- ggplot(list$err_rate_dat, aes(ntree, err_rate)) +
    geom_point(shape = 18, color = "dodgerblue3", size = 5) +
    xlab("Number of trees in the forest") +
    ylab("AUC") +
    theme_minimal()
  return(list)
}

############ Using the function to plot AUC for forests with 25, 50,..., 1000 trees
set.seed(420)
pred_plot2000420 <- pred_ntree_plot_cv_AUC(25, 1000, data, 5)

pred_plot2000420$plot

############### Survival curves of subjects
evt_time <- RSF$time.interest
surv_oob <- RSF$survival.oob %>%
  as.data.frame() %>%
  dplyr::select(V286) %>% 
  cbind(RSF$survival.oob) %>%
  as.data.frame() %>%
  arrange(desc(V286)) %>%
  dplyr::select(!V286)

surv_ind <- data.frame(x = 1:nrow(surv_oob), z = surv_oob)

surv_pivot <- surv_ind %>% 
  pivot_longer(cols = colnames(surv_ind)[-1]) %>% 
  mutate(evt_time = rep(evt_time, nrow(surv_oob)))

mean <- data.frame(x = rep(length(evt_time)+1, length(evt_time)), evt_time = evt_time, value = colMeans(RSF$survival.oob))

ggplot(surv_pivot, aes(evt_time, value, group = x, color = x)) + 
  geom_line() +
  scale_colour_gradientn(name = "x", colours = viridis(8)[-(1:2)]) +
  xlab("Time since diagnosis (years)") +
  ylab("Survival probability") +
  geom_line(data = mean, color = "black", size = 1.5) +
  theme_minimal() +
  theme(legend.position = "none")

##################### Get 1 tree from the forest and plot it

set.seed(2)
RSF_imp <- rfsrc(Surv(OS, status) ~ . - OS - status, data = data.imp[,-(49:52)])

tree1_RSF_imp <- get.tree(RSF_imp, tree.id = 1, time = 5, surv.type = "surv")

plot(tree1_RSF_imp)

############### Predictions on test
set.seed(123)
pred_RSF <- predict(RSF, newdata = test, na.action = "na.impute", importance = "TRUE", block.size = 1000)

## Predicted survival curves
pred_evt_time <- pred_RSF$time.interest
pred_surv_oob <- pred_RSF$survival %>%
  as.data.frame() %>%
  dplyr::select(V286) %>% 
  cbind(pred_RSF$survival) %>%
  as.data.frame() %>%
  arrange(desc(V286)) %>%
  dplyr::select(!V286)

pred_surv_ind <- data.frame(x = 1:nrow(pred_surv_oob), z = pred_surv_oob)

pred_surv_pivot <- pred_surv_ind %>% 
  pivot_longer(cols = colnames(pred_surv_ind)[-1]) %>% 
  mutate(pred_evt_time = rep(pred_evt_time, nrow(pred_surv_oob)))

pred_mean <- data.frame(x = rep(length(pred_evt_time)+1, length(pred_evt_time)), pred_evt_time = pred_evt_time,
                   value = colMeans(pred_RSF$survival))

ggplot(pred_surv_pivot, aes(pred_evt_time, value, group = x, color = x)) + 
  geom_line() +
  scale_colour_gradientn(name = "x", colours = viridis(8)[-(1:2)]) +
  xlab("Time since diagnosis (years)") +
  ylab("Survival probability") +
  geom_line(data = pred_mean, color = "black", size = 1.5) +
  theme_minimal() +
  theme(legend.position = "none")
