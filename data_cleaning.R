##############################
# Data cleaning + imputation #
##############################

library(tidyverse)
library(missForest)

## Load non-cleaned training and test data as "data" and "test"

# data <- ...
# test <- ...


# Scaling OS such that it is in years instead of days
data$OS <- data$OS / 365.24
test$OS <- test$OS / 365.24


# Coding factors

data[,22:ncol(data)] <- -1*(data[,22:ncol(data)] - 2)
for (i in 22:ncol(data)){
  data[,i] <- as.factor(data[,i])
}

data$ENODAL <- as.factor(data$ENODAL)
data$Stadium..Ann.Arbor. <- as.factor(data$Stadium..Ann.Arbor.)
data$ECOG <- as.factor(data$ECOG)

test[,22:ncol(test)] <- -1*(test[,22:ncol(test)] - 2)
for (i in 22:ncol(test)){
  test[,i] <- as.factor(test[,i])
}

test$ENODAL <- as.factor(test$ENODAL)
test$Stadium..Ann.Arbor. <- as.factor(test$Stadium..Ann.Arbor.)
test$ECOG <- as.factor(test$ECOG)


## Imputing all covariates using tree based algorithm

set.seed(2)
data.imp <- missForest(data[, c(1,2,5:ncol(data))], maxiter = 20, ntree = 1000)
test.imp <- missForest(test[, c(1,2,5:ncol(test))], maxiter = 20, ntree = 1000)


## Joining and separating training and test data, ensuring same categories in all variables
n_train <- nrow(data)
n_test <- nrow(test)

full_data <- rbind(data, test)

data <- full_data[(1:n_train),]
test <- full_data[((n_train+1):(n_train+n_test)),]

data.imp <- cbind(data$OS, data$status, data.imp$ximp)
test.imp <- cbind(test$OS, test$status, test.imp$ximp)

colnames(data.imp)[1:2] <- c("OS", "status")
colnames(test.imp)[1:2] <- c("OS", "status")



# Making IPS variable
# https://www.msdmanuals.com/professional/multimedia/clinical-calculator/international-prognostic-score-in-hodgkin-lymphoma
data.imp[, "IPS"] <- (data.imp[, "Albumin"] < 40) + (data.imp[, "Hemoglobin"] < 105) + (data.imp[, "sex"]=="male") +
  (as.numeric(data.imp[, "Stadium..Ann.Arbor."]) > 3) + (data.imp[, "Alder"] >= 45) + (data.imp[, "Leukocytter.mia.L."] >= 15) + 
  (data.imp[, "Lymfocytter.mia.L."] < 0.6 | (data.imp[, "Lymfocytter.mia.L."] / data.imp[, "Leukocytter.mia.L."] < 0.08))

# We regard these as cont. covariates (even though treated as cat. when imputing):
data.imp[, "organ.inv"]           <- rowSums(sapply(data.imp[,22:(ncol(data.imp)-1)], function(x) as.numeric(as.character(x))))
data.imp[, "ECOG"]                <- as.numeric(data.imp[, "ECOG"]) - 1
data.imp[, "Stadium..Ann.Arbor."] <- as.numeric(data.imp[, "Stadium..Ann.Arbor."]) - 1
data.imp[, "ENODAL"]              <- as.numeric(data.imp[, "ENODAL"]) - 1

# ---------
# The same is done for test.imp
# ---------

# Making IPS variable

# https://www.msdmanuals.com/professional/multimedia/clinical-calculator/international-prognostic-score-in-hodgkin-lymphoma
test.imp[, "IPS"] <- (test.imp[, "Albumin"] < 40) + (test.imp[, "Hemoglobin"] < 105) + (test.imp[, "sex"]=="male") +
  (as.numeric(test.imp[, "Stadium..Ann.Arbor."]) > 3) + (test.imp[, "Alder"] >= 45) + (test.imp[, "Leukocytter.mia.L."] >= 15) + 
  (test.imp[, "Lymfocytter.mia.L."] < 0.6 | (test.imp[, "Lymfocytter.mia.L."] / test.imp[, "Leukocytter.mia.L."] < 0.08))

# We regard these as cont. covariates (even though treated as cat. when imputing):
test.imp[, "organ.inv"]           <- rowSums(sapply(test.imp[,22:(ncol(test.imp)-1)], function(x) as.numeric(as.character(x))))
test.imp[, "ECOG"]                <- as.numeric(test.imp[, "ECOG"]) - 1
test.imp[, "Stadium..Ann.Arbor."] <- as.numeric(test.imp[, "Stadium..Ann.Arbor."]) - 1
test.imp[, "ENODAL"]              <- as.numeric(test.imp[, "ENODAL"]) - 1

# ---------


## Censoring covariate (only for imputed data)
data.imp$cens <- 1 - data.imp$status
test.imp$cens <- 1 - test.imp$status

## Covariate for collapsing IPS levels
data.imp$IPS2 <- case_when(data.imp$IPS %in% 5:7 ~ 5L, 
                           TRUE ~ data.imp$IPS)
test.imp$IPS2 <- case_when(test.imp$IPS %in% 5:7 ~ 5L, 
                           TRUE ~ test.imp$IPS)

