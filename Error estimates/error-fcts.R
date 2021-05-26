
##############################
# Error estimation functions #
##############################

source("data_cleaning.R")
source("fcts-fit-pred-loss.R")

library(survival)
library(splines)
library(glmnet)


### Training error function

'%notin%' <- Negate('%in%')

training_error <- function(data, data.RSF, w_inv, delta, models, type = "AUC", t = 5, alpha = 0.5, covars.std = covars.std, covars.std.test = covars.std.test){
  out <- rep(NA, length(models))
  
  for (j in 1:length(models)){
    
    # Estimation of t-year survival of patients in data:
    if("survreg" %in% class(models[[j]])){
      pred <- fit_pred.weibull(data.test=data, model=models[[j]], t=t)
    }
    if("coxph" %in% class(models[[j]])){
      pred <- fit_pred.coxph(data.test=data, model=models[[j]])
    }
    if("iregnet" %in% class(models[[j]])){
      pred <- fit_pred.weibull.lasso(data.fit=data, covars.std.fit = covars.std, covars.std.test = covars.std.test, fit=F,  model=models[[j]], t)
      
    }
    if("coxnet" %in% class(models[[j]])){
      pred <- fit_pred.coxph.lasso(data, covars.std.fit = covars.std, covars.std.test = covars.std.test, fit=F,  model=models[[j]])
    }
    
    
    if("rfsrc" %in% class(models[[j]])){
      pred <- fit_pred.RSF(data.test=data.RSF,fit=F, model=models[[j]])
    }
    
    if("survfit" %in% class(models[[j]])){
      pred <- fit_pred.IPS(data.test=data, model=models[[j]])
    }
    
    # Calculating training error based on type:
    if(type=="AUC"  & "rfsrc" %notin% class(models[[j]])){
      out[j] <- AUC(pred, w_inv, delta, data, t)
    }
    if(type=="mis_error"  & "rfsrc" %notin% class(models[[j]])){
      out[j] <- mis_error(pred, data, w_inv, delta, t, alpha)
    }
    if(type=="AUC" & "rfsrc" %in% class(models[[j]])){
      out[j] <- AUC(pred, w_inv, delta, data, t)
    }
    
    if(type=="mis_error" & "rfsrc" %in% class(models[[j]])){
      out[j] <- mis_error(pred, data, w_inv, delta, t, alpha)
    }
    
  }
  
  out
}



### CV_K function

# If only RSF model input, also give data argument (imputed data)
CV_K <- function(data = NULL, data.RSF = NULL, nfold, delta, models = list(), type = "AUC", t = 5, alpha = 0.5){
  
  sample.data <- sample(rep(1:nfold, l = nrow(data))) # sampling nfold training- and validation sets
  out <- rep(NA, length(models))
  
  for (j in 1:length(models)){
    
    loss <- rep(NA, nfold)
    
    for (i in 1:nfold){
      
      test.index <- which(sample.data == i, arr.ind = TRUE)
      data.test  <- data[test.index, ]
      data.fit   <- data[-test.index, ]
      pred       <- rep(NA, length(test.index))
      
      # Estimating t-year survival based on model class (models are fitted to data.fit):
      if("survreg" %in% class(models[[j]])){
        form <- paste(attr(models[[j]]$terms, "term.labels"), collapse = " + ")
        pred <- fit_pred.weibull(data.fit = data.fit, data.test = data.test, formula = form, t = t)
      }
      if("coxph" %in% class(models[[j]])){
        form <- paste(attr(models[[j]]$terms, "term.labels"), collapse = " + ")
        pred <- fit_pred.coxph(data.fit = data.fit, data.test = data.test, formula = form)
      }
      if("iregnet" %in% class(models[[j]])){
        ## standardising training data; using same transformations on test data.
        # covars.std.fit: standardising covariates used in the j'th model
        covars         <- rownames(models[[j]]$beta)[-1]
        covars.std.fit <- data.matrix(data.fit[, covars])
        covars.std.fit[, c("sex", "bsym", "LDH")] <- covars.std.fit[, c("sex", "bsym", "LDH")] - 1
        
        centers <- attr(scale(covars.std.fit[, !(colnames(covars.std.fit) %in% c("sex", "bsym", "LDH"))]), "scaled:center")
        scales  <- attr(scale(covars.std.fit[, !(colnames(covars.std.fit) %in% c("sex", "bsym", "LDH"))]), "scaled:scale")
        
        covars.std.fit[, !(colnames(covars.std.fit) %in% c("sex", "bsym", "LDH"))] <- scale(covars.std.fit[, !(colnames(covars.std.fit) %in% c("sex", "bsym", "LDH"))])
        
        # using same transformations on test data
        covars.std.test <- data.matrix(data.test[ ,covars])
        covars.std.test[, c("sex", "bsym", "LDH")] <- covars.std.test[, c("sex", "bsym", "LDH")] - 1
        for (k in names(scales)){
          covars.std.test[, k] <- sapply(covars.std.test[, k], function(x) (x - centers[k]) / scales[k])
        }
        
        pred <- fit_pred.weibull.lasso(data.fit = data.fit, covars.std.fit = covars.std.fit, covars.std.test = covars.std.test, fit = T, model = models[[j]], t = t)
      }
      if("coxnet" %in% class(models[[j]])){
        ## standardising training data; using same transformations on test data.
        # covars.std.fit: standardising covariates used in the j'th model
        covars         <- rownames(models[[j]]$beta)[-1]
        covars.std.fit <- data.matrix(data.fit[, covars])
        covars.std.fit[, c("sex", "bsym", "LDH")] <- covars.std.fit[, c("sex", "bsym", "LDH")] - 1
        
        centers <- attr(scale(covars.std.fit[, !(colnames(covars.std.fit) %in% c("sex", "bsym", "LDH"))]), "scaled:center")
        scales  <- attr(scale(covars.std.fit[, !(colnames(covars.std.fit) %in% c("sex", "bsym", "LDH"))]), "scaled:scale")
        
        covars.std.fit[, !(colnames(covars.std.fit) %in% c("sex", "bsym", "LDH"))] <- scale(covars.std.fit[, !(colnames(covars.std.fit) %in% c("sex", "bsym", "LDH"))])
        
        # using same transformations on test data
        covars.std.test <- data.matrix(data.test[ ,covars])
        covars.std.test[, c("sex", "bsym", "LDH")] <- covars.std.test[, c("sex", "bsym", "LDH")] - 1
        for (k in names(scales)){
          covars.std.test[, k] <- sapply(covars.std.test[, k], function(x) (x - centers[k]) / scales[k])
        }
        
        pred <- fit_pred.coxph.lasso(data.fit = data.fit, covars.std.fit = covars.std.fit, covars.std.test = covars.std.test, fit = T,  model = models[[j]], t = t)
      }
      if("rfsrc" %in% class(models[[j]])){
        data.test.RSF  <- data.RSF[test.index, ] # data needs to be RSF data in this case
        data.fit.RSF   <- data.RSF[-test.index, ]
        pred           <- fit_pred.RSF(data.fit = data.fit.RSF, data.test = data.test.RSF, fit = T)
      }
      if("survfit" %in% class(models[[j]])){
        pred <- fit_pred.IPS(data.fit = data.fit, data.test = data.test, formula = "IPS2")
      }
      
      
      # Calculating weights for censoring distribution on test set:
      cens.mod <- survfit(Surv(OS, cens) ~ 1, data = data.test)
      
      S_t   <- summary(cens.mod, t)$surv
      S_y   <- sapply(data.test[,"OS"], function(x) summary(cens.mod, x)$surv)
      S_y[data.test[, "status"] == 0] <- 1 # to avoid dividing by 0. these S_y values are not used.
      
      w_inv   <- (data.test[, "OS"] > t) * 1/S_t + (data.test[, "OS"] <= t & data.test[, "status"] == 1) * 1/S_y
      delta_i <- delta[test.index]
      
      
      # Calculating loss fct. for the i'th test data set with corresponding survival estimates
      if (type == "AUC"){
        loss[i] <- AUC(pred, w_inv, delta_i, data.test, t)
      }
      if (type == "mis_err"){
        # nrow(data.test) is multiplied to counteract for the mis_error function; mean is taken over N before outputting result
        loss[i] <- nrow(data.test) * mis_error(pred, test.data = data.test, w_inv, delta_i, t, alpha)
      }
    
    } # i loop closed
    
    # Calculating the average loss fct.
    if (type == "AUC"){
      out[j] <- mean(loss)
    }
    if (type == "mis_err"){
      out[j] <- 1/nrow(data) * sum(loss)
    }
    
  } # j loop closed
  
  out
}



### Err^(1) function

Err_1 <- function(data = NULL, data.RSF = NULL, w_inv, delta, models = list(), type = "mis_err", t = 5, alpha = 0.5, B = 50){
  
  # BS[, b] contains the i's chosen for the b'th bootstrap sample
  BS <- matrix(sample(1:nrow(data), nrow(data)*B, replace = T), ncol = B)
  
  # test.index[, b] contains the rows not contained in the b'th bootstrap sample (T/F)
  test.index <- apply(matrix(1:B), 1, function(b) sapply(1:nrow(data), function(i) !(any(i %in% BS[, b]))))
  
  # C[i] contains the number of bootstrap samples where i'th observation is not included
  C <- apply(test.index, 1, sum)
  
  out <- rep(NA, length(models))
  
  for (j in 1:length(models)){
    
    # Estimation of t-year survival of patients in data.
    # pred[i,b] contains estimated survival of i'th subject based on b'th bootstrap model
    if("survreg" %in% class(models[[j]])){
      form <- paste(attr(models[[j]]$terms, "term.labels"), collapse = " + ")
      pred <- apply(matrix(1:B), 1, function(b) fit_pred.weibull(data.fit = data[BS[,b], ], data.test = data, formula = form, t = t))
    }
    if("coxph" %in% class(models[[j]])){
      form <- paste(attr(models[[j]]$terms, "term.labels"), collapse = " + ")
      pred <- apply(matrix(1:B), 1, function(b) fit_pred.coxph(data.fit = data[BS[,b], ], data.test = data, formula = form))
    }
    if("iregnet" %in% class(models[[j]])){
      ## standardising training data; using same transformations on test data.
      # covars.std.fit[[b]]: standardising covariates used in the j'th model for the b'th bootstrap sample
      covars         <- rownames(models[[j]]$beta)[-1]
      
      covars.std.fit  <- list()
      covars.std.test <- list()
      for (b in 1:B){
        covars.std.fit[[b]] <- data.matrix(data[BS[, b], covars])
        covars.std.fit[[b]][, c("sex", "bsym", "LDH")] <- covars.std.fit[[b]][, c("sex", "bsym", "LDH")] - 1
        
        centers <- attr(scale(covars.std.fit[[b]][, !(colnames(covars.std.fit[[b]]) %in% c("sex", "bsym", "LDH"))]), "scaled:center")
        scales  <- attr(scale(covars.std.fit[[b]][, !(colnames(covars.std.fit[[b]]) %in% c("sex", "bsym", "LDH"))]), "scaled:scale")
        
        covars.std.fit[[b]][, !(colnames(covars.std.fit[[b]]) %in% c("sex", "bsym", "LDH"))] <- scale(covars.std.fit[[b]][, !(colnames(covars.std.fit[[b]]) %in% c("sex", "bsym", "LDH"))])
        
        # using same transformations on test data
        # covars.std.test[[b]]: standardising covariates used in the j'th model for the b'th bootstrap sample
        covars.std.test[[b]] <- data.matrix(data[, covars])
        covars.std.test[[b]][, c("sex", "bsym", "LDH")] <- covars.std.test[[b]][, c("sex", "bsym", "LDH")] - 1
        for (i in names(scales)){
          covars.std.test[[b]][, i] <- sapply(covars.std.test[[b]][, i], function(x) (x - centers[i]) / scales[i])
        }
      }
      
      pred <- apply(matrix(1:B), 1, function(b) fit_pred.weibull.lasso(data.fit = data[BS[,b], ], covars.std.fit = covars.std.fit[[b]], covars.std.test = covars.std.test[[b]], fit = T, model = models[[j]], t = t))
      pred <- matrix(unlist(pred), ncol = B, byrow = F)
    }
    if("coxnet" %in% class(models[[j]])){
      ## standardising training data; using same transformations on test data.
      # covars.std.fit[[b]]: standardising covariates used in the j'th model for the b'th bootstrap sample
      covars         <- rownames(models[[j]]$beta)[-1]
      
      covars.std.fit  <- list()
      covars.std.test <- list()
      for (b in 1:B){
        covars.std.fit[[b]] <- data.matrix(data[BS[, b], covars])
        covars.std.fit[[b]][, c("sex", "bsym", "LDH")] <- covars.std.fit[[b]][, c("sex", "bsym", "LDH")] - 1
        
        centers <- attr(scale(covars.std.fit[[b]][, !(colnames(covars.std.fit[[b]]) %in% c("sex", "bsym", "LDH"))]), "scaled:center")
        scales  <- attr(scale(covars.std.fit[[b]][, !(colnames(covars.std.fit[[b]]) %in% c("sex", "bsym", "LDH"))]), "scaled:scale")
        
        covars.std.fit[[b]][, !(colnames(covars.std.fit[[b]]) %in% c("sex", "bsym", "LDH"))] <- scale(covars.std.fit[[b]][, !(colnames(covars.std.fit[[b]]) %in% c("sex", "bsym", "LDH"))])
        
        # using same transformations on test data
        # covars.std.test[[b]]: standardising covariates used in the j'th model for the b'th bootstrap sample
        covars.std.test[[b]] <- data.matrix(data[, covars])
        covars.std.test[[b]][, c("sex", "bsym", "LDH")] <- covars.std.test[[b]][, c("sex", "bsym", "LDH")] - 1
        for (i in names(scales)){
          covars.std.test[[b]][, i] <- sapply(covars.std.test[[b]][, i], function(x) (x - centers[i]) / scales[i])
        }
      }
      
      pred <- apply(matrix(1:B), 1, function(b) fit_pred.coxph.lasso(data.fit = data[BS[,b], ], covars.std.fit = covars.std.fit[[b]], covars.std.test = covars.std.test[[b]], fit = T, model = models[[j]], t = t))
      pred <- matrix(unlist(pred), ncol = B, byrow = F)
    }
    if("rfsrc" %in% class(models[[j]])){
      pred <- apply(matrix(1:B), 1, function(b) fit_pred.RSF(data.fit = data.RSF[BS[,b], ], data.test = data.RSF, fit = T))
    }
    if("survfit" %in% class(models[[j]])){
      pred <- apply(matrix(1:B), 1, function(b) fit_pred.IPS(data.fit = data[BS[,b], ], data.test = data, formula = "IPS2"))
    }
    
    # converting predictions to T/F (T: death predicted)
    pred <- apply(pred, 2, function(i) i < alpha)
    
    # Calculating IPCW misclassification error
    true_val  <- data[, "OS"] < t & data[, "status"] == 1
    err_terms <- sapply(1:nrow(data), function(i) (1 / C[i]) * w_inv[i] * delta[i] * sum(true_val[i] != pred[i, test.index[i, ]]))
    out[j]    <- mean(err_terms, na.rm = T) # na.rm = T ensures that cases of C[i] = 0 are not included
  }
  
  out
}


### gamma function
# NB: Type argument not used (only implemented for misclassification error)
gamma <- function(data = NULL, data.RSF = NULL, w_inv = w_inv.data, delta = delta.data, models = list(), type = "mis_error", t = 5, alpha = 0.5){
  
  n   <- nrow(data)
  out <- rep(NA, length(models))
  
  for (j in 1:length(models)){
    
    # Estimation of t-year survival of patients in data:
    if("survreg" %in% class(models[[j]])){
      pred <- fit_pred.weibull(data.fit = data, data.test = data, model = models[[j]], t = t)
    }
    if("coxph" %in% class(models[[j]])){
      pred <- fit_pred.coxph(data.fit = data, data.test = data, model = models[[j]])
    }
    if("iregnet" %in% class(models[[j]])){
      pred <- fit_pred.weibull.lasso(data.fit = data, covars.std.fit = covars.std, covars.std.test = covars.std, fit = F,  model = models[[j]], t = t)
    }
    if("coxnet" %in% class(models[[j]])){
      pred <- fit_pred.coxph.lasso(data.fit = data, covars.std.fit = covars.std, covars.std.test = covars.std, fit = F,  model = models[[j]], t = t)
    }
    if("rfsrc" %in% class(models[[j]])){
      pred <- fit_pred.RSF(data.fit = data.RSF, data.test = data.RSF, fit = F, model = models[[j]])
    }
    if("survfit" %in% class(models[[j]])){
      pred <- fit_pred.IPS(data.fit = data, data.test = data, model = models[[j]])
    }
    
    # Calculating gamma by shifting prediction vector recursively and summing:
    out[j] <- sum(sapply(2:length(pred), function(x) mis_error(pred[c(x:length(pred), 1:(x-1))], data, w_inv, delta, t, alpha)))
    out[j] <- out[j] + mis_error(pred, data, w_inv, delta, t, alpha)
    out[j] <- 1/n * out[j] # 1/n^2 from definition, but 1/n* already applied in miss_error
  }
  
  out
}




