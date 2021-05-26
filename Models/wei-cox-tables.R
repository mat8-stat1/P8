
#################################################
# Coefficient tables for Weibull and Cox models #
#################################################

library(xtable)


#### Sourcing necessary files (loading all model objects)
source("data_cleaning.R")
source("IPS.R")
source("weibull.R")
source("coxph.R")
source("RSF.R")


### Weibull coefficient table:
mod.names.wei <- c("weibull.aic", "weibull.aic.spl", "weibull.bic", "weibull.bic.spl", "weibull.lasso")

wei.tab           <- as.data.frame(matrix(nrow = length(weibull.lasso$beta) + 3, ncol = 5))
colnames(wei.tab) <- mod.names.wei
rownames(wei.tab) <- c(rownames(weibull.lasso$beta)[1:2], "ns(Alder, df = 3)1", "ns(Alder, df = 3)2", "ns(Alder, df = 3)3", rownames(weibull.lasso$beta)[3:21])

wei.tab[names(weibull.aic$coefficients), 1]     <- weibull.aic$coefficients %>% round(4)
wei.tab[names(weibull.aic.spl$coefficients), 2] <- weibull.aic.spl$coefficients %>% round(4)
wei.tab[names(weibull.bic$coefficients), 3]     <- weibull.bic$coefficients %>% round(4)
wei.tab[names(weibull.bic.spl$coefficients), 4] <- weibull.bic.spl$coefficients %>% round(4)
wei.tab[c(1:2, 6:24), 5]                        <- weibull.lasso$beta %>% round(4)

wei.tab["LDH", "weibull.aic"]     <- wei.tab["LDHElevated", "weibull.aic"]
wei.tab["sex", "weibull.aic.spl"] <- wei.tab["sexmale", "weibull.aic.spl"]

wei.tab[wei.tab == 0] <- NA
wei.tab <- wei.tab[1:24,]

wei.tab[] <- lapply(wei.tab, function(x) format(x, scientific = FALSE))
wei.tab[wei.tab=="     NA"] <- "\\cellcolor{Gray!50}"

wei.tab <- xtable(wei.tab)
align(wei.tab) <- c("l", rep("r", 5))
caption(wei.tab) <- "Coefficient estimates of the selected Weibull models. Grey cells represent non-included covariates for the respective models. The ns() notation refers to the coefficient estimates of a natural cubic spline."
label(wei.tab) <- "tab:wei-coef"




### Coxph coefficient table:
mod.names.cox <- c("coxph.aic", "coxph.aic.spl", "coxph.bic", "coxph.bic.spl", "coxph.lasso")

cox.tab           <- as.data.frame(matrix(nrow = length(coxph.lasso$beta) + 3, ncol = 5))
colnames(cox.tab) <- mod.names.cox
rownames(cox.tab) <- c(rownames(coxph.lasso$beta)[1], "ns(Alder, df = 3)1", "ns(Alder, df = 3)2", "ns(Alder, df = 3)3", rownames(coxph.lasso$beta)[2:20])

cox.tab[names(coxph.aic$coefficients), 1]     <- coxph.aic$coefficients %>% round(4)
cox.tab[names(coxph.aic.spl$coefficients), 2] <- coxph.aic.spl$coefficients %>% round(4)
cox.tab[names(coxph.bic$coefficients), 3]     <- coxph.bic$coefficients %>% round(4)
cox.tab[names(coxph.bic.spl$coefficients), 4] <- coxph.bic.spl$coefficients %>% round(4)
cox.tab[c(1, 5:23), 5]                        <- coxph.lasso$beta %>% round(4)

cox.tab["bsym", "coxph.aic"] <- cox.tab["bsymJa", "coxph.aic"]
cox.tab["bsym", "coxph.aic.spl"] <- cox.tab["bsymJa", "coxph.aic.spl"]
cox.tab["LDH", "coxph.aic"]     <- cox.tab["LDHElevated", "coxph.aic"]
cox.tab["LDH", "coxph.aic.spl"]     <- cox.tab["LDHElevated", "coxph.aic.spl"]

cox.tab[cox.tab == 0] <- NA
cox.tab <- cox.tab[1:23,]

cox.tab[] <- lapply(cox.tab, function(x) format(x, scientific = FALSE))
cox.tab[cox.tab=="     NA"] <- "\\cellcolor{Gray!50}"

cox.tab <- xtable(cox.tab)
align(cox.tab) <- c("l", rep("r", 5))
caption(cox.tab) <- "Coefficient estimates of the selected Cox proportional hazards models. Grey cells represent non-included covariates for the respective models. The ns() notation refers to the coefficient estimates of a natural cubic spline."
label(cox.tab) <- "tab:coxph-coef"








