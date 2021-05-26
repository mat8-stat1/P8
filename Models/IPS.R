
###################################################################
# Fitting IPS model, plotting curves and table of 5-year survival #
###################################################################

library(survival)
library(survminer) #ggsurvplot
library(grDevices) #pdf
library(xtable)

source("data_cleaning.R")

data.imp_IPS <- data.imp %>% mutate(IPS_groups = case_when(IPS %in% 5:7 ~ 5L, TRUE ~ IPS))

IPS <- survfit(Surv(OS,status)
               ~ IPS_groups, data = data.imp_IPS)

ggsurvplot(IPS, data = data.imp_IPS,
           conf.int = TRUE,
           linetype = "solid",
           legend = "top",
           legend.title = "IPS groups",
           legend.labs = 0:5,
           xlab ="Time since diagnosis (years)")

sum_IPS <- summary(IPS, times = 5)

sum_IPS_dat <- data.frame(0:5, sum_IPS$n.risk, sum_IPS$n.event, paste0(round(sum_IPS$surv,3), " (", round(sum_IPS$std.err,3), ")"),
           paste0("[", round(sum_IPS$lower,3), ", ", round(sum_IPS $upper,3), "]"))
colnames(sum_IPS_dat) <- c("IPS group", "At risk", "Cumulative number of events", "Survival probability (sd)",
                           "95% confidence interval")

sum_IPS_dat
