
#########################################
# Plotting histograms of all covariates #
#########################################

source("data_cleaning.R")
source("RSF.R")
library(ggplot2)
library(tidyverse)

################### Plot Alder
Alder_mean <- data %>% 
  dplyr::select(Alder) %>%
  summarize(mean = mean(Alder)) %>% 
  round(1)

ggplot(data, aes(x = Alder)) +
  geom_histogram(aes(y = ..density..), # the histogram will display "density" on its y-axis
                 bins = 26, colour = "black", fill = "steelblue1") +
  geom_density(alpha = .35, fill="chartreuse3") + #overlay with a transparent (alpha value) density plot
  geom_vline(aes(xintercept=mean(Alder),
                 color="Mean"), size=1.7, show.legend = FALSE) +
  geom_vline(aes(xintercept=median(Alder),
                 color="Median"), linetype = "dashed", size=1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(Alder, 0.25), color = "1. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(Alder, 0.75), color = "3. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  scale_color_manual(name = "Statistics", values = c(Median = "chartreuse", `1. quantile` = "chartreuse",
                                                     `3. quantile` = "chartreuse", Mean = "gold"),
                     limits = c("Mean", "1. quantile", "Median", "3. quantile")) +
  geom_text(data = Alder_mean, aes(x = mean, label = mean), 
            y = 0.01, x=55, color = "gold", size = 6) +
  xlab("") +
  ylab("Density") +
  theme_minimal()

################### Plot sex
df <- data %>% 
  dplyr::select(sex) %>% 
  group_by(sex) %>%
  summarise(counts = n())

ggplot(df, aes(x = sex, y = counts)) +
  geom_bar(fill = "steelblue1", stat = "identity") +
  geom_text(aes(label = counts), vjust = -0.3) +
  xlab("") +
  theme_minimal() +
  ylim(c(0,800))

################### Plot stadium
df <- data %>% 
  dplyr::select(Stadium..Ann.Arbor.) %>% 
  group_by(Stadium..Ann.Arbor.) %>%
  summarise(counts = n())

ggplot(df, aes(x = Stadium..Ann.Arbor., y = counts)) +
  geom_bar(fill = "steelblue1", stat = "identity") +
  geom_text(aes(label = counts), vjust = -0.3) +
  xlab("") +
  ylim(c(0,530)) +
  theme_minimal()

################## Plot ECOG
df <- data %>% 
  dplyr::select(ECOG) %>%
  group_by(ECOG) %>%
  summarise(counts = n())

ggplot(df, aes(x = ECOG, y = counts)) +
  geom_bar(fill = "steelblue1", stat = "identity") +
  geom_text(aes(label = counts), vjust = -0.3) +
  xlab("") +
  ylim(c(0, 870)) +
  theme_minimal()

################## Plot bsym
df <- data %>% 
  dplyr::select(bsym) %>%
  group_by(bsym) %>%
  summarise(counts = n())

ggplot(df, aes(x = bsym, y = counts)) +
  geom_bar(fill = "steelblue1", stat = "identity") +
  geom_text(aes(label = counts), vjust = -0.3) +
  xlab("") +
  ylim(c(0,680)) +
  theme_minimal()

################# Plot nodal
df <- data %>% 
  dplyr::select(nodal) %>%
  group_by(nodal) %>%
  summarise(counts = n())

ggplot(df, aes(x = factor(nodal), y = counts)) +
  geom_bar(fill = "steelblue1", stat = "identity") +
  geom_text(aes(label = counts), vjust = -0.3) +
  xlab("") +
  theme_minimal() +
  ylim(c(0,350))

################## Plot ENODAL
df <- data %>% 
  dplyr::select(ENODAL) %>%
  group_by(ENODAL) %>%
  summarise(counts = n())

ggplot(df, aes(x = factor(ENODAL), y = counts)) +
  geom_bar(fill = "steelblue1", stat = "identity") +
  geom_text(aes(label = counts), vjust = -0.3) +
  xlab("") +
  theme_minimal() +
  ylim(c(0,1000))

################## Plot tumordiameter
tumordiameter_mean <- data %>% 
  dplyr::select(tumordiameter) %>%
  drop_na %>% 
  summarize(mean = mean(tumordiameter)) %>% 
  round(1)

data_tumordiameterNA <- data %>% dplyr::select(tumordiameter) %>% drop_na

ggplot(data_tumordiameterNA, aes(x = tumordiameter)) +
  geom_histogram(aes(y = ..density..), # the histogram will display "density" on its y-axis
                 bins = 25, colour = "black", fill = "steelblue1") +
  geom_density(alpha = .35, fill="chartreuse3") + #overlay with a transparent (alpha value) density plot
  geom_vline(aes(xintercept=mean(tumordiameter),
                 color="Mean"), size=1.7, show.legend = FALSE) +
  geom_vline(aes(xintercept=median(tumordiameter),
                 color="Median"), linetype = "dashed", size=1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(tumordiameter, 0.25), color = "1. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(tumordiameter, 0.75), color = "3. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  scale_color_manual(name = "Statistics", values = c(Median = "chartreuse", `1. quantile` = "chartreuse",
                                                     `3. quantile` = "chartreuse", Mean = "gold"),
                     limits = c("Mean", "1. quantile", "Median", "3. quantile")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(2,2,2,2), shape=c(NA,NA,NA,NA)))) +
  geom_text(data = tumordiameter_mean, aes(x = mean, label = mean), 
            y = .04, x=7, color = "gold", size = 6) +
  xlab("") +
  ylab("Density") +
  theme_minimal()

################## Plot Hemoglobin
Hemoglobin_mean <- data %>% 
  dplyr::select(Hemoglobin) %>%
  drop_na %>% 
  summarize(mean = mean(Hemoglobin)) %>% 
  round(1)

data_HemoglobinNA <- data %>% dplyr::select(Hemoglobin) %>% drop_na

ggplot(data_HemoglobinNA, aes(x = Hemoglobin)) +
  geom_histogram(aes(y = ..density..), # the histogram will display "density" on its y-axis
                 bins = 25, colour = "black", fill = "steelblue1") +
  geom_density(alpha = .35, fill="chartreuse3") + #overlay with a transparent (alpha value) density plot
  geom_vline(aes(xintercept=mean(Hemoglobin),
                 color="Mean"), size=1.7, show.legend = FALSE) +
  geom_vline(aes(xintercept=median(Hemoglobin),
                 color="Median"), linetype = "dashed", size=1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(Hemoglobin, 0.25), color = "1. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(Hemoglobin, 0.75), color = "3. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  scale_color_manual(name = "Statistics", values = c(Median = "chartreuse", `1. quantile` = "chartreuse",
                                                     `3. quantile` = "chartreuse", Mean = "gold"),
                     limits = c("Mean", "1. quantile", "Median", "3. quantile")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(2,2,2,2), shape=c(NA,NA,NA,NA)))) +
  geom_text(data = Hemoglobin_mean, aes(x = mean, label = mean), 
            y = 0.01, x=110, color = "gold", size = 6) +
  xlab("") +
  ylab("Density") +
  theme_minimal()

################## Plot LDH
df <- data %>% 
  dplyr::select(LDH) %>%
  mutate(LDH = as.character(LDH)) %>% 
  mutate(LDH = case_when(is.na(LDH) ~ "Not Available", T ~ LDH)) %>% 
  mutate(LDH = as.factor(LDH)) %>% 
  group_by(LDH) %>%
  summarise(counts = n())

ggplot(df, aes(x = factor(LDH), y = counts)) +
  geom_bar(fill = "steelblue1", stat = "identity") +
  geom_text(aes(label = counts), vjust = -0.3) +
  xlab("") +
  theme_minimal() +
  ylim(c(0,910))

################## Plot Albumin
Albumin_mean <- data %>% 
  dplyr::select(Albumin) %>%
  drop_na %>% 
  summarize(mean = mean(Albumin)) %>% 
  round(1)

data_AlbuminNA <- data %>% dplyr::select(Albumin) %>% drop_na

ggplot(data_AlbuminNA, aes(x = Albumin)) +
  geom_histogram(aes(y = ..density..), # the histogram will display "density" on its y-axis
                 bins = 25, colour = "black", fill = "steelblue1") +
  geom_density(alpha = .35, fill="chartreuse3") + #overlay with a transparent (alpha value) density plot
  geom_vline(aes(xintercept=mean(Albumin),
                 color="Mean"), size=1.7, show.legend = FALSE) +
  geom_vline(aes(xintercept=median(Albumin),
                 color="Median"), linetype = "dashed", size=1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(Albumin, 0.25), color = "1. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(Albumin, 0.75), color = "3. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  scale_color_manual(name = "Statistics", values = c(Median = "chartreuse", `1. quantile` = "chartreuse",
                                                     `3. quantile` = "chartreuse", Mean = "gold"),
                     limits = c("Mean", "1. quantile", "Median", "3. quantile")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(2,2,2,2), shape=c(NA,NA,NA,NA)))) +
  geom_text(data = Albumin_mean, aes(x = mean, label = mean), 
            y = .03, x=33, color = "gold", size = 6) +
  xlab("") +
  ylab("Density") +
  theme_minimal()

################## Plot Lymfocytter.mia.L.
Lymfocytter.mia.L._mean <- data %>% 
  dplyr::select(Lymfocytter.mia.L.) %>%
  drop_na %>% 
  summarize(mean = mean(Lymfocytter.mia.L.)) %>% 
  round(1)

data_Lymfocytter.mia.L.NA <- data %>% dplyr::select(Lymfocytter.mia.L.) %>% drop_na

ggplot(data_Lymfocytter.mia.L.NA, aes(x = Lymfocytter.mia.L.)) +
  geom_histogram(aes(y = ..density..), # the histogram will display "density" on its y-axis
                 bins = 25, colour = "black", fill = "steelblue1") +
  geom_density(alpha = .35, fill="chartreuse3") + #overlay with a transparent (alpha value) density plot
  geom_vline(aes(xintercept=mean(Lymfocytter.mia.L.),
                 color="Mean"), size=1.7, show.legend = FALSE) +
  geom_vline(aes(xintercept=median(Lymfocytter.mia.L.),
                 color="Median"), linetype = "dashed", size=1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(Lymfocytter.mia.L., 0.25), color = "1. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(Lymfocytter.mia.L., 0.75), color = "3. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  scale_color_manual(name = "Statistics", values = c(Median = "chartreuse", `1. quantile` = "chartreuse",
                                                     `3. quantile` = "chartreuse", Mean = "gold"),
                     limits = c("Mean", "1. quantile", "Median", "3. quantile")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(2,2,2,2), shape=c(NA,NA,NA,NA)))) +
  geom_text(data = Lymfocytter.mia.L._mean, aes(x = mean, label = mean), 
            y = .3, x=1.1, color = "gold", size = 6) +
  xlab("") +
  ylab("Density") +
  xlim(c(-.5, 10)) +
  theme_minimal()

################## Plot Leukocytter.mia.L.
Leukocytter.mia.L._mean <- data %>% 
  dplyr::select(Leukocytter.mia.L.) %>%
  drop_na %>% 
  summarize(mean = mean(Leukocytter.mia.L.)) %>% 
  round(1)

data_Leukocytter.mia.L.NA <- data %>% dplyr::select(Leukocytter.mia.L.) %>% drop_na

ggplot(data_Leukocytter.mia.L.NA, aes(x = Leukocytter.mia.L.)) +
  geom_histogram(aes(y = ..density..), # the histogram will display "density" on its y-axis
                 bins = 25, colour = "black", fill = "steelblue1") +
  geom_density(alpha = .35, fill="chartreuse3") + #overlay with a transparent (alpha value) density plot
  geom_vline(aes(xintercept=mean(Leukocytter.mia.L.),
                 color="Mean"), size=1.7, show.legend = FALSE) +
  geom_vline(aes(xintercept=median(Leukocytter.mia.L.),
                 color="Median"), linetype = "dashed", size=1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(Leukocytter.mia.L., 0.25), color = "1. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(Leukocytter.mia.L., 0.75), color = "3. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  scale_color_manual(name = "Statistics", values = c(Median = "chartreuse", `1. quantile` = "chartreuse",
                                                     `3. quantile` = "chartreuse", Mean = "gold"),
                     limits = c("Mean", "1. quantile", "Median", "3. quantile")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(2,2,2,2), shape=c(NA,NA,NA,NA)))) +
  geom_text(data = Leukocytter.mia.L._mean, aes(x = mean, label = mean), 
            y = .05, x=7.6, color = "gold", size = 6) +
  xlab("") +
  ylab("Density") +
  xlim(c(0,40)) +
  theme_minimal()

################## Plot Beta2Microglobulin
Beta2Microglobulin_mean <- data %>% 
  dplyr::select(Beta2Microglobulin) %>%
  drop_na %>% 
  summarize(mean = mean(Beta2Microglobulin)) %>% 
  round(3)

data_Beta2MicroglobulinNA <- data %>% dplyr::select(Beta2Microglobulin) %>% drop_na

ggplot(data_Beta2MicroglobulinNA, aes(x = Beta2Microglobulin)) +
  geom_histogram(aes(y = ..density..), # the histogram will display "density" on its y-axis
                 bins = 25, colour = "black", fill = "steelblue1") +
  geom_density(alpha = .35, fill="chartreuse3") + #overlay with a transparent (alpha value) density plot
  geom_vline(aes(xintercept=mean(Beta2Microglobulin),
                 color="Mean"), size=1.7, show.legend = FALSE) +
  geom_vline(aes(xintercept=median(Beta2Microglobulin),
                 color="Median"), linetype = "dashed", size=1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(Beta2Microglobulin, 0.25), color = "1. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(Beta2Microglobulin, 0.75), color = "3. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  scale_color_manual(name = "Statistics", values = c(Median = "chartreuse", `1. quantile` = "chartreuse",
                                                     `3. quantile` = "chartreuse", Mean = "gold"),
                     limits = c("Mean", "1. quantile", "Median", "3. quantile")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(2,2,2,2), shape=c(NA,NA,NA,NA)))) +
  geom_text(data = Beta2Microglobulin_mean, aes(x = mean, label = mean), 
            y = 30, x=0.0047, color = "gold", size = 6) +
  xlab("") +
  ylab("Density") +
  xlim(c(0,0.015)) +
  theme_minimal()

################## Plot IgG
IgG_mean <- data %>% 
  dplyr::select(IgG) %>%
  drop_na %>% 
  summarize(mean = mean(IgG)) %>% 
  round(1)

data_IgGNA <- data %>% dplyr::select(IgG) %>% drop_na

ggplot(data_IgGNA, aes(x = IgG)) +
  geom_histogram(aes(y = ..density..), # the histogram will display "density" on its y-axis
                 bins = 25, colour = "black", fill = "steelblue1") +
  geom_density(alpha = .35, fill="chartreuse3") + #overlay with a transparent (alpha value) density plot
  geom_vline(aes(xintercept=mean(IgG),
                 color="Mean"), size=1.7, show.legend = FALSE) +
  geom_vline(aes(xintercept=median(IgG),
                 color="Median"), linetype = "dashed", size=1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(IgG, 0.25), color = "1. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(IgG, 0.75), color = "3. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  scale_color_manual(name = "Statistics", values = c(Median = "chartreuse", `1. quantile` = "chartreuse",
                                                     `3. quantile` = "chartreuse", Mean = "gold"),
                     limits = c("Mean", "1. quantile", "Median", "3. quantile")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(2,2,2,2), shape=c(NA,NA,NA,NA)))) +
  geom_text(data = IgG_mean, aes(x = mean, label = mean), 
            y = 0.012, x=17.9, color = "gold", size = 6) +
  xlab("") +
  ylab("Density") +
  xlim(c(-0.15, 50)) +
  theme_minimal()

################## Plot IgM
IgM_mean <- data %>% 
  dplyr::select(IgM) %>%
  drop_na %>% 
  summarize(mean = mean(IgM)) %>% 
  round(1)

data_IgMNA <- data %>% dplyr::select(IgM) %>% drop_na

ggplot(data_IgMNA, aes(x = IgM)) +
  geom_histogram(aes(y = ..density..), # the histogram will display "density" on its y-axis
                 bins = 25, colour = "black", fill = "steelblue1") +
  geom_density(alpha = .35, fill="chartreuse3") + #overlay with a transparent (alpha value) density plot
  geom_vline(aes(xintercept=mean(IgM),
                 color="Mean"), size=1.7, show.legend = FALSE) +
  geom_vline(aes(xintercept=median(IgM),
                 color="Median"), linetype = "dashed", size=1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(IgM, 0.25), color = "1. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(IgM, 0.75), color = "3. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  scale_color_manual(name = "Statistics", values = c(Median = "chartreuse", `1. quantile` = "chartreuse",
                                                     `3. quantile` = "chartreuse", Mean = "gold"),
                     limits = c("Mean", "1. quantile", "Median", "3. quantile")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(2,2,2,2), shape=c(NA,NA,NA,NA)))) +
  geom_text(data = IgM_mean, aes(x = mean, label = mean), 
            y = 0.1, x=1.7, color = "gold", size = 6) +
  xlab("") +
  ylab("Density") +
  xlim(c(-0.1, 7)) +
  theme_minimal()

################## Plot IgA
IgA_mean <- data %>% 
  dplyr::select(IgA) %>%
  drop_na %>% 
  summarize(mean = mean(IgA)) %>% 
  round(1)

data_IgANA <- data %>% dplyr::select(IgA) %>% drop_na

ggplot(data_IgANA, aes(x = IgA)) +
  geom_histogram(aes(y = ..density..), # the histogram will display "density" on its y-axis
                 bins = 25, colour = "black", fill = "steelblue1") +
  geom_density(alpha = .35, fill="chartreuse3") + #overlay with a transparent (alpha value) density plot
  geom_vline(aes(xintercept=mean(IgA),
                 color="Mean"), size=1.7, show.legend = FALSE) +
  geom_vline(aes(xintercept=median(IgA),
                 color="Median"), linetype = "dashed", size=1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(IgA, 0.25), color = "1. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(IgA, 0.75), color = "3. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  scale_color_manual(name = "Statistics", values = c(Median = "chartreuse", `1. quantile` = "chartreuse",
                                                     `3. quantile` = "chartreuse", Mean = "gold"),
                     limits = c("Mean", "1. quantile", "Median", "3. quantile")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(2,2,2,2), shape=c(NA,NA,NA,NA)))) +
  geom_text(data = IgA_mean, aes(x = mean, label = mean), 
            y = 0.07, x=3.8, color = "gold", size = 6) +
  xlab("") +
  ylab("Density") +
  xlim(c(-0.1, 12)) +
  theme_minimal()

################## Plot CreatineMikroMol
CreatineMikroMol_mean <- data %>% 
  dplyr::select(CreatineMikroMol) %>%
  drop_na %>% 
  summarize(mean = mean(CreatineMikroMol)) %>% 
  round(1)

data_CreatineMikroMolNA <- data %>% dplyr::select(CreatineMikroMol) %>% drop_na

ggplot(data_CreatineMikroMolNA, aes(x = CreatineMikroMol)) +
  geom_histogram(aes(y = ..density..), # the histogram will display "density" on its y-axis
                 bins = 25, colour = "black", fill = "steelblue1") +
  geom_density(alpha = .35, fill="chartreuse3") + #overlay with a transparent (alpha value) density plot
  geom_vline(aes(xintercept=mean(CreatineMikroMol),
                 color="Mean"), size=1.7, show.legend = FALSE) +
  geom_vline(aes(xintercept=median(CreatineMikroMol),
                 color="Median"), linetype = "dashed", size=1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(CreatineMikroMol, 0.25), color = "1. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(CreatineMikroMol, 0.75), color = "3. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  scale_color_manual(name = "Statistics", values = c(Median = "chartreuse", `1. quantile` = "chartreuse",
                                                     `3. quantile` = "chartreuse", Mean = "gold"),
                     limits = c("Mean", "1. quantile", "Median", "3. quantile")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(2,2,2,2), shape=c(NA,NA,NA,NA)))) +
  geom_text(data = CreatineMikroMol_mean, aes(x = mean, label = mean), 
            y = 0.005, x=89, color = "gold", size = 6) +
  xlab("") +
  ylab("Density") +
  xlim(c(20, 230)) +
  theme_minimal()

################## Plot Thrombocytter.mia.L.
Thrombocytter.mia.L._mean <- data %>% 
  dplyr::select(Thrombocytter.mia.L.) %>%
  drop_na %>% 
  summarize(mean = mean(Thrombocytter.mia.L.)) %>% 
  round(1)

data_Thrombocytter.mia.L.NA <- data %>% dplyr::select(Thrombocytter.mia.L.) %>% drop_na

ggplot(data_Thrombocytter.mia.L.NA, aes(x = Thrombocytter.mia.L.)) +
  geom_histogram(aes(y = ..density..), # the histogram will display "density" on its y-axis
                 bins = 25, colour = "black", fill = "steelblue1") +
  geom_density(alpha = .35, fill="chartreuse3") + #overlay with a transparent (alpha value) density plot
  geom_vline(aes(xintercept=mean(Thrombocytter.mia.L.),
                 color="Mean"), size=1.7, show.legend = FALSE) +
  geom_vline(aes(xintercept=median(Thrombocytter.mia.L.),
                 color="Median"), linetype = "dashed", size=1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(Thrombocytter.mia.L., 0.25), color = "1. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  geom_vline(aes(xintercept = quantile(Thrombocytter.mia.L., 0.75), color = "3. quantile"), linetype = "dashed", size = 1.2, show.legend = FALSE) +
  scale_color_manual(name = "Statistics", values = c(Median = "chartreuse", `1. quantile` = "chartreuse",
                                                     `3. quantile` = "chartreuse", Mean = "gold"),
                     limits = c("Mean", "1. quantile", "Median", "3. quantile")) +
  guides(colour = guide_legend(override.aes = list(linetype=c(2,2,2,2), shape=c(NA,NA,NA,NA)))) +
  geom_text(data = Thrombocytter.mia.L._mean, aes(x = mean, label = mean), 
            y = 0.001, x=440, color = "gold", size = 6) +
  xlab("") +
  ylab("Density") +
  theme_minimal()

################## Plot organ.inv categorical
df <- data.imp %>% 
  dplyr::select(organ.inv) %>%
  group_by(organ.inv) %>%
  summarise(counts = n())

ggplot(df, aes(x = factor(organ.inv), y = counts)) +
  geom_bar(fill = "steelblue1", stat = "identity") +
  geom_text(aes(label = counts), vjust = -0.3) +
  xlab("") +
  theme_minimal() +
  ylim(c(0,320))
