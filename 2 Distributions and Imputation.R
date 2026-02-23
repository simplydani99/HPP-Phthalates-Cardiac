# Programmer: Danielle Stevens
# Purpose of syntax: HPP analysis examining phthalates, alternatives, & fetal cardiac outcomes
# Last updated on: 01/26/2026


#libraries
library(pheatmap)
library(dichromat)
library(haven)
library(knitr)
library(stringr)
library(tidyverse)
library(plyr)
library(tidyr)
library(dplyr) 
library(reshape2)
library(geepack)
library(ggplot2)
library(gridExtra)
library(gtable)
library(grid)
library(itsadug)
library(lattice)
library(scales)
library(mgcv)
library(ggbreak) 
library(patchwork)
library(mice)
library(qgcomp)
library(qgcompint)
library(patchwork)
library(cowplot)
library(RColorBrewer)
library(janitor)


#####################
## Read in dataset ##
#####################

file <- read_sas("J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Derived_data\\hpp_cardiac_analytic.sas7bdat")


###########################
## Set working directory ##
###########################

setwd("J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results")


##################
## Define lists ##
##################

# Covariates
unadj.covlist <- c("spl1", "spl2")
covlist <- c('spl1', "spl2", 'drugs_new', 'alcohol_new', 'age', 'bmi', 'smoking','gender_nb', 'highest_edu', 'race_new', 'any_prior_preg')
covlist.nospl1 <- c('drugs_new', 'alcohol_new', 'age', 'bmi', 'smoking','gender_nb', 'highest_edu', 'race_new', 'any_prior_preg')

# Exposures
long_explist_pht <- c('mep', 'dnbp', 'dibp', 'mbzp', 'mcpp', 'dehp', 'dinp', 'mcnp')
long_explist_lmw <- c("mep", "dnbp", "dibp")
long_explist_hmw <- c('mbzp', 'mcpp', 'dehp', 'dinp', 'mcnp')
long_explist_alt <- c('dehtp','dinch')
long_explist <- c(long_explist_pht, long_explist_alt)
long_explist_log <- paste0('ln_', long_explist)
long_explist_iqr <- paste0('iqr_', long_explist)
long_explist_iqr_pht <- paste0('iqr_', long_explist_pht)
long_explist_iqr_lmw <- paste0('iqr_', long_explist_lmw)
long_explist_iqr_hmw <- paste0('iqr_', long_explist_hmw)
long_explist_iqr_alt <- paste0('iqr_', long_explist_alt)

explist_pht <- c('preg_mep', 'preg_dnbp', 'preg_dibp', 'preg_mbzp', 'preg_mcpp', 'preg_dehp', 'preg_dinp', 'preg_mcnp')
explist_lmw <- c("preg_mep", "preg_dnbp", "preg_dibp")
explist_hmw <- c('preg_mbzp', 'preg_mcpp', 'preg_dehp', 'preg_dinp', 'preg_mcnp')
explist_alt <- c('preg_dehtp','preg_dinch')
explist <- c(explist_pht, explist_alt)
explist_log <- paste0('ln_', explist)
explist_iqr <- paste0('iqr_', explist)
explist_iqr_pht <- paste0('iqr_', explist_pht)
explist_iqr_lmw <- paste0('iqr_', explist_lmw)
explist_iqr_hmw <- paste0('iqr_', explist_hmw)
explist_iqr_alt <- paste0('iqr_', explist_alt)

t1_explist_pht <- c('tri1_mep', 'tri1_dnbp', 'tri1_dibp', 'tri1_mbzp', 'tri1_mcpp', 'tri1_dehp', 'tri1_dinp', 'tri1_mcnp')
t1_explist_alt <- c('tri1_dehtp','tri1_dinch')
t1_explist_lmw <- c("tri1_mep", "tri1_dnbp", "tri1_dibp")
t1_explist_hmw <- c('tri1_mbzp', 'tri1_mcpp', 'tri1_dehp', 'tri1_dinp', 'tri1_mcnp')
t1_explist <- c(t1_explist_pht, t1_explist_alt)
t1_explist_log <- paste0('ln_', t1_explist)
t1_explist_iqr <- paste0('iqr_', t1_explist)
t1_explist_iqr_pht <- paste0('iqr_', t1_explist_pht)
t1_explist_iqr_lmw <- paste0('iqr_', t1_explist_lmw)
t1_explist_iqr_hmw <- paste0('iqr_', t1_explist_hmw)
t1_explist_iqr_alt <- paste0('iqr_', t1_explist_alt)

t2_explist_pht <- c('tri2_mep', 'tri2_dnbp', 'tri2_dibp', 'tri2_mbzp', 'tri2_mcpp', 'tri2_dehp', 'tri2_dinp', 'tri2_mcnp')
t2_explist_alt <- c('tri2_dehtp','tri2_dinch')
t2_explist_lmw <- c("tri2_mep", "tri2_dnbp", "tri2_dibp")
t2_explist_hmw <- c('tri2_mbzp', 'tri2_mcpp', 'tri2_dehp', 'tri2_dinp', 'tri2_mcnp')
t2_explist <- c(t2_explist_pht, t2_explist_alt)
t2_explist_log <- paste0('ln_', t2_explist)
t2_explist_iqr <- paste0('iqr_', t2_explist)
t2_explist_iqr_pht <- paste0('iqr_', t2_explist_pht)
t2_explist_iqr_lmw <- paste0('iqr_', t2_explist_lmw)
t2_explist_iqr_hmw <- paste0('iqr_', t2_explist_hmw)
t2_explist_iqr_alt <- paste0('iqr_', t2_explist_alt)

t3_explist_pht <- c('tri3_mep', 'tri3_dnbp', 'tri3_dibp', 'tri3_mbzp', 'tri3_mcpp', 'tri3_dehp', 'tri3_dinp', 'tri3_mcnp')
t3_explist_alt <- c('tri3_dehtp','tri3_dinch')
t3_explist_lmw <- c("tri3_mep", "tri3_dnbp", "tri3_dibp")
t3_explist_hmw <- c('tri3_mbzp', 'tri3_mcpp', 'tri3_dehp', 'tri3_dinp', 'tri3_mcnp')
t3_explist <- c(t3_explist_pht, t3_explist_alt)
t3_explist_log <- paste0('ln_', t3_explist)
t3_explist_iqr <- paste0('iqr_', t3_explist)
t3_explist_iqr_pht <- paste0('iqr_', t3_explist_pht)
t3_explist_iqr_lmw <- paste0('iqr_', t3_explist_lmw)
t3_explist_iqr_hmw <- paste0('iqr_', t3_explist_hmw)
t3_explist_iqr_alt <- paste0('iqr_', t3_explist_alt)

abbrevs.exp <- c('MEP', 'ΣDnBP', 'ΣDiBP',  'MBzP', 'MCPP', 'ΣDEHP', 'ΣDiNP', 'MCNP', 'ΣDEHTP','ΣDINCH')

outlist <- c('chest_area', 'heart_area', 'z_heart_area', 'z_chest_area', 'cardiothoracic_ratio',
  'rv_inlet_length', 'lv_inlet_length', 'z_rv_inlet_length', 'z_lv_inlet_length',
  'mmode_rv_tapse', 'mmode_lv_mapse', 'z_rv_tapse', 'z_lv_mapse',
  'mmode_rv_fs', 'mmode_lv_fs',
  'std_rv_tei', 'std_lvent_tei',
  'tricuspid_e', 'tricuspid_a', 'tricuspid_earatio',
  'mitral_e', 'mitral_a', 'mitral_earatio')


##############################################
##### Get variables set up for analysis ######
##############################################

#Lowercase file names
names(file) <- tolower(names(file))

#Replace blank with missing
file[file == ""] <- NA

#Dataset structure: Update covariate names to conform with old names prior to Kate dataset updates in 2025
file <- file %>% mutate(record = as.factor(record),
                        visit = as.factor(visit),
                        smoking = as.factor(smoking),
                        gender_nb = as.factor(sex_nb),
                        highest_edu = as.factor(highest_edu_r),
                        race_new = as.factor(race_eth),
                        any_prior_preg = as.factor(any_prior_preg),
                        clinicsite = as.factor(clinicsite),
                        alcohol_new = as.factor(alcohol_new),
                        drugs_new = as.factor(drugs_new),
                        insurance = as.factor(insurance_r),
                        marital_status = as.factor(marital_status_r),
                        employment = as.factor(employment),
                        income = as.factor(income))


#######################
## Plots of outcomes ##
#######################

heart_area <- file[, c('gaweeks','heart_area',"cat_heart_area")]
heart_area <- heart_area[complete.cases(heart_area),]
heart_areaPlot <- ggplot(data = heart_area, aes(x = gaweeks, y = heart_area, color=as.factor(cat_heart_area), shape = as.factor(cat_heart_area))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nHeart area (mm2)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c('black', "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
heart_areaPlot

z_heart_area <- file[, c('gaweeks','z_heart_area',"cat_z_heart_area")]
z_heart_area <- z_heart_area[complete.cases(z_heart_area),]
z_heart_areaPlot <- ggplot(data = z_heart_area, aes(x = gaweeks, y = z_heart_area, color=as.factor(cat_z_heart_area), shape = as.factor(cat_z_heart_area))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nHeart area (z-score)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c('black', "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
z_heart_areaPlot

chest_area <- file[, c('gaweeks','chest_area',"cat_chest_area")]
chest_area <- chest_area[complete.cases(chest_area),]
chest_areaPlot <- ggplot(data = chest_area, aes(x = gaweeks, y = chest_area, color=as.factor(cat_chest_area), shape = as.factor(cat_chest_area))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nChest area (mm2)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c('black', "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
chest_areaPlot

z_chest_area <- file[, c('gaweeks','z_chest_area',"cat_z_chest_area")]
z_chest_area <- z_chest_area[complete.cases(z_chest_area),]
z_chest_areaPlot <- ggplot(data = z_chest_area, aes(x = gaweeks, y = z_chest_area, color=as.factor(cat_z_chest_area), shape = as.factor(cat_z_chest_area))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nChest area (z-score)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c('black', "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
z_chest_areaPlot

crt <- file[, c('gaweeks','cardiothoracic_ratio',"cat_cardiothoracic_ratio")]
crt <- crt[complete.cases(crt),]
crtPlot <- ggplot(data = crt, aes(x = gaweeks, y = cardiothoracic_ratio, color=as.factor(cat_cardiothoracic_ratio), shape = as.factor(cat_cardiothoracic_ratio))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nCardiothoracic ratio (%)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c('black', "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
crtPlot



rv_inlet_length <- file[, c('gaweeks','rv_inlet_length',"cat_rv_inlet_length")]
rv_inlet_length <- rv_inlet_length[complete.cases(rv_inlet_length),]
rv_inlet_lengthPlot <- ggplot(data = rv_inlet_length, aes(x = gaweeks, y = rv_inlet_length, color=as.factor(cat_rv_inlet_length), shape = as.factor(cat_rv_inlet_length))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nRV inlet length (mm)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c('black', "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
rv_inlet_lengthPlot

z_rv_inlet_length <- file[, c('gaweeks','z_rv_inlet_length',"cat_z_rv_inlet_length")]
z_rv_inlet_length <- z_rv_inlet_length[complete.cases(z_rv_inlet_length),]
z_rv_inlet_lengthPlot <- ggplot(data = z_rv_inlet_length, aes(x = gaweeks, y = z_rv_inlet_length, color=as.factor(cat_z_rv_inlet_length), shape = as.factor(cat_z_rv_inlet_length))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nRV inlet length (z-score)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c('black', "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
z_rv_inlet_lengthPlot

lv_inlet_length <- file[, c('gaweeks','lv_inlet_length',"cat_lv_inlet_length")]
lv_inlet_length <- lv_inlet_length[complete.cases(lv_inlet_length),]
lv_inlet_lengthPlot <- ggplot(data = lv_inlet_length, aes(x = gaweeks, y = lv_inlet_length, color=as.factor(cat_lv_inlet_length), shape = as.factor(cat_lv_inlet_length))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nLV inlet length (mm)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c('black', "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
lv_inlet_lengthPlot

z_lv_inlet_length <- file[, c('gaweeks','z_lv_inlet_length',"cat_z_lv_inlet_length")]
z_lv_inlet_length <- z_lv_inlet_length[complete.cases(z_lv_inlet_length),]
z_lv_inlet_lengthPlot <- ggplot(data = z_lv_inlet_length, aes(x = gaweeks, y = z_lv_inlet_length, color=as.factor(cat_z_lv_inlet_length), shape = as.factor(cat_z_lv_inlet_length))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nLV inlet length (z-score)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c('black', "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
z_lv_inlet_lengthPlot




mmode_rv_fs <- file[, c('gaweeks','mmode_rv_fs',"cat_mmode_rv_fs")]
mmode_rv_fs <- mmode_rv_fs[complete.cases(mmode_rv_fs),]
mmode_rv_fsPlot <- ggplot(data = mmode_rv_fs, aes(x = gaweeks, y = mmode_rv_fs, color=as.factor(cat_mmode_rv_fs), shape = as.factor(cat_mmode_rv_fs))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nRV fractional shortening (%)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c("black","grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  theme(text = element_text(size=14), 
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
mmode_rv_fsPlot

mmode_lv_fs <- file[, c('gaweeks','mmode_lv_fs',"cat_mmode_lv_fs")]
mmode_lv_fs <- mmode_lv_fs[complete.cases(mmode_lv_fs),]
mmode_lv_fsPlot <- ggplot(data = mmode_lv_fs, aes(x = gaweeks, y = mmode_lv_fs, color=as.factor(cat_mmode_lv_fs), shape = as.factor(cat_mmode_lv_fs))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nLV fractional shortening (%)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c("black","grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  theme(text = element_text(size=14), 
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
mmode_lv_fsPlot

std_rv_tei <- file[, c('gaweeks','std_rv_tei',"cat_std_rv_tei")]
std_rv_tei <- std_rv_tei[complete.cases(std_rv_tei),]
std_rv_teiPlot <- ggplot(data = std_rv_tei, aes(x = gaweeks, y = std_rv_tei, color=as.factor(cat_std_rv_tei), shape = as.factor(cat_std_rv_tei))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nRV MPI") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c("black","grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  theme(text = element_text(size=14), 
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
std_rv_teiPlot

std_lvent_tei <- file[, c('gaweeks','std_lvent_tei',"cat_std_lvent_tei")]
std_lvent_tei <- std_lvent_tei[complete.cases(std_lvent_tei),]
std_lvent_teiPlot <- ggplot(data = std_lvent_tei, aes(x = gaweeks, y = std_lvent_tei, color=as.factor(cat_std_lvent_tei), shape = as.factor(cat_std_lvent_tei))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nLV MPI") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c("black","grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  theme(text = element_text(size=14), 
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
std_lvent_teiPlot



tricuspid_e <- file[, c('gaweeks','tricuspid_e',"cat_tricuspid_e")]
tricuspid_e <- tricuspid_e[complete.cases(tricuspid_e),]
tricuspid_ePlot <- ggplot(data = tricuspid_e, aes(x = gaweeks, y = tricuspid_e, color=as.factor(cat_tricuspid_e), shape = as.factor(cat_tricuspid_e))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nTricuspid E (cm/s)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c("black", "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
tricuspid_ePlot

tricuspid_a <- file[, c('gaweeks','tricuspid_a',"cat_tricuspid_a")]
tricuspid_a <- tricuspid_a[complete.cases(tricuspid_a),]
tricuspid_aPlot <- ggplot(data = tricuspid_a, aes(x = gaweeks, y = tricuspid_a, color=as.factor(cat_tricuspid_a), shape = as.factor(cat_tricuspid_a))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nTricuspid A (cm/s)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c("black", "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
tricuspid_aPlot

tricuspid_earatio <- file[, c('gaweeks','tricuspid_earatio',"cat_tricuspid_earatio")]
tricuspid_earatio <- tricuspid_earatio[complete.cases(tricuspid_earatio),]
tricuspid_earatioPlot <- ggplot(data = tricuspid_earatio, aes(x = gaweeks, y = tricuspid_earatio, color=as.factor(cat_tricuspid_earatio), shape = as.factor(cat_tricuspid_earatio))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nTricuspid E/A") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c("black", "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
tricuspid_earatioPlot




mitral_e <- file[, c('gaweeks','mitral_e',"cat_mitral_e")]
mitral_e <- mitral_e[complete.cases(mitral_e),]
mitral_ePlot <- ggplot(data = mitral_e, aes(x = gaweeks, y = mitral_e, color=as.factor(cat_mitral_e), shape = as.factor(cat_mitral_e))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nMitral E (cm/s)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c("black", "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
mitral_ePlot

mitral_a <- file[, c('gaweeks','mitral_a',"cat_mitral_a")]
mitral_a <- mitral_a[complete.cases(mitral_a),]
mitral_aPlot <- ggplot(data = mitral_a, aes(x = gaweeks, y = mitral_a, color=as.factor(cat_mitral_a), shape = as.factor(cat_mitral_a))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nMitral A (cm/s)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c("black", "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
mitral_aPlot

mitral_earatio <- file[, c('gaweeks','mitral_earatio',"cat_mitral_earatio")]
mitral_earatio <- mitral_earatio[complete.cases(mitral_earatio),]
mitral_earatioPlot <- ggplot(data = mitral_earatio, aes(x = gaweeks, y = mitral_earatio, color=as.factor(cat_mitral_earatio), shape = as.factor(cat_mitral_earatio))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nMitral E/A") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c("black", "grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  scale_y_continuous(labels=label_number(accuracy = 1))+
  theme(text = element_text(size=14),
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
mitral_earatioPlot



mmode_rv_tapse <- file[, c('gaweeks','mmode_rv_tapse',"cat_mmode_rv_tapse")]
mmode_rv_tapse <- mmode_rv_tapse[complete.cases(mmode_rv_tapse),]
mmode_rv_tapsePlot <- ggplot(data = mmode_rv_tapse, aes(x = gaweeks, y = mmode_rv_tapse, color=as.factor(cat_mmode_rv_tapse), shape = as.factor(cat_mmode_rv_tapse))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nRV annulus displacement (mm)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c("black","grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  theme(text = element_text(size=14), 
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
mmode_rv_tapsePlot

mmode_lv_mapse <- file[, c('gaweeks','mmode_lv_mapse',"cat_mmode_lv_mapse")]
mmode_lv_mapse <- mmode_lv_mapse[complete.cases(mmode_lv_mapse),]
mmode_lv_mapsePlot <- ggplot(data = mmode_lv_mapse, aes(x = gaweeks, y = mmode_lv_mapse, color=as.factor(cat_mmode_lv_mapse), shape = as.factor(cat_mmode_lv_mapse))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nLV annulus displacement (mm)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c("black","grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  theme(text = element_text(size=14), 
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
mmode_lv_mapsePlot

z_rv_tapse <- file[, c('gaweeks','z_rv_tapse',"cat_z_rv_tapse")]
z_rv_tapse <- z_rv_tapse[complete.cases(z_rv_tapse),]
z_rv_tapsePlot <- ggplot(data = z_rv_tapse, aes(x = gaweeks, y = z_rv_tapse, color=as.factor(cat_z_rv_tapse), shape = as.factor(cat_z_rv_tapse))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nRV annulus displacement (z-score)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c("black","grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  theme(text = element_text(size=14), 
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
z_rv_tapsePlot

z_lv_mapse <- file[, c('gaweeks','z_lv_mapse',"cat_z_lv_mapse")]
z_lv_mapse <- z_lv_mapse[complete.cases(z_lv_mapse),]
z_lv_mapsePlot <- ggplot(data = z_lv_mapse, aes(x = gaweeks, y = z_lv_mapse, color=as.factor(cat_z_lv_mapse), shape = as.factor(cat_z_lv_mapse))) + 
  geom_point(size=2) +  xlab("Gestational age (weeks)\n") + ylab("\nLV annulus displacement (z-score)") + 
  scale_colour_manual(name="No. of measures \nby visit", values=c("black","grey37","grey69")) +
  scale_shape_discrete(name="No. of measures \nby visit") + theme_bw() +
  theme(text = element_text(size=14), 
        legend.position = "right", legend.text=element_text(size=11), legend.title = element_text(size=12)) 
z_lv_mapsePlot


comboPlot1 = grid.arrange(heart_areaPlot, z_heart_areaPlot, 
                          chest_areaPlot, z_chest_areaPlot,
                          ncol=2, widths = c(6,6))
comboPlot1

comboPlot2 = grid.arrange(crtPlot,
                          comboPlot1,
                          ncol=2, widths = c(6,12))
comboPlot2



comboPlot3 = grid.arrange(lv_inlet_lengthPlot, z_lv_inlet_lengthPlot,
                          rv_inlet_lengthPlot, z_rv_inlet_lengthPlot,
                          ncol=2, widths = c(6,6))
comboPlot3


ggsave(plot = comboPlot2, filename = "Distribution plot for cardiac and chest areas.png", width = 18, height = 6, dpi=800)
ggsave(plot = comboPlot3, filename = "Distribution plot for inlet lengths.png", width = 12, height = 6, dpi=800)



comboPlot1 = grid.arrange(mmode_lv_mapsePlot, z_lv_mapsePlot,
                          mmode_rv_tapsePlot, z_rv_tapsePlot,
                          ncol=2, widths = c(6,6))
comboPlot1

comboPlot2 = grid.arrange(mitral_earatioPlot,
                          mitral_ePlot,
                          mitral_aPlot,
                          ncol=3, widths = c(6,6,6))
comboPlot2

comboPlot3 = grid.arrange(tricuspid_earatioPlot,
                          tricuspid_ePlot,
                          tricuspid_aPlot,
                          ncol=3, widths = c(6,6,6))
comboPlot3

comboPlot4 = grid.arrange(comboPlot2,comboPlot3,
                          nrow=2)


comboPlot5 = grid.arrange(mmode_lv_fsPlot,
                          mmode_rv_fsPlot,
                          ncol=2, widths = c(6,6))
comboPlot5

comboPlot6 = grid.arrange(std_lvent_teiPlot,
                          std_rv_teiPlot,
                          ncol=2, widths = c(6,6))
comboPlot6

ggsave(plot = comboPlot1, filename = "Distribution plot for annulus displacement.png", width = 12, height = 8, dpi=800)
ggsave(plot = comboPlot4, filename = "Distribution plot for tricuspid and mitral e and a.png", width = 18, height = 6, dpi=800)
ggsave(plot = comboPlot5, filename = "Distribution plot for fractional shortening.png", width = 12, height = 6, dpi=800)
ggsave(plot = comboPlot6, filename = "Distribution plot for MPI.png", width = 12, height = 6, dpi=800)




## What is the pattern for an individual over time?
repplotfunct <- function(yvar){
  ggplot(data=file, aes(x = gaweeks, y = yvar, color=as.factor(record))) + 
    geom_line(size=1) + xlab("Gestational age (weeks)\n") +  theme_bw() +
    scale_y_continuous(labels=label_number(accuracy = 1))+
    theme(text = element_text(size=14), legend.position = "none") 
}

repplotfunct(file$z_chest_area) 
repplotfunct(file$chest_area) 
repplotfunct(file$z_heart_area) 
repplotfunct(file$heart_area) 
repplotfunct(file$cardiothoracic_ratio) 

repplotfunct(file$z_rv_tapse)
repplotfunct(file$mmode_rv_tapse)
repplotfunct(file$z_lv_mapse) 
repplotfunct(file$mmode_lv_mapse) 

repplotfunct(file$mmode_rv_fs)
repplotfunct(file$std_rv_tei)
repplotfunct(file$mmode_lv_fs) 
repplotfunct(file$std_lvent_tei) 

repplotfunct(file$tricuspid_e)
repplotfunct(file$tricuspid_a)
repplotfunct(file$tricuspid_earatio)
repplotfunct(file$mitral_e) 
repplotfunct(file$mitral_a) 
repplotfunct(file$mitral_earatio) 



### What is the pattern over time for the population?
### Run by BW status

## Read in sas dataset of predicted values
cardiac.pred.1 <- read_sas("J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Derived_data\\hpp_cardiac_pred.sas7bdat")
names(cardiac.pred.1) <- tolower(names(cardiac.pred.1))
cardiac.pred.1[cardiac.pred.1 == ""] <- NA
cardiac.pred.1$gaweeks <- round(cardiac.pred.1$ga_us)
cardiac.pred.1[1:100, c("record","gaweeks","fit_std_rv_tei")]


plot.function <- function(fit.input, lcl.input, ucl.input, ylabel) {
  
  mean.predict.out <- data.frame(
    fit = aggregate(fit.input, list(cardiac.pred.1$gaweeks, cardiac.pred.1$bwga_zcat), FUN=mean, na.rm=TRUE),
    lcl = aggregate(lcl.input, list(cardiac.pred.1$gaweeks, cardiac.pred.1$bwga_zcat), FUN=mean, na.rm=TRUE),
    ucl = aggregate(ucl.input, list(cardiac.pred.1$gaweeks, cardiac.pred.1$bwga_zcat), FUN=mean, na.rm=TRUE)
  )
  
  
  # Plot averaged predicted data
  output_plot <- ggplot(mean.predict.out, aes(x = fit.Group.1, y = fit.x, ymax = ucl.x, ymin = lcl.x, col=fit.Group.2)) +
    geom_line() +
    scale_colour_brewer(palette='Dark2', breaks = c("LGA","AGA","SGA")) +
    scale_x_continuous(breaks = seq(12, 38, by = 2), labels = seq(12, 38, by = 2)) +
    labs(
      x = "Gestational age (weeks)\n ",
      y = ylabel
    ) +
    theme_classic() +
    theme(
      text = element_text(size = 19),
      legend.position = "right",
      plot.title.position = "plot", legend.title = element_blank())
  output_plot}

z_heart_area.plot <- plot.function(cardiac.pred.1$fit_z_heart_area, cardiac.pred.1$lcl_z_heart_area, cardiac.pred.1$ucl_z_heart_area, "Heart area  (z-score)")
ggsave(filename = paste0("Plot of ", "Heart area (z-score)", ".png"), plot = z_heart_area.plot, width = 6, height = 6, dpi=800)
z_chest_area.plot <- plot.function(cardiac.pred.1$fit_z_chest_area, cardiac.pred.1$lcl_z_chest_area, cardiac.pred.1$ucl_z_chest_area, "chest area  (z-score)")
ggsave(filename = paste0("Plot of ", "Chest area (z-score)", ".png"), plot = z_chest_area.plot, width = 6, height = 6, dpi=800)
cardiothoracic_ratio.plot <- plot.function(cardiac.pred.1$fit_cardiothoracic_ratio, cardiac.pred.1$lcl_cardiothoracic_ratio, cardiac.pred.1$ucl_cardiothoracic_ratio, "Cardiothoracic Ratio")
ggsave(filename = paste0("Plot of ", "Cardiothoracic ratio", ".png"), plot = cardiothoracic_ratio.plot, width = 6, height = 6, dpi=800)

std_rv_tei.plot <- plot.function(cardiac.pred.1$fit_std_rv_tei, cardiac.pred.1$lcl_std_rv_tei, cardiac.pred.1$ucl_std_rv_tei, "RV MPI")
ggsave(filename = paste0("Plot of ", "RV MPI", ".png"), plot = std_rv_tei.plot, width = 6, height = 6, dpi=800)
std_lvent_tei.plot <- plot.function(cardiac.pred.1$fit_std_lvent_tei, cardiac.pred.1$lcl_std_lvent_tei, cardiac.pred.1$ucl_std_lvent_tei, "LV MPI")
ggsave(filename = paste0("Plot of ", "LV MPI", ".png"), plot = std_lvent_tei.plot, width = 6, height = 6, dpi=800)

z_rv_tapse.plot <- plot.function(cardiac.pred.1$fit_z_rv_tapse, cardiac.pred.1$lcl_z_rv_tapse, cardiac.pred.1$ucl_z_rv_tapse, "RV AD (z-score)")
ggsave(filename = paste0("Plot of ", "RV AD (z-score)", ".png"), plot = z_rv_tapse.plot, width = 6, height = 6, dpi=800)
z_lv_mapse.plot <- plot.function(cardiac.pred.1$fit_z_lv_mapse, cardiac.pred.1$lcl_z_lv_mapse, cardiac.pred.1$ucl_z_lv_mapse, "LV AD (z-score)")
ggsave(filename = paste0("Plot of ", "LV AD (z-score)", ".png"), plot = z_lv_mapse.plot, width = 6, height = 6, dpi=800)

mmode_rv_fs.plot <- plot.function(cardiac.pred.1$fit_mmode_rv_fs, cardiac.pred.1$lcl_mmode_rv_fs, cardiac.pred.1$ucl_mmode_rv_fs, "RV FS")
ggsave(filename = paste0("Plot of ", "RV FS", ".png"), plot = mmode_rv_fs.plot, width = 6, height = 6, dpi=800)
mmode_lv_fs.plot <- plot.function(cardiac.pred.1$fit_mmode_lv_fs, cardiac.pred.1$lcl_mmode_lv_fs, cardiac.pred.1$ucl_mmode_lv_fs, "LV FS")
ggsave(filename = paste0("Plot of ", "LV FS", ".png"), plot = mmode_lv_fs.plot, width = 6, height = 6, dpi=800)

mitral_e.plot <- plot.function(cardiac.pred.1$fit_mitral_e, cardiac.pred.1$lcl_mitral_e, cardiac.pred.1$ucl_mitral_e, "Mitral E")
ggsave(filename = paste0("Plot of ", "Mitral E", ".png"), plot = mitral_e.plot, width = 6, height = 6, dpi=800)
mitral_a.plot <- plot.function(cardiac.pred.1$fit_mitral_a, cardiac.pred.1$lcl_mitral_a, cardiac.pred.1$ucl_mitral_a, "Mitral A")
ggsave(filename = paste0("Plot of ", "Mitral A", ".png"), plot = mitral_a.plot, width = 6, height = 6, dpi=800)
mitral_earatio.plot <- plot.function(cardiac.pred.1$fit_mitral_earatio, cardiac.pred.1$lcl_mitral_earatio, cardiac.pred.1$ucl_mitral_earatio, "Mitral E/A")
ggsave(filename = paste0("Plot of ", "Mitral EA", ".png"), plot = mitral_earatio.plot, width = 6, height = 6, dpi=800)

tricuspid_e.plot <- plot.function(cardiac.pred.1$fit_tricuspid_e, cardiac.pred.1$lcl_tricuspid_e, cardiac.pred.1$ucl_tricuspid_e, "Tricuspid E")
ggsave(filename = paste0("Plot of ", "Tricuspid E", ".png"), plot = tricuspid_e.plot, width = 6, height = 6, dpi=800)
tricuspid_a.plot <- plot.function(cardiac.pred.1$fit_tricuspid_a, cardiac.pred.1$lcl_tricuspid_a, cardiac.pred.1$ucl_tricuspid_a, "Tricuspid A")
ggsave(filename = paste0("Plot of ", "Tricuspid A", ".png"), plot = tricuspid_a.plot, width = 6, height = 6, dpi=800)
tricuspid_earatio.plot <- plot.function(cardiac.pred.1$fit_tricuspid_earatio, cardiac.pred.1$lcl_tricuspid_earatio, cardiac.pred.1$ucl_tricuspid_earatio, "Tricuspid E/A")
ggsave(filename = paste0("Plot of ", "Tricuspid EA", ".png"), plot = tricuspid_earatio.plot, width = 6, height = 6, dpi=800)





###########################
## Exposure Distribution ##
###########################
med.iqr <- data.frame()
funct.exp.dist <- function(inputvars, counttype, vartype){
  dat1 <- subset(file[,inputvars],)
  med.iqr <- data.frame(t(data.frame(lapply(dat1, quantile, probs=c(0.5, 0.25, 0.75), na.rm=TRUE))))
  colnames(med.iqr) = c("Median", "Lower", "Upper")
  med.iqr$MedIQR1 <- paste0(round(med.iqr$Median,1)," (",round(med.iqr$Lower,1),", ",round(med.iqr$Upper,1),")")
  med.iqr$MedIQR2 <- paste0(round(med.iqr$Median,2)," (",round(med.iqr$Lower,2),", ",round(med.iqr$Upper,2),")")
  med.iqr$Variable <- vartype
  
  dat2 <- subset(file, complete.cases(file[ ,counttype]) & file[ ,counttype]>0)
  dat2 <- dat2[!duplicated(dat2$record),]
  med.iqr$N <- length(unique(dat2$record))
  med.iqr$count.samples <- colSums(dat2[,counttype])
  
  med.iqr
}

t1.dist <- funct.exp.dist(t1_explist, 'count_exp_tri1', "Early")
t2.dist <- funct.exp.dist(t2_explist, 'count_exp_tri2', "Mid")
t3.dist <- funct.exp.dist(t3_explist, 'count_exp_tri3', "Late")
pa.dist <- funct.exp.dist(explist, 'count_exp_preg', "Pregnancy-average")

exp.dist <- rbind(t1.dist,t2.dist,t3.dist,pa.dist)
exp.dist$Label <- rep(abbrevs.exp,4)
exp.dist.wide = exp.dist %>% pivot_wider(id_cols = c("Label"),
                                         names_from = c("Variable"),
                                         values_from = c("MedIQR2"))

write.csv(exp.dist, "Exposure Distribution Table in Analytic Sample including ns.csv")
write.csv(exp.dist.wide, "Exposure Distribution Table in Analytic Sample.csv")


## Median and IQR for number of samples? ##
dat2 <- subset(file, complete.cases(file[ ,'count_exp_preg']) & file[ ,'count_exp_preg']>0)
dat2 <- dat2[!duplicated(dat2$record),]
quantile(dat2$count_exp_preg, prob=c(0.5, 0.25, 0.75))




###########################
## Exposure correlations ##
###########################

#Correlations between exposures

# Limit to just one visit so that you don't have those with > 1 visit weighted more in correlation matrix
corfile <- file[match(unique(file$record), file$record), ]

#Correlations between exposures
cormat <- cor(corfile[,explist], use = "complete.obs", method="spearman")
colnames(cormat) <- abbrevs.exp
rownames(cormat) <- abbrevs.exp
min(cormat)
mean(cormat)

#plot
num_labels <- matrix("", nrow=nrow(cormat), ncol=ncol(cormat))
num_labels[upper.tri((cormat))] <- round(cormat[upper.tri(cormat)], digits = 2)
cormat[upper.tri(cormat)] <- NA
p <- pheatmap(cormat,
              color=colorRampPalette(c('navy', "white", "firebrick3"))(100),
              breaks = seq(-1,1, length.out = 101),
              cluster_rows=F, 
              cluster_cols=F,
              legend=T,
              legend_breaks=-1:1,
              legend_labels=c(-1,0,1),
              cellwidth=60, 
              cellheight=60,
              fontsize=25,
              display_numbers=num_labels,
              number_color="black",
              na_col=)


png(height=1800, width=1800, file="HPP3D Cardiac Exposures Correlation Matrix.png", type = "cairo")
p
dev.off()



###################################
## Impute missing covariate data ##
###################################

## Just covariates - this is for a clinical paper and we don't want to do a long'l prediction so let's stick with covariates 
covariates <- data.frame(unique(file[, c('record', 'gender_nb', 'age', 'bmi', 'any_prior_preg', 'race_new', 'highest_edu', 'smoking', 'alcohol_new', 'drugs_new', 'clinicsite', 'income','employment','insurance','marital_status')]))

## Run mice on dataset
sapply(covariates, function(x) sum(is.na(x)))
summary(covariates)

## Initialize mice
init <- mice(covariates, maxit=0, remove_collinear = FALSE, remove.constant  = FALSE)
meth <- init$method
predM = init$predictorMatrix
predM[1,] = 0
predM[,1] = 0
meth

## Run
mids.original <- mice(covariates, 
                      method=meth, 
                      predictorMatrix=predM,
                      m = 10, # Num. of imputed datasets
                      maxit = 20, # Num. of iterations per imputed dataset 
                      seed=505,
                      remove.collinear=FALSE)

## Examine logged errors
head(mids.original$loggedEvents) 
warnings()

## Check methods
mids.original$meth

## Double check the polyreg and logreg to be sure missing didn't become it's own category
table(mice::complete(mids.original)$highest_edu, useNA = "always")
table(mice::complete(mids.original)$gender_nb, useNA = "always")
#All of them look good

## Post-imputation diagnostics
my.plot <- plot(mids.original)
my.plot
# Save as pdf
myheight <- 10
mywidth <- 10
pdf("Covariate imputation trace plots.pdf", # initiate pdf file
    height = myheight, width = mywidth, 
    paper = "letter") 
my.plot
dev.off() #stop writing pdf

## Make into dataframe for editing
imp.cov = complete(mids.original, "long", include = T)


#####################################
## Merge with outcomes & exposures ##
#####################################

expout <- file[, -which(names(file) %in% c('gender_nb', 'age', 'bmi', 'any_prior_preg', 'race_new', 'highest_edu', 'smoking', 'alcohol_new', 'drugs_new', 'clinicsite', 'income','employment','insurance','marital_status'))]
expcovout <- merge(expout, imp.cov, by.x = "record")

8734/794 #Size of imputed dataset divided by size of original dataset
str(expcovout) #Correct # records'


###############################################
## Log-transform & IQR standardize exposures ##
###############################################


for(i in 1:length(explist)){ 
  ## Log-transform
  
  # Run for pregnancy-averaged 
  expcovout[explist_log[i]] = log(expcovout[explist[i]])
  # Run for T1
  expcovout[t1_explist_log[i]] = log(expcovout[t1_explist[i]])
  # Run for T2
  expcovout[t2_explist_log[i]] = log(expcovout[t2_explist[i]])
  # Run for T3
  expcovout[t3_explist_log[i]] = log(expcovout[t3_explist[i]])
  # Run for longform exposure
  expcovout[long_explist_log[i]] = log(expcovout[long_explist[i]])
  
  
  ## IQR-standardize
  
  # Determine IQR - only do for preg avg 
  q = quantile(expcovout[explist_log[i]],  quantiles = c(0.25, 0.75), na.rm=T)
  iqr = q[[4]] - q[[2]] # IQR for pregnancy-averaged exposure #
  # Run for pregnancy-averaged 
  expcovout[explist_iqr[i]] = expcovout[explist_log[i]]/iqr
  # Run for T1
  expcovout[t1_explist_iqr[i]] = expcovout[t1_explist_log[i]]/iqr
  # Run for T2
  expcovout[t2_explist_iqr[i]] = expcovout[t2_explist_log[i]]/iqr
  # Run for T3
  expcovout[t3_explist_iqr[i]] = expcovout[t3_explist_log[i]]/iqr
  # Run for longform exposure
  expcovout[long_explist_iqr[i]] = expcovout[long_explist_log[[i]]]/iqr
}

quantile(expcovout[explist_log[1]],  quantiles = c(0.25, 0.75), na.rm=T)
4.56-3.22
1.34 ## This should be IQR for MEP

## Fix the few inf and NaNs that this creates in the dataframe
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))
expcovout <- replace(expcovout, is.nan(expcovout), NA)
expcovout[expcovout==Inf] = NA
expcovout[expcovout==-Inf] = NA


## Dataset
mymids <- subset(expcovout, expcovout$.imp > 0)
str(mymids) #Correct # records
sum(is.na(mymids$highest_edu)) #NO missing - this is the correct subset

## Order
mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]

## Checks
mymids[1:20, c(".imp","record","visit","preg_mep","ln_preg_mep","iqr_preg_mep")]
log(5.5)
1.697/1.334 ## This is correct for first participant
mymids[1:50,c('.imp', 'record', 'visit', 'ln_tri1_mep', 'ln_tri2_mep', 'ln_tri3_mep', 'ln_preg_mep', 'iqr_preg_mep', 'mep', 'ln_mep','iqr_mep')]
1.4480/1.34 ## This is correct for first participant



## Save this imputed exposure, outcome, covariate dataset ##
write.csv(expcovout, "J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Derived_data\\Imputed HPP Cardiac Dataset.csv")

