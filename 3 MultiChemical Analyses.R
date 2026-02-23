# Programmer: Danielle Stevens
# Purpose of syntax: HPP analysis examining phthalates, alternatives, & fetal cardiac outcomes
# Last updated on: 12/15/2025


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



###########################
## Set working directory ##
###########################


setwd("J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results")


##################
## Define lists ##
##################

# Covariates
unadj.covlist <- c("spl1", "spl2")
covlist <- c('spl1', "spl2", 'drugs_new', 'alcohol_new', 'age', 'bmi', 'smoking','gender_nb', 'highest_edu', 'race_new', 'any_prior_preg', 'clinicsite')

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

outcomes.list <- list(chest_area = "Chest area",
                      heart_area = "Heart area",
                      z_heart_area = "Heart area (z-scored)",
                      z_chest_area = "Chest area (z-scored)",
                      cardiothoracic_ratio = "Cardiothoracic ratio",
                      rv_inlet_length = "RV inlet length",
                      lv_inlet_length = "LV inlet length",
                      z_rv_inlet_length = "RV inlet length (z-scored)",
                      z_lv_inlet_length = "LV inlet length (z-scored)",
                      mmode_rv_tapse = "RV annulus displacement",
                      mmode_lv_mapse = "LV annulus displacement",
                      z_rv_tapse = "RV annulus displacement (z-scored)",
                      z_lv_mapse= "LV annulus displacement (z-scored)",
                      mmode_rv_fs = "RV fractional shortening",
                      mmode_lv_fs = "LV fractional shortening",
                      std_rv_tei = "RV MPI",
                      std_lvent_tei = "LV MPI",
                      tricuspid_e = "Tricuspid E",
                      tricuspid_a = "Tricuspid A",
                      tricuspid_earatio = "Tricuspid EA",
                      mitral_e = "Mitral E",
                      mitral_a = "Mitral A",
                      mitral_earatio = "Mitral EA")


##########################
## Read in imputed data ##
##########################

expcovout <- read.csv("J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Derived_data\\Imputed HPP Cardiac Dataset.csv")
expcovout$clinicsite <- as.factor(expcovout$clinicsite)

############################
## Quantile g-computation ##
############################

### Setting parameters for model
samp.size <- 295 #Sample size
imp.length <- 10  #No. of imputations 

### First we want to test the model form for our qgcomp models

## Formula
qgformula <- as.formula(paste('cardiothoracic_ratio', "~", paste(c(explist_iqr, 'spl1' ,'spl2'), collapse = "+"))) #formula for model
qgformula

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0)
mymids <- subset(mymids, is.na(mymids$cardiothoracic_ratio) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Analysis - includes id for clustering over time
## We run the noboot here to get out the weights
qg.fit.imp <- list(
  analyses = lapply(1:imp.length, 
                    function(x) qgcomp.noboot(f = qgformula,
                                              expnms = explist_iqr, 
                                              family = gaussian(), 
                                              id="record",
                                              q = NULL, 
                                              bayes = TRUE, 
                                              data = mymids[which(mymids$.imp==x),])))

## Weights 
wts = as.data.frame(t(sapply(qg.fit.imp$analyses, function(x) c(-x$neg.weights, x$pos.weights)[explist_iqr])))
eachwt = do.call(c, wts)
expwts = data.frame(Plotlabel = rep(abbrevs.exp, each=nrow(wts)), Exposure = rep(explist_iqr, each=nrow(wts)), Weight=eachwt)
expwts$Plotlabel <- factor(expwts$Plotlabel, levels = abbrevs.exp)
qg.wts.plot <- ggplot(data=expwts)+ theme_classic() +
  geom_point(aes(x=Plotlabel, y=Weight)) + 
  scale_x_discrete(labels = abbrevs.exp, limits = abbrevs.exp) + 
  geom_hline(aes(yintercept=0))
qg.wts.plot

# Pool across imputations using mitools package (on my computer, pool() function won't work)
# We want to convert the estimate & cis to estimate percent difference for doubling in exposures
library(mitools)
MIcombine(qg.fit.imp$analyses)
MIcombinevals <- function(MIcombineRes,digits,outcome) {
  est <- MIcombineRes$coefficients
  se <- sqrt(diag(MIcombineRes$variance))
  x <- 1 - 0.05/2
  tStat <- MIcombineRes$coefficients/sqrt(diag(MIcombineRes$variance))
  pval <- round(2*pt(-abs(tStat),df=MIcombineRes$df),digits)
  lcl <- (MIcombineRes$coefficients) - qt(x,df=MIcombineRes$df)*(se)
  ucl <- (MIcombineRes$coefficients) + qt(x,df=MIcombineRes$df)*(se)
  ci <- str_c(round(lcl,2), ", ", round(ucl,2))
  estci <- str_c(round(est,2), " (", round(lcl,2), ", ", round(ucl,2), ")")
  label <- outcome
  results <- cbind(outcome, est, ci, estci, lcl, ucl, pval)
  results[2,]
} #Calculate the lcl and ucl based on https://bookdown.org/mwheymans/bookmi/rubins-rules.html#degrees-of-freedom-and-p-values
#Pval calculation taken from existing f'n at https://thestatsgeek.com/2020/11/05/p-values-after-multiple-imputation-using-mitools-in-r/ 
qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses),3,'outcome_mod')
qg.psi
pool(qg.fit.imp$analyses) # Works now - replicates these results perfectly so will continue to use original code




##########################
## cardiothoracic_ratio ##
##########################

### Set up function
### This function uses no bootstrapping and runs to accommodate adjustments for other chemicals (which may or may not be appropriate - see correlation matrix)

myqgcompmod.adj <- function(out, covs, exps, mod.description){
  ## Formula
  qgformula <- as.formula(paste(out, "~", paste(c(explist_iqr, covs), collapse = "+"))) #formula for model
  
  ## Analysis
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, 
                      function(x) qgcomp.noboot(f = qgformula,
                                                expnms = exps, 
                                                family = gaussian(), 
                                                id="record",
                                                q = NULL, 
                                                bayes = TRUE, 
                                                #seed = 125,
                                                data = mymids[which(mymids$.imp==x),])))
  
  qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses),3, mod.description)
  qg.psi
  
  ## Weights 
  wts = as.data.frame(t(sapply(qg.fit.imp$analyses, function(x) c(-x$neg.weights, x$pos.weights)[exps])))
  eachwt = do.call(c, wts)
  expwts = data.frame(Plotlabel = rep(exps, each=nrow(wts)), Exposure = rep(exps, each=nrow(wts)), Weight=eachwt)
  expwts$Plotlabel <- factor(expwts$Plotlabel, levels = exps)
  qg.wts.plot <- ggplot(data=expwts)+ theme_classic() +
    geom_point(aes(x=Plotlabel, y=Weight)) + 
    scale_x_discrete(labels = exps, limits = exps) + 
    geom_hline(aes(yintercept=0))
  qg.wts.plot
  
  qg.results <- list()
  qg.results[[1]] <- qg.psi
  qg.results[[2]] <- qg.wts.plot
  names(qg.results) <- c("Betas","Wts_plot")
  
  qg.results
}



## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0)
mymids <- subset(mymids, is.na(mymids$cardiothoracic_ratio) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run without bootstrapping ##
unadj.all <- myqgcompmod.adj('cardiothoracic_ratio', unadj.covlist, explist_iqr, "Unadjusted model for all, Cardiothoracic Ratio")
adj.all <- myqgcompmod.adj('cardiothoracic_ratio', covlist, explist_iqr, "Adjusted model for all, Cardiothoracic Ratio")
unadj.phth <- myqgcompmod.adj('cardiothoracic_ratio', unadj.covlist, explist_iqr_pht, "Unadjusted model for phthalates, Cardiothoracic Ratio")
adj.phth <- myqgcompmod.adj('cardiothoracic_ratio', covlist, explist_iqr_pht, "Adjusted model for phthalates, Cardiothoracic Ratio")
unadj.lmw <- myqgcompmod.adj('cardiothoracic_ratio', unadj.covlist, explist_iqr_lmw, "Unadjusted model for LMWP, Cardiothoracic Ratio")
adj.lmw <- myqgcompmod.adj('cardiothoracic_ratio', covlist, explist_iqr_lmw, "Adjusted model for LMWP, Cardiothoracic Ratio")
unadj.hmw <- myqgcompmod.adj('cardiothoracic_ratio', unadj.covlist, explist_iqr_hmw, "Unadjusted model for HMWP, Cardiothoracic Ratio")
adj.hmw <- myqgcompmod.adj('cardiothoracic_ratio', covlist, explist_iqr_hmw, "Adjusted model for HMWP, Cardiothoracic Ratio")
unadj.alt <- myqgcompmod.adj('cardiothoracic_ratio', unadj.covlist, explist_iqr_alt, "Unadjusted model for alternatives, Cardiothoracic Ratio")
adj.alt <- myqgcompmod.adj('cardiothoracic_ratio', covlist, explist_iqr_alt, "Adjusted model for alternatives, Cardiothoracic Ratio")

combine.results <- rbind(n, num_vis, unadj.all$Betas, adj.all$Betas, unadj.phth$Betas, adj.phth$Betas, unadj.lmw$Betas, adj.lmw$Betas, unadj.hmw$Betas, adj.hmw$Betas, unadj.alt$Betas, adj.alt$Betas)
write.csv(combine.results, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Preg Average Cardiothoracic Ratio From QGCOMP no boot.csv")


### This function uses bootstrapping and runs without adjustments for other chemicals (which may or may not be appropriate - see correlation matrix)
myqgcompmod <- function(out, covs, exps, mod.description){
  ## Formula
  qgformula <- as.formula(paste(out, "~", paste(c(exps, covs), collapse = "+"))) #formula for model
  
  ## Analysis
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, 
                      function(x) qgcomp.boot(f = qgformula,
                                              expnms = exps, 
                                              family = gaussian(), 
                                              id="record",
                                              q = NULL, 
                                              bayes = TRUE, 
                                              seed = 125,
                                              data = mymids[which(mymids$.imp==x),])))
  
  qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses),3, mod.description)
  qg.psi
}


## Run models
unadj.all <- myqgcompmod('cardiothoracic_ratio', unadj.covlist, explist_iqr, "Unadjusted model for all, Cardiothoracic Ratio")
adj.all <- myqgcompmod('cardiothoracic_ratio', covlist, explist_iqr, "Adjusted model for all, Cardiothoracic Ratio")
unadj.phth <- myqgcompmod('cardiothoracic_ratio', unadj.covlist, explist_iqr_pht, "Unadjusted model for phthalates, Cardiothoracic Ratio")
adj.phth <- myqgcompmod('cardiothoracic_ratio', covlist, explist_iqr_pht, "Adjusted model for phthalates, Cardiothoracic Ratio")
unadj.lmw <- myqgcompmod('cardiothoracic_ratio', unadj.covlist, explist_iqr_lmw, "Unadjusted model for LMWP, Cardiothoracic Ratio")
adj.lmw <- myqgcompmod('cardiothoracic_ratio', covlist, explist_iqr_lmw, "Adjusted model for LMWP, Cardiothoracic Ratio")
unadj.hmw <- myqgcompmod('cardiothoracic_ratio', unadj.covlist, explist_iqr_hmw, "Unadjusted model for HMWP, Cardiothoracic Ratio")
adj.hmw <- myqgcompmod('cardiothoracic_ratio', covlist, explist_iqr_hmw, "Adjusted model for HMWP, Cardiothoracic Ratio")
unadj.alt <- myqgcompmod('cardiothoracic_ratio', unadj.covlist, explist_iqr_alt, "Unadjusted model for alternatives, Cardiothoracic Ratio")
adj.alt <- myqgcompmod('cardiothoracic_ratio', covlist, explist_iqr_alt, "Adjusted model for alternatives, Cardiothoracic Ratio")

combine.results <- rbind(n, num_vis, unadj.all, adj.all, unadj.phth, adj.phth, unadj.lmw, adj.lmw, unadj.hmw, adj.hmw, unadj.alt, adj.alt)
write.csv(combine.results, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Preg Average Cardiothoracic Ratio From QGCOMP.csv")


### This function uses no bootstrapping and runs without adjustments for other chemicals (which may or may not be appropriate - see correlation matrix)
myqgcompmod <- function(out, covs, exps, mod.description){
  ## Formula
  qgformula <- as.formula(paste("`", out, "` ~ ", paste(c(exps, covs), collapse = "+"), sep=""))
  
  ## Analysis
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, 
                      function(x) qgcomp.noboot(f = qgformula,
                                                expnms = exps, 
                                                family = gaussian(), 
                                                id="record",
                                                q = NULL, 
                                                bayes = TRUE, 
                                                data = mymids[which(mymids$.imp==x),])))
  
  qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses),3, mod.description)
  qg.psi
}


## Run models
unadj.all <- myqgcompmod('cardiothoracic_ratio', unadj.covlist, explist_iqr, "Unadjusted model for all, Cardiothoracic Ratio")
adj.all <- myqgcompmod('cardiothoracic_ratio', covlist, explist_iqr, "Adjusted model for all, Cardiothoracic Ratio")
unadj.phth <- myqgcompmod('cardiothoracic_ratio', unadj.covlist, explist_iqr_pht, "Unadjusted model for phthalates, Cardiothoracic Ratio")
adj.phth <- myqgcompmod('cardiothoracic_ratio', covlist, explist_iqr_pht, "Adjusted model for phthalates, Cardiothoracic Ratio")
unadj.lmw <- myqgcompmod('cardiothoracic_ratio', unadj.covlist, explist_iqr_lmw, "Unadjusted model for LMWP, Cardiothoracic Ratio")
adj.lmw <- myqgcompmod('cardiothoracic_ratio', covlist, explist_iqr_lmw, "Adjusted model for LMWP, Cardiothoracic Ratio")
unadj.hmw <- myqgcompmod('cardiothoracic_ratio', unadj.covlist, explist_iqr_hmw, "Unadjusted model for HMWP, Cardiothoracic Ratio")
adj.hmw <- myqgcompmod('cardiothoracic_ratio', covlist, explist_iqr_hmw, "Adjusted model for HMWP, Cardiothoracic Ratio")
unadj.alt <- myqgcompmod('cardiothoracic_ratio', unadj.covlist, explist_iqr_alt, "Unadjusted model for alternatives, Cardiothoracic Ratio")
adj.alt <- myqgcompmod('cardiothoracic_ratio', covlist, explist_iqr_alt, "Adjusted model for alternatives, Cardiothoracic Ratio")

n.numvis <- str_c(n, ", ",num_vis)
combine.results <- cbind(n.numvis,
                         rbind(unadj.all, adj.all, 
                         unadj.phth, adj.phth, unadj.lmw, adj.lmw, 
                         unadj.hmw, adj.hmw, unadj.alt, adj.alt))

combine.results <- rbind(n, num_vis, unadj.all, adj.all, unadj.phth, adj.phth, unadj.lmw, adj.lmw, unadj.hmw, adj.hmw, unadj.alt, adj.alt)
write.csv(combine.results, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Preg Average Cardiothoracic Ratio From QGCOMP no boot no adjustments.csv")


## Limit & order - check models with and without co-chemical adjustments for LV MPI
mymids <- subset(expcovout, expcovout$.imp > 0)
mymids <- subset(mymids, is.na(mymids$std_lvent_tei) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]
adj.phth.unadj <- myqgcompmod('std_lvent_tei', covlist, explist_iqr_pht, "Adjusted model for phthalates, LV MPI")
adj.phth.adj <- myqgcompmod.adj('std_lvent_tei', covlist, explist_iqr_pht, "Adjusted model for phthalates, LV MPI")

## Limit & order - check models with and without co-chemical adjustments for RV MPI
mymids <- subset(expcovout, expcovout$.imp > 0)
mymids <- subset(mymids, is.na(mymids$std_rv_tei) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]
adj.phth.unadj <- myqgcompmod('std_rv_tei', covlist, explist_iqr_pht, "Adjusted model for phthalates, RV MPI")
adj.phth.adj <- myqgcompmod.adj('std_rv_tei', covlist, explist_iqr_pht, "Adjusted model for phthalates, RV MPI")



####
#### For CRT: The results for phthalates with adjustment for replacements varies wildly and is null once we take Replacements out of the model
####           This null association for phthalates makes more sense for CRT given the single-chemical effect estimates
#### There is little difference for MPI outcomes
#### Note that the SEs appear stable across models with and without co-chemical adjustments
#### Based on this and discussions with co-authors, we will run models without co-chemical adjustments
#### Run with bootstrapping, which is technically the correct way to run these models, though we note strong similarities between effect estimates between bootstrapped and non-bootstrapped qgcomp
#### 





####
#### Function for models
#### 

myqgcompmod <- function(out, covs, exps, mod.description, data){
  qgformula <- reformulate(c(exps, covs), response = out)
  
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, function(x) {
      dat <- data[data$.imp==x, ]
      qgcomp.boot(f = qgformula,
                  expnms = exps,
                  family = gaussian(),
                  id="record",
                  q = NULL,
                  bayes = TRUE,
                  seed = 125,
                  data = dat)
    })
  )
  
  qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses), 3, mod.description)
  qg.psi
}

####
#### Wrapper to run across exposures and outcomes
####

run_outcome_models <- function(outcome_name, outcome_label){
  
  ## Limit & order
  mymids <- subset(expcovout, expcovout$.imp > 0)
  
  # Outcome & exposure filtering
  mymids <- subset(mymids, !is.na(mymids[[outcome_name]]))
  mymids <- mymids[complete.cases(mymids[explist_iqr]), ]
  mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]
  
  ## Relevant info for analysis
  num_vis <- nrow(mymids)/imp.length
  n <- length(unique(mymids$record))
  n.numvis <- str_c(n, ", ",num_vis)
 
  ## Run models - Preg Average exposure
  unadj.all <- myqgcompmod(outcome_name, unadj.covlist, explist_iqr, paste("Unadjusted model for all,", outcome_label), mymids)
  adj.all   <- myqgcompmod(outcome_name, covlist, explist_iqr, paste("Adjusted model for all,", outcome_label), mymids)
  
  unadj.phth <- myqgcompmod(outcome_name, unadj.covlist, explist_iqr_pht, paste("Unadjusted model for phthalates,", outcome_label), mymids)
  adj.phth   <- myqgcompmod(outcome_name, covlist, explist_iqr_pht, paste("Adjusted model for phthalates,", outcome_label), mymids)
  
  unadj.lmw <- myqgcompmod(outcome_name, unadj.covlist, explist_iqr_lmw, paste("Unadjusted model for LMWP,", outcome_label), mymids)
  adj.lmw   <- myqgcompmod(outcome_name, covlist, explist_iqr_lmw, paste("Adjusted model for LMWP,", outcome_label), mymids)
  
  unadj.hmw <- myqgcompmod(outcome_name, unadj.covlist, explist_iqr_hmw, paste("Unadjusted model for HMWP,", outcome_label), mymids)
  adj.hmw   <- myqgcompmod(outcome_name, covlist, explist_iqr_hmw, paste("Adjusted model for HMWP,", outcome_label), mymids)
  
  unadj.alt <- myqgcompmod(outcome_name, unadj.covlist, explist_iqr_alt, paste("Unadjusted model for alternatives,", outcome_label), mymids)
  adj.alt   <- myqgcompmod(outcome_name, covlist, explist_iqr_alt, paste("Adjusted model for alternatives,", outcome_label), mymids)
  
  ## Combine results
  combine.results <- cbind(n.numvis,
                           rbind(unadj.all, adj.all, 
                                 unadj.phth, adj.phth, unadj.lmw, adj.lmw, 
                                 unadj.hmw, adj.hmw, unadj.alt, adj.alt))
  
  ## Save to CSV
  out_file <- paste0("J:/StevensLab/The HPP3D Study/HPP3D Phthalates and Cardiac_DRS/Analysis/Results/Preg Average ", 
                     gsub(" ", "_", outcome_label), "_From_QGCOMP.csv")
  write.csv(combine.results, file = out_file)
}

# Loop through outcomes in outcomes.list
for (out in names(outcomes.list)) {
  run_outcome_models(out, outcomes.list[[out]])
}





################################################
## Secondary analyses: By trimester exposures ##
################################################

########################
## Interaction Models ##
########################

expcovout$visit2 <- as.factor(ifelse(expcovout$visit == 1, 1,
                           ifelse(expcovout$visit == 4, 2,
                           ifelse(expcovout$visit == 7, 3, NA))))
summary(expcovout$visit2)

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0)
# Outcome
mymids <- subset(mymids, is.na(mymids$cardiothoracic_ratio) == FALSE) #rows = 7330
length(unique(mymids$record)) # n=293
# Exposure
mymids <- mymids[complete.cases(mymids[long_explist_iqr]), ] #6180
length(unique(mymids$record)) # n=293
mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]
# Formula
qgformula <- as.formula(paste('cardiothoracic_ratio', "~", paste(c(long_explist_iqr), collapse = "+"), "+ visit2")) #formula for model

## Analysis
qg.fit.imp <- list(
  analyses = lapply(1:imp.length, 
                    function(x) qgcomp.emm.glm.ee(f = qgformula,
                                                  expnms = long_explist_iqr,
                                                  emmvar = "visit2",
                                                  family = gaussian(), 
                                                  id="record",
                                                  q = NULL, 
                                                  bayes = TRUE, 
                                                  data = mymids[which(mymids$.imp==x),])))
# Identify interaction terms
int_terms <- grep("^mixture:visit2", 
                  names(coef(qg.fit.imp$analyses[[1]])), 
                  value = TRUE)

extract_int <- function(fit, int_terms){
  b <- coef(fit)[int_terms]
  V <- vcov(fit)[int_terms, int_terms, drop = FALSE]
  list(b = b, V = V)
}

ests <- lapply(qg.fit.imp$analyses, extract_int, int_terms = int_terms)

# Pool across imputations to get global Wald test p-value
MI_Wald_test <- function(ests){
  m <- length(ests)
  p <- length(ests[[1]]$b)
  
  # Stack estimates
  Bmat <- do.call(cbind, lapply(ests, function(x) x$b))
  Vlist <- lapply(ests, function(x) x$V)
  
  # Mean coefficient
  bbar <- rowMeans(Bmat)
  
  # Within-imputation variance
  W <- Reduce("+", Vlist) / m
  
  # Between-imputation variance
  B <- (Bmat - bbar) %*% t(Bmat - bbar) / (m - 1)
  
  # Total variance
  Tmat <- W + (1 + 1/m) * B
  
  # Wald statistic
  Wstat <- t(bbar) %*% solve(Tmat) %*% bbar
  
  # df (Barnard & Rubin)
  lambda <- (1 + 1/m) * sum(diag(B %*% solve(W)))
  df <- (m - 1) * (1 + 1/lambda)^2
  
  pval <- 1 - pchisq(Wstat, df = p)
  
  list(W = Wstat, df = df, p = pval, bbar = bbar, Tmat = Tmat)
}
global_test <- MI_Wald_test(ests)$p


#Get strata-specific effects: same visit, different imputations
getstrateffects(qg.fit.imp$analyses[[1]], emmval = 1)
getstrateffects(qg.fit.imp$analyses[[2]], emmval = 1)
getstrateffects(qg.fit.imp$analyses[[3]], emmval = 1)

#Different visits, same imputation
getstrateffects(qg.fit.imp$analyses[[3]], emmval = 1)
getstrateffects(qg.fit.imp$analyses[[3]], emmval = 2)
getstrateffects(qg.fit.imp$analyses[[3]], emmval = 3)


## Create new covariate lists
#timeint.unadj.covlist <- c('visit2')
timeint.covlist <- c('visit2', 'drugs_new', 'alcohol_new', 'age', 'bmi', 'smoking','gender_nb', 'highest_edu', 'race_new', 'any_prior_preg', "clinicsite")


## Convert to a function ##
myqgcompmod <- function(out, covs, exps, mod.description, data){
  qgformula <- reformulate(c(exps, covs), response = out)
  
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, function(x) {
      dat <- data[data$.imp==x, ]
      qgcomp.emm.glm.ee(f = qgformula,
                  expnms = exps,
                  emmvar = "visit2",
                  family = gaussian(),
                  id="record",
                  q = NULL,
                  bayes = TRUE,
                  seed = 125,
                  data = dat)
    })
  )
int_terms <- grep("^mixture:visit2", 
                    names(coef(qg.fit.imp$analyses[[1]])), 
                    value = TRUE)
ests <- lapply(qg.fit.imp$analyses, extract_int, int_terms = int_terms)

global_test <- MI_Wald_test(ests)$p

qg.psi <- c(global_test, mod.description)
qg.psi 
}

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0)

# Outcome & exposure filtering
mymids <- subset(mymids, !is.na(mymids[['mitral_earatio']]))
mymids <- mymids[complete.cases(mymids[long_explist_iqr]), ]
mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))
n.numvis <- str_c(n, ", ",num_vis)
myqgcompmod('mitral_earatio', timeint.covlist, long_explist_iqr_pht, paste("Adjusted model for all,", 'outcome_label'), mymids)


####
#### Wrapper to run across exposures and outcomes
####

run_outcome_models.long <- function(outcome_name, outcome_label){
  
  ## Limit & order
  mymids <- subset(expcovout, expcovout$.imp > 0)
  
  # Outcome & exposure filtering
  mymids <- subset(mymids, !is.na(mymids[[outcome_name]]))
  mymids <- mymids[complete.cases(mymids[long_explist_iqr]), ]
  mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]
  
  ## Relevant info for analysis
  num_vis <- nrow(mymids)/imp.length
  n <- length(unique(mymids$record))
  n.numvis <- str_c(n, ", ",num_vis)
  
  ## Run models - Preg Average exposure
#  unadj.all <- myqgcompmod(outcome_name, timeint.unadj.covlist, long_explist_iqr, paste("Unadjusted model for all,", outcome_label), mymids)
  adj.all   <- myqgcompmod(outcome_name, timeint.covlist, long_explist_iqr, paste("Adjusted model for all,", outcome_label), mymids)
  
#  unadj.phth <- myqgcompmod(outcome_name, timeint.unadj.covlist, long_explist_iqr_pht, paste("Unadjusted model for phthalates,", outcome_label), mymids)
  adj.phth   <- myqgcompmod(outcome_name, timeint.covlist, long_explist_iqr_pht, paste("Adjusted model for phthalates,", outcome_label), mymids)
  
#  unadj.lmw <- myqgcompmod(outcome_name,  timeint.unadj.covlist, long_explist_iqr_lmw, paste("Unadjusted model for LMWP,", outcome_label), mymids)
  adj.lmw   <- myqgcompmod(outcome_name, timeint.covlist, long_explist_iqr_lmw, paste("Adjusted model for LMWP,", outcome_label), mymids)
  
#  unadj.hmw <- myqgcompmod(outcome_name,  timeint.unadj.covlist, long_explist_iqr_hmw, paste("Unadjusted model for HMWP,", outcome_label), mymids)
  adj.hmw   <- myqgcompmod(outcome_name, timeint.covlist, long_explist_iqr_hmw, paste("Adjusted model for HMWP,", outcome_label), mymids)
  
#  unadj.alt <- myqgcompmod(outcome_name,  timeint.unadj.covlist, long_explist_iqr_alt, paste("Unadjusted model for alternatives,", outcome_label), mymids)
  adj.alt   <- myqgcompmod(outcome_name, timeint.covlist, long_explist_iqr_alt, paste("Adjusted model for alternatives,", outcome_label), mymids)
  
  ## Combine results
  combine.results <- cbind(n.numvis,
                           rbind(adj.all, 
                                 adj.phth, adj.lmw, 
                                 adj.hmw, adj.alt))
  
  
  ## Save to CSV
  out_file <- paste0("J:/StevensLab/The HPP3D Study/HPP3D Phthalates and Cardiac_DRS/Analysis/Results/Categorical Time Interaction for ", 
                     gsub(" ", "_", outcome_label), "_From_QGCOMP.csv")
  write.csv(combine.results, file = out_file)
}

# Loop through outcomes in outcomes.list
for (out in names(outcomes.list)) {
  run_outcome_models.long(out, outcomes.list[[out]])
}



## For only E/A based on significance of the p-values, we will now re-run stratified for exposures at t1, t2, t3
## Limit to adjusted models and LMW, HMW, Replacement chemical mixtures
myqgcompmod <- function(out, covs, exps, mod.description, data){
  qgformula <- reformulate(c(exps, covs), response = out)
  
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, function(x) {
      dat <- data[data$.imp==x, ]
      qgcomp.boot(f = qgformula,
                  expnms = exps,
                  family = gaussian(),
                  id="record",
                  q = NULL,
                  bayes = TRUE,
                  seed = 125,
                  data = dat)
    })
  )
  
  qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses), 3, mod.description)
  qg.psi
}

####
#### Wrapper to run across T1 exposures and outcomes
####

run_outcome_models <- function(outcome_name, outcome_label){
  
  ## Limit & order
  mymids <- subset(expcovout, expcovout$.imp > 0)
  
  # Outcome & exposure filtering
  mymids <- subset(mymids, !is.na(mymids[[outcome_name]]))
  mymids <- mymids[complete.cases(mymids[t1_explist_iqr]), ]
  mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]
  
  ## Relevant info for analysis
  num_vis <- nrow(mymids)/imp.length
  n <- length(unique(mymids$record))
  n.numvis <- str_c(n, ", ",num_vis)
  
  ## Run models - Preg Average exposure
  adj.phth   <- myqgcompmod(outcome_name, covlist, t1_explist_iqr_pht, paste("Adjusted model for phthalates,", outcome_label), mymids)
  adj.lmw   <- myqgcompmod(outcome_name, covlist, t1_explist_iqr_lmw, paste("Adjusted model for LMWP,", outcome_label), mymids)
  adj.hmw   <- myqgcompmod(outcome_name, covlist, t1_explist_iqr_hmw, paste("Adjusted model for HMWP,", outcome_label), mymids)
  adj.alt   <- myqgcompmod(outcome_name, covlist, t1_explist_iqr_alt, paste("Adjusted model for alternatives,", outcome_label), mymids)
  
  ## Combine results
  combine.results <- cbind(n.numvis, rbind(adj.phth, adj.lmw, adj.hmw, adj.alt))
  
  ## Save to CSV
  out_file <- paste0("J:/StevensLab/The HPP3D Study/HPP3D Phthalates and Cardiac_DRS/Analysis/Results/T1 ", 
                     gsub(" ", "_", outcome_label), "_From_QGCOMP.csv")
  write.csv(combine.results, file = out_file)
}

run_outcome_models('tricuspid_e', "Tricuspid E")
run_outcome_models('tricuspid_a', "Tricuspid A")
run_outcome_models('tricuspid_earatio', "Tricuspid EA")
run_outcome_models('mitral_e', "Mitral E")
run_outcome_models('mitral_a', "Mitral A")
run_outcome_models('mitral_earatio', "Mitral EA")





####
#### Wrapper to run across T2 exposures and outcomes
####

run_outcome_models <- function(outcome_name, outcome_label){
  
  ## Limit & order
  mymids <- subset(expcovout, expcovout$.imp > 0)
  
  # Outcome & exposure filtering
  mymids <- subset(mymids, !is.na(mymids[[outcome_name]]))
  mymids <- mymids[complete.cases(mymids[t2_explist_iqr]), ]
  mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]
  
  ## Relevant info for analysis
  num_vis <- nrow(mymids)/imp.length
  n <- length(unique(mymids$record))
  n.numvis <- str_c(n, ", ",num_vis)
  
  ## Run models - Preg Average exposure
  adj.phth   <- myqgcompmod(outcome_name, covlist, t2_explist_iqr_pht, paste("Adjusted model for phthalates,", outcome_label), mymids)
  adj.lmw   <- myqgcompmod(outcome_name, covlist, t2_explist_iqr_lmw, paste("Adjusted model for LMWP,", outcome_label), mymids)
  adj.hmw   <- myqgcompmod(outcome_name, covlist, t2_explist_iqr_hmw, paste("Adjusted model for HMWP,", outcome_label), mymids)
  adj.alt   <- myqgcompmod(outcome_name, covlist, t2_explist_iqr_alt, paste("Adjusted model for alternatives,", outcome_label), mymids)
  
  ## Combine results
  combine.results <- cbind(n.numvis, rbind(adj.phth, adj.lmw, adj.hmw, adj.alt))
  
  ## Save to CSV
  out_file <- paste0("J:/StevensLab/The HPP3D Study/HPP3D Phthalates and Cardiac_DRS/Analysis/Results/T2 ", 
                     gsub(" ", "_", outcome_label), "_From_QGCOMP.csv")
  write.csv(combine.results, file = out_file)
}

run_outcome_models('tricuspid_e', "Tricuspid E")
run_outcome_models('tricuspid_a', "Tricuspid A")
run_outcome_models('tricuspid_earatio', "Tricuspid EA")
run_outcome_models('mitral_e', "Mitral E")
run_outcome_models('mitral_a', "Mitral A")
run_outcome_models('mitral_earatio', "Mitral EA")




####
#### Wrapper to run across T3 exposures and outcomes
####

run_outcome_models <- function(outcome_name, outcome_label){
  
  ## Limit & order
  mymids <- subset(expcovout, expcovout$.imp > 0)
  
  # Outcome & exposure filtering
  mymids <- subset(mymids, !is.na(mymids[[outcome_name]]))
  mymids <- mymids[complete.cases(mymids[t3_explist_iqr]), ]
  mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]
  
  ## Relevant info for analysis
  num_vis <- nrow(mymids)/imp.length
  n <- length(unique(mymids$record))
  n.numvis <- str_c(n, ", ",num_vis)
  
  ## Run models - Preg Average exposure
  adj.phth   <- myqgcompmod(outcome_name, covlist, t3_explist_iqr_pht, paste("Adjusted model for phthalates,", outcome_label), mymids)
  adj.lmw   <- myqgcompmod(outcome_name, covlist, t3_explist_iqr_lmw, paste("Adjusted model for LMWP,", outcome_label), mymids)
  adj.hmw   <- myqgcompmod(outcome_name, covlist, t3_explist_iqr_hmw, paste("Adjusted model for HMWP,", outcome_label), mymids)
  adj.alt   <- myqgcompmod(outcome_name, covlist, t3_explist_iqr_alt, paste("Adjusted model for alternatives,", outcome_label), mymids)
  
  ## Combine results
  combine.results <- cbind(n.numvis, rbind(adj.phth, adj.lmw, adj.hmw, adj.alt))
  
  ## Save to CSV
  out_file <- paste0("J:/StevensLab/The HPP3D Study/HPP3D Phthalates and Cardiac_DRS/Analysis/Results/T3 ", 
                     gsub(" ", "_", outcome_label), "_From_QGCOMP.csv")
  write.csv(combine.results, file = out_file)
}

run_outcome_models('tricuspid_e', "Tricuspid E")
run_outcome_models('tricuspid_a', "Tricuspid A")
run_outcome_models('tricuspid_earatio', "Tricuspid EA")
run_outcome_models('mitral_e', "Mitral E")
run_outcome_models('mitral_a', "Mitral A")
run_outcome_models('mitral_earatio', "Mitral EA")








##########################################
## Secondary analyses: Sex interactions ##
##########################################

########################
## Interaction Models ##
########################

## Ensure coding
mymids$gender_nb <- as.factor(mymids$gender_nb)
table(mymids$gender_nb)

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0)
# Outcome
mymids <- subset(mymids, is.na(mymids$cardiothoracic_ratio) == FALSE) #rows = 7330
length(unique(mymids$record)) # n=293
# Exposure
mymids <- mymids[complete.cases(mymids[explist_iqr]), ] #6180
length(unique(mymids$record)) # n=293
mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]

# Formula
qgformula <- as.formula(paste('cardiothoracic_ratio', "~", paste(c(explist_iqr, unadj.covlist), collapse = "+"), "+ gender_nb")) #formula for model
## Analysis
qg.fit.imp <- list(
  analyses = lapply(1:imp.length, 
                    function(x) qgcomp.emm.glm.ee(f = qgformula,
                                                  expnms = explist_iqr,
                                                  emmvar = "gender_nb",
                                                  family = gaussian(), 
                                                  id="record",
                                                  q = NULL, 
                                                  bayes = TRUE, 
                                                  data = mymids[which(mymids$.imp==x),])))
# Pool across imputations using mitools package 
# We want to get out the interaction term p-value
MIcombine(qg.fit.imp$analyses)
MIcombine.pval <- function(MIcombineRes,digits,outcome) {
  est <- MIcombineRes$coefficients
  se <- sqrt(diag(MIcombineRes$variance))
  x <- 1 - 0.05/2
  tStat <- MIcombineRes$coefficients/sqrt(diag(MIcombineRes$variance))
  pval <- round(2*pt(-abs(tStat),df=MIcombineRes$df),digits)
  label <- outcome
  results <- cbind(outcome, pval["mixture:gender_nb"])
} 
#Pval calculation taken from existing f'n at https://thestatsgeek.com/2020/11/05/p-values-after-multiple-imputation-using-mitools-in-r/ 
qg.psi <- MIcombine.pval(MIcombine(qg.fit.imp$analyses),3,'outcome_mod')
qg.psi


## Create new covariate lists
sexint.unadj.covlist <- c('spl1', 'spl2', 'gender_nb')
covlist

## Convert to a function ##
myqgcompmod <- function(out, covs, exps, mod.description, data){
  qgformula <- reformulate(c(exps, covs), response = out)
  
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, function(x) {
      dat <- data[data$.imp==x, ]
      qgcomp.emm.glm.ee(f = qgformula,
                        expnms = exps,
                        emmvar = "gender_nb",
                        family = gaussian(),
                        id="record",
                        q = NULL,
                        bayes = TRUE,
                        seed = 125,
                        data = dat)
    })
  )
  
  qg.psi <- MIcombine.pval(MIcombine(qg.fit.imp$analyses),3,mod.description)
  qg.psi 
}


####
#### Wrapper to run across exposures and outcomes
####

run_outcome_models <- function(outcome_name, outcome_label){
  
  ## Limit & order
  mymids <- subset(expcovout, expcovout$.imp > 0)
  
  # Outcome & exposure filtering
  mymids <- subset(mymids, !is.na(mymids[[outcome_name]]))
  mymids <- mymids[complete.cases(mymids[explist_iqr]), ]
  mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]
  
  ## Relevant info for analysis
  num_vis <- nrow(mymids)/imp.length
  n <- length(unique(mymids$record))
  n.numvis <- str_c(n, ", ",num_vis)
  
  ## Run models - Preg Average exposure
  unadj.all <- myqgcompmod(outcome_name, sexint.unadj.covlist, explist_iqr, paste("Unadjusted model for all,", outcome_label), mymids)
  adj.all   <- myqgcompmod(outcome_name, covlist, explist_iqr, paste("Adjusted model for all,", outcome_label), mymids)
  
  unadj.phth <- myqgcompmod(outcome_name, sexint.unadj.covlist, explist_iqr_pht, paste("Unadjusted model for phthalates,", outcome_label), mymids)
  adj.phth   <- myqgcompmod(outcome_name, covlist, explist_iqr_pht, paste("Adjusted model for phthalates,", outcome_label), mymids)
  
  unadj.lmw <- myqgcompmod(outcome_name,  sexint.unadj.covlist, explist_iqr_lmw, paste("Unadjusted model for LMWP,", outcome_label), mymids)
  adj.lmw   <- myqgcompmod(outcome_name, covlist, explist_iqr_lmw, paste("Adjusted model for LMWP,", outcome_label), mymids)
  
  unadj.hmw <- myqgcompmod(outcome_name,  sexint.unadj.covlist, explist_iqr_hmw, paste("Unadjusted model for HMWP,", outcome_label), mymids)
  adj.hmw   <- myqgcompmod(outcome_name, covlist, explist_iqr_hmw, paste("Adjusted model for HMWP,", outcome_label), mymids)
  
  unadj.alt <- myqgcompmod(outcome_name,  sexint.unadj.covlist, explist_iqr_alt, paste("Unadjusted model for alternatives,", outcome_label), mymids)
  adj.alt   <- myqgcompmod(outcome_name, covlist, explist_iqr_alt, paste("Adjusted model for alternatives,", outcome_label), mymids)
  
  ## Combine results
  combine.results <- cbind(n.numvis,
                           rbind(unadj.all, adj.all, 
                                 unadj.phth, adj.phth, unadj.lmw, adj.lmw, 
                                 unadj.hmw, adj.hmw, unadj.alt, adj.alt))
  
  
  ## Save to CSV
  out_file <- paste0("J:/StevensLab/The HPP3D Study/HPP3D Phthalates and Cardiac_DRS/Analysis/Results/Sex Interaction for ", 
                     gsub(" ", "_", outcome_label), "_From_QGCOMP.csv")
  write.csv(combine.results, file = out_file)
}

# Loop through outcomes in outcomes.list
for (out in names(outcomes.list)) {
  run_outcome_models(out, outcomes.list[[out]])
}



## For only inlet lengths based on significance of the p-values, we will now re-run for male and female sex stratified models
## Limit to adjusted models and LMW, HMW, Replacement chemical mixtures
myqgcompmod <- function(out, covs, exps, mod.description, data){
  qgformula <- reformulate(c(exps, covs), response = out)
  
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, function(x) {
      dat <- data[data$.imp==x, ]
      qgcomp.boot(f = qgformula,
                  expnms = exps,
                  family = gaussian(),
                  id="record",
                  q = NULL,
                  bayes = TRUE,
                  seed = 125,
                  data = dat)
    })
  )
  
  qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses), 3, mod.description)
  qg.psi
}

# Covariates
covlist.sexstrat <- c('spl1', "spl2", 'drugs_new', 'alcohol_new', 'age', 'bmi', 'smoking', 'highest_edu', 'race_new', 'any_prior_preg', "clinicsite")

####
#### Wrapper to run across outcomes
####

run_outcome_models <- function(outcome_name, outcome_label, sex, sex_label){
  
  ## Limit & order
  mymids <- subset(expcovout, expcovout$.imp > 0)

  ## Outcome & exposure filtering
  mymids <- subset(mymids, !is.na(mymids[[outcome_name]]))
  mymids <- mymids[complete.cases(mymids[explist_iqr]), ]
  mymids <- mymids[order(mymids$.imp, mymids$record, mymids$ga_us),]

  ## Limit by infant sex
  mymids <- subset(mymids, mymids$gender_nb == sex)  

  ## Relevant info for analysis
  num_vis <- nrow(mymids)/imp.length  # Might not be whole number given imputation of gender_nb variable
  n <- length(unique(mymids$record))
  n.numvis <- str_c(n, ", ",num_vis)
  
  ## Run models - Preg Average exposure
  adj.phth   <- myqgcompmod(outcome_name, covlist.sexstrat, explist_iqr_pht, paste("Adjusted model for phthalates,", outcome_label), mymids)
  adj.lmw   <- myqgcompmod(outcome_name, covlist.sexstrat, explist_iqr_lmw, paste("Adjusted model for LMWP,", outcome_label), mymids)
  adj.hmw   <- myqgcompmod(outcome_name, covlist.sexstrat, explist_iqr_hmw, paste("Adjusted model for HMWP,", outcome_label), mymids)
  adj.alt   <- myqgcompmod(outcome_name, covlist.sexstrat, explist_iqr_alt, paste("Adjusted model for alternatives,", outcome_label), mymids)
  
  ## Combine results
  combine.results <- cbind(sex_label, n.numvis, rbind(adj.phth, adj.lmw, adj.hmw, adj.alt))
  
  ## Save to CSV
  out_file <- paste0("J:/StevensLab/The HPP3D Study/HPP3D Phthalates and Cardiac_DRS/Analysis/Results/Sex Stratified ", 
                     gsub(" ", "_", outcome_label), "_", sex_label, "_From_QGCOMP.csv")
  write.csv(combine.results, file = out_file)
}

run_outcome_models('z_lv_inlet_length', "LV inlet length (z-scored)", 1, "Male")
run_outcome_models('z_rv_inlet_length', "RV inlet length (z-scored)", 1, "Male")
run_outcome_models('z_lv_inlet_length', "LV inlet length (z-scored)", 2, "Female")
run_outcome_models('z_rv_inlet_length', "RV inlet length (z-scored)", 2, "Female")





####################################################################
## Secondary analyses: Estimate effects on change over time/slope ##
####################################################################

####
#### Set up function 
#### 


### TEST ###
## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_cardiothoracic_ratio) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Formula
qgformula <- as.formula(paste('v_cardiothoracic_ratio', "~", paste(c(explist_iqr, 'first_cardiothoracic_ratio'), collapse = "+"))) #formula for model

## Analysis
qg.fit.imp <- list(
  analyses = lapply(1:imp.length, 
                    function(x) qgcomp.boot(f = qgformula,
                                            expnms = explist_iqr, 
                                            family = gaussian(), 
                                            q = 4, 
                                            bayes = TRUE, 
                                            seed = 125,
                                            data = mymids[which(mymids$.imp==x),])))

## Plot results from first imputation - no good way to pool as of Feb 2025 - results are identical for the unadjusted analyses
linepred <- subset(modelbound.boot(qg.fit.imp$analyses[[1]], alpha = 0.05, pwonly = FALSE), quantile == 0 | quantile == 3)
linepred$quartile <- as.character(ifelse(linepred$quantile==0, 'First', "Fourth"))
linepred$label <- 'outcome'

linepred.plot <- ggplot(linepred, aes(x=quartile, y=linpred, ymin=ll.simul, ymax=ul.simul)) + 
  geom_point() + 
  geom_errorbar(width=0.2, position=position_dodge(width=0.5))+
  labs(x = "Quartile",
       y = "Slope of cardiothoracic ratio\n 13-30 weeks",
       title = "A)") +
  theme_classic() + theme(text = element_text(size=19), 
                          legend.position = "none",
                          plot.title.position = "plot") 
linepred.plot
## This appears to work --> scale up to function




myqgcompmod.unadj <- function(out, first.out, exps, mod.description, plotlabel){
  ## Formula
  qgformula <- as.formula(paste(out, "~", paste(c(explist_iqr, first.out), collapse = "+"))) #formula for model
  
  ## Analysis - IQR for beta estimates
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, 
                      function(x) qgcomp.boot(f = qgformula,
                                              expnms = exps, 
                                              family = gaussian(), 
                                              q = NULL, 
                                              bayes = TRUE, 
                                              seed = 125,
                                              data = mymids[which(mymids$.imp==x),])))
  
  qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses),3, mod.description)
  
  
  ## Analysis - quartiles for predicted values
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, 
                      function(x) qgcomp.boot(f = qgformula,
                                              expnms = exps, 
                                              family = gaussian(), 
                                              q = 4, 
                                              bayes = TRUE, 
                                              seed = 125,
                                              data = mymids[which(mymids$.imp==x),])))
  
  ## Plot results from first imputation - no good way to pool as of Feb 2025 
  linepred <- subset(modelbound.boot(qg.fit.imp$analyses[[1]], alpha = 0.05, pwonly = FALSE), quantile == 0 | quantile == 3)
  linepred$quartile <- as.character(ifelse(linepred$quantile==0, 'First', "Fourth"))
  linepred$label <- mod.description
  
  linepred.plot <- ggplot(linepred, aes(x=quartile, y=linpred, ymin=ll.simul, ymax=ul.simul)) + 
    geom_point() + 
    geom_errorbar(width=0.2, position=position_dodge(width=0.5))+
    labs(x = "Quartile",
         y = plotlabel,
         title = "A)") +
    theme_classic() + theme(text = element_text(size=19), 
                            legend.position = "none",
                            plot.title.position = "plot") 
  
  ### Combine results
  qg.results <- list()
  qg.results[[1]] <- qg.psi
  qg.results[[2]] <- linepred
  qg.results[[3]] <- linepred.plot
  names(qg.results) <- c("Estimates","Predictions","Predictions Plot")
  
  qg.results
}

myqgcompmod <- function(out, first.out, covs, exps, mod.description, plotlabel){
  ## Formula
  qgformula <- as.formula(paste(out, "~", paste(c(explist_iqr, first.out, covs), collapse = "+"))) #formula for model
  
  ## Analysis - IQR for beta estimates
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, 
                      function(x) qgcomp.boot(f = qgformula,
                                              expnms = exps, 
                                              family = gaussian(), 
                                              q = NULL, 
                                              bayes = TRUE, 
                                              seed = 125,
                                              data = mymids[which(mymids$.imp==x),])))
  
  qg.psi <- MIcombinevals(MIcombine(qg.fit.imp$analyses),3, mod.description)
  
  
  ## Analysis - quartiles for predicted values
  qg.fit.imp <- list(
    analyses = lapply(1:imp.length, 
                      function(x) qgcomp.boot(f = qgformula,
                                              expnms = exps, 
                                              family = gaussian(), 
                                              q = 4, 
                                              bayes = TRUE, 
                                              seed = 125,
                                              data = mymids[which(mymids$.imp==x),])))
  
  ## Plot results from first imputation - no good way to pool as of Feb 2025 
  linepred <- subset(modelbound.boot(qg.fit.imp$analyses[[1]], alpha = 0.05, pwonly = FALSE), quantile == 0 | quantile == 3)
  linepred$quartile <- as.character(ifelse(linepred$quantile==0, 'First', "Fourth"))
  linepred$label <- mod.description
  
  linepred.plot <- ggplot(linepred, aes(x=quartile, y=linpred, ymin=ll.simul, ymax=ul.simul)) + 
    geom_point() + 
    geom_errorbar(width=0.2, position=position_dodge(width=0.5))+
    labs(x = "Quartile",
         y = plotlabel,
         title = "A)") +
    theme_classic() + theme(text = element_text(size=19), 
                            legend.position = "none",
                            plot.title.position = "plot") 
  
  ### Combine results
  qg.results <- list()
  qg.results[[1]] <- qg.psi
  qg.results[[2]] <- linepred
  qg.results[[3]] <- linepred.plot
  names(qg.results) <- c("Estimates","Predictions","Predictions Plot")
  
  qg.results
}

# Covariates without spline terms
covlist.nospl1 <- c('drugs_new', 'alcohol_new', 'age', 'bmi', 'smoking','gender_nb', 'highest_edu', 'race_new', 'any_prior_preg', 'clinicsite')
  
  

##########################
## z_heart_area ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_z_heart_area) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_z_heart_area', 'first_z_heart_area', explist_iqr, "Unadjusted model for all, Heart Area", "Slope of Heart Area\n 13-30 weeks")
adj.all <- myqgcompmod('v_z_heart_area', 'first_z_heart_area', covlist.nospl1, explist_iqr, "Adjusted model for all, Heart Area", "Slope of Heart Area\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_z_heart_area', 'first_z_heart_area', explist_iqr_pht, "Unadjusted model for phthalates, Heart Area", "Slope of Heart Area\n 13-30 weeks")
adj.phth <- myqgcompmod('v_z_heart_area', 'first_z_heart_area', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, Heart Area", "Slope of Heart Area\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_z_heart_area', 'first_z_heart_area', explist_iqr_lmw, "Unadjusted model for LMWP, Heart Area", "Slope of Heart Area\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_z_heart_area', 'first_z_heart_area', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, Heart Area", "Slope of Heart Area\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_z_heart_area', 'first_z_heart_area', explist_iqr_hmw, "Unadjusted model for HMWP, Heart Area", "Slope of Heart Area\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_z_heart_area', 'first_z_heart_area', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, Heart Area", "Slope of Heart Area\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_z_heart_area', 'first_z_heart_area', explist_iqr_alt, "Unadjusted model for alternatives, Heart Area", "Slope of Heart Area\n 13-30 weeks")
adj.alt <- myqgcompmod('v_z_heart_area', 'first_z_heart_area', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, Heart Area", "Slope of Heart Area\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Heart Area (Z-score) From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Heart Area (Z-score) From QGCOMP Predictions.csv")


##########################
## z_chest_area ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_z_chest_area) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_z_chest_area', 'first_z_chest_area', explist_iqr, "Unadjusted model for all, Chest Area", "Slope of Chest Area\n 13-30 weeks")
adj.all <- myqgcompmod('v_z_chest_area', 'first_z_chest_area', covlist.nospl1, explist_iqr, "Adjusted model for all, Chest Area", "Slope of Chest Area\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_z_chest_area', 'first_z_chest_area', explist_iqr_pht, "Unadjusted model for phthalates, Chest Area", "Slope of Chest Area\n 13-30 weeks")
adj.phth <- myqgcompmod('v_z_chest_area', 'first_z_chest_area', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, Chest Area", "Slope of Chest Area\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_z_chest_area', 'first_z_chest_area', explist_iqr_lmw, "Unadjusted model for LMWP, Chest Area", "Slope of Chest Area\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_z_chest_area', 'first_z_chest_area', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, Chest Area", "Slope of Chest Area\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_z_chest_area', 'first_z_chest_area', explist_iqr_hmw, "Unadjusted model for HMWP, Chest Area", "Slope of Chest Area\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_z_chest_area', 'first_z_chest_area', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, Chest Area", "Slope of Chest Area\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_z_chest_area', 'first_z_chest_area', explist_iqr_alt, "Unadjusted model for alternatives, Chest Area", "Slope of Chest Area\n 13-30 weeks")
adj.alt <- myqgcompmod('v_z_chest_area', 'first_z_chest_area', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, Chest Area", "Slope of Chest Area\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Chest Area (Z-score) From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Chest Area (Z-score) From QGCOMP Predictions.csv")


##########################
## cardiothoracic_ratio ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_cardiothoracic_ratio) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_cardiothoracic_ratio', 'first_cardiothoracic_ratio', explist_iqr, "Unadjusted model for all, Cardiothoracic Ratio", "Slope of cardiothoracic ratio\n 13-30 weeks")
adj.all <- myqgcompmod('v_cardiothoracic_ratio', 'first_cardiothoracic_ratio', covlist.nospl1, explist_iqr, "Adjusted model for all, Cardiothoracic Ratio", "Slope of cardiothoracic ratio\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_cardiothoracic_ratio', 'first_cardiothoracic_ratio', explist_iqr_pht, "Unadjusted model for phthalates, Cardiothoracic Ratio", "Slope of cardiothoracic ratio\n 13-30 weeks")
adj.phth <- myqgcompmod('v_cardiothoracic_ratio', 'first_cardiothoracic_ratio', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, Cardiothoracic Ratio", "Slope of cardiothoracic ratio\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_cardiothoracic_ratio', 'first_cardiothoracic_ratio', explist_iqr_lmw, "Unadjusted model for LMWP, Cardiothoracic Ratio", "Slope of cardiothoracic ratio\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_cardiothoracic_ratio', 'first_cardiothoracic_ratio', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, Cardiothoracic Ratio", "Slope of cardiothoracic ratio\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_cardiothoracic_ratio', 'first_cardiothoracic_ratio', explist_iqr_hmw, "Unadjusted model for HMWP, Cardiothoracic Ratio", "Slope of cardiothoracic ratio\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_cardiothoracic_ratio', 'first_cardiothoracic_ratio', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, Cardiothoracic Ratio", "Slope of cardiothoracic ratio\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_cardiothoracic_ratio', 'first_cardiothoracic_ratio', explist_iqr_alt, "Unadjusted model for alternatives, Cardiothoracic Ratio", "Slope of cardiothoracic ratio\n 13-30 weeks")
adj.alt <- myqgcompmod('v_cardiothoracic_ratio', 'first_cardiothoracic_ratio', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, Cardiothoracic Ratio", "Slope of cardiothoracic ratio\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Cardiothoracic Ratio From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Cardiothoracic Ratio From QGCOMP Predictions.csv")



##########################
## z_lv_inlet_length ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_z_lv_inlet_length) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_z_lv_inlet_length', 'first_z_lv_inlet_length', explist_iqr, "Unadjusted model for all, LV Inlet Length", "Slope of LV Inlet Length\n 13-30 weeks")
adj.all <- myqgcompmod('v_z_lv_inlet_length', 'first_z_lv_inlet_length', covlist.nospl1, explist_iqr, "Adjusted model for all, LV Inlet Length", "Slope of LV Inlet Length\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_z_lv_inlet_length', 'first_z_lv_inlet_length', explist_iqr_pht, "Unadjusted model for phthalates, LV Inlet Length", "Slope of LV Inlet Length\n 13-30 weeks")
adj.phth <- myqgcompmod('v_z_lv_inlet_length', 'first_z_lv_inlet_length', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, LV Inlet Length", "Slope of LV Inlet Length\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_z_lv_inlet_length', 'first_z_lv_inlet_length', explist_iqr_lmw, "Unadjusted model for LMWP, LV Inlet Length", "Slope of LV Inlet Length\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_z_lv_inlet_length', 'first_z_lv_inlet_length', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, LV Inlet Length", "Slope of LV Inlet Length\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_z_lv_inlet_length', 'first_z_lv_inlet_length', explist_iqr_hmw, "Unadjusted model for HMWP, LV Inlet Length", "Slope of LV Inlet Length\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_z_lv_inlet_length', 'first_z_lv_inlet_length', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, LV Inlet Length", "Slope of LV Inlet Length\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_z_lv_inlet_length', 'first_z_lv_inlet_length', explist_iqr_alt, "Unadjusted model for alternatives, LV Inlet Length", "Slope of LV Inlet Length\n 13-30 weeks")
adj.alt <- myqgcompmod('v_z_lv_inlet_length', 'first_z_lv_inlet_length', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, LV Inlet Length", "Slope of LV Inlet Length\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of LV Inlet Length (Z-score) From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of LV Inlet Length (Z-score) From QGCOMP Predictions.csv")



##########################
## z_rv_inlet_length ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_z_rv_inlet_length) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_z_rv_inlet_length', 'first_z_rv_inlet_length', explist_iqr, "Unadjusted model for all, RV Inlet Length", "Slope of RV Inlet Length\n 13-30 weeks")
adj.all <- myqgcompmod('v_z_rv_inlet_length', 'first_z_rv_inlet_length', covlist.nospl1, explist_iqr, "Adjusted model for all, RV Inlet Length", "Slope of RV Inlet Length\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_z_rv_inlet_length', 'first_z_rv_inlet_length', explist_iqr_pht, "Unadjusted model for phthalates, RV Inlet Length", "Slope of RV Inlet Length\n 13-30 weeks")
adj.phth <- myqgcompmod('v_z_rv_inlet_length', 'first_z_rv_inlet_length', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, RV Inlet Length", "Slope of RV Inlet Length\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_z_rv_inlet_length', 'first_z_rv_inlet_length', explist_iqr_lmw, "Unadjusted model for LMWP, RV Inlet Length", "Slope of RV Inlet Length\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_z_rv_inlet_length', 'first_z_rv_inlet_length', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, RV Inlet Length", "Slope of RV Inlet Length\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_z_rv_inlet_length', 'first_z_rv_inlet_length', explist_iqr_hmw, "Unadjusted model for HMWP, RV Inlet Length", "Slope of RV Inlet Length\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_z_rv_inlet_length', 'first_z_rv_inlet_length', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, RV Inlet Length", "Slope of RV Inlet Length\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_z_rv_inlet_length', 'first_z_rv_inlet_length', explist_iqr_alt, "Unadjusted model for alternatives, RV Inlet Length", "Slope of RV Inlet Length\n 13-30 weeks")
adj.alt <- myqgcompmod('v_z_rv_inlet_length', 'first_z_rv_inlet_length', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, RV Inlet Length", "Slope of RV Inlet Length\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of RV Inlet Length (Z-score) From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of RV Inlet Length (Z-score) From QGCOMP Predictions.csv")



##########################
## mmode_rv_fs ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_mmode_rv_fs) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_mmode_rv_fs', 'first_mmode_rv_fs', explist_iqr, "Unadjusted model for all, RV Fractional Shortening", "Slope of RV Fractional Shortening\n 13-30 weeks")
adj.all <- myqgcompmod('v_mmode_rv_fs', 'first_mmode_rv_fs', covlist.nospl1, explist_iqr, "Adjusted model for all, RV Fractional Shortening", "Slope of RV Fractional Shortening\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_mmode_rv_fs', 'first_mmode_rv_fs', explist_iqr_pht, "Unadjusted model for phthalates, RV Fractional Shortening", "Slope of RV Fractional Shortening\n 13-30 weeks")
adj.phth <- myqgcompmod('v_mmode_rv_fs', 'first_mmode_rv_fs', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, RV Fractional Shortening", "Slope of RV Fractional Shortening\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_mmode_rv_fs', 'first_mmode_rv_fs', explist_iqr_lmw, "Unadjusted model for LMWP, RV Fractional Shortening", "Slope of RV Fractional Shortening\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_mmode_rv_fs', 'first_mmode_rv_fs', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, RV Fractional Shortening", "Slope of RV Fractional Shortening\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_mmode_rv_fs', 'first_mmode_rv_fs', explist_iqr_hmw, "Unadjusted model for HMWP, RV Fractional Shortening", "Slope of RV Fractional Shortening\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_mmode_rv_fs', 'first_mmode_rv_fs', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, RV Fractional Shortening", "Slope of RV Fractional Shortening\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_mmode_rv_fs', 'first_mmode_rv_fs', explist_iqr_alt, "Unadjusted model for alternatives, RV Fractional Shortening", "Slope of RV Fractional Shortening\n 13-30 weeks")
adj.alt <- myqgcompmod('v_mmode_rv_fs', 'first_mmode_rv_fs', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, RV Fractional Shortening", "Slope of RV Fractional Shortening\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of RV Fractional Shortening From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of RV Fractional Shortening From QGCOMP Predictions.csv")



##########################
## mmode_lv_fs ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_mmode_lv_fs) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))


## Run models
unadj.all <- myqgcompmod.unadj('v_mmode_lv_fs', 'first_mmode_lv_fs', explist_iqr, "Unadjusted model for all, LV Fractional Shortening", "Slope of LV Fractional Shortening\n 13-30 weeks")
adj.all <- myqgcompmod('v_mmode_lv_fs', 'first_mmode_lv_fs', covlist.nospl1, explist_iqr, "Adjusted model for all, LV Fractional Shortening", "Slope of LV Fractional Shortening\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_mmode_lv_fs', 'first_mmode_lv_fs', explist_iqr_pht, "Unadjusted model for phthalates, LV Fractional Shortening", "Slope of LV Fractional Shortening\n 13-30 weeks")
adj.phth <- myqgcompmod('v_mmode_lv_fs', 'first_mmode_lv_fs', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, LV Fractional Shortening", "Slope of LV Fractional Shortening\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_mmode_lv_fs', 'first_mmode_lv_fs', explist_iqr_lmw, "Unadjusted model for LMWP, LV Fractional Shortening", "Slope of LV Fractional Shortening\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_mmode_lv_fs', 'first_mmode_lv_fs', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, LV Fractional Shortening", "Slope of LV Fractional Shortening\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_mmode_lv_fs', 'first_mmode_lv_fs', explist_iqr_hmw, "Unadjusted model for HMWP, LV Fractional Shortening", "Slope of LV Fractional Shortening\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_mmode_lv_fs', 'first_mmode_lv_fs', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, LV Fractional Shortening", "Slope of LV Fractional Shortening\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_mmode_lv_fs', 'first_mmode_lv_fs', explist_iqr_alt, "Unadjusted model for alternatives, LV Fractional Shortening", "Slope of LV Fractional Shortening\n 13-30 weeks")
adj.alt <- myqgcompmod('v_mmode_lv_fs', 'first_mmode_lv_fs', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, LV Fractional Shortening", "Slope of LV Fractional Shortening\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of LV Fractional Shortening From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of LV Fractional Shortening From QGCOMP Predictions.csv")



##########################
## std_rv_tei ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_std_rv_tei) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))


## Run models
unadj.all <- myqgcompmod.unadj('v_std_rv_tei', 'first_std_rv_tei', explist_iqr, "Unadjusted model for all, RV MPI", "Slope of RV MPI\n 13-30 weeks")
adj.all <- myqgcompmod('v_std_rv_tei', 'first_std_rv_tei', covlist.nospl1, explist_iqr, "Adjusted model for all, RV MPI", "Slope of RV MPI\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_std_rv_tei', 'first_std_rv_tei', explist_iqr_pht, "Unadjusted model for phthalates, RV MPI", "Slope of RV MPI\n 13-30 weeks")
adj.phth <- myqgcompmod('v_std_rv_tei', 'first_std_rv_tei', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, RV MPI", "Slope of RV MPI\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_std_rv_tei', 'first_std_rv_tei', explist_iqr_lmw, "Unadjusted model for LMWP, RV MPI", "Slope of RV MPI\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_std_rv_tei', 'first_std_rv_tei', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, RV MPI", "Slope of RV MPI\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_std_rv_tei', 'first_std_rv_tei', explist_iqr_hmw, "Unadjusted model for HMWP, RV MPI", "Slope of RV MPI\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_std_rv_tei', 'first_std_rv_tei', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, RV MPI", "Slope of RV MPI\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_std_rv_tei', 'first_std_rv_tei', explist_iqr_alt, "Unadjusted model for alternatives, RV MPI", "Slope of RV MPI\n 13-30 weeks")
adj.alt <- myqgcompmod('v_std_rv_tei', 'first_std_rv_tei', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, RV MPI", "Slope of RV MPI\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of RV MPI From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of RV MPI From QGCOMP Predictions.csv")



##########################
## std_lvent_tei ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_std_lvent_tei) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_std_lvent_tei', 'first_std_lvent_tei', explist_iqr, "Unadjusted model for all, LV MPI", "Slope of LV MPI\n 13-30 weeks")
adj.all <- myqgcompmod('v_std_lvent_tei', 'first_std_lvent_tei', covlist.nospl1, explist_iqr, "Adjusted model for all, LV MPI", "Slope of LV MPI\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_std_lvent_tei', 'first_std_lvent_tei', explist_iqr_pht, "Unadjusted model for phthalates, LV MPI", "Slope of LV MPI\n 13-30 weeks")
adj.phth <- myqgcompmod('v_std_lvent_tei', 'first_std_lvent_tei', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, LV MPI", "Slope of LV MPI\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_std_lvent_tei', 'first_std_lvent_tei', explist_iqr_lmw, "Unadjusted model for LMWP, LV MPI", "Slope of LV MPI\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_std_lvent_tei', 'first_std_lvent_tei', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, LV MPI", "Slope of LV MPI\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_std_lvent_tei', 'first_std_lvent_tei', explist_iqr_hmw, "Unadjusted model for HMWP, LV MPI", "Slope of LV MPI\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_std_lvent_tei', 'first_std_lvent_tei', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, LV MPI", "Slope of LV MPI\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_std_lvent_tei', 'first_std_lvent_tei', explist_iqr_alt, "Unadjusted model for alternatives, LV MPI", "Slope of LV MPI\n 13-30 weeks")
adj.alt <- myqgcompmod('v_std_lvent_tei', 'first_std_lvent_tei', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, LV MPI", "Slope of LV MPI\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of LV MPI From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of LV MPI From QGCOMP Predictions.csv")



##########################
## tricuspid_a ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_tricuspid_a) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_tricuspid_a', 'first_tricuspid_a', explist_iqr, "Unadjusted model for all, Tricuspid A", "Slope of Tricuspid A\n 13-30 weeks")
adj.all <- myqgcompmod('v_tricuspid_a', 'first_tricuspid_a', covlist.nospl1, explist_iqr, "Adjusted model for all, Tricuspid A", "Slope of Tricuspid A\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_tricuspid_a', 'first_tricuspid_a', explist_iqr_pht, "Unadjusted model for phthalates, Tricuspid A", "Slope of Tricuspid A\n 13-30 weeks")
adj.phth <- myqgcompmod('v_tricuspid_a', 'first_tricuspid_a', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, Tricuspid A", "Slope of Tricuspid A\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_tricuspid_a', 'first_tricuspid_a', explist_iqr_lmw, "Unadjusted model for LMWP, Tricuspid A", "Slope of Tricuspid A\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_tricuspid_a', 'first_tricuspid_a', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, Tricuspid A", "Slope of Tricuspid A\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_tricuspid_a', 'first_tricuspid_a', explist_iqr_hmw, "Unadjusted model for HMWP, Tricuspid A", "Slope of Tricuspid A\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_tricuspid_a', 'first_tricuspid_a', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, Tricuspid A", "Slope of Tricuspid A\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_tricuspid_a', 'first_tricuspid_a', explist_iqr_alt, "Unadjusted model for alternatives, Tricuspid A", "Slope of Tricuspid A\n 13-30 weeks")
adj.alt <- myqgcompmod('v_tricuspid_a', 'first_tricuspid_a', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, Tricuspid A", "Slope of Tricuspid A\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Tricuspid A From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Tricuspid A From QGCOMP Predictions.csv")



##########################
## tricuspid_e ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_tricuspid_e) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_tricuspid_e', 'first_tricuspid_e', explist_iqr, "Unadjusted model for all, Tricuspid E", "Slope of Tricuspid E\n 13-30 weeks")
adj.all <- myqgcompmod('v_tricuspid_e', 'first_tricuspid_e', covlist.nospl1, explist_iqr, "Adjusted model for all, Tricuspid E", "Slope of Tricuspid E\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_tricuspid_e', 'first_tricuspid_e', explist_iqr_pht, "Unadjusted model for phthalates, Tricuspid E", "Slope of Tricuspid E\n 13-30 weeks")
adj.phth <- myqgcompmod('v_tricuspid_e', 'first_tricuspid_e', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, Tricuspid E", "Slope of Tricuspid E\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_tricuspid_e', 'first_tricuspid_e', explist_iqr_lmw, "Unadjusted model for LMWP, Tricuspid E", "Slope of Tricuspid E\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_tricuspid_e', 'first_tricuspid_e', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, Tricuspid E", "Slope of Tricuspid E\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_tricuspid_e', 'first_tricuspid_e', explist_iqr_hmw, "Unadjusted model for HMWP, Tricuspid E", "Slope of Tricuspid E\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_tricuspid_e', 'first_tricuspid_e', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, Tricuspid E", "Slope of Tricuspid E\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_tricuspid_e', 'first_tricuspid_e', explist_iqr_alt, "Unadjusted model for alternatives, Tricuspid E", "Slope of Tricuspid E\n 13-30 weeks")
adj.alt <- myqgcompmod('v_tricuspid_e', 'first_tricuspid_e', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, Tricuspid E", "Slope of Tricuspid E\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Tricuspid E From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Tricuspid E From QGCOMP Predictions.csv")


##########################
## tricuspid_earatio ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_tricuspid_earatio) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_tricuspid_earatio', 'first_tricuspid_earatio', explist_iqr, "Unadjusted model for all, Tricuspid EA", "Slope of Tricuspid E/A\n 13-30 weeks")
adj.all <- myqgcompmod('v_tricuspid_earatio', 'first_tricuspid_earatio', covlist.nospl1, explist_iqr, "Adjusted model for all, Tricuspid EA", "Slope of Tricuspid E/A\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_tricuspid_earatio', 'first_tricuspid_earatio', explist_iqr_pht, "Unadjusted model for phthalates, Tricuspid EA", "Slope of Tricuspid E/A\n 13-30 weeks")
adj.phth <- myqgcompmod('v_tricuspid_earatio', 'first_tricuspid_earatio', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, Tricuspid EA", "Slope of Tricuspid E/A\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_tricuspid_earatio', 'first_tricuspid_earatio', explist_iqr_lmw, "Unadjusted model for LMWP, Tricuspid EA", "Slope of Tricuspid E/A\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_tricuspid_earatio', 'first_tricuspid_earatio', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, Tricuspid EA", "Slope of Tricuspid E/A\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_tricuspid_earatio', 'first_tricuspid_earatio', explist_iqr_hmw, "Unadjusted model for HMWP, Tricuspid EA", "Slope of Tricuspid E/A\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_tricuspid_earatio', 'first_tricuspid_earatio', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, Tricuspid EA", "Slope of Tricuspid E/A\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_tricuspid_earatio', 'first_tricuspid_earatio', explist_iqr_alt, "Unadjusted model for alternatives, Tricuspid EA", "Slope of Tricuspid E/A\n 13-30 weeks")
adj.alt <- myqgcompmod('v_tricuspid_earatio', 'first_tricuspid_earatio', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, Tricuspid EA", "Slope of Tricuspid E/A\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Tricuspid EA From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Tricuspid EA From QGCOMP Predictions.csv")



##########################
## mitral_a ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_mitral_a) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_mitral_a', 'first_mitral_a', explist_iqr, "Unadjusted model for all, Mitral A", "Slope of Mitral A\n 13-30 weeks")
adj.all <- myqgcompmod('v_mitral_a', 'first_mitral_a', covlist.nospl1, explist_iqr, "Adjusted model for all, Mitral A", "Slope of Mitral A\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_mitral_a', 'first_mitral_a', explist_iqr_pht, "Unadjusted model for phthalates, Mitral A", "Slope of Mitral A\n 13-30 weeks")
adj.phth <- myqgcompmod('v_mitral_a', 'first_mitral_a', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, Mitral A", "Slope of Mitral A\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_mitral_a', 'first_mitral_a', explist_iqr_lmw, "Unadjusted model for LMWP, Mitral A", "Slope of Mitral A\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_mitral_a', 'first_mitral_a', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, Mitral A", "Slope of Mitral A\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_mitral_a', 'first_mitral_a', explist_iqr_hmw, "Unadjusted model for HMWP, Mitral A", "Slope of Mitral A\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_mitral_a', 'first_mitral_a', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, Mitral A", "Slope of Mitral A\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_mitral_a', 'first_mitral_a', explist_iqr_alt, "Unadjusted model for alternatives, Mitral A", "Slope of Mitral A\n 13-30 weeks")
adj.alt <- myqgcompmod('v_mitral_a', 'first_mitral_a', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, Mitral A", "Slope of Mitral A\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Mitral A From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Mitral A From QGCOMP Predictions.csv")



##########################
## mitral_e ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_mitral_e) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_mitral_e', 'first_mitral_e', explist_iqr, "Unadjusted model for all, Mitral E", "Slope of Mitral E\n 13-30 weeks")
adj.all <- myqgcompmod('v_mitral_e', 'first_mitral_e', covlist.nospl1, explist_iqr, "Adjusted model for all, Mitral E", "Slope of Mitral E\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_mitral_e', 'first_mitral_e', explist_iqr_pht, "Unadjusted model for phthalates, Mitral E", "Slope of Mitral E\n 13-30 weeks")
adj.phth <- myqgcompmod('v_mitral_e', 'first_mitral_e', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, Mitral E", "Slope of Mitral E\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_mitral_e', 'first_mitral_e', explist_iqr_lmw, "Unadjusted model for LMWP, Mitral E", "Slope of Mitral E\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_mitral_e', 'first_mitral_e', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, Mitral E", "Slope of Mitral E\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_mitral_e', 'first_mitral_e', explist_iqr_hmw, "Unadjusted model for HMWP, Mitral E", "Slope of Mitral E\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_mitral_e', 'first_mitral_e', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, Mitral E", "Slope of Mitral E\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_mitral_e', 'first_mitral_e', explist_iqr_alt, "Unadjusted model for alternatives, Mitral E", "Slope of Mitral E\n 13-30 weeks")
adj.alt <- myqgcompmod('v_mitral_e', 'first_mitral_e', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, Mitral E", "Slope of Mitral E\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Mitral E From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Mitral E From QGCOMP Predictions.csv")


##########################
## mitral_earatio ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_mitral_earatio) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_mitral_earatio', 'first_mitral_earatio', explist_iqr, "Unadjusted model for all, Mitral EA", "Slope of Mitral E/A\n 13-30 weeks")
adj.all <- myqgcompmod('v_mitral_earatio', 'first_mitral_earatio', covlist.nospl1, explist_iqr, "Adjusted model for all, Mitral EA", "Slope of Mitral E/A\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_mitral_earatio', 'first_mitral_earatio', explist_iqr_pht, "Unadjusted model for phthalates, Mitral EA", "Slope of Mitral E/A\n 13-30 weeks")
adj.phth <- myqgcompmod('v_mitral_earatio', 'first_mitral_earatio', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, Mitral EA", "Slope of Mitral E/A\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_mitral_earatio', 'first_mitral_earatio', explist_iqr_lmw, "Unadjusted model for LMWP, Mitral EA", "Slope of Mitral E/A\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_mitral_earatio', 'first_mitral_earatio', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, Mitral EA", "Slope of Mitral E/A\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_mitral_earatio', 'first_mitral_earatio', explist_iqr_hmw, "Unadjusted model for HMWP, Mitral EA", "Slope of Mitral E/A\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_mitral_earatio', 'first_mitral_earatio', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, Mitral EA", "Slope of Mitral E/A\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_mitral_earatio', 'first_mitral_earatio', explist_iqr_alt, "Unadjusted model for alternatives, Mitral EA", "Slope of Mitral E/A\n 13-30 weeks")
adj.alt <- myqgcompmod('v_mitral_earatio', 'first_mitral_earatio', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, Mitral EA", "Slope of Mitral E/A\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Mitral EA From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of Mitral EA From QGCOMP Predictions.csv")


##########################
## z_lv_mapse ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_z_lv_mapse) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_z_lv_mapse', 'first_z_lv_mapse', explist_iqr, "Unadjusted model for all, LV Annulus Displacement", "Slope of LV Annulus Displacement\n 13-30 weeks")
adj.all <- myqgcompmod('v_z_lv_mapse', 'first_z_lv_mapse', covlist.nospl1, explist_iqr, "Adjusted model for all, LV Annulus Displacement", "Slope of LV Annulus Displacement\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_z_lv_mapse', 'first_z_lv_mapse', explist_iqr_pht, "Unadjusted model for phthalates, LV Annulus Displacement", "Slope of LV Annulus Displacement\n 13-30 weeks")
adj.phth <- myqgcompmod('v_z_lv_mapse', 'first_z_lv_mapse', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, LV Annulus Displacement", "Slope of LV Annulus Displacement\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_z_lv_mapse', 'first_z_lv_mapse', explist_iqr_lmw, "Unadjusted model for LMWP, LV Annulus Displacement", "Slope of LV Annulus Displacement\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_z_lv_mapse', 'first_z_lv_mapse', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, LV Annulus Displacement", "Slope of LV Annulus Displacement\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_z_lv_mapse', 'first_z_lv_mapse', explist_iqr_hmw, "Unadjusted model for HMWP, LV Annulus Displacement", "Slope of LV Annulus Displacement\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_z_lv_mapse', 'first_z_lv_mapse', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, LV Annulus Displacement", "Slope of LV Annulus Displacement\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_z_lv_mapse', 'first_z_lv_mapse', explist_iqr_alt, "Unadjusted model for alternatives, LV Annulus Displacement", "Slope of LV Annulus Displacement\n 13-30 weeks")
adj.alt <- myqgcompmod('v_z_lv_mapse', 'first_z_lv_mapse', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, LV Annulus Displacement", "Slope of LV Annulus Displacement\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of LV Annulus Displacement From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of LV Annulus Displacement From QGCOMP Predictions.csv")


##########################
## z_rv_tapse ##
##########################

## Limit & order
mymids <- subset(expcovout, expcovout$.imp > 0 & is.na(expcovout$v_z_rv_tapse) == FALSE)
mymids <- mymids[order(mymids$.imp, mymids$record),]
mymids <- mymids %>% group_by(.imp, record) %>% filter(row_number() == 1) %>% ungroup()

## Relevant info for analysis
num_vis <- nrow(mymids)/imp.length
n <- length(unique(mymids$record))

## Run models
unadj.all <- myqgcompmod.unadj('v_z_rv_tapse', 'first_z_rv_tapse', explist_iqr, "Unadjusted model for all, RV Annulus Displacement", "Slope of RV Annulus Displacement\n 13-30 weeks")
adj.all <- myqgcompmod('v_z_rv_tapse', 'first_z_rv_tapse', covlist.nospl1, explist_iqr, "Adjusted model for all, RV Annulus Displacement", "Slope of RV Annulus Displacement\n 13-30 weeks")
unadj.phth <- myqgcompmod.unadj('v_z_rv_tapse', 'first_z_rv_tapse', explist_iqr_pht, "Unadjusted model for phthalates, RV Annulus Displacement", "Slope of RV Annulus Displacement\n 13-30 weeks")
adj.phth <- myqgcompmod('v_z_rv_tapse', 'first_z_rv_tapse', covlist.nospl1, explist_iqr_pht, "Adjusted model for phthalates, RV Annulus Displacement", "Slope of RV Annulus Displacement\n 13-30 weeks")
unadj.lmw <- myqgcompmod.unadj('v_z_rv_tapse', 'first_z_rv_tapse', explist_iqr_lmw, "Unadjusted model for LMWP, RV Annulus Displacement", "Slope of RV Annulus Displacement\n 13-30 weeks")
adj.lmw <- myqgcompmod('v_z_rv_tapse', 'first_z_rv_tapse', covlist.nospl1, explist_iqr_lmw, "Adjusted model for LMWP, RV Annulus Displacement", "Slope of RV Annulus Displacement\n 13-30 weeks")
unadj.hmw <- myqgcompmod.unadj('v_z_rv_tapse', 'first_z_rv_tapse', explist_iqr_hmw, "Unadjusted model for HMWP, RV Annulus Displacement", "Slope of RV Annulus Displacement\n 13-30 weeks")
adj.hmw <- myqgcompmod('v_z_rv_tapse', 'first_z_rv_tapse', covlist.nospl1, explist_iqr_hmw, "Adjusted model for HMWP, RV Annulus Displacement", "Slope of RV Annulus Displacement\n 13-30 weeks")
unadj.alt <- myqgcompmod.unadj('v_z_rv_tapse', 'first_z_rv_tapse', explist_iqr_alt, "Unadjusted model for alternatives, RV Annulus Displacement", "Slope of RV Annulus Displacement\n 13-30 weeks")
adj.alt <- myqgcompmod('v_z_rv_tapse', 'first_z_rv_tapse', covlist.nospl1, explist_iqr_alt, "Adjusted model for alternatives, RV Annulus Displacement", "Slope of RV Annulus Displacement\n 13-30 weeks")

combine.results.1 <- rbind(n, num_vis, unadj.all$Estimates, adj.all$Estimates, unadj.phth$Estimates, adj.phth$Estimates, unadj.lmw$Estimates, adj.lmw$Estimates, unadj.hmw$Estimates, adj.hmw$Estimates, unadj.alt$Estimates, adj.alt$Estimates)
combine.results.2 <- rbind(unadj.all$Predictions, adj.all$Predictions, unadj.phth$Predictions, adj.phth$Predictions, unadj.lmw$Predictions, adj.lmw$Predictions, unadj.hmw$Predictions, adj.hmw$Predictions, unadj.alt$Predictions, adj.alt$Predictions)
write.csv(combine.results.1, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of RV Annulus Displacement From QGCOMP Estimates.csv")
write.csv(combine.results.2, file="J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Results\\Secondary Analysis Preg Average Slope of RV Annulus Displacement From QGCOMP Predictions.csv")





############################
## Sensitivity analysis   ##
## Complete case analysis ##
############################


####
#### Function for models
#### 


myqgcompmod <- function(out, covs, exps, mod.description){
  
  ## Limit & order
  clean.dat <- subset(expcovout[,c('record', exps, covs, out)], expcovout$.imp == 0)
  clean.dat <- clean.dat[complete.cases(clean.dat), ]
  
  ## Relevant info for analysis
  num_vis <- nrow(clean.dat)
  n <- length(unique(clean.dat$record))
  
  ## Formula
  qgformula <- reformulate(c(exps, covs), response = out)
  
  ## Model
  qg.fit.imp <- qgcomp.boot(f = qgformula,
                            expnms = exps, 
                            family = gaussian(), 
                            id="record",
                            q = NULL, 
                            bayes = TRUE, 
                            seed = 125,
                            data = clean.dat)
  est<-qg.fit.imp$psi
  lcl<-qg.fit.imp$ci[1]
  ucl<-qg.fit.imp$ci[2]
  pval<-qg.fit.imp$pval[2]
  estci <- str_c(round(est,2), " (", round(lcl,2), ", ", round(ucl,2), ")")
  results <- cbind(mod.description, n, num_vis, est, estci, pval)
  results
}

####
#### Wrapper to run across exposures and outcomes
####

run_outcome_models <- function(outcome_name, outcome_label){

  ## Run models - Preg Average exposure
  unadj.all <- myqgcompmod(outcome_name, unadj.covlist, explist_iqr, paste("Unadjusted model for all,", outcome_label))
  adj.all   <- myqgcompmod(outcome_name, covlist, explist_iqr, paste("Adjusted model for all,", outcome_label))
  
  unadj.phth <- myqgcompmod(outcome_name, unadj.covlist, explist_iqr_pht, paste("Unadjusted model for phthalates,", outcome_label))
  adj.phth   <- myqgcompmod(outcome_name, covlist, explist_iqr_pht, paste("Adjusted model for phthalates,", outcome_label))
  
  unadj.lmw <- myqgcompmod(outcome_name, unadj.covlist, explist_iqr_lmw, paste("Unadjusted model for LMWP,", outcome_label))
  adj.lmw   <- myqgcompmod(outcome_name, covlist, explist_iqr_lmw, paste("Adjusted model for LMWP,", outcome_label))
  
  unadj.hmw <- myqgcompmod(outcome_name, unadj.covlist, explist_iqr_hmw, paste("Unadjusted model for HMWP,", outcome_label))
  adj.hmw   <- myqgcompmod(outcome_name, covlist, explist_iqr_hmw, paste("Adjusted model for HMWP,", outcome_label))
  
  unadj.alt <- myqgcompmod(outcome_name, unadj.covlist, explist_iqr_alt, paste("Unadjusted model for alternatives,", outcome_label))
  adj.alt   <- myqgcompmod(outcome_name, covlist, explist_iqr_alt, paste("Adjusted model for alternatives,", outcome_label))
  
  ## Combine results
  combine.results <- rbind(unadj.all, adj.all, unadj.phth, adj.phth, unadj.lmw, adj.lmw, unadj.hmw, adj.hmw, unadj.alt, adj.alt)
  
  ## Save to CSV
  out_file <- paste0("J:/StevensLab/The HPP3D Study/HPP3D Phthalates and Cardiac_DRS/Analysis/Results/Sensitivity Analysis CC Preg Average ", 
                     gsub(" ", "_", outcome_label), "_From_QGCOMP.csv")
  write.csv(combine.results, file = out_file)
}

# Loop through outcomes in outcomes.list
for (out in names(outcomes.list)) {
  run_outcome_models(out, outcomes.list[[out]])
}

