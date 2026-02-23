# Programmer: Danielle Stevens
# Purpose of syntax: Tables and Figures for HPP3D analysis examining gestational biomarkers of phthalates & replacement chemicals, & fetal cardiac outcomes
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
#library(hrbrthemes)
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

outcomes.list <- list(cardiothoracic_ratio = "Cardiothoracic ratio",
                      z_heart_area = "Heart area (z-scored)",
                      z_chest_area = "Chest area (z-scored)",
                      z_lv_inlet_length = "LV inlet length (z-scored)",
                      z_rv_inlet_length = "RV inlet length (z-scored)",
                      z_lv_mapse= "LV annulus displacement (z-scored)",
                      z_rv_tapse = "RV annulus displacement (z-scored)",
                      mmode_lv_fs = "LV fractional shortening",
                      mmode_rv_fs = "RV fractional shortening",
                      std_lvent_tei = "LV MPI",
                      std_rv_tei = "RV MPI",
                      mitral_earatio = "Mitral EA",
                      mitral_e = "Mitral E",
                      mitral_a = "Mitral A",
                      tricuspid_earatio = "Tricuspid EA",
                      tricuspid_e = "Tricuspid E",
                      tricuspid_a = "Tricuspid A",
                      heart_area = "Heart area",
                      chest_area = "Chest area",
                      lv_inlet_length = "LV inlet length",
                      rv_inlet_length = "RV inlet length",
                      mmode_lv_mapse = "LV annulus displacement",
                      mmode_rv_tapse = "RV annulus displacement")



###############################
## Pregnancy-average results ##
## This is primary analysis  ##
###############################

readindata <- function(unadj.adj, outcome_label){
  mchem <- read.csv(paste0("Preg Average ", gsub(" ", "_", outcome_label), "_From_QGCOMP.csv"))
  
  adj.mchem <- mchem[grep(unadj.adj, mchem$outcome), ]
  adj.mchem <- tidyr::separate(data = adj.mchem, col = outcome, into = c("Description", "Outcome"), sep="\\,")
  
  adj.mchem$name <- ifelse(grepl("all", adj.mchem$Description), "Overall Mixture",
                    ifelse(grepl("phthalates", adj.mchem$Description), "Phthalate Mixture", 
                    ifelse(grepl("LMW", adj.mchem$Description), "LMW Mixture",
                    ifelse(grepl("HMW", adj.mchem$Description), "HMW Mixture", "Replacement Mixture"))))
  
  estdat <- t(adj.mchem[, c('name', 'estci')])
  estdat <- janitor::row_to_names(estdat, row_number=1)
  
  outname <- data.frame('n.numvis' = adj.mchem[1,'n.numvis'],
                        estdat)
  return(outname)
}


## Unadjusted models ##
merged_results <- dplyr::bind_rows(
  lapply(names(outcomes.list), function(out) {
    df <- readindata("Unadjusted", outcomes.list[[out]])
    df$outcome <- outcomes.list[[out]]  
    return(df)
  })
)
head(merged_results)   
write.csv(merged_results, "Unadjusted Results for Primary Analyses.csv")


## Adjusted models ##
merged_results <- dplyr::bind_rows(
  lapply(names(outcomes.list), function(out) {
    df <- readindata("Adjusted", outcomes.list[[out]])
    df$outcome <- outcomes.list[[out]]  
    return(df)
  })
)
head(merged_results)   
write.csv(merged_results, "Adjusted Results for Primary Analyses.csv")




#################################################
##  P-values for visit interaction models      ##
#################################################

## Read & format data function
readindata.pint <- function(unadj.adj, outcome_label){
  adj.mchem <- read.csv(paste0("Categorical Time Interaction for ", gsub(" ", "_", outcome_label), "_From_QGCOMP.csv"))
  adj.mchem <- separate(data = adj.mchem, col = X.2, into = c("Description", "Outcome"), sep="\\,")
  adj.mchem$name <- ifelse(grepl("all", adj.mchem$Description), "Overall Mixture",
                           ifelse(grepl("phthalates", adj.mchem$Description), "Phthalate Mixture", 
                                  ifelse(grepl("LMW", adj.mchem$Description), "LMW Mixture",
                                         ifelse(grepl("HMW", adj.mchem$Description), "HMW Mixture", "Replacement Mixture"))))  
  outname <- t(adj.mchem[, c('name', 'X.1')])
  outname <- row_to_names(outname, row_number=1)
  outname <- data.frame(cbind(adj.mchem$Outcome[[1]], outname))
  return(outname)
}

## Adjusted models ##
merged_results <- dplyr::bind_rows(
  lapply(names(outcomes.list), function(out) {
    df <- readindata.pint("Adjusted", outcomes.list[[out]])
    df$outcome <- outcomes.list[[out]]  
    return(df)
  })
)
head(merged_results)   
write.csv(merged_results, "Adjusted P-Int for Time Interaction.csv")





#################################################
##  P-values for Sex interaction models      ##
#################################################

## Read & format data function
readindata.pint <- function(unadj.adj, outcome_label){
  mchem <- read.csv(paste0("Sex Interaction for ", gsub(" ", "_", outcome_label), "_From_QGCOMP.csv"))
  adj.mchem <- mchem[grep(unadj.adj, mchem$outcome), ]
  adj.mchem <- separate(data = adj.mchem, col = outcome, into = c("Description", "Outcome"), sep="\\,")
  adj.mchem$name <- ifelse(grepl("all", adj.mchem$Description), "Overall Mixture",
                           ifelse(grepl("phthalates", adj.mchem$Description), "Phthalate Mixture", 
                                  ifelse(grepl("LMW", adj.mchem$Description), "LMW Mixture",
                                         ifelse(grepl("HMW", adj.mchem$Description), "HMW Mixture", "Replacement Mixture"))))  
  outname <- t(adj.mchem[, c('name', 'X.1')])
  outname <- row_to_names(outname, row_number=1)
  outname <- data.frame(cbind(adj.mchem$Outcome[[1]], outname))
  return(outname)
}

## Adjusted models ##
merged_results <- dplyr::bind_rows(
  lapply(names(outcomes.list), function(out) {
    df <- readindata.pint("Adjusted", outcomes.list[[out]])
    df$outcome <- outcomes.list[[out]]  
    return(df)
  })
)
head(merged_results)   
write.csv(merged_results, "Adjusted P-Int for Sex Interaction.csv")



############################################
### Trimester-specific associations      ###
### Limit to sign p-values/patterns only ###
############################################

outcomes.list.2 <- list(mitral_earatio = "Mitral EA",
                        mitral_e = "Mitral E",
                        mitral_a = "Mitral A",
                        tricuspid_earatio = "Tricuspid EA",
                        tricuspid_e = "Tricuspid E",
                        tricuspid_a = "Tricuspid A")


readindata <- function(trimester, outcome_label){
  mchem <- read.csv(paste0(trimester, " ", gsub(" ", "_", outcome_label), "_From_QGCOMP.csv"))
  
  adj.mchem <- mchem[grep('Adjusted', mchem$outcome), ]
  adj.mchem <- tidyr::separate(data = adj.mchem, col = outcome, into = c("Description", "Outcome"), sep="\\,")
  
  adj.mchem$name <- ifelse(grepl("phthalates", adj.mchem$Description), "Phthalate Mixture", 
                                  ifelse(grepl("LMW", adj.mchem$Description), "LMW Mixture",
                                         ifelse(grepl("HMW", adj.mchem$Description), "HMW Mixture", "Replacement Mixture")))
  
  estdat <- t(adj.mchem[, c('name', 'estci')])
  estdat <- janitor::row_to_names(estdat, row_number=1)
  
  outname <- data.frame('n.numvis' = adj.mchem[1,'n.numvis'],
                        'Trimester' = trimester,
                        estdat)
  return(outname)
}

readindata("T1", "Mitral EA")

t1_results <- dplyr::bind_rows(
  lapply(names(outcomes.list.2), function(out) {
    df <- readindata("T1", outcomes.list.2[[out]])
    df$outcome <- outcomes.list.2[[out]]  
    return(df)
  })
)
t2_results <- dplyr::bind_rows(
  lapply(names(outcomes.list.2), function(out) {
    df <- readindata("T2", outcomes.list.2[[out]])
    df$outcome <- outcomes.list.2[[out]]  
    return(df)
  })
)
t3_results <- dplyr::bind_rows(
  lapply(names(outcomes.list.2), function(out) {
    df <- readindata("T3", outcomes.list.2[[out]])
    df$outcome <- outcomes.list.2[[out]]  
    return(df)
  })
)

merged_results <- rbind(t1_results, t2_results, t3_results)
write.csv(merged_results, "Adjusted Trimester Specific Results for Secondary Analyses.csv")



################################
##  Predicted change results  ##
## This is secondary analysis ##
################################

readindata.est <- function(dataname, outcome_label){
mchem <- read.csv(dataname)
adj.mchem <- mchem[grep('Adjusted', mchem$outcome), ]

adj.mchem <- separate(data = adj.mchem, col = outcome, into = c("Description", "Outcome"), sep="\\,")
adj.mchem$name <- ifelse(grepl("all", adj.mchem$Description), "Overall Mixture",
                         ifelse(grepl("phthalates", adj.mchem$Description), "Phthalate Mixture", 
                                ifelse(grepl("LMW", adj.mchem$Description), "LMW Mixture",
                                       ifelse(grepl("HMW", adj.mchem$Description), "HMW Mixture", "Replacement Mixture"))))  
outname <- data.frame(t(adj.mchem[, c('name', 'estci')]))
outname <- row_to_names(outname, row_number=1)
outname <- data.frame(Outcome=outcome_label, 
                      outname)
outname
}


cardiothoracic_ratio.est <- readindata.est("Secondary Analysis Preg Average Slope of Cardiothoracic Ratio From QGCOMP Estimates.csv", "Cardiothoracic ratio")
z_chest_area.est <- readindata.est("Secondary Analysis Preg Average Slope of Chest Area (Z-score) From QGCOMP Estimates.csv", "Chest area (z-score)")
z_heart_area.est <- readindata.est("Secondary Analysis Preg Average Slope of Heart Area (Z-score) From QGCOMP Estimates.csv", "Heart area (z-score)")
z_lv_inlet_length.est <- readindata.est("Secondary Analysis Preg Average Slope of LV Inlet Length (Z-score) From QGCOMP Estimates.csv", "LV Inlet length (z-score)")
z_rv_inlet_length.est <- readindata.est("Secondary Analysis Preg Average Slope of RV Inlet Length (Z-score) From QGCOMP Estimates.csv", "RV Inlet length (z-score)")
z_lv_mapse.est <- readindata.est("Secondary Analysis Preg Average Slope of LV Annulus Displacement From QGCOMP Estimates.csv", "LV annulus displacement (z-score)")
z_rv_tapse.est <- readindata.est("Secondary Analysis Preg Average Slope of RV Annulus Displacement From QGCOMP Estimates.csv", "RV annulus displacement (z-score)")
mmode_lv_fs.est <- readindata.est("Secondary Analysis Preg Average Slope of LV fractional shortening From QGCOMP Estimates.csv", "LV fractional shortening")
mmode_rv_fs.est <- readindata.est("Secondary Analysis Preg Average Slope of RV Fractional Shortening From QGCOMP Estimates.csv", "RV fractional shortening")
std_lvent_tei.est <- readindata.est("Secondary Analysis Preg Average Slope of LV MPI From QGCOMP Estimates.csv", "LV MPI")
std_rv_tei.est <- readindata.est("Secondary Analysis Preg Average Slope of RV MPI From QGCOMP Estimates.csv", "RV MPI")
tricuspid_earatio.est <- readindata.est("Secondary Analysis Preg Average Slope of Tricuspid EA From QGCOMP Estimates.csv", "Tricuspid EA")
tricuspid_e.est <- readindata.est("Secondary Analysis Preg Average Slope of Tricuspid E From QGCOMP Estimates.csv", "Tricuspid E")
tricuspid_a.est <- readindata.est("Secondary Analysis Preg Average Slope of Tricuspid A From QGCOMP Estimates.csv", "Tricuspid A")
mitral_earatio.est <- readindata.est("Secondary Analysis Preg Average Slope of Mitral EA From QGCOMP Estimates.csv", "Mitral EA")
mitral_e.est <- readindata.est("Secondary Analysis Preg Average Slope of Mitral E From QGCOMP Estimates.csv", "Mitral E")
mitral_a.est <- readindata.est("Secondary Analysis Preg Average Slope of Mitral A From QGCOMP Estimates.csv", "Mitral A")

out.tab.est <- rbind(cardiothoracic_ratio.est, z_heart_area.est, z_chest_area.est, z_lv_inlet_length.est, z_rv_inlet_length.est,
                 z_lv_mapse.est, z_rv_tapse.est, 
                 mmode_lv_fs.est, mmode_rv_fs.est, 
                 std_lvent_tei.est, std_rv_tei.est,
                 mitral_earatio.est, mitral_e.est, mitral_a.est,
                 tricuspid_earatio.est, tricuspid_e.est, tricuspid_a.est)
out.tab.est$order <- 3
write.csv(out.tab.est, "Adjusted Results for Secondary Analyses for Slope Estimates.csv")



## Predictions
readindata.pred <- function(dataname, outlab){
  mchem <- read.csv(dataname)
  adj.mchem <- mchem[grep("Adjusted", mchem$label), ]
  adj.mchem <- separate(data = adj.mchem, col = label, into = c("Description", "Outcome"), sep="\\,")
  adj.mchem$name <- ifelse(grepl("all", adj.mchem$Description), "Overall Mixture",
                           ifelse(grepl("phthalates", adj.mchem$Description), "Phthalate Mixture", 
                                  ifelse(grepl("LMW", adj.mchem$Description), "LMW Mixture",
                                         ifelse(grepl("HMW", adj.mchem$Description), "HMW Mixture", "Replacement Mixture"))))  
  adj.mchem$estci <- str_c(round(adj.mchem$linpred,2),  " (", round(adj.mchem$ll.pw,2), ", ", round(adj.mchem$ul.pw,2), ")")
  first <- as.data.frame( t(subset(adj.mchem[, c('Outcome', 'name', "estci")], adj.mchem$quantile == 0)) )
  first <- row_to_names(first, row_number=2)
  first$quartile <- "First"
  
  fourth <- as.data.frame( t(subset(adj.mchem[, c('Outcome', 'name', "estci")], adj.mchem$quantile == 3)) )
  fourth <- row_to_names(fourth, row_number=2)
  fourth$quartile <- "Fourth"
  
  outname <- rbind(first, fourth)
  outname$Outcome <- outlab
  outname
}

cardiothoracic_ratio.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Cardiothoracic Ratio From QGCOMP Predictions.csv", "Cardiothoracic ratio")
z_chest_area.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Chest Area (Z-score) From QGCOMP Predictions.csv", "Chest area (z-score)")
z_heart_area.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Heart Area (Z-score) From QGCOMP Predictions.csv", "Heart area (z-score)")
z_lv_inlet_length.pred <- readindata.pred("Secondary Analysis Preg Average Slope of LV Inlet Length (Z-score) From QGCOMP Predictions.csv", "LV Inlet length (z-score)")
z_rv_inlet_length.pred <- readindata.pred("Secondary Analysis Preg Average Slope of RV Inlet Length (Z-score) From QGCOMP Predictions.csv", "RV Inlet length (z-score)")
z_lv_mapse.pred <- readindata.pred("Secondary Analysis Preg Average Slope of LV Annulus Displacement From QGCOMP Predictions.csv", "LV annulus displacement (z-score)")
z_rv_tapse.pred <- readindata.pred("Secondary Analysis Preg Average Slope of RV Annulus Displacement From QGCOMP Predictions.csv", "RV annulus displacement (z-score)")
mmode_lv_fs.pred <- readindata.pred("Secondary Analysis Preg Average Slope of LV fractional shortening From QGCOMP Predictions.csv", "LV fractional shortening")
mmode_rv_fs.pred <- readindata.pred("Secondary Analysis Preg Average Slope of RV Fractional Shortening From QGCOMP Predictions.csv", "RV fractional shortening")
std_lvent_tei.pred <- readindata.pred("Secondary Analysis Preg Average Slope of LV MPI From QGCOMP Predictions.csv", "LV MPI")
std_rv_tei.pred <- readindata.pred("Secondary Analysis Preg Average Slope of RV MPI From QGCOMP Predictions.csv", "RV MPI")
tricuspid_earatio.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Tricuspid EA From QGCOMP Predictions.csv", "Tricuspid EA")
tricuspid_e.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Tricuspid E From QGCOMP Predictions.csv", "Tricuspid E")
tricuspid_a.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Tricuspid A From QGCOMP Predictions.csv", "Tricuspid A")
mitral_earatio.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Mitral EA From QGCOMP Predictions.csv", "Mitral EA")
mitral_e.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Mitral E From QGCOMP Predictions.csv", "Mitral E")
mitral_a.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Mitral A From QGCOMP Predictions.csv", "Mitral A")

out.tab.pred <- rbind(cardiothoracic_ratio.pred, z_heart_area.pred, z_chest_area.pred, z_lv_inlet_length.pred, z_rv_inlet_length.pred,
                 z_lv_mapse.pred, z_rv_tapse.pred, 
                 mmode_lv_fs.pred, mmode_rv_fs.pred, 
                 std_lvent_tei.pred, std_rv_tei.pred,
                 mitral_earatio.pred, mitral_e.pred, mitral_a.pred,
                 tricuspid_earatio.pred, tricuspid_e.pred, tricuspid_a.pred)
out.tab.pred$order <- ifelse(out.tab.pred$quartile == "First", 1, 2)
write.csv(out.tab.pred, "Adjusted Results for Secondary Analyses for Slope Predictions.csv")


##################################
##  Complete Case analysis      ##
## This is sensitivity analysis ##
##################################

readindata <- function(unadj.adj, outcome_label){
  mchem <- read.csv(paste0("Sensitivity Analysis CC Preg Average ", gsub(" ", "_", outcome_label), "_From_QGCOMP.csv"))
  mchem$n.numvis <- str_c(round(as.numeric(mchem$n),2), ", ", round(as.numeric(mchem$num_vis),2))
  
  adj.mchem <- mchem[grep(unadj.adj, mchem$mod.description), ]
  adj.mchem <- tidyr::separate(data = adj.mchem, col = mod.description, into = c("Description", "Outcome"), sep="\\,")
  
  adj.mchem$name <- ifelse(grepl("all", adj.mchem$Description), "Overall Mixture",
                           ifelse(grepl("phthalates", adj.mchem$Description), "Phthalate Mixture", 
                                  ifelse(grepl("LMW", adj.mchem$Description), "LMW Mixture",
                                         ifelse(grepl("HMW", adj.mchem$Description), "HMW Mixture", "Replacement Mixture"))))
  
  estdat <- t(adj.mchem[, c('name', 'estci')])
  estdat <- janitor::row_to_names(estdat, row_number=1)
  
  outname <- data.frame('n.numvis' = adj.mchem[1,'n.numvis'],
                        estdat)
  return(outname)
}


merged_results <- dplyr::bind_rows(
  lapply(names(outcomes.list), function(out) {
    df <- readindata("Adjusted", outcomes.list[[out]])
    df$outcome <- outcomes.list[[out]]  
    return(df)
  })
)
head(merged_results)   
write.csv(merged_results, "Adjusted Results Sensitivity Analysis CC Analysis.csv")








#############################
## Create Plots of Results ##
#############################

## Read & format function
readindata <- function(outcome_label){
  mchem <- read.csv(paste0("Preg Average ", gsub(" ", "_", outcome_label), "_From_QGCOMP.csv"))
  adj.mchem <- mchem[grep("Adjusted", mchem$outcome), ]
  
  adj.mchem <- separate(data = adj.mchem, col = outcome, into = c("Description", "Outcome"), sep="\\,")
  adj.mchem$name <- ifelse(grepl("all", adj.mchem$Description), "Overall Mixture",
                           ifelse(grepl("phthalates", adj.mchem$Description), "Phthalate Mixture", 
                                  ifelse(grepl("LMW", adj.mchem$Description), "LMW Mixture",
                                         ifelse(grepl("HMW", adj.mchem$Description), "HMW Mixture", "Replacement Mixture"))))
  adj.mchem[, c('Outcome', 'name', 'est', 'lcl', 'ucl')]
}

merged_results <- dplyr::bind_rows(
  lapply(names(outcomes.list), function(out) {
    df <- readindata(outcomes.list[[out]])
    df$outcome <- outcomes.list[[out]]  
    return(df)
  })
)
head(merged_results)   


## Drop some mixtures
merged_results = merged_results[which(merged_results$name != "Overall Mixture"), ]
merged_results = merged_results[which(merged_results$name != "Phthalate Mixture"), ]

## Order
merged_results = merged_results %>% mutate(Outcome = factor(outcome, 
                                              levels = unlist(outcomes.list, use.names=FALSE)),
                             name = factor(name,
                                           levels = c("LMW Mixture",
                                                      "HMW Mixture",
                                                      "Replacement Mixture")))



## Separate columns
merged_results <- separate(data = merged_results, col = outcome, into = c("ventricle", "rest"), sep=" ", extra='merge')
merged_results$Ventricle <- ifelse(merged_results$ventricle == "LV", "Left Ventricle", 
                            ifelse(merged_results$ventricle == "RV", "Right Ventricle", merged_results$ventricle))


####
#### Figure 1a: Plot for CRT 
####
cardiothoracic_ratio = subset(merged_results, merged_results$Outcome == "Cardiothoracic ratio")
crt.adj.beta.plot <- ggplot(cardiothoracic_ratio, aes(x=as.numeric(est), xmin=as.numeric(lcl), xmax=as.numeric(ucl), y=name, col="black")) +
  geom_point(col="black") +
  geom_errorbar(width=0.2, position=position_dodge(width=0.5), col="black") +
  geom_vline(xintercept=(0), linetype="dashed", color="grey57") +
  xlab("Average difference (95% CI) \nfor IQR increase \nin concentrations") +  ylab(" ") +
  scale_y_discrete(labels = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture')),
                   limits = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture'))) +
  theme(text = element_text(size=19, color="black"),
        legend.position = "none")
crt.adj.beta.plot


####
#### Figure 1b: Plot for Chest and Heart Areas
####
chestheart = subset(merged_results, merged_results$Outcome == "Chest area (z-scored)" | merged_results$Outcome == "Heart area (z-scored)" )
chestheart.adj.beta.plot <- ggplot(chestheart, aes(x=as.numeric(est), xmin=as.numeric(lcl), xmax=as.numeric(ucl), y=name, shape=Outcome, col="black")) +
  geom_point(size=3, col="black", position=position_dodge(width=0.5)) +
  geom_errorbar(width=0.2, col="black", position=position_dodge(width=0.5)) +
  geom_vline(xintercept=(0), linetype="dashed", color="grey57") +
  xlab("Average difference (95% CI) \nfor IQR increase \nin concentrations") +  ylab(" ") +
  scale_y_discrete(labels = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture')),
                   limits = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture'))) +
  scale_shape_manual(values=c(18,8), name='Outcome', guide = guide_legend(title=" ", reverse = TRUE)) +
  theme(text = element_text(size=19, color="black"),
        legend.position = "right")
chestheart.adj.beta.plot






####
#### Plot for inlet length
####
il = subset(merged_results, merged_results$rest == "inlet length (z-scored)")
il$Label = "Inlet Length (z-scored)"
il.adj.beta.plot <-  ggplot(il, aes(x=as.numeric(est), xmin=as.numeric(lcl), xmax=as.numeric(ucl), y=name, col="black", shape=Ventricle)) +
  facet_wrap(~Label) + 
  geom_point(size=3, col="black", position=position_dodge(width=0.5)) +
  geom_errorbar(width=0.2, col="black", position=position_dodge(width=0.5)) +
  geom_vline(xintercept=(0), linetype="dashed", color="grey57") +
  xlab(" ") +  ylab(" ") +
  scale_y_discrete(labels = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture')),
                   limits = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture'))) +
  scale_shape_manual(values=c(17,15),name='Ventricle', guide = guide_legend(title=" ", reverse = TRUE)) +
  theme(text = element_text(size=19, color="black"),
        legend.position = "none")
il.adj.beta.plot


####
#### Plot for Annulus Displacement
####
ad = subset(merged_results, merged_results$rest == "annulus displacement")
ad$Label = "Annulus Displacement (z-score)"
ad.adj.beta.plot <-  ggplot(ad, aes(x=as.numeric(est), xmin=as.numeric(lcl), xmax=as.numeric(ucl), y=name, col="black", shape=Ventricle)) +
  facet_wrap(~Label) + 
  geom_point(size=3, col="black", position=position_dodge(width=0.5)) +
  geom_errorbar(width=0.2, col="black", position=position_dodge(width=0.5)) +
  geom_vline(xintercept=(0), linetype="dashed", color="grey57") +
  xlab(" ") +  ylab(" ") +
  scale_y_discrete(labels = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture')),
                   limits = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture'))) +
  scale_shape_manual(values=c(17,15),name='Ventricle', guide = guide_legend(title=" ", reverse = TRUE)) +
  theme(text = element_text(size=19, color="black"),
        legend.position = "none")
ad.adj.beta.plot

####
#### Plot for Fractional Shortening
####
fs = subset(merged_results, merged_results$rest == "fractional shortening")
fs$Label = "Fractional Shortening"
fs.adj.beta.plot <-  ggplot(fs, aes(x=as.numeric(est), xmin=as.numeric(lcl), xmax=as.numeric(ucl), y=name, col="black", shape=Ventricle)) +
  facet_wrap(~Label) + 
  geom_point(size=3, col="black", position=position_dodge(width=0.5)) +
  geom_errorbar(width=0.2, col="black", position=position_dodge(width=0.5)) +
  geom_vline(xintercept=(0), linetype="dashed", color="grey57") +
  xlab(" ") +  ylab(" ") +
  scale_y_discrete(labels = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture')),
                   limits = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture'))) +
  scale_shape_manual(values=c(17,15),name='Ventricle', guide = guide_legend(title=" ", reverse = TRUE)) +
  theme(text = element_text(size=19, color="black"),
        legend.position = "none")
fs.adj.beta.plot


####
#### Plot for MPI
####
mpi = subset(merged_results, merged_results$rest == "MPI")
mpi.adj.beta.plot <-  ggplot(mpi, aes(x=as.numeric(est), xmin=as.numeric(lcl), xmax=as.numeric(ucl), y=name, col="black", shape=Ventricle)) +
  facet_wrap(~rest) + 
  geom_point(size=3, col="black", position=position_dodge(width=0.5)) +
  geom_errorbar(width=0.2, col="black", position=position_dodge(width=0.5)) +
  geom_vline(xintercept=(0), linetype="dashed", color="grey57") +
  xlab(" ") +  ylab(" ") +
  scale_y_discrete(labels = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture')),
                   limits = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture'))) +
  scale_shape_manual(values=c(17,15),name='Ventricle', guide = guide_legend(title=" ", reverse = TRUE)) +
  theme(text = element_text(size=19, color="black"),
        legend.position = "none")
mpi.adj.beta.plot


####
#### Plot for E/A 
####
ea = subset(merged_results, merged_results$rest == "EA")
ea$Label = "E/A"
ea$Ventricle = ifelse(ea$ventricle == "Tricuspid", "Right Ventricle/Tricuspid", "Left Ventricle/Mitral")
ea.adj.beta.plot <-  ggplot(ea, aes(x=as.numeric(est), xmin=as.numeric(lcl), xmax=as.numeric(ucl), y=name, col="black", shape=Ventricle)) +
  facet_wrap(~Label) + 
  geom_point(size=3, col="black", position=position_dodge(width=0.5)) +
  geom_errorbar(width=0.2, col="black", position=position_dodge(width=0.5)) +
  geom_vline(xintercept=(0), linetype="dashed", color="grey57") +
  xlab(" ") +  ylab(" ") +
  scale_y_discrete(labels = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture')),
                   limits = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture'))) +
  scale_shape_manual(values=c(17,15),name='Ventricle', guide = guide_legend(title=" ", reverse = TRUE)) +
  theme(text = element_text(size=19, color="black"),
        legend.position = "right")
ea.adj.beta.plot




####
#### Primary results for CRT and its derivatives 
####
combo.plot <- plot_grid(crt.adj.beta.plot + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"), axis.title.x=element_blank()),
                        chestheart.adj.beta.plot + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), , axis.title.x=element_blank()),
                        labels = c('A','B'),
                        rel_widths = c(0.29,0.31),
                        ncol=2,
                        label_size = 20)
combo.plot <- ggdraw(add_sub(combo.plot, "Average difference (95% CI) \nfor IQR increase in concentrations", size=20,
                             vpadding=grid::unit(0.5,"lines")))
ggsave(combo.plot, filename="Primary results for CRT and its derivatives.jpg", width = 12, height = 6, units = "in", dpi=900)


####
#### Primary results for all functional measures 
####
combo.plot <- plot_grid(ad.adj.beta.plot + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm")),
                        fs.adj.beta.plot + theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()),
                        mpi.adj.beta.plot + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm"), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()),
                        ea.adj.beta.plot + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()),
                        labels = c('A','B','C','D'),
                        rel_widths = c(0.29,0.2,0.2,0.31),
                        ncol=4,
                        label_size = 20)
combo.plot <- ggdraw(add_sub(combo.plot, "Average difference (95% CI) \nfor IQR increase \nin concentrations", size=20,
                             vpadding=grid::unit(0.5,"lines")))

ggsave(combo.plot, filename="Primary results for all functional measures.jpg", width = 24, height = 6, units = "in", dpi=900)


####
#### Primary results for left and right ventricle measures 
####
####
#### Plot for MPI
####
mpi = subset(merged_results, merged_results$rest == "MPI")
mpi.adj.beta.plot <-  ggplot(mpi, aes(x=as.numeric(est), xmin=as.numeric(lcl), xmax=as.numeric(ucl), y=name, col="black", shape=Ventricle)) +
  facet_wrap(~rest) + 
  geom_point(size=3, col="black", position=position_dodge(width=0.5)) +
  geom_errorbar(width=0.2, col="black", position=position_dodge(width=0.5)) +
  geom_vline(xintercept=(0), linetype="dashed", color="grey57") +
  xlab(" ") +  ylab(" ") +
  scale_y_discrete(labels = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture')),
                   limits = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture'))) +
  scale_shape_manual(values=c(17,15),name='Ventricle', guide = guide_legend(title=" ", reverse = TRUE)) +
  theme(text = element_text(size=19, color="black"),
        legend.position = "right")
mpi.adj.beta.plot

combo.plot <- plot_grid(il.adj.beta.plot + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm")),
                        ad.adj.beta.plot + theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()),
                        fs.adj.beta.plot + theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()),
                        mpi.adj.beta.plot + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()),
                        labels = c('A','B','C','D'),
                        rel_widths = c(0.29,0.2,0.2,0.29),
                        ncol=4,
                        label_size = 20)
combo.plot <- ggdraw(add_sub(combo.plot, "Average difference (95% CI) \nfor IQR increase in concentrations", size=20,
                             vpadding=grid::unit(0.5,"lines")))

ggsave(combo.plot, filename="Primary results for left and right ventricle measures.jpg", width = 24, height = 6, units = "in", dpi=900)








####
#### Primary results for EA and its derivatives 
####


####
#### Plot for E/A 
####
ea = subset(merged_results, merged_results$rest == "EA")
ea$Label = "E/A"
ea.adj.beta.plot <-  ggplot(ea, aes(x=as.numeric(est), xmin=as.numeric(lcl), xmax=as.numeric(ucl), y=name, col="black", shape=ventricle)) +
  facet_wrap(~Label) + 
  geom_point(size=3, col="black", position=position_dodge(width=0.5)) +
  geom_errorbar(width=0.2, col="black", position=position_dodge(width=0.5)) +
  geom_vline(xintercept=(0), linetype="dashed", color="grey57") +
  xlab(" ") +  ylab(" ") +
  scale_y_discrete(labels = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture')),
                   limits = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture'))) +
  scale_shape_manual(values=c(17,15),name='ventricle', guide = guide_legend(title=" ", reverse = TRUE)) +
  theme(text = element_text(size=19, color="black"),
        legend.position = "none")
ea.adj.beta.plot


####
#### Plot for E
####
ea = subset(merged_results, merged_results$rest == "E")
ea$Label = "E"
e.adj.beta.plot <-  ggplot(ea, aes(x=as.numeric(est), xmin=as.numeric(lcl), xmax=as.numeric(ucl), y=name, col="black", shape=ventricle)) +
  facet_wrap(~Label) + 
  geom_point(size=3, col="black", position=position_dodge(width=0.5)) +
  geom_errorbar(width=0.2, col="black", position=position_dodge(width=0.5)) +
  geom_vline(xintercept=(0), linetype="dashed", color="grey57") +
  xlab(" ") +  ylab(" ") +
  scale_y_discrete(labels = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture')),
                   limits = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture'))) +
  scale_shape_manual(values=c(17,15),name='ventricle', guide = guide_legend(title=" ", reverse = TRUE)) +
  theme(text = element_text(size=19, color="black"),
        legend.position = "none")
e.adj.beta.plot


####
#### Plot for A 
####
ea = subset(merged_results, merged_results$rest == "A")
ea$Label = "A"
a.adj.beta.plot <-  ggplot(ea, aes(x=as.numeric(est), xmin=as.numeric(lcl), xmax=as.numeric(ucl), y=name, col="black", shape=ventricle)) +
  facet_wrap(~Label) + 
  geom_point(size=3, col="black", position=position_dodge(width=0.5)) +
  geom_errorbar(width=0.2, col="black", position=position_dodge(width=0.5)) +
  geom_vline(xintercept=(0), linetype="dashed", color="grey57") +
  xlab(" ") +  ylab(" ") +
  scale_y_discrete(labels = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture')),
                   limits = rev(c('LMW Mixture',
                                  'HMW Mixture',
                                  'Replacement Mixture'))) +
  scale_shape_manual(values=c(17,15),name='ventricle', guide = guide_legend(title=" ", reverse = TRUE)) +
  theme(text = element_text(size=19, color="black"),
        legend.position = "right")
a.adj.beta.plot


####
#### Primary results for EA and its derivatives 
####
combo.plot <- plot_grid(ea.adj.beta.plot + theme(plot.margin = unit(c(0.5, 1, 0.5, 0.5), "cm")),
                        e.adj.beta.plot + theme(plot.margin = unit(c(0.5, 1, 0.5, 1), "cm"), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()),
                        a.adj.beta.plot + theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"), axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()),
                        labels = c('A','B','C'),
                        rel_widths = c(0.31,0.24,0.28),
                        ncol=3,
                        label_size = 20)
combo.plot <- ggdraw(add_sub(combo.plot, "Average difference (95% CI) \nfor IQR increase in concentrations", size=20,
                             vpadding=grid::unit(0.5,"lines")))

ggsave(combo.plot, filename="Primary results for EA and its derivatives.jpg", width = 20, height = 6, units = "in", dpi=900)




################################
## Secondary analyses: Slopes ##
################################

## Predictions
readindata.pred <- function(dataname, outlab){
  mchem <- read.csv(dataname)
  adj.mchem <- mchem[grep("Adjusted", mchem$label), ]
  adj.mchem <- separate(data = adj.mchem, col = label, into = c("Description", "Outcome"), sep="\\,")
  adj.mchem$name <- ifelse(grepl("all", adj.mchem$Description), "Overall Mixture",
                           ifelse(grepl("phthalates", adj.mchem$Description), "Phthalate Mixture", 
                                  ifelse(grepl("LMW", adj.mchem$Description), "LMW Mixture",
                                         ifelse(grepl("HMW", adj.mchem$Description), "HMW Mixture", "Replacement Mixture"))))  

  outname <- adj.mchem[, c('Outcome', 'name', 'quartile', 'linpred', 'll.simul', 'ul.simul')]
  return(outname)
  }

cardiothoracic_ratio.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Cardiothoracic Ratio From QGCOMP Predictions.csv", "Cardiothoracic ratio")
z_chest_area.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Chest Area (Z-score) From QGCOMP Predictions.csv", "Chest area (z-score)")
z_heart_area.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Heart Area (Z-score) From QGCOMP Predictions.csv", "Heart area (z-score)")
z_lv_inlet_length.pred <- readindata.pred("Secondary Analysis Preg Average Slope of LV Inlet Length (Z-score) From QGCOMP Predictions.csv", "LV Inlet length (z-score)")
z_rv_inlet_length.pred <- readindata.pred("Secondary Analysis Preg Average Slope of RV Inlet Length (Z-score) From QGCOMP Predictions.csv", "LV Inlet length (z-score)")
z_lv_mapse.pred <- readindata.pred("Secondary Analysis Preg Average Slope of LV Annulus Displacement From QGCOMP Predictions.csv", "LV annulus displacement (z-score)")
z_rv_tapse.pred <- readindata.pred("Secondary Analysis Preg Average Slope of RV Annulus Displacement From QGCOMP Predictions.csv", "RV annulus displacement (z-score)")
mmode_lv_fs.pred <- readindata.pred("Secondary Analysis Preg Average Slope of LV fractional shortening From QGCOMP Predictions.csv", "LV fractional shortening")
mmode_rv_fs.pred <- readindata.pred("Secondary Analysis Preg Average Slope of RV Fractional Shortening From QGCOMP Predictions.csv", "RV fractional shortening")
std_lvent_tei.pred <- readindata.pred("Secondary Analysis Preg Average Slope of LV MPI From QGCOMP Predictions.csv", "LV MPI")
std_rv_tei.pred <- readindata.pred("Secondary Analysis Preg Average Slope of RV MPI From QGCOMP Predictions.csv", "RV MPI")
tricuspid_earatio.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Tricuspid EA From QGCOMP Predictions.csv", "Tricuspid EA")
tricuspid_e.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Tricuspid E From QGCOMP Predictions.csv", "Tricuspid E")
tricuspid_a.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Tricuspid A From QGCOMP Predictions.csv", "Tricuspid A")
mitral_earatio.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Mitral EA From QGCOMP Predictions.csv", "Mitral EA")
mitral_e.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Mitral E From QGCOMP Predictions.csv", "Mitral E")
mitral_a.pred <- readindata.pred("Secondary Analysis Preg Average Slope of Mitral A From QGCOMP Predictions.csv", "Mitral A")

out.tab.pred <- rbind(cardiothoracic_ratio.pred, z_chest_area.pred, z_heart_area.pred, z_lv_inlet_length.pred, z_rv_inlet_length.pred,
                      z_lv_mapse.pred, z_rv_tapse.pred, 
                      mmode_lv_fs.pred, mmode_rv_fs.pred, 
                      std_lvent_tei.pred, std_rv_tei.pred,
                      mitral_earatio.pred, mitral_e.pred, mitral_a.pred,
                      tricuspid_earatio.pred, tricuspid_e.pred, tricuspid_a.pred)


## Drop some mixtures
out.tab.pred = out.tab.pred[which(out.tab.pred$name != "Overall Mixture"), ]
out.tab.pred = out.tab.pred[which(out.tab.pred$name != "Phthalate Mixture"), ]


new.out.list <- unique(out.tab.pred$Outcome)

## Order
out.tab.pred = out.tab.pred %>% mutate(Outcome = factor(Outcome, 
                                                        levels = new.out.list),
                                       name = factor(name,
                                                     levels = c("LMW Mixture",
                                                                "HMW Mixture",
                                                                "Replacement Mixture")),
                                       quartile = factor(quartile,
                                                         levels = c("First",
                                                                    "Fourth")))

adj.pred.plot <- ggplot(out.tab.pred, aes(x=name, ymin=as.numeric(ll.simul), ymax=as.numeric(ul.simul), y=as.numeric(linpred ), 
                                                       group=quartile, shape=quartile, col="black")) +
  facet_wrap(~Outcome, scales="free") + 
  geom_point(position=position_dodge(width=0.5), size=3, col="black") +
  geom_errorbar(width=0.2, position=position_dodge(width=0.5), col="black") +
  labs(x="",
       y="Predicted slope (95% CI) \nsetting concentrations to\n first and fourth quartiles",
       title="") +
  scale_x_discrete(labels=label_wrap(10))+
  scale_shape_manual(values=c(1,19), guide = guide_legend(title="Quartile", reverse = TRUE)) +
  theme(text = element_text(size=19, color="black"),
        #        axis.text.x = element_text(angle=45),
        plot.title.position = "plot")
adj.pred.plot
## This was done for visualization of the different patterns that were difficult to see in the supplemental table




### Figure 1b: Do for CRT. Include a line at the median slope for the sample?
## Calculate median

expcovout <- read.csv("J:\\StevensLab\\The HPP3D Study\\HPP3D Phthalates and Cardiac_DRS\\Analysis\\Derived_data\\Imputed HPP Cardiac Dataset.csv")
clean.dat <- subset(expcovout, expcovout$.imp == 0)
clean.dat <- clean.dat %>% group_by(record) %>% filter(row_number() == 1) %>% ungroup()
median(clean.dat$v_cardiothoracic_ratio, na.rm=TRUE)
# 8.066223

## Plot
crt <- subset(out.tab.pred, out.tab.pred$Outcome == " Cardiothoracic Ratio")
crt.adj.pred.plot <- ggplot(crt, aes(x=name, ymin=as.numeric(ll.simul), ymax=as.numeric(ul.simul), y=as.numeric(linpred ), 
                                          group=quartile, shape=quartile, col="black")) +
  geom_point(position=position_dodge(width=0.5), size=3, col="black") +
  geom_errorbar(width=0.2, position=position_dodge(width=0.5), col="black") +
  geom_hline(yintercept=(8.066223), linetype="dashed", color="grey57") +
    labs(x="",
       y="Predicted slope (95% CI) \nsetting concentrations to\n first and fourth quartiles",
       title="") +
  scale_x_discrete(labels=label_wrap(10))+
  scale_shape_manual(values=c(1,19), guide = guide_legend(title="Quartile", reverse = TRUE)) +
  theme(text = element_text(size=19, color="black"),
        #        axis.text.x = element_text(angle=45),
        plot.title.position = "plot")
crt.adj.pred.plot

## Combine with 1a
combo.plot <- plot_grid(crt.adj.beta.plot + theme(plot.margin = unit(c(1, 1, 0.5, 0.5), "cm")),
                        crt.adj.pred.plot + theme(plot.margin = unit(c(0.5, 0.5, 1, 1), "cm")),
                        labels = c('A', 'B'),
                        ncol=2,
                        label_size = 20)
ggsave(combo.plot, filename="CRT Results for primary and secondary analyses of slopes.jpg", width = 13, height = 6, units = "in", dpi=900)




