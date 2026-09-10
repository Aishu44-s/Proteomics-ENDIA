# Author: Helena Oakey
# Date: Feb 23, 2024

#============
# LIBRARIES
#============


##################
#install.packages("R2admb")
#install.packages("glmmADMB", 
#    repos=c("http://glmmadmb.r-forge.r-project.org/repos",
#            getOption("repos")),
#    type="source")

library(glmmADMB)
#https://stats.oarc.ucla.edu/r/faq/random-coefficient-poisson-models/

####################################################
setwd("C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/NEW ANALYSIS/Helena_code")

# Confirm catch.R is now found:
file.exists("catch.R")               # must return TRUE
file.exists("LiandJiBH_function.R")  # must return TRUE

####################################################

library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
#library(openxlsx)
library(car)# for Anova function

source(file = 'catch.R')


#==========
# Models
#==========

#---------------


#NOTE: response changed from log(btmean) to btmean.
#btmean is copied directly from a VSN-normalised column, which is already on
#a generalised-log2 scale - wrapping it in log() again was a double transform
#(compressing an already-log-scale value further), not a crash risk (VSN output
#here is strictly positive, ~5.9-23.2) but a distortion of effect-size
#interpretability and possibly residual normality. Removed 2024- see chat log /
#methods note for the residual-diagnostic comparison that motivated this.
model = list(
  
  
  #UNADJUSTED 1
  #nest and mother id matching factors no maternal_T1D interaction
  form_unadj1_rm="btmean~casecat*corrected_visit_sample+(1|riskset)+(1|mother_id)",#1
  form_unadj1_m="btmean~casecat*corrected_visit_sample+(1|mother_id)",#2
  form_unadj1_r="btmean~casecat* corrected_visit_sample+(1|riskset)",#3
  
  
  #UNADJUSTED 2
  #add other matching factors
  #----
  #glm adjusted random mother random risk set + matching
  form_unadj2_rm="btmean~casecat*corrected_visit_sample+infant_sex+timefromonsetmonths+(1|riskset)+(1|mother_id)",#4
  form_unadj2_m="btmean~casecat*corrected_visit_sample+infant_sex+timefromonsetmonths+(1|mother_id)",#5
  form_unadj2_r="btmean~casecat*corrected_visit_sample+infant_sex+timefromonsetmonths+(1|riskset)",#6
  
  
  #ADJUSTED model-----
  #add matching factors+confounders 
  form_adj_rm="btmean~casecat* corrected_visit_sample+maternal_T1D+infant_sex+timefromonsetmonths+ maternal_HLA+maternal_age_at_conception+maternal_conception_bmi+maternal_smoking+(1|riskset)+(1|mother_id)",#7
  form_adj_m="btmean~casecat* corrected_visit_sample+maternal_T1D+infant_sex+timefromonsetmonths+ maternal_HLA+maternal_age_at_conception+maternal_conception_bmi+maternal_smoking+(1|mother_id)",#8
  form_adj_r="btmean~casecat* corrected_visit_sample+maternal_T1D+infant_sex+timefromonsetmonths+ maternal_HLA+maternal_age_at_conception+maternal_conception_bmi+maternal_smoking+(1|riskset)",#9
  
  
  
  #UNADJUSTED 1
  #nest and mother id matching factors maternal_T1D interaction
  #----
  form_unadj1_int_rm="btmean~casecat*corrected_visit_sample*maternal_T1D+(1|riskset)+(1|mother_id)",#10
  form_unadj1_int_m="btmean~casecat*corrected_visit_sample*maternal_T1D+(1|mother_id)", #11
  form_unadj1_int_r="btmean~casecat*corrected_visit_sample*maternal_T1D+(1|riskset)",#12
  
  #UNADJUSTED 2
  #add other matching factors maternal_T1D interaction
  #----
  #glm adjusted random mother random risk set + matching
  form_unadj2_int_rm="btmean~casecat*corrected_visit_sample*maternal_T1D+infant_sex+timefromonsetmonths+(1|riskset)+(1|mother_id)",#13
  form_unadj2_int_m="btmean~casecat*corrected_visit_sample*maternal_T1D+infant_sex+timefromonsetmonths+(1|mother_id)",#14
  form_unadj2_int_r="btmean~casecat*corrected_visit_sample*maternal_T1D+infant_sex+timefromonsetmonths+(1|riskset)",#15
  
  
  #ADJUSTED model-----
  #add matching factors+confounders +maternal_T1D interaction
  form_adj_int_rm="btmean~casecat* corrected_visit_sample*maternal_T1D+infant_sex+timefromonsetmonths+ maternal_HLA+maternal_age_at_conception+maternal_conception_bmi+maternal_smoking+(1|riskset)+(1|mother_id)",#16
  form_adj_int_m="btmean~casecat* corrected_visit_sample*maternal_T1D+infant_sex+timefromonsetmonths+ maternal_HLA+maternal_age_at_conception+maternal_conception_bmi+maternal_smoking+(1|mother_id)",#17
  form_adj_int_r="btmean~casecat* corrected_visit_sample*maternal_T1D+infant_sex+timefromonsetmonths+ maternal_HLA+maternal_age_at_conception+maternal_conception_bmi+maternal_smoking+(1|riskset)"#18
  
  
  
)


#==============
#SET UP
#==============
# Create empty lists to store the parameter, se, jackstrap se, convergence and warnings 
# latter 4 are not needed for the cl model


names(mother)
#response/outcome variables
resp


#https://www.ssc.wisc.edu/sscc/pubs/MM/MM_TestEffects.html
#comparing fixed using lrt fit models with REML=FALSE

#ANALYSIS
#create a list with two levels to collect model summaries
analysis.ls<-replicate(length(resp),rep(list(NULL),length(model)),simplify=F)
class(analysis.ls)
#names for the first level
names(analysis.ls)<-resp
analysis.ls

#names for second level
analysis.ls<-setNames(lapply(seq_along(analysis.ls), function(x) { 
  setNames(analysis.ls[[x]],names(model)) })
  ,names(analysis.ls))
analysis.ls

names(analysis.ls)

#WARNINGS 1
#create a list with two levels to collect warnings
warn1.ls<-replicate(length(resp),rep(list(NULL),length(model)),simplify=F)
class(warn1.ls)
#names for the first level
names(warn1.ls)<-resp
warn1.ls

#names for second level
warn1.ls<-setNames(lapply(seq_along(warn1.ls), function(x) { 
  setNames(warn1.ls[[x]],names(model)) })
  ,names(warn1.ls))
warn1.ls

#CONVERGENCE
#create a list with two levels to collect warnings
conv.ls<-replicate(length(resp),rep(list(NULL),length(model)),simplify=F)
class(conv.ls)
#names for the first level
names(conv.ls)<-resp
conv.ls

#names for second level
conv.ls<-setNames(lapply(seq_along(conv.ls), function(x) { 
  setNames(conv.ls[[x]],names(model)) })
  ,names(conv.ls))
conv.ls


#LOG LIKELIHOOD
loglik.ls<-replicate(length(resp),rep(list(NULL),length(model)),simplify=F)
class(loglik.ls)
#names for the first level
names(loglik.ls)<-resp
loglik.ls

#names for second level
loglik.ls<-setNames(lapply(seq_along(loglik.ls), function(x) { 
  setNames(loglik.ls[[x]],names(model)) })
  ,names(loglik.ls))
loglik.ls

#ANOVA
#create a list with two levels to collect anova
anova.ls<-replicate(length(resp),rep(list(NULL),length(model)),simplify=F)
class(anova.ls)
#names for the first level
names(anova.ls)<-resp
anova.ls

#names for second level
anova.ls<-setNames(lapply(seq_along(anova.ls), function(x) { 
  setNames(anova.ls[[x]],names(model)) })
  ,names(anova.ls))
anova.ls



#AIC/BIC
#create a list with two levels to collect IC
AIC.ls<-replicate(length(resp),rep(list(NULL),length(model)),simplify=F)
class(AIC.ls)
#names for the first level
names(AIC.ls)<-resp
AIC.ls

#names for second level
AIC.ls<-setNames(lapply(seq_along(AIC.ls), function(x) { 
  setNames(AIC.ls[[x]],names(model)) })
  ,names(AIC.ls))
AIC.ls



#PARAMETER ESTIMATES FIXED
parameter_estimates_f.ls<-replicate(length(resp),rep(list(NULL),length(model)),simplify=F)
class(parameter_estimates_f.ls)
#names for the first level
names(parameter_estimates_f.ls)<-resp
parameter_estimates_f.ls

#names for second level
parameter_estimates_f.ls<-setNames(lapply(seq_along(parameter_estimates_f.ls), function(x) { 
  setNames(parameter_estimates_f.ls[[x]],names(model)) })
  ,names(parameter_estimates_f.ls))
parameter_estimates_f.ls



#PARAMETER ESTIMATES RANDOM
parameter_estimates_r.ls<-replicate(length(resp),rep(list(NULL),length(model)),simplify=F)
class(parameter_estimates_r.ls)
#names for the first level
names(parameter_estimates_r.ls)<-resp
parameter_estimates_r.ls

#names for second level
parameter_estimates_r.ls<-setNames(lapply(seq_along(parameter_estimates_r.ls), function(x) { 
  setNames(parameter_estimates_r.ls[[x]],names(model)) })
  ,names(parameter_estimates_r.ls))
parameter_estimates_r.ls


#============================
# Modelling
#============================

#https://www.jtrive.com/posts/spline-smoothing-r/
#https://stats.stackexchange.com/questions/546829/how-to-interpret-lm-coefficients-when-using-bs-splines-in-r

# mother_id: Identifier for each individual
# casecat: Group membership 2 levels ( control vs. case INFANT IA)
# caset1dcat: Group membership 4 levels ( t1d vs non-t1d MOTHER |control vs. case INFANT IA)
# corrected_visit_sample: Time points
# resp: Outcome variable of interest
#------------
#INITIAL
#------------
i=1
j=1
j=4
j=10
j=13
j=16


form<-NULL
warn1<-NULL
mod<-NULL
conv<-NULL

#for each response
for( i in 1:length(resp)){
  cat("i=",i,"\n")
  print(resp[i])
  #for each model
  for( j in 1:length(names(model))){
    cat("i=",i,"\n")
    cat("j=",j,"\n")
    rm(form,warn1,mod,conv)
    
    #name model
    k<-names(model)[j]
    print(k)
    
    # Get data
    tmp<-mother
    w<-which(names(tmp)==resp[i]);w
    names(tmp)[w]<-"btmean"
    tmp<-tmp[!is.na(tmp$btmean),]
    
    
    #Model formula
    form<-as.formula(model[[j]])
    
    #Fit model
    warn1<-ch(
      mod<-lmer(formula=form,weight=weightTEDDY,data=tmp)
    )
    conv<-mod@optinfo$conv
    
    #Warnings and convergence capture
    warn1.ls[[i]][[j]]<-warn1
    conv.ls[[i]][[j]]<-conv
    
    #Model parameters
    analysis.ls[[i]][[j]]<-summary(mod)
    
    
    #Wald p-values
    anova.ls[[i]][[j]]<-Anova(mod)
    
    
    
    #AIC/BIC
    AIC.ls[[i]][[j]]<-summary(mod)$AICtab
    
    
    #Parameter estimates fixed and random
    parameter_estimates_f.ls[[i]][[j]]<-coef(summary(mod))
    parameter_estimates_r.ls[[i]][[j]]<-as.data.frame(VarCorr(mod))
    names(summary(mod))
    
    
    #Log-likelihood test
    loglik.ls[[i]][[j]]<-summary(mod)$logLik
    
    
    
  }
}


#--------------------------------
#warnings/convergence
#--------------------------------
warn1.ls[[1]]
class(warn1.ls)




warn_e.mat <- sapply(warn1.ls, function(x) sapply(x, function(y) !is.null(y$e)))
rownames(warn_e.mat) <- names(model)
colnames(warn_e.mat) <- resp

conv.ls
conv_opt.mat<-mapply(function(x) (sapply(x,function(y)(y[["opt"]] ))),conv.ls)
conv_opt.mat

conv_code.mat<-mapply(function(x) (sapply(x,function(y)(y[[2]][["code"]] ))),conv.ls)
conv_code.mat

conv_messages.mat<-mapply(function(x) (sapply(x,function(y)(y[[2]][["messages"]] ))),conv.ls)
conv_messages.mat

#MAIN ONE
warn_w.mat <- sapply(warn1.ls, function(x) sapply(x, function(y) !is.null(y$w)))
rownames(warn_w.mat) <- names(model)
colnames(warn_w.mat) <- resp

#save(conv.ls,warn1.ls, file="warn.RData")

#need to rerun some models 

#-------------------------------------------------------------------------------
# Checking models and proteins that failed to converge
#-------------------------------------------------------------------------------

had_warning <- sapply(warn1.ls, function(protein_results) {
  sapply(protein_results, function(one_model) !is.null(one_model$w))
})
# had_warning is now a clean logical matrix: 18 rows (models) x 451 cols (proteins)
# TRUE = this fit actually threw a warning, FALSE = it didn't

failed_idx <- which(had_warning, arr.ind = TRUE)
colnames(failed_idx) <- c("model_j", "protein_i")

nrow(failed_idx)   # sanity check — should be a small number, not 8118

had_error <- sapply(warn1.ls, function(protein_results) {
  sapply(protein_results, function(one_model) !is.null(one_model$e))
})
error_idx <- which(had_error, arr.ind = TRUE)
colnames(error_idx) <- c("model_j", "protein_i")
nrow(error_idx)

table(names(model)[failed_idx[,"model_j"]])   # which formulas are failing most
table(resp[failed_idx[,"protein_i"]])          # which proteins are failing most


#------------
#RERUN MODELS 
#------------
#rerun models which failed convergence (29 pairs, identified via warn1.ls / had_warning check)

nrow(failed_idx)

for (idx in 1:nrow(failed_idx)) {
  
  i <- failed_idx[idx, "protein_i"]
  j <- failed_idx[idx, "model_j"]
  
  cat("i=", i, "\n")
  cat("j=", j, "\n")
  print(resp[i])
  
  rm(form, warn1, mod, conv)
  
  #name model
  k <- names(model)[j]
  
  # Get data
  tmp <- mother
  w <- which(names(tmp) == resp[i]); w
  names(tmp)[w] <- "btmean"
  tmp <- tmp[!is.na(tmp$btmean), ]
  
  #Model formula
  form <- as.formula(model[[j]])
  
  #Fit model
  warn1 <- ch(
    mod <- lmer(formula = form, weight = weightTEDDY, data = tmp,
                control = lmerControl(optimizer = "bobyqa",
                                      optCtrl = list(maxfun = 2e5)))
  )
  
  conv <- mod@optinfo$conv
  
  #Warnings and convergence capture
  warn1.ls[[i]][[j]] <- warn1
  conv.ls[[i]][[j]]  <- conv
  
  #Model parameters
  analysis.ls[[i]][[j]] <- summary(mod)
  
  #Wald p-values
  anova.ls[[i]][[j]] <- Anova(mod)
  
  #AIC/BIC
  AIC.ls[[i]][[j]] <- summary(mod)$AICtab
  
  #Parameter estimates fixed and random
  parameter_estimates_f.ls[[i]][[j]] <- coef(summary(mod))
  parameter_estimates_r.ls[[i]][[j]] <- as.data.frame(VarCorr(mod))
  
  #Log-likelihood test
  loglik.ls[[i]][[j]] <- summary(mod)$logLik
  
}

# ============================================================================
# MASTER TABLE: all fixed-effect results for every protein x model x term
# Add after the re-run loop, before "Checking which models got fixed"
# ============================================================================

master_results <- dplyr::bind_rows(
  
  lapply(seq_along(resp), function(i) {
    
    dplyr::bind_rows(
      
      lapply(seq_along(model), function(j) {
        
        coef_tab <- parameter_estimates_f.ls[[i]][[j]]
        
        # Skip a protein-model pair if no coefficient table was created
        if (is.null(coef_tab)) return(NULL)
        
        out <- as.data.frame(coef_tab)
        out$Term <- rownames(out)
        rownames(out) <- NULL
        
        # Identify the protein and fitted model
        out$Protein <- resp[i]
        out$ModelName <- names(model)[j]
        out$Formula <- unname(model[[j]])
        
        # Model-quality flags
        out$HadWarning <- !is.null(warn1.ls[[i]][[j]]$w)
        out$HadError <- !is.null(warn1.ls[[i]][[j]]$e)
        
        conv <- conv.ls[[i]][[j]]
        out$OptimizerCode <- if (!is.null(conv$opt)) conv$opt else NA
        
        out$ConvergenceMessage <- if (!is.null(conv$lme4$messages)) {
          paste(unlist(conv$lme4$messages), collapse = " | ")
        } else {
          NA_character_
        }
        
        # Raw coefficient-level significance only
        p_col <- grep("^Pr\\(", names(out), value = TRUE)
        if (length(p_col) == 1) {
          out$RawPValue <- out[[p_col]]
          out$RawP_lt_0.05 <- out$RawPValue < 0.05
        } else {
          out$RawPValue <- NA_real_
          out$RawP_lt_0.05 <- NA
        }
        
        out
        
      })
    )
  })
)

# Put identifying columns first
master_results <- master_results %>%
  dplyr::select(
    Protein, ModelName, Formula, Term,
    Estimate, `Std. Error`, df, `t value`, RawPValue, RawP_lt_0.05,
    HadWarning, HadError, OptimizerCode, ConvergenceMessage,
    dplyr::everything()
  )

write.csv(
  master_results,
  "master_results_proteomics_raw.csv",
  row.names = FALSE
)

message("Saved: master_results_proteomics_raw.csv")

#---------------------------------------------------------------------------
# Checking which models got fixed
#---------------------------------------------------------------------------
still_failing <- sapply(1:nrow(failed_idx), function(idx) {
  i <- failed_idx[idx, "protein_i"]
  j <- failed_idx[idx, "model_j"]
  !is.null(warn1.ls[[i]][[j]]$w)
})

sum(still_failing)          # how many still have a warning after bobyqa
failed_idx[still_failing, ] # which specific ones remain
