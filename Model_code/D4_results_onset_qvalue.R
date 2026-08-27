#https://www.nature.com/articles/s41598-020-73883-7
#https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#what-are-the-p-values-listed-by-summaryglmerfit-etc.-are-they-reliable
#https://stackoverflow.com/questions/77067138/predict-function-produces-different-values-than-hand-calculation-in-a-glmer
#https://bbolker.github.io/mixedmodels-misc/glmmFAQ.html
#https://entnemdept.ufl.edu/Hahn/generalized-linear-mixed-models.html#poisson-and-nb-glmms
#https://stats.stackexchange.com/questions/264437/how-to-determine-degrees-of-freedom-in-linear-mixed-effect-regression
#http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
#======================
#LIBRARY
#======================
#if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install("qvalue")

library(qvalue)
library(stats) #p.adjust
library(ggplot2)
#======================
#FUNCTIONS
#=======================

source("qvaluetruncfunction.R")

#--------------------------------------------------------
#q-value wrapper (replaces p.adjust_LiJiBH.mat)
#--------------------------------------------------------

q.adjust.mat <- function(pvalue.mat) {
  qvalue.mat <- pvalue.mat
  for (r in 1:nrow(pvalue.mat)) {
    p <- as.numeric(pvalue.mat[r, ])
    if (all(is.na(p))) {
      qvalue.mat[r, ] <- NA
    } else {
      qobj <- qvalue_truncp(p)
      qvalue.mat[r, ] <- qobj$qvalues
    }
  }
  return(qvalue.mat)
}

#function to calculate loglikelihood ratio test
#CAUTION ONLY SUITABLE FOR ADDING 1 RANDOM EFFECT as df=1
#x is nested within y
lrt<-function(x,y){
  teststat<- -2*(x-y)
  p.val <- pchisq(teststat, df = 1, lower.tail = FALSE)
  return(p.val)
}

#function to calculate Walds test
waldtest<-function(est,se,df){
  teststat<-(est)^2 / unlist(se)^2
  p.val <- 1 - pchisq(teststat, df = df)
  return(p.val)
}


resp
resp3
cbind(resp,resp3)

#====================
#CORRELATION MATRIX
#====================
#NOTE: this matrix was previously required by the Li and Ji BH method.
#Now that adjustment uses q-values (q.adjust.mat / qvalue_truncp) instead,
#c.mat is only used below for the exploratory correlation plots (corrplot) -
#it no longer feeds into any p-value adjustment.
names(mother)


c.mat<-cor(as.matrix(mother[, resp]),use="pairwise.complete.obs")
c.mat
dim(c.mat)
colnames(c.mat)<-resp3
rownames(c.mat)<-resp3


library(corrplot)
#corrplot(c.mat)
#corrplot(c.mat,method='color',order='AOE',diag=F)
#corrplot.mixed(c.mat,order='AOE')

# Basic
jpeg("corrplot_basic.jpeg", width = 5000, height = 5000, res = 300)
corrplot(c.mat,
         tl.cex = 0.2,
         tl.col = "black")
dev.off()

# Colour method with AOE ordering
jpeg("corrplot_color_AOE.jpeg", width = 5000, height = 5000, res = 300)
corrplot(c.mat,
         method = 'color',
         order  = 'AOE',
         diag   = FALSE,
         tl.cex = 0.2,
         tl.col = "black")
dev.off()

#================
#LMER
#=================

#---------------
#loglik extraction
#---------------

#---extract  Log likelihood values into matrix
loglik.mat<-mapply(function(x) (sapply(x,function(y)(y[[1]] ))),loglik.ls)
loglik.mat;class(loglik.mat)
write.csv(loglik.mat, "loglik_mat.csv")
#mother id important risk set not


#RANDOM EFFECTS -test significance of
#riskset
lrt(unlist(loglik.mat[2,]),unlist(loglik.mat[1,]))
lrt(unlist(loglik.mat[2,]),unlist(loglik.mat[1,])) <0.05
#mother_id
lrt(unlist(loglik.mat[3,]),unlist(loglik.mat[1,]))
lrt(unlist(loglik.mat[3,]),unlist(loglik.mat[1,]))< 0.05

#riskset
lrt(unlist(loglik.mat[5,]),unlist(loglik.mat[4,]))
lrt(unlist(loglik.mat[5,]),unlist(loglik.mat[4,])) <0.05
#mother_id
lrt(unlist(loglik.mat[6,]),unlist(loglik.mat[4,]))
lrt(unlist(loglik.mat[6,]),unlist(loglik.mat[4,]))<0.05

#riskset
lrt(unlist(loglik.mat[8,]),unlist(loglik.mat[7,]))
lrt(unlist(loglik.mat[8,]),unlist(loglik.mat[7,]))<0.05
#mother_id
lrt(unlist(loglik.mat[9,]),unlist(loglik.mat[7,]))
lrt(unlist(loglik.mat[9,]),unlist(loglik.mat[7,]))<0.05

#riskset
#riskset important IFN-G, 
lrt(unlist(loglik.mat[11]),unlist(loglik.mat[10,]))
lrt(unlist(loglik.mat[11]),unlist(loglik.mat[10,]))<0.05
#mother_id
lrt(unlist(loglik.mat[12,]),unlist(loglik.mat[10,]))
lrt(unlist(loglik.mat[12,]),unlist(loglik.mat[10,]))<0.05


#riskset
lrt(unlist(loglik.mat[14,]),unlist(loglik.mat[13,]))
lrt(unlist(loglik.mat[14,]),unlist(loglik.mat[13,]))<0.05
#mother_id
lrt(unlist(loglik.mat[15,]),unlist(loglik.mat[13,]))
lrt(unlist(loglik.mat[15,]),unlist(loglik.mat[13,]))<0.05

#riskset
lrt(unlist(loglik.mat[17,]),unlist(loglik.mat[16,]))
lrt(unlist(loglik.mat[17,]),unlist(loglik.mat[16,]))<0.05
#mother_id
lrt(unlist(loglik.mat[18,]),unlist(loglik.mat[16,]))
lrt(unlist(loglik.mat[18,]),unlist(loglik.mat[16,]))<0.05




#------------------------
#testing of fixed effects
#------------------------
#EXTRACT P-VALUES 

#-----------
#investigate output of anova---------

anova.ls[1]
sink(file="anova.txt")
print(anova.ls)
sink()

#one response all models
#chisquare of 3rd model
sapply(anova.ls,function(x) x[[3]])

#all models have btmean 
#form_adj_des_rm
#btmean p-value
as.data.frame(anova.ls[[1]][[4]])[1,3]

#btmean:t1d m13-15
#interaction pvalue unadjusted
#form_unadj_t1d_int_rm
as.data.frame(anova.ls[[1]][[13]])[3,3]

#btmean:t1d m16-18
#interaction pvalue
#form_adj_des_t1d_int_rm
as.data.frame(anova.ls[[1]][[16]])[6,3]


#this is first row of anova table for all responses, model 1, third position 
#i.e pvalue of btmean
mapply(function(x) x[[1]][[3]], anova.ls)


#this is 3rd row of anova table for all responses, model 13
#form_unadj_rm
mapply(function(x) x[[13]][[3,3]], anova.ls)


#this is the btmean chi-square estimate
tmp<-mapply(function(x) (sapply(x,function(y)(y[[1]]))),anova.ls)
tmp

#this is the btmean p-value estimate
tmp<-mapply(function(x) (sapply(x,function(y)(y[[3]]))),anova.ls)
tmp



#------
#this is the case v control p-value------------------------
#models1-18
pvalue_line1<-mapply(function(x) (sapply(x,function(y)(y[[3]][1]))),anova.ls)
pvalue_line1
pvalue_cc<-pvalue_line1
#pvalue_cc[19:30,]<-NA
pvalue_cc

write.csv(pvalue_cc,"pvalue_cc_cal.csv")

#q-values
qvalue_cc<-q.adjust.mat(pvalue_cc)
qvalue_cc<-as.data.frame(qvalue_cc)
qvalue_cc

write.csv(qvalue_cc,"qvalue_cc_cal.csv")




#------
#this is the visit p-value------------------------
pvalue_line2<-mapply(function(x) (sapply(x,function(y)(y[[3]][2]))),anova.ls)
pvalue_visit<-pvalue_line2
write.csv(pvalue_visit,"pvalue_visit_cal.csv")

#q-values
qvalue_visit<-q.adjust.mat(pvalue_visit)
View(qvalue_visit)
write.csv(qvalue_visit,"qvalue_visit_cal.csv")


#------
#this is the cc &visit interaction p-value--------------
#only applicable to models 13-18, 
#column 3 (pvalue) 3 row (interaction)
pvalue_line3<-mapply(function(x) (sapply(x,function(y)(y[[3]][3]))),anova.ls)
#column 3 (pvalue) 4 row (interaction)
pvalue_line4<-mapply(function(x) (sapply(x,function(y)(y[[3]][4]))),anova.ls)
#column 3 (pvalue) 5 row (interaction)
pvalue_line5<-mapply(function(x) (sapply(x,function(y)(y[[3]][5]))),anova.ls)
#column 3 (pvalue) 6 row (interaction)
pvalue_line6<-mapply(function(x) (sapply(x,function(y)(y[[3]][6]))),anova.ls)
#column 3 (pvalue) 10 row (interaction)
pvalue_line10<-mapply(function(x) (sapply(x,function(y)(y[[3]][10]))),anova.ls)



#models 1-3 
pvalue_ccvisit<-pvalue_line3
pvalue_ccvisit[4:6,]<-pvalue_line5[4:6,]
pvalue_ccvisit[7:9,]<-pvalue_line10[7:9,]
pvalue_ccvisit[10:12,]<-pvalue_line4[10:12,]
pvalue_ccvisit[13:15,]<-pvalue_line6[13:15,]
pvalue_ccvisit[16:18,]<-pvalue_line10[16:18,]

write.csv(pvalue_ccvisit,"pvalue_ccvisit_cal.csv")

#q-values
qvalue_ccvisit<-q.adjust.mat(pvalue_ccvisit)
qvalue_ccvisit<-as.data.frame(qvalue_ccvisit)
View(qvalue_ccvisit)


write.csv(qvalue_ccvisit,"qvalue_ccvisit_cal.csv")



#------
#this is the cc &visit &t1d interaction p-value--------------

#column 3 (pvalue) 7 row (interaction)
pvalue_line7<-mapply(function(x) (sapply(x,function(y)(y[[3]][7]))),anova.ls)
#column 3 (pvalue) 9 row (interaction)
pvalue_line9<-mapply(function(x) (sapply(x,function(y)(y[[3]][9]))),anova.ls)

#column 3 (pvalue) 13 row (interaction)
pvalue_line13<-mapply(function(x) (sapply(x,function(y)(y[[3]][13]))),anova.ls)

pvalue_t1dccvisit<-pvalue_line13
pvalue_t1dccvisit[10:12,]<-pvalue_line7[10:12,]
pvalue_t1dccvisit[13:15,]<-pvalue_line9[13:15,]

write.csv(pvalue_t1dccvisit,"pvalue_t1dccvisit_cal.csv")

#q-values
qvalue_t1dccvisit<-q.adjust.mat(pvalue_t1dccvisit)
qvalue_t1dccvisit<-as.data.frame(qvalue_t1dccvisit)
View(qvalue_t1dccvisit)


write.csv(qvalue_t1dccvisit,"qvalue_t1dccvisit_cal.csv")


#------
#this is the cc &t1d interaction p-value--------------


#column 3 (pvalue) 11 row (interaction)
pvalue_line11<-mapply(function(x) (sapply(x,function(y)(y[[3]][11]))),anova.ls)

pvalue_t1dcc<-pvalue_line5
pvalue_t1dcc[1:9,]<-NA

pvalue_t1dcc[13:15,]<-pvalue_line7[13:15,]
pvalue_t1dcc[16:18,]<-pvalue_line11[16:18,]

write.csv(pvalue_t1dcc,"pvalue_t1dcc_cal.csv")

#q-values
qvalue_t1dcc<-q.adjust.mat(pvalue_t1dcc)
qvalue_t1dcc<-as.data.frame(qvalue_t1dcc)
View(qvalue_t1dcc)

write.csv(qvalue_t1dcc,"qvalue_t1dcc_cal.csv")



#------
#this is the visit &t1d interaction p-value--------------

#column 3 (pvalue) 8 row (interaction)
pvalue_line8<-mapply(function(x) (sapply(x,function(y)(y[[3]][8]))),anova.ls)

#column 3 (pvalue) 12 row (interaction)
pvalue_line12<-mapply(function(x) (sapply(x,function(y)(y[[3]][12]))),anova.ls)

pvalue_t1dvisit<-pvalue_line6
pvalue_t1dvisit[1:9,]<-NA

pvalue_t1dvisit[13:15,]<-pvalue_line8[13:15,]
pvalue_t1dvisit[16:18,]<-pvalue_line12[16:18,]


write.csv(pvalue_t1dvisit,"pvalue_t1dvisit_cal.csv")



#q-values
qvalue_t1dvisit<-q.adjust.mat(pvalue_t1dvisit)
qvalue_t1dvisit<-as.data.frame(qvalue_t1dvisit)
View(qvalue_t1dvisit)

write.csv(qvalue_t1dvisit,"qvalue_t1dvisit_cal.csv")


#---------------
#AIC extraction
#---------------

#---extract AIC values into matrix
AIC.mat<-mapply(function(x) (sapply(x,function(y)(y[[1]] ))),AIC.ls)
AIC.mat;class(AIC.mat)
write.csv(AIC.mat, "AIC_mat.csv")


#-------------------------
#PARAMETER ESTIMATES FIXED
#-------------------------

#----===========
#BTMEAN
#----===========


#parameter_estimates_f.ls[[1]]

#this is the btmean parameter estimate------------------------
#parameter_estimates_btmean<-mapply(function(x) (sapply(x,function(y)as.data.frame(y[[2]][1]))),parameter_estimates_f.ls)
#parameter_estimates_btmean<-parameter_estimates_btmean[1:18,]
#head(parameter_estimates_btmean)
#parameter_estimates_btmean<- data.frame(apply(parameter_estimates_btmean,2,function(x) as.numeric(x)))
#row is model , column is response

#note divide by 4 as dilution factor
#parameter_estimates_btmean<-parameter_estimates_btmean/4




#----===========
#STANDARD DEVIATION
#----===========


#sapply(parameter_estimates_f.ls,function(y)as.data.frame(y[[2]])[2,2])
#parameter_estimates_btmean_se<-mapply(function(x) (sapply(x,function(y)as.data.frame(y)[2,2])),parameter_estimates_f.ls)
#parameter_estimates_btmean_se
#parameter_estimates_btmean_se<-parameter_estimates_btmean_se[1:18,]
#head(parameter_estimates_btmean_se)
#parameter_estimates_btmean_se<- data.frame(apply(parameter_estimates_btmean_se,2,function(x) as.numeric(x)))

#parameter_estimates_btmean_se<-parameter_estimates_btmean_se/4
