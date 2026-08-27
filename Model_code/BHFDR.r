https://www.r-bloggers.com/2023/07/the-benjamini-hochberg-procedure-fdr-and-p-value-adjusted-explained/

p_values_log_a<-cbind(pval_visit_lme_log_a,pval_t1d_lme_log_a,pval_case_lme_log_a,
pval_visit_t1d_lme_log_a,pval_visit_case_lme_log_a,pval_t1d_case_lme_log_a,
pval_visit_t1d_case_lme_log_a)

p_values_sqrt_a<-cbind(pval_visit_lme_sqrt_a,pval_t1d_lme_sqrt_a,pval_case_lme_sqrt_a,
pval_visit_t1d_lme_sqrt_a,pval_visit_case_lme_sqrt_a,pval_t1d_case_lme_sqrt_a,
pval_visit_t1d_case_lme_sqrt_a)

final_p_values_a<-p_values_log_a

resp
final_p_values_a[9,]<-p_values_sqrt_a[9,]

#-------------------------
#FDR

#i=3
#pvalues<-final_p_values_a[,i]
pvalues<-c(final_p_values_a)

#m=length(pvalues)
m=length(resp)
ranks<-rank(pvalues, ties.method = "last")

p_m_over_k<-pvalues*m/ranks

for (r in length(pvalues):1) {
  print(p_m_over_k[ranks>=r])
}
pvalues_adj<-c()
for (i in 1:length(pvalues)) {
  
  # find the rank
  tmp_rank<-ranks[i]
  
  # get all the p_m_over_k that are greater or equal to this rank
  # and get the min value
  pvalues_adj<-c(pvalues_adj, min(1,min(p_m_over_k[ranks>=tmp_rank])))
}
print(pvalues_adj)

d1<-dim(final_p_values_a)[1]
d2<-dim(final_p_values_a)[2]

pvalues_adj<-matrix(pvalues_adj,nrow=d1,ncol=d2,byrow=F)
#check read in the best way
#tmp<-matrix(pvalues,nrow=d1,ncol=d2,byrow=F)

#---------------------------------------- 
fdrs<-p.adjust(pvalues, method="fdr")
print(fdrs)
fdrs<-matrix(fdrs,nrow=d1,ncol=d2,byrow=F)
fdrs
final_p_values_a

#----------------------------------
library(FDRestimation)
i=3
p.fdr(pvalues=final_p_values_a,threshold="0.10")


#-----------------------------------
#https://github.com/StoreyLab/qvalue
#install.packages("devtools")
library("devtools")
install_github("jdstorey/qvalue")

library(qvalue)


#i=3
#pvalues<-final_p_values_a[,i]
pvalues<-c(final_p_values_a)


qvalues<-qvalue(p=pvalues,fdr.level =0.1)

fdrs<-matrix(qvalues$qvalues,nrow=d1,ncol=d2,byrow=F)
fdrs


