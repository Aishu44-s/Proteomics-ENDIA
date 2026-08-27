
#==============================
#CALCULATE FDR ADJUSTED PVALUES
#==============================
library(qvalue)
unadjpvalue<-pvalue_cc

d1<-dim(unadjpvalue)[1];d1
d2<-dim(unadjpvalue)[2];d2

adjpvalue<-matrix(NA, nrow=d1,ncol=(d2))
colnames(adjpvalue)<-names(unadjpvalue)
rownames(adjpvalue)<-rownames(unadjpvalue)

source("qvaluetruncfunction.R")

#qvalue
for(i in 1:d1){
p<-unlist(unadjpvalue[i,]);p

q2<-qvalue_truncp(p,fdr.level=0.1);q2
q2<-qvalue(p,fdr.level=0.1);q2
adjpvalue[i,]<-q2$qvalue

}

