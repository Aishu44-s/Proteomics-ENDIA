##BASED ON FILE FROM DOMAGOJ
#Helena Oakey

#correlations matrix on relevant responses
#c.mat<-cor(as.matrix(onset[,4:17]),use="pairwise.complete.obs")
#c.mat


#p.raw is the pvalues and cmat is the correlations matrix
#p.raw is matrix of response in columns, models in rows
#will pick one model so want to each model separately
#for testing
#p.mat<-pvalue_btmean
#p.mat<-p.mat[1:5,1:5]
#p.df<-as.data.frame(p.mat)


p.adjust_LiJiBH.mat<-function(p.df,c.mat){

egn.val<-eigen(x=c.mat,symmetric="TRUE")$values
#get rid of any negative values
#note if neg values implies matrix not positive definite, may want to bend matrix to be positive def.
egn.val<-egn.val[egn.val>0]


#Li and Ji 2005
#doi:10.1038/sj.hdy.6800717
#Eqn 5
#Effective number of independent tests given values are correlated
M.eff<-sum(ifelse(egn.val>=1,1,0)+ (egn.val-floor(egn.val)))

#create matrix to contain results
result.mat<-matrix(NA,nrow=nrow(p.df),ncol=ncol(p.df))
rownames(result.mat)<-rownames(p.df)
colnames(result.mat)<-colnames(p.df)


#for each model in matrix
for(j in 1:nrow(p.df)){
#row corresponds to same model different responses
p.raw<-unlist(p.df[j,])

#rank raw values
p.rank=rank(x=p.raw,ties.method="max")

M<-length(p.raw)

#Eqn 11
#Adjust p-values
pmdivj = p.raw * M.eff *( M - 1)/(M + p.rank * M.eff - p.rank - M.eff)


#BH q-values
p.adjusted<-c()
for (i in 1:length(p.rank)){
r<-p.rank[i];r
p.adjusted[i]<-min(1,min(pmdivj[p.rank>=r]))
}#for i
result.mat[j,]<-p.adjusted
}#for j

l<-list(M.eff,result.mat)
names(l)<-c("M.eff","LiJi_BH_adjusted_p.values")
return(l)

}


