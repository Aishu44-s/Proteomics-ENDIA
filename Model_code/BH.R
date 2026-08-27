https://www.r-bloggers.com/2023/07/the-benjamini-hochberg-procedure-fdr-and-p-value-adjusted-explained/

1.Sort the p-values in ascending order.
2.Assign a rank or position (denoted as “i”) to each p-value based on its position in the sorted list.
3.Calculate the adjusted p-value for each p-value using the formula: 
Adjusted p-value padjusted(i)=min{1,minj≥i{mp(j)/j}}
 where:
padjusted(i) represents the adjusted p-value for the i-th ranked p-value.
p(j) represents the j-th ranked p-value.
m represents the total number of hypotheses being tested. 
In this formula, for each p-value at rank i, the adjusted p-value is calculated as the minimum value between 1
and the minimum ratio obtained by dividing m (the total number of hypotheses) by j (the rank) 
for all j greater than or equal to i. 
The purpose of taking the minimum ratio is to find the smallest value that controls the FDR 
while considering all hypotheses ranked at or above the current p-value. 
Taking the minimum ensures that the adjusted p-value is conservative, 
accounting for the multiple-testing scenario.
Report the adjusted p-values alongside the original p-values. 
The adjusted p-values provide a corrected measure of statistical significance that 
accounts for multiple testing.




p.raw<-p_jack_unadj_m #p
#p.raw<-c(0.01,0.001, 0.05, 0.20, 0.15, 0.15)

fdrs<-p.adjust(p.raw, method="BH");fdrs

p.rank<-rank(x=p.raw,ties.method="last")#j
length(p.rank)#m

pmdivj<-p.raw*length(p.raw)/p.rank



p.adjusted<-c()
for (i in 1:length(p.rank)){
r<-p.rank[i];r
p.adjusted[i]<-min(1,min(pmdivj[p.rank>=r]))
}
p.adjusted


