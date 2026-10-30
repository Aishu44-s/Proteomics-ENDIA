
#===================================
#LIBRARIES and FUNCTIONS
#===================================
#rm(list=ls())

"%notin%"<-Negate("%in%")

library(readxl)

#######################
#TH17 DATA
#######################
#includes standards, controls and blank
#mother<-read.csv("/Volumes/ENDIA22-Q5126/ENDIA_Metadata/20240820_Mother_code/Maternal example code-selected/cyto_mother_cal.csv",header=T)
mother<-read.csv("C:/Users/Dell/Documents/BIOX7011/ENDIA-Proteomics/mother_vsn_merged.csv",header=T)
names(mother)
dim(mother)


#===========================
# SET UP DATA FOR ANALYSIS
#===========================
head(mother)
tail(mother)
str(mother)

mother$corrected_visit_sample<-factor(mother$corrected_visit_sample)
table(mother$corrected_visit_sample,exclude=NULL)

mother$mother_id<-factor(mother$mother_id)

summary(mother$maternal_age_at_conception)

#mother$riskset<-factor(mother$nest_singletime)

mother$riskset<-factor(mother$deidentified_nest_id_new)

mother$maternal_T1D <-as.factor(mother$maternal_T1D )
table(mother$"maternal_T1D" )

mother$caset1dcat<-paste("IA",mother$case,"T1D",mother$maternal_T1D,sep="_")
mother$caset1dcat<-factor(mother$caset1dcat)
table(mother$caset1dcat)

mother$casecat<-factor(mother$case)
table(mother$casecat)


head(mother$maternal_HLA)
table(mother$maternal_HLA,exclude=NULL)
which(is.na(mother$maternal_HLA))
mother$maternal_HLA<-factor(mother$maternal_HLA)

table(mother$maternal_smoking )
#mother$maternal_smoking<-factor(mother$maternal_smoking,levels=c(0,1),labels=c("No","Yes"))
mother$maternal_smoking<-factor(mother$maternal_smoking)


#mother$age_sample_collection_month<-mother$age_at_sample_collection_day/365*12
# I do not have age at sample collection date in my meta data

mother$age_sample_collection_month <- mother$maternal_age_at_conception * 12


mother$infant_sex<-factor(mother$infant_sex)



names(mother)[40:490] # my protein columns
resp <- names(mother)[grepl("^(sp|tr)\\.", names(mother))]
cat("Number of proteins:", length(resp), "\n")


#names of response variables

#CALIBRATED RESPONSES
#resp<-names(mother)[36:40]
#resp


#resp3<-c(
#"CD40 Ligand", 
#"IFN-\u03b3",
#"IL-1\u03b2",
#"IL-4",
#"IL-6")

#cbind(resp,resp3)


#checking before running model:
cat("\nKey variables:\n")
cat("weightTEDDY NAs:", sum(is.na(mother$weightTEDDY)), "\n")
cat("timefromonsetmonths NAs:", sum(is.na(mother$timefromonsetmonths)), "\n")
cat("casecat levels:", levels(mother$casecat), "\n")
cat("corrected_visit_sample levels:", levels(mother$corrected_visit_sample), "\n")
cat("maternal_T1D levels:", levels(mother$maternal_T1D), "\n")
cat("riskset levels (n):", length(levels(mother$riskset)), "\n")
cat("mother_id levels (n):", length(levels(mother$mother_id)), "\n")
cat("maternal_smoking levels:", levels(mother$maternal_smoking), "\n")
cat("maternal_HLA NAs:", sum(is.na(mother$maternal_HLA)), "\n")
cat("maternal_conception_bmi NAs:", sum(is.na(mother$maternal_conception_bmi)), "\n")
table(mother$caset1dcat)