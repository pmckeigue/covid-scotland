#########################################################################################################################################################################
########### COVID-19 Case Control Analysis   ################################################################################################
#########################################################################################################################################################################

###############################################################################
#  Created on 14th Apr 2020
#  ------------------------------------------------------------
#  Author:  Amanda Weir

# rm(list=ls()) 


library(tidyverse)
library(janitor)
library(scales)
library(foreign)
library(Hmisc)
library(knitr)
library(rmarkdown)
library(lubridate)
library(readr)
library(dplyr)
library(survival)
library(epiDisplay)
library(data.table)
library(MASS)
library(caret)

options(tibble.print_max = Inf)
options(tibble.width = Inf)

# Read in data

cov <- read.csv(paste0(westdata_path, "/Case_control/data/case_control_linked_20200416.csv"))
cov[is.na(cov)]<-""			#removes NAs and replaces them with blanks
head(cov)

#Recodes
cov$COVID<-cov$CASE
cov$SPECDATE2<-format(as.Date(cov$SPECDATE,"%m/%d/%Y"), "%d/%m/%Y")
cov$wk<-isoweek(dmy(cov$SPECDATE2))
cov$Date.Death2<-sub(" ", "", cov$Date.Death)
cov[cov$Date.Death!="",]$Date.Death2<-format(as.Date(cov[cov$Date.Death!="",]$Date.Death,"%m/%d/%Y"), "%d/%m/%Y")
cov[is.na(cov)]<-""			#removes NAs and replaces them with blanks
cov$TimeDeath<-rep(9999, length(cov[,1]))
cov[cov$Date.Death2!="",]$TimeDeath<-as.Date(cov[cov$Date.Death2!="",]$Date.Death2, "%d/%m/%Y")-as.Date(cov[cov$Date.Death2!="",]$SPECDATE2, "%d/%m/%Y")
 
cov$dead28b<-rep(9999,length(cov[,1]))  # -1 died before positive test (controls only) / 0 alive / 1 died within 28 days of positive test / 2 died after 28 days positive test
cov[cov$TimeDeath<0,]$dead28b<--1
cov[cov$TimeDeath>=0&cov$TimeDeath<=28,]$dead28b<-1
cov[cov$TimeDeath>28,]$dead28b<-2
cov[cov$TimeDeath==9999,]$dead28b<-0
table(cov$CASE,cov$dead28b)

# Remove controls who have died on or before the test data of matched control
cov<-cov[!cov$CHI%in%cov[cov$COVID==0&cov$TimeDeath<=0,]$CHI,]

# Recode age groups
table(cov$agegp)
cov$agegp2<-cov$AGE
table(cov$AGE)

cov$agegp2<-replace(rep(1,length(cov[,1])),cov$AGE>=50,2) #three broad age groups as per protocol
cov$agegp2<-replace(cov$agegp2,cov$AGE>70,3)
table(cov$agegp2)
cov$agegp2<-as.factor(cov$agegp2)
levels(cov$agegp2)<-c("<50","50-70","70+")
table(cov$agegp2)

cov$agegp3<-replace(rep(1,length(cov[,1])),cov$AGE>=5,2)  # five year age bands for Paul
cov$agegp3<-replace(cov$agegp3,cov$AGE>=10,3)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=15,4)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=20,5)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=25,6)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=30,7)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=35,8)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=40,9)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=45,10)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=50,11)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=55,12)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=60,13)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=65,14)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=70,15)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=75,16)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=80,17)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=85,18)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=90,19)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=95,20)
cov$agegp3<-replace(cov$agegp3,cov$AGE>=100,21)
table(cov$agegp3)
cov$agegp3<-as.factor(cov$agegp3)
levels(cov$agegp3)<-c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44",
  "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95-99","100-104")
table(cov$agegp3)

cov$agegp4<-replace(rep(1,length(cov[,1])),cov$AGE>=35,2) #five year age bands but <35 and 95+ groups merged to avoid deductive disclosure
cov$agegp4<-replace(cov$agegp4,cov$AGE>=40,3)
cov$agegp4<-replace(cov$agegp4,cov$AGE>=45,4)
cov$agegp4<-replace(cov$agegp4,cov$AGE>=50,5)
cov$agegp4<-replace(cov$agegp4,cov$AGE>=55,6)
cov$agegp4<-replace(cov$agegp4,cov$AGE>=60,7)
cov$agegp4<-replace(cov$agegp4,cov$AGE>=65,8)
cov$agegp4<-replace(cov$agegp4,cov$AGE>=70,9)
cov$agegp4<-replace(cov$agegp4,cov$AGE>=75,10)
cov$agegp4<-replace(cov$agegp4,cov$AGE>=80,11)
cov$agegp4<-replace(cov$agegp4,cov$AGE>=85,12)
cov$agegp4<-replace(cov$agegp4,cov$AGE>=90,13)
cov$agegp4<-replace(cov$agegp4,cov$AGE>=95,24)
table(cov$agegp4)
cov$agegp4<-as.factor(cov$agegp4)
levels(cov$agegp4)<-c("<35","35-39","40-44",
                      "45-49","50-54","55-59","60-64","65-69","70-74","75-79","80-84","85-89","90-94","95+")
table(cov$agegp4)

# Recode SIMD
cov$simd<-cov$simd2020_sc_quintile
cov[cov$simd=="",]$simd<-"NK"
cov$simd<-as.factor(cov$simd)

#healthboard

table(cov$HBRES)
cov <- cov %>%
  mutate(HBRES2 = hb_recode(HBRES)) #create new variable HB2014 using lookup table hb_recode (Megans code)
cov[cov$HBRES2=="",]$HBRES2<-"NK"
table(cov$HBRES2)
cov[cov$HBRES2=="Orkney",]$HBRES2<-"Island Boards"
cov[cov$HBRES2=="Shetland",]$HBRES2<-"Island Boards"
cov[cov$HBRES2=="Western Isles",]$HBRES2<-"Island Boards"
table(cov$HBRES2)

#  CCI group
cov <- within(cov,{
  CCIgrp <- NA
  CCIgrp [is.na(CCIscore)] <- "[0]No_stay"
  CCIgrp [CCIscore ==0] <- "[0]No_score"
  CCIgrp [CCIscore >=1 & CCIscore <=2] <- "[1-2]Mild"
  CCIgrp [CCIscore >=3 & CCIscore <=5] <- "[3-5]Moderate"
  CCIgrp [CCIscore >=6] <- "[6+]Severe"  })
cov$CCIgrp <- as.factor(cov$CCIgrp)
table(cov$CCIgrp)

## Categorise controls to same group as cases i.e. A, B, C
caselu<-cov[cov$CASE==1,][,c("CASE_NO","group")]
colnames(caselu)[2]<-"casegroup"
cov<-merge(cov, caselu, by=c("CASE_NO"),all.x=T)
table(cov$CASE,cov$casegroup)


# Split Resp and Asthma up until get updated data
cov$resp_x<-cov$resp_c
cov[cov$asthma_c==1,]$resp_x<-0

# Sum up co-morbidities
cov$comorb<-cov$resp_x+cov$lung_c+cov$heart_c+cov$asthma_c+cov$diab.merge+cov$neuro_c+cov$immuno_b+
                    cov$renal+cov$liver+
                    cov$heart.ace.bnf+cov$heart.nonace.bnf+cov$angio.bnf+cov$hiv.bnf
table(cov$comorb)
cov$comorbgp<-cov$comorb
cov[cov$comorb>=7,]$comorbgp<-7
table(cov$comorbgp)

# Compare ACE, ARB, ACE/ARB, Non ACE/ARB, each Non ACE/ARB seperately

cov$acearb<-rep(0,length(cov[,1]))
cov[cov$heart.ace.bnf==1,]$acearb<-1
cov[cov$angio.bnf==1,]$acearb<-1

#factorise

cols <- c("SEX", "agegp2", "simd","HBRES2", "wk" ,              #demographics
          "emerg","icu.hdu.ccu","inpat","PIS",           #other health care exposure variables in lookback upto 2 weeks prior to positive test
         "CCIgrp","comorbgp","resp_c","lung_c","heart_c","asthma_c","diab_c","neuro_c",    #comorbidities
          "immuno_b","renal","liver","hiv.bnf",             #comorbidities
          "heart.ace.bnf","heart.nonace.bnf","angio.bnf", "acearb"  ,         #drugs/estimate of hypertension
          "icu","Level3Day","Renal_Replacement_Therapy", "Advanced_Ventilation")     #severe disease
cov[cols] <- lapply(cov[cols], factor)

# ### Set reference category for non-ordered factors

levels(cov$simd)
cov <- within(cov, simd <- relevel(simd, ref="1"))

#####
#### Frequency tables - function for looping through multiple variables --------------------------------------

### Demographics, ECOSS variables, week number of test -----------------------

# variables for looping through
vars <- c(#"group",
          "agegp2", "SEX", "simd", "HBRES2","wk",   #ethnicity
          #"emerg","icu.hdu.ccu","inpat","PIS","SMR01",
          #"CCIgrp","comorbgp","resp_c", "lung_c", "heart_c", "asthma_c", "diab.merge", "neuro_c", "immuno_b", 
          #"renal", "liver","hiv.bnf","ms.bnf",
          #"heart.ace.bnf","heart.nonace.bnf","angio.bnf","acearb", 
          "icu", "Level3Day","Renal_Replacement_Therapy", "Advanced_Ventilation")

#table of cases by grouping and age
#vars2<-c("agegp2")
#vars2<-c("agegp3")
#vars2<-c("agegp4","comorbgp")

func1 <- lapply(vars, function(x) {
  
  Case_A_N <- table(subset(cov, (cov$group=="A" & cov$COVID==1))[[x]])
  Case_A_ColPerc <- round_half_up(prop.table(Case_A_N)*100, 2)
  tb1 <- cbind(Case_A_N, Case_A_ColPerc)
  
  Case_B_N <- table(subset(cov, (cov$group=="B" & cov$COVID==1))[[x]])
  Case_B_ColPerc <- round_half_up(prop.table(Case_B_N)*100, 2)
  tb2 <- cbind(tb1, Case_B_N, Case_B_ColPerc)
  tb2 <- as.data.frame(tb2)
  
  Case_C_N <- table(subset(cov, (cov$group=="C" & cov$COVID==1))[[x]])
  Case_C_ColPerc <- round_half_up(prop.table(Case_C_N)*100, 2)
  tb3 <- cbind(Case_C_N, Case_C_ColPerc)
  tb3 <- cbind(tb2,tb3)
  tb3 <- as.data.frame(tb3)
  
  Case_Total <- table(subset(cov, (cov$COVID==1))[[x]])
  Case_ColPerc <- round_half_up(prop.table(Case_Total)*100, 2)
  tb4 <- cbind(Case_Total, Case_ColPerc)  
  tb4 <- cbind(tb3,tb4)
  tb4 <- as.data.frame(tb4)
  
  tb4$Variable <- x
  tb4$Level <- rownames(tb4)
  tb4 <- tb4[c("Variable", "Level", "Case_A_N", "Case_A_ColPerc", "Case_B_N", "Case_B_ColPerc", "Case_C_N", "Case_C_ColPerc",
               "Case_Total", "Case_ColPerc")]
    print(tb4)
})
cases_group_age<-as.data.frame(do.call(rbind, func1))

totals<-c("Total","",dim(cov[cov$group=="A"&cov$COVID==1,,])[1],"", dim(cov[cov$group=="B"&cov$COVID==1,])[1],"", 
                          dim(cov[cov$group=="C"&cov$COVID==1,])[1],"", dim(cov[cov$COVID==1,])[1],"")
cases_group_age<-rbind(cases_group_age,totals)
cases_group_age<-as.data.frame(cases_group_age,row.names = NULL)
write.csv(cases_group_age, paste0(westdata_path, "/Case_control/output/CC_Cases_Demog.csv.csv"))

### Comoborbidites -----------------------

# variables for looping through
vars <- c("agegp2", "SEX", "simd", #"ethnicity"
          "CCIgrp","comorbgp",
          "resp_x", "lung_c", "heart_c", "asthma_c", "diab12", "neuro_c", "immuno_b", 
          "renal", "liver","hiv.bnf","ms.bnf",
          "heart.ace.bnf","heart.nonace.bnf","angio.bnf","acearb",
         # "icu", "Level3Day","Renal_Replacement_Therapy", "Advanced_Ventilation",    #used to assign severe disease
          "emerg","icu.hdu.ccu","inpat")

# function for looping through multiple vars to give total, col %, COVID, row %
func1 <- lapply(vars, function(x) {
  
  N_Total <- table(cov[[x]]) 
  
  N_CASE <- table(subset(cov, cov$COVID==1)[[x]])
  Col_Percent_CASE <- round_half_up(prop.table(N_CASE)*100, 2)
  tb1 <- cbind(N_CASE, Col_Percent_CASE)#, Row_Percent_COVID)

  N_CONTROL <- table(subset(cov, cov$COVID==0)[[x]]) # add as.vector() to prevent non-conformable arrays error
  Col_Percent_CONTROL <- round_half_up(prop.table(N_CONTROL)*100, 2)
  tb2 <- cbind(tb1, N_CONTROL, Col_Percent_CONTROL)
  tb2 <- as.data.frame(tb2)
  
  N_Total <- table(cov[[x]])
  Col_Percent_Total <- round_half_up(prop.table(N_Total)*100, 2)
  tb3 <- cbind(N_Total, Col_Percent_Total)
  tb3 <- cbind(tb2, tb3)
  tb3 <- as.data.frame(tb3)  
  
  tb3$Variable <- x
  tb3$Level <- rownames(tb3)
  tb3 <- tb3[c("Variable", "Level", "N_CASE", "Col_Percent_CASE",#"Row_Percent_CASE", 
               "N_CONTROL", "Col_Percent_CONTROL",#"Row_Percent_CONTROL",
              "N_Total", "Col_Percent_Total")]
  print(tb3)
})


# combine output into dataframe for saving
cctab_char.df <- as.data.frame(do.call(rbind, func1))
#cctab_charA.df <- as.data.frame(do.call(rbind, func1))
#cctab_charB.df <- as.data.frame(do.call(rbind, func1))
#cctab_charC.df <- as.data.frame(do.call(rbind, func1))

totals<-c("Total","",dim(cov[cov$COVID==1,,])[1],"", dim(cov[cov$COVID==0,])[1],"", 
          dim(cov)[1],"")
cctab_char.df<-rbind(cctab_char.df,totals)
cctab_char.df<-as.data.frame(cctab_char.df,row.names = NULL)


cc_char.table<-data.table(cctab_char.df)
write.csv(cctab_char.df, paste0(westdata_path, "/Case_control/output/CC_Cases_Controls_Char.csv"))


###----------------------------------------------------------------------------------------------------
### Conditional logistic regression - unadjusted OR across All Cases/Controls / Hospital & Severe / Severe  -----------------------------------------------------------
### controls matched on age, sex and GP practice so exclude age and sex from model. The effect of SIMD will be affected by matchin on GP


### DEMOGRAPHIC VARIABLES ###---------------------------------------

# unadjusted odds (OR - clogit)

## Dataset    #use all data for socio demographic data
#cov                                                      #All cases/controls
#cov[cov$casegroup%in%c("A","B"),]          #Hospital/Severe disease cases only with matched controls
#cov[cov$casegroup=="A",]                                 #Severe cases only with matched controls

# function to loop through list of vars
#can't add age or sex as matched on
#can't add ICU data as these only apply to severe cases

vars <- c(#"AGE","SEX",
          "simd")#, #"ethnicity",  #"HBRES2",
          #"CCIgrp", #"resp_x", "lung_c", 
          #"heart_c", #"asthma_c", 
          #"diab12")#, #"neuro_c", "immuno_b", 
          #"renal", "liver","hiv.bnf",#,"ms.bnf",
          #"heart.ace.bnf","heart.nonace.bnf","angio.bnf","acearb")
          #"icu", "Level3Day", "Renal_Replacement_Therapy", 
          #"Advanced_Ventilation")
          #"emerg","CCIgrp","icu.hdu.ccu","inpat")

formula <- lapply(vars, function(i) {paste("COVID ~", i, "+ strata(CASE_NO)")})

result <- lapply(formula, function(x) {
  z1 <- clogit(as.formula(x), data=cov[cov$casegroup%in%c("A"),])  #data=cov   #data=cov[cov$casegroup%in%c("A","B"),]   #data=cov[cov$casegroup%in%c("A"),]
  out <- summary(z1)
  #  names(out)
  outdf <- as.data.frame(out$coefficients)
  #  dimnames(outdf)[[1]]
  outdf <- cbind(Factor=dimnames(outdf)[[1]], outdf)
  # print(names(outdf))
  outdf$lower <- exp(outdf$coef - (1.96*outdf$`se(coef)`))
  outdf$upper <- exp(outdf$coef + (1.96*outdf$`se(coef)`))
  #  print(outdf)
  
  tb <- outdf[c("Factor", "exp(coef)", "lower", "upper", "Pr(>|z|)")]
  names(tb) <- c("Factor", "OR", "LCL", "UCL", "P")
  tb$OR <- round(tb$OR, 2)
  tb$LCL <- round(tb$LCL, 2)
  tb$UCL <- round(tb$UCL, 2)
  tb$P <- ifelse(tb$P < 0.001, "<0.001", round(tb$P, 3))
  print(tb)
})

cc_clogit_unadj_ALL.table <- as.data.frame(do.call(rbind, result))
cc_clogit_unadj_HOSPSEV.table <- as.data.frame(do.call(rbind, result))
cc_clogit_unadj_SEV.table <- as.data.frame(do.call(rbind, result))

write.csv(cc_clogit.table, paste0(westdata_path, "/Case_control/output/CC_clogit_unadj_ALL.csv"))
write.csv(cc_clogit_unadj_HOSPSEV.table, paste0(westdata_path, "/Case_control/output/CC_clogit_unadj_HOSPSEV.csv"))
write.csv(cc_clogit_unadj_SEV.table, paste0(westdata_path, "/Case_control/output/CC_clogit_unadj_SEV.csv"))


# ### cond.logistic reg - adjusted ----------------------------------------

z <- clogit(COVID  ~ simd + 
              CCIgrp +
              heart_c + 
              diab12 + 
              strata(CASE_NO),
            data=cov[cov$casegroup%in%c("A"),])  #data=cov[cov$casegroup%in%c("A","B","C"),]   #data=cov[cov$casegroup%in%c("A","B"),]    #data=cov[cov$casegroup%in%c("A"),]
summary(z)


#anova(z, z1, test="Chisq")   #compare model fitting

# output table of full and final model
out <- summary(z)
outdf <- as.data.frame(out$coefficients)
outdf <- cbind(Factor=dimnames(outdf)[[1]], outdf)
outdf$lower <- exp(outdf$coef - (1.96*outdf$`se(coef)`))
outdf$upper <- exp(outdf$coef + (1.96*outdf$`se(coef)`))
tb <- outdf[c("Factor", "exp(coef)", "lower", "upper", "Pr(>|z|)")]
tb$"exp(coef)" <- round(tb$"exp(coef)", 2)
tb$lower <- round(tb$lower, 2)
tb$upper <- round(tb$upper, 2)
tb$`Pr(>|z|)` <- ifelse(tb$`Pr(>|z|)` < 0.001, "<0.001", round(tb$`Pr(>|z|)`, 3))
names(tb) <- c("Factor", "OR", "LCL", "UCL", "P")
tb




##---------------------------------------------------------------------------------------
## Assign training/test data set for model building and validation of models  -------------
##---------------------------------------------------------------------------------------


## Conditional logistic regression - multivariate model - SEVERE disease ------------------------------------------------------
# controls matched on age and sex so exclude from model

#cov[cov$casegroup%in%c("A"),]
#sev<-sev[sev$AGE<75,]  #Exclude older ages
#sev<-sev[sev$AGE<70,]
#sev<-sev[sev$AGE<60,]

#x <- data.table::data.table(sev[sev$COVID==1,])#choose training data from all cases with an equal split across A, B, C
x <- data.table::data.table(cov[cov$COVID==1,])#
x[, id := 1:.N]
x[, my.training.flag := 
    #ifelse(id %in% x[createDataPartition(x$CASE_NO, 
    ifelse(id %in% x[createDataPartition(c(x$group), 
                  p = 0.7, list = F)]$id, 1, 0)]
with(x, table(COVID, my.training.flag))
with(x, table(COVID, my.training.flag, group))
cov$training<-rep(0,length(cov[,1]))
cov$test<-rep(0,length(cov[,1]))
cov[cov$CASE_NO%in%x[x$my.training.flag==1,]$CASE_NO,]$training<-1
cov[cov$training==0,]$test<-1
table(cov$training)
table(cov$test)


# function to loop through list of vars
#can't add age or sex as matched on
#can't add ICU data as these only apply to severe cases

vars <- c("resp_x", "lung_c", 
          "heart_c", "asthma_c", 
          "diab12", "neuro_c", "immuno_b", 
          "renal", "liver")#,
          #"hiv.bnf")#,"ms.bnf",  #no HIV, MS cases
          #"emerg",#"CCIgrp",
          #"icu.hdu.ccu","inpat")

formula <- lapply(vars, function(i) {paste("COVID ~", i, "+ strata(CASE_NO)")})

result <- lapply(formula, function(x) {
  z1 <- clogit(as.formula(x), data=cov[cov$training==1&cov$casegroup%in%c("A"),])

    out <- summary(z1)
  #  names(out)
    outdf <- as.data.frame(out$coefficients)
  #  dimnames(outdf)[[1]]
  outdf <- cbind(Factor=dimnames(outdf)[[1]], outdf)
  # print(names(outdf))
  outdf$lower <- exp(outdf$coef - (1.96*outdf$`se(coef)`))
  outdf$upper <- exp(outdf$coef + (1.96*outdf$`se(coef)`))
  #  print(outdf)
  
  tb <- outdf[c("Factor", "exp(coef)", "lower", "upper", "Pr(>|z|)")]
  names(tb) <- c("Factor", "OR", "LCL", "UCL", "P")
  tb$OR <- round(tb$OR, 2)
  tb$LCL <- round(tb$LCL, 2)
  tb$UCL <- round(tb$UCL, 2)
  tb$P <- ifelse(tb$P < 0.001, "<0.001", round(tb$P, 3))
  print(tb)
})

severe.clogit.uni <- as.data.frame(do.call(rbind, result),row.names = NULL)
cc_severe_clogit.univ.table<-data.table(severe.clogit.uni)
write.csv(cc_severe_clogit.univ.table, paste0(westdata_path, "/Case_control/output/COVID_CC_severe_clogit_univ.csv"))


# ### check for confounding between new vars -----------------------
# change in OR <10% in most categories = no strong evidence of confounding
#sev<-cov[cov$casegroup=="A",]
  
with(cov[cov$training==1&cov$casegroup%in%c("A"),], table(asthma_c, SEX, COVID))
z <- clogit(COVID ~ resp_x + heart_c + asthma_c + 
                + renal + diab12 + neuro_c + immuno_b + emerg + inpat + strata(CASE_NO), data=cov[cov$training==1&cov$casegroup%in%c("A"),])
z1 <- clogit(COVID ~ resp_x + heart_c + asthma_c + 
              + renal + diab12 + neuro_c + immuno_b + emerg + inpat +
              immuno_b*AGE + strata(CASE_NO), data=cov[cov$training==1&cov$casegroup%in%c("A"),])

z.out <- round(exp(cbind(coef(z),confint(z))),3)
z1.out <- round(exp(cbind(coef(z1),confint(z1))),3)
z.output <- merge(z.out, z1.out, by="row.names")
z.output$pc <- round((z.output$V1.y - z.output$V1.x)/z.output$V1.y*100, 1)
names(z.output) <- c("Risk factor", "Crude OR", "Crude LCL", "Crude UCL",
                     "Adjusted OR", "Adjusted LCL", "Adjusted UCL", "Percent change")
z.output

summary(z)
summary(z1)
anova(z, z1, test="Chisq")

## Conditional logistic regression - Forward/Backward selection ------------------------------------------------------------------------------
# Forwards and Bacwards selection adding covariates

#covred<-cov[c("COVID","training","casegroup","CASE_NO",
#              "AGE", "SEX",
#              "resp_x", "lung_c", 
#              "heart_c", "asthma_c", 
#              "diab12", "neuro_c", "immuno_b", 
#              "renal", "liver")]    #"hiv.bnf",#,"ms.bnf",
#              #"emerg","CCIgrp","icu.hdu.ccu","inpat"))]     #reduce dataset to key variables for logistic regression

#remove.vars <- c("COVID", "training", "casegroup","CASE_NO")

# Our base model with minimal predictors

#null = clogit(COVID ~ resp_x + strata(CASE_NO), data = covred[covred$training==1&covred$casegroup=="A",])
#summary(null)
null = clogit(COVID ~ resp_x + strata(CASE_NO), data = cov[cov$training==1&cov$casegroup=="A",])
summary(null)

# If interaction terms are enabled

#full = clogit(COVID ~ strata(CASE_NO)+ 
#                resp_x + lung_c + asthma_c + heart_c + diab12 + neuro_c + immuno_b +
#                renal + liver + asthma_c*SEX, data = covred[covred$training==1&covred$casegroup=="A",])
full = clogit(COVID ~ strata(CASE_NO)+ 
                resp_x + lung_c + asthma_c + heart_c + diab12 + neuro_c + immuno_b +
                renal + liver +#,
                emerg + inpat + icu.hdu.ccu, 
                data = cov[cov$training==1&cov$casegroup=="A",])
summary(full)

stepf <- stepAIC(null, scope = list(upper = full, lower = null), family="binomial",
                 direction = "forward",trace=T)
stepf$anova


stepb <- stepAIC(full, scope = list(upper = full, lower = null), family="binomial",
                 direction = "backward",trace=T)
stepb$anova


#modf <- clogit(COVID ~ paste(as.list(stepf$formula)[[3]])[2], data = cov[cov$training==1&cov$casegroup=="A",])
summary(stepf)
summary(stepb)
final.f<-clogit(COVID ~ strata(CASE_NO)+ 
                  resp_x + heart_c + +asthma_c + renal+ neuro_c + diab12 + immuno_b +
                  emerg + inpat ,#emerg, icu.hdu.ccu,inpat,
                  data = cov[cov$training==1&cov$casegroup=="A",])
out <- summary(final.f)
outdf <- as.data.frame(out$coefficients)
outdf <- cbind(Factor=dimnames(outdf)[[1]], outdf)
outdf$lower <- exp(outdf$coef - (1.96*outdf$`se(coef)`))
outdf$upper <- exp(outdf$coef + (1.96*outdf$`se(coef)`))
tb <- outdf[c("Factor", "exp(coef)", "lower", "upper", "Pr(>|z|)")]
tb$"exp(coef)" <- round(tb$"exp(coef)", 2)
tb$lower <- round(tb$lower, 2)
tb$upper <- round(tb$upper, 2)
tb$`Pr(>|z|)` <- ifelse(tb$`Pr(>|z|)` < 0.001, "<0.001", round(tb$`Pr(>|z|)`, 3))
names(tb) <- c("Factor", "OR", "LCL", "UCL", "P")
tb

# Remove cohort specific data from saved model for export
modfs <- data.table::copy(modf)
modfs[c("residuals", "weights", "y", "fitted.values", "data", "prior.weights", "na.action", 
        "effects", "offset", "predictors", "linear.predictors", "model")] <- NULL
modfs$y = 0
summary(modfs)

save(modfs, file = paste0(lfp$data, "scottish.forward.model.RData"))







