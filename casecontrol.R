## analysis script for case-control study
if(exists("cl")) {
    # parallel::stopCluster(cl)
    showConnections(all = TRUE)
    closeAllConnections()
}
rm(list=ls())
gc()

## required system packages: libcurl4-openssl-dev, pandoc, pandoc-citeproc
##   libssl-dev libxml2-dev
## texlive-full or at least texlive-latex-extra, texlive-luatex

## install all required R packages
list.of.packages = c("car", 
                     "survival", 
                     "MASS", 
                     "wevid", 
                     "rmarkdown",
                     "bookdown",
                     "rticles",
                     "kableExtra",
                     "pander", 
                     "ggplot2", 
                     "doParallel", 
                     "readxl", 
                     "DescTools", 
                     "icd", 
                     "gam", 
                     "dplyr",
                     "data.table")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages) > 0) {install.packages(new.packages)}
lapply(list.of.packages, require, character.only=T)

library(car)
library(survival)
library(MASS)
library(wevid)
library(rmarkdown)
library(pander)
library(ggplot2)
library(doParallel)
library(readxl)
library(DescTools)
library(icd)
library(gam)
library(data.table)

Rprof(tmp <- tempfile())

source("helperfunctions.R")

old <- TRUE
old <- FALSE # uncomment to use latest linkage

#stepwise <- TRUE   
stepwise <- FALSE ## uncomment to save time if old version still valid

fatal.predict <- TRUE
fatal.predict <- FALSE

recode.dmtype <- function(x) {
    r <- as.factor(car::recode(x, 
                        "c(0, 10, 17, 97)='Not diabetic';
                         c(1, 101, 102)='Type 1 diabetes';
                         c(2, 202, 203)='Type 2 diabetes'; 
                         3:9='Other/unknown type';
                         11:16='Other/unknown type';
                         18:96='Other/unknown type';
                         98:100='Other/unknown type'"))
    r <- factor(r, levels=levels(r)[c(1, 3, 4, 2)])
    return(r)
}

icdToInt <- function(x) {
    N <- length(x)
    int3 <- integer(N)
    lastchar <- substr(x, 3, 3)
    lastchar.asc <- DescTools::CharToAsc(lastchar)
    lastchar.digit <- lastchar.asc >= 48 & lastchar.asc <= 57
    int1 <- as.integer(DescTools::CharToAsc(substr(x, 1, 1)))
    int2 <- as.integer(substr(x, 2, 2))
    int3[lastchar.digit] <- as.integer(lastchar[lastchar.digit]) + 10
                                        # integer range 10 to 19
    int3[!lastchar.digit] <- lastchar.asc[!lastchar.digit]
                                        # integer range 65 to 90
    1000 * int1 + 100 * int2 + int3
}

##################################################################

if(old) { # case-control dataset 18 July
    datadir <- "./data/"
    cc.filename <- paste0(datadir, "CC_linked_ANON_2020-06-18.rds")
    diagnoses.filename <- paste0(datadir, "CC_SMR01_ICD10_x25_2020-06-18.rds")
    procedures.filename <- paste0(datadir, "CC_SMR01_OPCS4_MAIN.x25_ANON_2020-06-18.rds")
    rapid.filename <- paste0(datadir, "CC_RAPID_ANON_2020-06-18.rds")
    chistatus.filename <- paste0(datadir, "CC_EXTENDED_STATUS_ANON_2020-06-18.rds")
    sicsag.filename <- paste0(datadir, "CC_SICSAG_ANON_2020-06-18.rds")
    scrips.filename <- paste0(datadir, "CC_PIS_x15_ANON_2020-06-18.rds")
    onomap <- as.data.table(readRDS(paste0(datadir, "ONOMAP_ANON_2020-06-18.rds")), key="ANON_ID")
} else { # case-control dataset 23 July
   datadir <- "./data/2020-07-09/"
   cc.filename <- paste0(datadir, "CC_linked_ANON_2020-07-23.rds")
   diagnoses.filename <- paste0(datadir, "CC_SMR01_ICD10_x25_ANON_2020-07-23.rds")
   procedures.filename <-  paste0(datadir, "CC_SMR01_OPCS4_MAIN.25_ANON_2020-07-23.rds")
   rapid.filename <- paste0(datadir, "CC_RAPID_ANON_2020-07-23.rds")
   chistatus.filename <- paste0(datadir, "CHI_Extended_status.rds")
   chistatus <- read.csv(paste0(datadir, "CHI_Extended_status.csv"))
   saveRDS(chistatus, file=chistatus.filename)
   sicsag.filename <- paste0(datadir, "CC_SICSAG_ANON_2020-07-23.rds")
   scrips.filename <- paste0(datadir, "CC_PIS_x15_ANON_2020-07-23.rds")

   smr04.filename <- paste0(datadir, "CC_SMR04__ANON_2020-07-28.rds")
   smr06.filename <- paste0(datadir, "CC_SMR06_ICD10_ANON_2020-07-23.rds")
#    "CC_PIS_15_ANON_2020-07-23.rds"
#    "CC_SMR01_ICD10_25_ANON_2020-07-23.rds"
    shielded <-
        as.data.table(read.csv("./data/update31Aug20/CC_shielding_patients_anon.csv"),
                      key="ANON_ID")
    setnames(shielded, "group", "shield.group")
    setnames(shielded, "batch", "shield.batch")
    setnames(shielded, "DuplicateFlag", "shield.dupflag")
    setnames(shielded, "Removal", "shield.removal")
   
}

batches <- read.table("data/BatchReference.csv", sep="\t", header=TRUE)
batches$Batch <- as.integer(gsub("Batch", "", batches$Batch))
batches$Date.Sent <- gsub("Friday ", "", batches$Date.Sent)
batches$Date.Sent <- gsub("^([[:digit:]]+)([stndrdth]{2} )", "\\1 ", batches$Date.Sent)
batches$Date.Sent <- as.Date(batches$Date.Sent, "%d %B %Y")
shielded <- merge(shielded, batches, all.x=TRUE, by.x="shield.batch", by.y="Batch")
setkey(shielded, ANON_ID)

cc.all <- as.data.table(readRDS(cc.filename), key="ANON_ID") # 203987 records

if(!old) { # include cases ascertained through discharge diagnosis
    ## case ascertainment uses three sources in order:
    ##   test positives, discharge diagnoses, death certs
    ## is.case appears to encode case ascertainment through test positivity
    ## diag.case is a 0-1 variable encoding cases ascertained through discharge diagnosis
    ## cod.case encodes ascertainment through death cert
    ## 77 of these have death cert with COVID as underlying cause
    with(cc.all, table(is.case, cod.case)) # 142 have is.case==0 and diag.case=1
    with(cc.all, table(is.case, diag.case)) # 14 have diag.case==1 and is.case==0
    with(cc.all[diag.case==1 & cod.case==0],
         table(icu + hdu > 0, !is.na(DATE_OF_DEATH))) # 4 cases ascertained through discharge diagnosis were not in critical care and have a death certificate with no mention of COVID.
    ## Sharon & Jen have assigned these a dummy specimen date 7 days before date of admission
    ## we reconstruct the date of admission as 7 days after the specimen date,
    ## then code as severe if the date of death was within 28 days of admission.
    cc.all[, death28daycovid := is.case==1 &
                                        # covid_cod==0 & icu + hdu == 0 &
                 !is.na(DATE_OF_DEATH) &
                 DATE_OF_DEATH - SPECIMENDATE <= 35]

    cc.all[, CASE := as.integer(is.case | diag.case | covid_ucod==1)]
    ## reconstruct ECOSS_POSITIVE as is.case excluding diag.case and cod.case
    cc.all[, ECOSS_POSITIVE := ifelse(is.case & diag.case==0 & cod.case==0,
                                      "Positive", "Negative/not done")]
}


#######################################################################
## cases-only table for diabetes paper -- this should be a separate script

casesonly <- as.data.table(readRDS("./data/2020-09-09/all_cases.rds"), key="ANON_ID")
                                        # latest specimen date 29 Sep 20
                                        # latest death date 31 July 20
sicsag.new <- as.data.table(readRDS("./data/2020-09-09/CC_SICSAG_ANON.rds"), key="ANON_ID") # latest hospital date 6 Sep 2020
casesonly$DATE_OF_DEATH <- as.Date(casesonly$DATE_OF_DEATH)
casesonly <- casesonly[SPECIMENDATE <= as.Date("2020-08-25") &
                       (is.na(DATE_OF_DEATH) |
                        DATE_OF_DEATH <= as.Date("2020-07-31")), ]
sicsag.new <- sicsag.new[AdmitUnit <= as.Date("2020-08-23"), ]
casesonly <- sicsag.new[casesonly]
setnames(casesonly, "AgeYear", "AGE")
setnames(casesonly, "SEX", "sex")
severe.cases <- casesonly[(!is.na(covidICUorHDU) & (covidICUorHDU == 1 | covidICUorHDU ==3)) |
                          dead28 == 1 |
                          covid_cod == 1, ]
fatal.cases <- casesonly[dead28 == 1 | covid_cod == 1, ]

case.freqs <- with(severe.cases[, .(AGE, sex, SPECIMENDATE)],
                   table(AGE, sex, exclude=NULL))
death.freqs <- with(fatal.cases[, .(AGE, sex, SPECIMENDATE)], table(AGE, sex, exclude=NULL))
severe.cases <- severe.cases[,  .(AGE, sex, SPECIMENDATE)]
fatal.cases <- fatal.cases[,  .(AGE, sex, SPECIMENDATE)]

save(severe.cases, fatal.cases, case.freqs, death.freqs,
     file="casefreqs31July.RData") ## national cumulative cases and deaths by sex and one year age group

################################################################
   
sicsag <- as.data.table(readRDS(sicsag.filename), key="ANON_ID")

## FIXME: discharge date missing from new diagnoses table

diagnoses <- as.data.table(readRDS(diagnoses.filename), key="ANON_ID")
diagnoses <- cc.all[, .(ANON_ID, SPECIMENDATE)][diagnoses] 

procedures <- readRDS(procedures.filename)
rapid <- as.data.table(readRDS(rapid.filename), key="ANON_ID")

cases.status <-
    as.data.table(readRDS(chistatus.filename), key="ANON_ID")
setnames(cases.status, "CHI_EXTENDED_STATUS", "EXTENDED_STATUS", skip_absent=TRUE)


scrips.firsttime <- TRUE
if(scrips.firsttime) { 
    scrips <- readRDS(scrips.filename)
    scrips[, daysbefore := as.integer(SPECIMENDATE - dispensed_date)]
    ## restrict to last 240 days before cutoff of (specimendate - 15) 
    scrips <- scrips[daysbefore - 15 <= 240]
    ## save this file so that other reads from it will have this restriction applied 
    saveRDS(scrips, paste0(datadir, "scrips.last240days.rds"))
}
scrips.filename <- paste0(datadir, "scrips.last240days.rds")
scrips <-
    as.data.table(readRDS(scrips.filename)[, c("ANON_ID",
                                               "dispensed_date",
                                               "approved_name",
                                               "bnf_item_code", 
                                               "bnf_paragraph_code",
                                               "bnf_paragraph_description",
                                               "daysbefore")], key="ANON_ID")
source("icdchapters.R")

###############################################

setnames(cc.all, "CASE_NO", "stratum", skip_absent=TRUE)
setnames(cc.all, "SEX", "sex", skip_absent=TRUE)
setnames(cc.all, "imumune", "immune", skip_absent=TRUE)
setnames(cc.all, "CAREHOME", "care.home", skip_absent=TRUE)
setnames(cc.all, "simd", "SIMD.quintile", skip_absent=TRUE)
setnames(cc.all, "DATE_OF_DEATH", "Date.Death", skip_absent=TRUE)
setnames(cc.all, "age", "AGE", skip_absent=TRUE)
setnames(cc.all, "AgeYear", "AGE", skip_absent=TRUE)
# setnames(cc.all, "AgeYear", "AGE") not sure why this was in twice

cc.all$CASE <- as.integer(cc.all$is.case)
cc.all$stratum <- as.integer(cc.all$stratum)

shielded <- shielded[ANON_ID %in% cc.all[, ANON_ID]]
cc.all <- shielded[cc.all]
setkey(cc.all, ANON_ID)
cc.all[is.na(shield.group), shield.group := 0]
cc.all[is.na(shield.batch), shield.batch := 0]
cc.all[, shield.group := as.factor(shield.group)]
cc.all[, shield.batch := as.factor(shield.batch)]
cc.all[, shield.any := as.factor(shield.group > 0)]
## pregnant with heart disease (code 6) grouped with Clinician identified
cc.all[, shield.group := car::recode(shield.group,
                                     "0='No shielding';
                                      1='Solid organ transplant';
                                      2='Specific cancers';
                                      3='Severe respiratory';
                                      4='Rare diseases';
                                      5='On immunosuppressants';
                                      6:7='Clinician identified'",
                                     as.factor=TRUE,
                                     levels=c(
                                         "No shielding",
                                         "Solid organ transplant",
                                         "Specific cancers",
                                         "Severe respiratory",
                                         "Rare diseases",
                                         "On immunosuppressants",
                                         "Clinician identified"
                                     ))]

## assign after.letter as 1 for those shielded if specimen date  > 10 days after letter 
cc.all[, after.letter := ifelse(Date.Sent + 10 < SPECIMENDATE, 1, 0)]
## assign after.letter as 0 for those not shielded
cc.all[is.na(Date.Sent), after.letter := 0]

###########################################

## merge SICSAG data and overwrite icu and hdu values incorrectly coded 0
sicsag <- sicsag[covidICUorHDU==1 | covidICUorHDU==3] ## what are codes 4 and 5?
sicsag <- cc.all[, .(ANON_ID, SPECIMENDATE)][sicsag]

## exclude entries to critical care more than 7 days before first positive test
sicsag <- sicsag[AdmitUnit - SPECIMENDATE >= -7]

## get last date of entry in table 
maxdate.sicsag <- max(sicsag$AdmitUnit)

## restrict to first SICSAG record for each ID
setkeyv(sicsag, c("ANON_ID", "AdmitUnit"))
sicsag <- sicsag[!duplicated(ANON_ID), c("ANON_ID", "AdmitUnit", "covidICUorHDU")]
setkey(sicsag, ANON_ID)
cc.all <- sicsag[cc.all]
cc.all[is.na(covidICUorHDU), covidICUorHDU := 0]

print(with(cc.all, table(covidICUorHDU, icu==1 | hdu==1)))
    ## overwrite icu or hdu field where incorrectly coded as 0
    ## this field covidICUorHDU is absent in new dataset 
cc.all[covidICUorHDU==1, icu := 1]
cc.all[covidICUorHDU==3, hdu := 1]

cc.all[, dispensing.days := as.integer(SPECIMENDATE - as.Date("2019-06-01"))]

## HAI is based on the ECDC definition of nosocomial infection

paste.colpercent(with(cc.all, table(ANON_ID %in% scrips$ANON_ID, CASE)))
paste.colpercent(with(cc.all, table(ANON_ID %in% diagnoses$ANON_ID, CASE)))

cc.all[, scrip.any := as.factor(as.integer(ANON_ID %in% scrips$ANON_ID))]
cc.all[, diag.any := as.factor(as.integer(ANON_ID %in% diagnoses$ANON_ID))]
cc.all[, scripordiag := as.factor(as.integer(diag.any=="1" | scrip.any=="1"))]

ids.noscrip.agegt75 <- subset(cc.all, CASE==0 & scrip.any==0 & AGE > 75)[["ANON_ID"]][1:20]
ids.nodiag.agegt75 <- subset(cc.all, CASE==0 & diag.any==0 & AGE > 75)[["ANON_ID"]][1:20]

save(ids.noscrip.agegt75, ids.nodiag.agegt75, file="anon_ids_notmatched.RData")

## exclude controls already dead on date of test of case they were matched to
controls.deceased <- with(cc.all, CASE==0 &
                                  !is.na(Date.Death) &
                                  Date.Death <= SPECIMENDATE)
cc.all <- cc.all[!controls.deceased]                         

## exclude cases classified as unobservable, for consistency with controls 
cc.all <- cases.status[cc.all]
if(old) {
    with(cc.all[CASE==1],
         table(Explanation, is.na(Date.Death)))
}
with(cc.all[CASE==1],
     table(EXTENDED_STATUS, is.na(Date.Death)))
cc.all <- subset(cc.all,
                 is.na(EXTENDED_STATUS) |
                 EXTENDED_STATUS == "C" |
                 (EXTENDED_STATUS == "D" & !is.na(Date.Death)))

## rows have been replicated within strata containing N cases so that each case and control appears N^2 times
## remove these replicated rows
cc.all <- distinct(.data=cc.all, ANON_ID, stratum, .keep_all = TRUE)
setkey(cc.all, ANON_ID)

numcases.strata <- tapply(cc.all$CASE, cc.all$stratum, sum) 
numctrls.strata <- tapply(1 - cc.all$CASE, cc.all$stratum, sum)
print(table(numctrls.strata, numcases.strata))

cc.all[, SIMD.quintile := as.factor(car::recode(SIMD.quintile, "'Unknown'=NA"))]

cc.all[, sex := car::recode(as.factor(sex), "1='Male'; 2='Female'")]
cc.all[, sex := factor(sex, levels=c("Female", "Male"))]

## generates a warning about taking a shallow copy
cc.all[, agegr20 := as.factor(car::recode(as.integer(AGE),
                              "0:39='0-39'; 40:59='40-59';
                               60:74='60-74'; 75:hi='75 or more'"))]
cc.all[, agegr3 :=
    as.factor(car::recode(AGE,
                          "0:59='0-60 years'; 60:74='60-74 years'; 75:hi='75+ years'"))]
cc.all[, agegr2 :=
    as.factor(car::recode(AGE,
                          "0:74='0-74 years'; 75:hi='75+ years'"))]

cc.all[NURSINGHOME==1, care.home := 1]
cc.all[, care.home := as.factor(car::recode(care.home, "0='Independent';
                                                        1='Care/nursing home'"))]
cc.all[, care.home := relevel(care.home, ref="Independent")]

with(cc.all[CASE==1], table(ECOSS_POSITIVE, covid_cod, exclude=NULL))
## 92 individuals are included as cases but do not have covid on death cert and their test result is negative or empty string

## relabel empty string as "No result"
cc.all[is.na(ECOSS_POSITIVE) | ECOSS_POSITIVE=="", ECOSS_POSITIVE := "No result"] 

## generate a factor variable with 3 levels
cc.all[, ecoss := as.factor(ECOSS_POSITIVE)]

## generate testpositive.case based on ecoss
cc.all[, testpositive.case := CASE==1 & ecoss=="Positive"]

## set testpositive case to TRUE for all cases with covid_cod=0
## they couldn't have been ascertained except through testing positive
cc.all[CASE==1 & covid_cod==0, testpositive.case := TRUE]

## now all cases without covid_cod are test-positive
with(cc.all[CASE==1], table(covid_cod, testpositive.case, exclude=NULL))
with(cc.all[CASE==1], table(ecoss, testpositive.case, exclude=NULL))

## all cases have nonmissing SPECIMENDATE
## controls should be assigned same SPECIMENDATE as the case they were matched to

## deathwithin28 is a logical variable that literally means death within 28 days of a positive test
## coded as FALSE for those who did not test positive
cc.all[, deathwithin28 := testpositive.case &
             !is.na(Date.Death) & 
             Date.Death - SPECIMENDATE >= 0 & Date.Death - SPECIMENDATE <= 28]

## Sharon's variable dead28 is assigned by this line
## Covid_CC_linkage_Part2_desktop.R:
## cc$dead28 <- ifelse(!is.na(cc$DATE_OF_DEATH) & cc$DATE_OF_DEATH >= cc$SPECIMENDATE & as.numeric(cc$DATE_OF_DEATH - cc$SPECIMENDATE) <=28, 1, 0)
## this evaluates to 1 for anyone classified as a case who dies within 28 days of specimen
## date even if this is a dummy specimen date.

with(cc.all[CASE==1], table(dead28, deathwithin28, exclude=NULL)) 
with(cc.all[CASE==1], table(icu, hdu, exclude=NULL)) 

cc.all[, criticalcare := icu==1 | hdu==1]

## adm28 should be 0 for cases who did not test positive
cc.all[, adm28 := adm28 & testpositive.case]

with(cc.all[CASE==1], table(inhosp, adm28, exclude=NULL))
## Covid_CC_linkage_Part2_desktop.R:rapid$days <- as.numeric(rapid$Admission.Date - rapid$SPECIMENDATE)
## Covid_CC_linkage_Part2_desktop.R:rapid$adm28 <- ifelse(rapid$days %in% c(0:28), 1, 0)

## new linkage does not have the variable nrs_covid_case
with(cc.all[CASE==1 & covid_ucod==1], table(deathwithin28, exclude=NULL))
## 3701 cases with covid_cod have deathwithin28==1
## all those with covid_ucod==1 have covid_cod==1 

with(cc.all[CASE==1 & testpositive.case], table(criticalcare, deathwithin28, exclude=NULL))
with(cc.all[CASE==1], table(criticalcare, testpositive.case, exclude=NULL))
with(cc.all[CASE==1], table(testpositive.case, deathwithin28, exclude=NULL))

## coding of case groups
cc.all[, group := "Unclassified"]

## define group A (severe cases)
## group A includes anyone certified with COVID as underlying cause
    cc.all[CASE==1 & (criticalcare | deathwithin28 | covid_ucod==1),
           group := "A"]

## assign remaining test-positive cases hospitalized within 28 days of positive test as group B
cc.all[testpositive.case & group=="Unclassified" & (adm28 > 0 | inhosp > 0), group := "B"]
## assign all remaining test-positive cases to group C
cc.all[testpositive.case & group=="Unclassified", group := "C"]
## assign remaining cases with mention on death cert to group D
cc.all[CASE==1 & group=="Unclassified" & covid_cod==1, group := "D"]

print(table(cc.all$CASE, cc.all$group, exclude=NULL))

## assign a logical variable for fatalcase, and a binary variable for fatal.casegroup
cc.all$fatalcase <- with(cc.all, CASE==1 & group=="A" & (deathwithin28 |
                                                         covid_ucod==1 |
                                                         death28daycovid))

cc.all[, fatal.casegroup := 0]
cc.all[CASE==1, fatal.casegroup := as.integer(fatalcase)]

## five categories of severe case
cc.all[, severe.casegr := 0]
cc.all[CASE==1 & group=="A" & criticalcare & !fatalcase,
       severe.casegr := 1]
cc.all[CASE==1 & group=="A" & testpositive.case & criticalcare & fatalcase, 
       severe.casegr := 2]
cc.all[CASE==1 & group=="A" & testpositive.case & !criticalcare & fatalcase,
       severe.casegr := 3]
cc.all[CASE==1 & group=="A" & !criticalcare & !testpositive.case & fatalcase,
       severe.casegr := 4]
cc.all[severe.casegr==0, severe.casegr := NA]
cc.all[, severe.casegr := factor(severe.casegr,
                                 labels=c("Critical care, non-fatal",
                                          "Critical care, fatal",
                                          "No critical care, test-positive, fatal",
                                          "No critical care, not test-positive, fatal"))]
print(table(cc.all$severe.casegr))

## assign controls in each stratum to same casegroup as case
## also assign controls in each stratum to same after.letter category as case
## for strata with 2 or more cases, assign all controls to highest group of case
casegroups <- subset(cc.all, subset=CASE==1, select=c("stratum", "group", "fatal.casegroup"))
colnames(casegroups)[2] <- "casegroup"
setkey(casegroups, casegroup) # orders by casegroup
casegroups <- casegroups[!duplicated(casegroups$stratum), ]
setkey(casegroups, stratum)
setkey(cc.all, stratum)
cc.all <- casegroups[cc.all]

## for cases, overwrite the casegroup field with the group assigned above
cc.all[CASE==1, casegroup := group]
table(cc.all$CASE, cc.all$casegroup, exclude=NULL)

## drop records with missing casegroup (controls in strata with no remaining classified case)
cc.all <- cc.all[!is.na(casegroup)]

with(cc.all[CASE==1], table(casegroup, deathwithin28, exclude=NULL))
print(paste.colpercent(with(cc.all[CASE==1], table(care.home, casegroup))))

### incidence and mortality using national population estimates #####
narrow.case <- cc.all$CASE==1 & cc.all$casegroup=="A"
narrow.fatalcase <- narrow.case & (cc.all$deathwithin28 | cc.all$covid_ucod==1)
table(narrow.case, narrow.fatalcase)

## broad.case adds in all those with mention of covid on death cert
broad.case <- cc.all$CASE==1 & (cc.all$casegroup=="A" | cc.all$casegroup=="D")
broad.fatalcase <- broad.case & (cc.all$deathwithin28 | cc.all$covid_cod==1)
table(broad.case, broad.fatalcase)

## for graphing, we want three categories
## severe cases as defined
## broad cases as used for diabetes report: severe, hospitalized or NRS
## fatal cases including NRS deaths
narrow.case.freqs <- with(cc.all[narrow.case], table(AGE, sex, exclude=NULL))
broad.case.freqs <- with(cc.all[broad.case], table(AGE, sex, exclude=NULL))
narrow.death.freqs <- with(cc.all[narrow.fatalcase], table(AGE, sex, exclude=NULL))
broad.death.freqs <- with(cc.all[broad.fatalcase], table(AGE, sex, exclude=NULL))

save(narrow.case.freqs, broad.case.freqs,
     narrow.death.freqs, broad.death.freqs, file="casefreqs.4cats.agesex.RData")

source("incidencemortality.R")

######## coding ethnicity ##############################

source("ethnic_assign.R")

cc.all[, ethnic5.smr := collapseto5.ethnicsmr(ETHNIC_SMR_LAST)]
## recode SMR ethnicity to 4 categories: White, Black, South Asian, Other
cc.all[, ethnic4.smr := as.factor(car::recode(ethnic5.smr,
                                              "'Chinese'='Other'"))]
cc.all[, ethnic4.smr := factor(ethnic4.smr,
                               levels=levels(ethnic4.smr)[c(4, 3, 1, 2)])]   

if(old) {
## ANON_ID may be duplicated in cc.all but not in onomap
    ## some IDs in onomap are not in cc.all
    onomap <- onomap[ANON_ID %in% cc.all$ANON_ID]
    setkey(onomap, ANON_ID)
    setkey(cc.all, ANON_ID)
    cc.all <- onomap[cc.all]
    
    cc.all[, group.onomap := group.onomap(OnolyticsType, GeographicalArea)]
    cc.all[, ethnic5.onomap := collapseto5.onomap.group(group.onomap)]
    
    ## tabulate ONOMAP ethnicity against SMR ethnicity
    table.ethnic <- table(cc.all$ethnic5.onomap, cc.all$ethnic5.smr, exclude=NULL)
    tn <- as.data.frame.matrix(table(cc.all$ethnic5.onomap, cc.all$ethnic5.smr))
    SouthAsian.sensitivity <- 100 * tn["South Asian", "South Asian"] / sum(tn[, "South Asian"])
    SA.trueneg <- sum(tn) - sum(tn[, "South Asian"])
    SA.trueneg.correct <- SA.trueneg - (sum(tn[, "South Asian"]) - tn["South Asian", "South Asian"]) 
    SouthAsian.specificity <- 100 * SA.trueneg.correct / SA.trueneg
    sum.xtabulate <- sum(tn)

    ## recode ONOMAP ethnicity to 4 categories: White, South Asian, Chinese, Other
    cc.all[, ethnic4.onomap := car::recode(ethnic5.onomap, "'Black'='Other'")]
    
    cc.all[, ethnic4.onomap := factor(ethnic4.onomap,
                                      levels=levels(ethnic4.onomap)[c(4, 3, 1, 2)])]
    
    ## recode to 3 categories: White, South Asian, Other
    cc.all[, ethnic3.onomap := car::recode(ethnic4.onomap, "'Chinese'='Other'")]
    cc.all[, ethnic3.onomap := factor(ethnic3.onomap,
                                      levels=levels(ethnic3.onomap)[c(3, 2, 1)])]
}

setkey(cc.all, ANON_ID)
######################################################
                 
source("bnfcodes.R")

source("drugs.R")

###########################################################################

###  diabetes based on combining Sci-Diabetes records with ICD codes and drugs 
## add in BNF codes 6.1 for diabetes drugs and
## E10 to E14 for diabetes

########## coding listed conditions ####################

ids.icd.diabetes <- unique(diagnoses$ANON_ID[grep("^E1[0-4]", diagnoses$ICD10)])
ids.bnf.diabetes <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 601])
ids.diabetes.extra <- unique(c(ids.icd.diabetes, ids.bnf.diabetes))

## missing recoded as zero
cc.all[is.na(dm.type), dm.type := 0]

## add in extra cases notified directly from SCI-Diabetes register, without assignment
## of diabetes type from SDRN database
cc.all[dm.type==0 & diab.reg==1, dm.type := 3]

#cat("Extra diabetes cases from SCI-Diabetes by diabetes type\n")
#print(table(cc.all$dm.type, cc.all$ANON_ID %in% ids.diabetes.extra))

## code diagnoses detected from discharges or BNF codes as unknown type
## we could classify those not on insulin as definite Type 2 but Helen says no
## REVISION: for consistency with the diabetes paper, we will not include the extra cases identified through diagnostic codes or drug codes - they may be transient/resolved

## cc.all$dm.type[cc.all$dm.type==0 & cc.all$ANON_ID %in% ids.diabetes.extra] <- 3

# recode diabetes type
cc.all[, dm.type := recode.dmtype(dm.type)]

## define an indicator variable for any diabetes
cc.all[, diabetes.any := as.integer(dm.type != "Not diabetic")]
cc.all[, diabetes.any := as.factor(car::recode(diabetes.any,
                                               "0='Not diabetic'; 1='Diabetic'"))]
cc.all[, diabetes.any := relevel(diabetes.any, ref="Not diabetic")]


 
objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

########## restrict to severe cases and matched controls ###################### 

cat("Restricting to severe cases and matched controls\n")
cc.severe <- cc.all[casegroup=="A"]

# if(!old)  .Internal(.invokeRestart(list(NULL, NULL), NULL))

#####################################################################

## use overlap join to assign chapter and subchapter to ICD diagnoses

diagnoses[, icdnum := icdToInt(ICD10)]
diagnoses[, icdnum2 := icdnum]
setkey(diagnoses, icdnum, icdnum2)
setkey(icd.subchapters, startnum, endnum)
diagnoses <- foverlaps(diagnoses, icd.subchapters[, .(chnum, subchnum, startnum, endnum)])

if(!old) { # hack to assign discharge date
    diagnoses[, DISCHARGE_DATE := as.Date("2019-12-31")]
}
    
diagnoses <- diagnoses[, .(ANON_ID, SPECIMENDATE, DISCHARGE_DATE,
                           ICD10, chnum, subchnum)]
    ids.icd.neoplasm.lastyear <-
    unique(diagnoses[as.integer(SPECIMENDATE - DISCHARGE_DATE) + 25 > 365 &
                     grepl("^C[0-9]|^D[0-4]", ICD10), ANON_ID])


########################################

## cast chapters in wide format and merge with cc.severe
icdchapters.wide <- data.table::dcast(diagnoses, ANON_ID ~ chnum, fun.aggregate=length,
                                      value.var="chnum")

## use Roman numerals for ICD chapters
chnums <- as.integer(colnames(icdchapters.wide)[-1])
colnames(icdchapters.wide)[-1] <-
    paste0("Ch.", as.roman(chnums), "_",
           icdchapters$shortname[match(chnums, icdchapters$chnum)])
icdchapters.wide <- as.data.table(icdchapters.wide, key="ANON_ID")

## drop rare chapters
cols.keep <- colSums(icdchapters.wide) >= 20
icdchapters.wide <- icdchapters.wide[, ..cols.keep]

cc.severe <- icdchapters.wide[cc.severe]
rm(icdchapters.wide)

icdcols <- grep("^Ch\\.", colnames(cc.severe), value=TRUE)

## recode as indicator variables: NA to 0, >1 to 1
cc.severe[, (icdcols) := lapply(.SD, recode.indicator), .SDcols = icdcols]

## cc.severe[, (icdcols)] returns the column names as a vector
## cc.severe[, .(icdcols)] returns the column names as a data.table
## cc.severe[, ..icdcols] returns a data.table selected by column

## this line returns a warning that "Both 'icdcols' and '..icdcols' exist in calling scope."
cc.severe[, num.icdchapters := rowSums(matrix(as.integer(as.matrix(cc.severe[, ..icdcols])),
                                              nrow=nrow(cc.severe)))]

cc.severe[, num.icdchapters.gr := as.factor(car::recode(num.icdchapters,
                          "0='No discharge records'; 1:2='1-2 ICD-10 chapters';
                           3:hi='3 or more chapters'"))]
cc.severe[, num.icdchapters.gr := relevel(num.icdchapters.gr, ref="No discharge records")]

## cast subchapters in wide format and merge with cc.severe
icdsubchapters.wide <- data.table::dcast(diagnoses, ANON_ID ~ subchnum, fun.aggregate=length,
                                         value.var="subchnum")
## column names are integers as character strings 
subchnums <- as.integer(colnames(icdsubchapters.wide)[-1])
colnames(icdsubchapters.wide)[-1] <-
    paste0("Ch_", 
           as.roman(icd.subchapters$chnum[match(subchnums, icd.subchapters$subchnum)]), ": ",
           icdsubchapters$start[match(subchnums, icd.subchapters$subchnum)], "-",
           icdsubchapters$end[match(subchnums, icd.subchapters$subchnum)], " ", 
           icdsubchapters$name[match(subchnums, icd.subchapters$subchnum)])
icdsubchapters.wide <- as.data.table(icdsubchapters.wide, key="ANON_ID")

## drop rare subchapters
cols.keep <- colSums(icdsubchapters.wide) >= 20
icdsubchapters.wide <- icdsubchapters.wide[, ..cols.keep]

cc.severe <- icdsubchapters.wide[cc.severe]
rm(icdsubchapters.wide)

## recode as indicators
icdcols <- grep("^Ch_", colnames(cc.severe), value=TRUE)
cc.severe[, (icdcols) := lapply(.SD, recode.indicator), .SDcols = icdcols]

###########################################

source("comorbidity.R")

## immune.any includes primary immunodeficiency and secondary immunosuppression

## 8 listed conditions designated by NHS
listed.conditions <- c("dm.type", "IHD.any", "heart.other.any", "oad.any",
                       "ckd.any", "neuro.any", "liver.any", "immune.any")

############ extract predefined disease categories #################
## as these are coded as factors, lowest level will be 1

cc.severe[, listed.any := 
    as.factor(as.integer(diabetes.any=="Diabetic" | IHD.any==1 |
                         heart.other.any==1 |
                         ckd.any==1 | oad.any==1 |
                         neuro.any==1 | liver.any==1 | immune.any==1))]

cc.severe[, diag.other := as.integer(listed.any==0 & diag.any==1)]
cc.severe[, diag.other := as.factor(car::recode(diag.other,
                                  "0='Listed condition or no admission';
                                   1='No listed condition, but other admission diagnosis'"))]

####################################################################

source("paper1tables.R")

######## stepwise regressions use saved version #####################

nfold <- 10
source("stepwise.R")

#####################################################


if(old) {
    source("pharmaco.R")
    rmarkdown::render("casecontrol.Rmd",
                      output_file="casecontrol_PMEDICINE_D20_024838June.pdf")
    rmarkdown::render("pharmaco.Rmd", output_file="pharmaco8June.pdf")
    saveRDS(cc.severe, file="data/cc.severe.rds")
} else {
    source("shielding.R")
}

## remove large objects from memory

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

#rmarkdown::render("casecontrol_script.Rmd", output_file="casecontrol_script.pdf")

Rprof()
print(summaryRprof(tmp)$by.total[1:20, ])


