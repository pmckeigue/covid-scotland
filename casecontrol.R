## analysis script for case-control study
if(exists("cl")) {
    # parallel::stopCluster(cl)
    showConnections(all = TRUE)
    closeAllConnections()
}
rm(list=ls())
gc()

linkdate <- "jun18"
linkdate <- "jan6"
linkdate <- "jan28"

controls <- TRUE
shielding <- TRUE
sicsag <- TRUE
pis <- TRUE

## stepwise <- TRUE   
stepwise <- FALSE ## uncomment to save time if old version still valid
fatal.predict <- TRUE
fatal.predict <- FALSE
numdays.letter <- 14

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


add.index <- function(dt)  { # why does this return NA values for ANON_ID? 
    dt$case.index <- 1:nrow(dt)
    dt$case.interval <- c(NA, as.integer(diff(dt$SPECIMENDATE)))
    return(dt[, .(ANON_ID, case.index, case.interval)])
}

year.to2020 <- function(x) {
    lubridate::year(x) <- 2020
    return(x)
}

year.to2021 <- function(x) {
    lubridate::year(x) <- 2021
    return(x)
}

cleanID <- function(filename) {
    namestring <- eval(parse(text=filename))
    if(grepl("\\.csv$", namestring)) {
        x <- fread(namestring)
    } else {
        x <- as.data.table(readRDS(namestring))
    }
    setnames(x, "anon_id", "ANON_ID", skip_absent=TRUE)
    setnames(x, "CC_ANON_ID", "ANON_ID", skip_absent=TRUE)
    x[, ANON_ID := as.integer(gsub(".* ", "", ANON_ID))] # remove date prefix
    namesstring <- gsub("\\.csv$", "\\.rds", namestring) # save with .rds extension
    saveRDS(x, file=namestring)
}

checkID <- function(rdsfilename) {
    namestring <- eval(parse(text=rdsfilename))
    if(!is.null(namestring)) {
        x <- readRDS(namestring)
        setnames(x, "anon_id", "ANON_ID", skip_absent=TRUE)
        setnames(x, "CC_ANON_ID", "ANON_ID", skip_absent=TRUE)
        cat(rdsfilename, x[["ANON_ID"]][1], "\n")
    }
}

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

##################################################################
batches <- as.data.table(read.table("data/BatchReference.csv", sep="\t", header=TRUE))
setnames(batches, "Batch", "shield.batch")
batches$shield.batch <- as.integer(gsub("Batch", "", batches$shield.batch))
setkey(batches, shield.batch)
batches$Date.Sent <- gsub("Friday ", "", batches$Date.Sent)
batches$Date.Sent <- gsub("^([[:digit:]]+)([stndrdth]{2} )", "\\1 ", batches$Date.Sent)
batches$Date.Sent <- as.Date(batches$Date.Sent, "%d %B %Y")
 
if(linkdate =="jun18") {
    sicsag <- FALSE
    datadir <- "./data/jun18"
    cc.filename <- paste0(datadir, "CC_linked_ANON_2020-06-18.rds")
    diagnoses.filename <- paste0(datadir, "CC_SMR01_ICD10_x25_2020-06-18.rds")
    procedures.filename <- paste0(datadir, "CC_SMR01_OPCS4_MAIN.x25_ANON_2020-06-18.rds")
    rapid.filename <- paste0(datadir, "CC_RAPID_ANON_2020-06-18.rds")
    chistatus.filename <- paste0(datadir, "CC_EXTENDED_STATUS_ANON_2020-06-18.rds")
    sicsag.filename <- paste0(datadir, "CC_SICSAG_ANON_2020-06-18.rds")
    scrips.filename <- paste0(datadir, "CC_PIS_x15_ANON_2020-06-18.rds")
    onomap <- RDStodt(paste0(datadir, "ONOMAP_ANON_2020-06-18.rds"), keyname="ANON_ID")
} else if(linkdate == "jan6") {
    cleanid <- FALSE
    checkid <- FALSE
    occup <- TRUE
    datadir <- "./data/2021-01-06/"
    ct.filename <-  paste0(datadir, "CT_Test_Results_2021-01-21.rds")
    ecoss.tests.filename <- paste0(datadir, "ecoss_tests_2021-01-14.rds")
    ecoss.cases.filename <- paste0(datadir, "ecoss_cases_2021-01-14.rds")
    teach.filename <- paste0(datadir, "teacher_hcw.rds")

    hcw <- RDStodt(teach.filename, keyname="ANON_ID")
    hcw[hcw_hhd=="teacher", occup := "Teacher"]
    hcw[hcw_hhd=="hcw", occup := role]
    hcw[is.na(occup), occup := "Other / undetermined"]
    hcw[, occup := car::recode(occup, "c('npf', 'undetermined')='Health care, not PF / undetermined'; 'pf_any'='Health care PF'")]

    rapid.filename <- paste0(datadir, "CC_RAPID_ANON_2021-01-07.rds")
    sicsag.filename <- paste0(datadir, "CC_SICSAG_ANON_2021-01-21.rds")
    cc.filename <- paste0(datadir, "CC_linked_ANON_2021-01-21.rds")
    diagnoses.filename <- paste0(datadir, "CC_SMR01_ICD10_x25_ANON_2021-01-07.rds")
    diagnoses.last25.filename <- paste0(datadir, "CC_SMR01_ICD10_25_ANON_2021-01-07.rds")
    shielding.full.filename <- paste0(datadir,
                                      "CC_shielding_patients_anon_20210106.rds")
    smr00.filename <- paste0(datadir, "CC_SMR00__ANON_2021-01-07.rds")
    smr04.filename <- paste0(datadir, "CC_SMR04__ANON_2021-01-07.rds")
    smr06.filename <- paste0(datadir, "CC_SMR06_ICD10_ANON_2021-01-07.rds")
    scrips.filename <- paste0(datadir, "CC_PIS_x15_ANON_2021-01-12.rds")
    scrips.last15.filename <- paste0(datadir, "CC_PIS_15_ANON_2021-01-12.rds")
    procedures.filename <- paste0(datadir, "CC_SMR01_OPCS4_MAIN.x25_ANON_2021-01-07.rds")
    procedures.last25.filename <-
        paste0(datadir, "CC_SMR01_OPCS4_MAIN.25_ANON_2021-01-07.rds")
    ## this is a character vector containing the names of the objects, not the filenames
    raw.filenames <- grep("\\.filename", objects(), value=TRUE)
    raw.filenames <- raw.filenames[raw.filenames != "scrips.filename"]
    if(cleanid) {
        for(f in raw.filenames) {
            cleanID(f)
        }
    } else if(checkid) {
        for(f in raw.filenames) {
            checkID(f)
        }
    }
    chistatus.filename <- NULL
} else if(linkdate == "jan28") {
    cleanid <- FALSE
    checkid <- FALSE
    controls <- TRUE
    shielding <- TRUE
    occup <- TRUE
    chistatus.filename <- NULL
    datadir <- "./data/2021-01-28/"
#    old.datadir <- "./data/old2021-01-28/"
    ## new ct file is missing the flag_symptomatic_covid field
    ct.filename <-  paste0(datadir, "CT_Test_Results_2021-01-28.rds")
    ecoss.tests.filename <- paste0(datadir, "ecoss_tests_2021-01-28.rds")
    if(occup) {
        hcw_static.fname <- "./data/HCW_data_cleaned/static_hcw_hhd.Rds"
        hcw.filename <- paste0(datadir, "hcw_anon_ids.rds")
        ## this has a field hcw_hhd_anon which is the ID used in the source HCW file
        hcw.static <- RDStodt(hcw_static.fname, keyname="anon_id")
        hcw.static <- hcw.static[!is.na(role)]
        hcw <- RDStodt(hcw.filename)
        hcw <- hcw[!is.na(hcw_hhd_anon) & hcw_hhd_anon != "missing_anon_id" & !is.na(ANON_ID)]
        setkey(hcw, hcw_hhd_anon)  ## have to set this in a separate step
        table(hcw.static$anon_id %in% hcw$hcw_hhd_anon)
        hcw <- hcw.static[hcw]
        rm(hcw.static)
        hcw <- hcw[, .(ANON_ID, role, role_sub)]
        hcw <- hcw[!is.na(role)]
        setkey(hcw, ANON_ID)
   }
    rapid.filename <- paste0(datadir, "CC_RAPID_ANON_2021-01-28.rds")
    sicsag.filename <- paste0(datadir, "CC_SICSAG_ANON_2021-01-28.rds")
    cc.filename <- paste0(datadir, "CC_linked_ANON_2021-01-28.rds")
    diagnoses.filename <- paste0(datadir, "CC_SMR01_ICD10_x25_ANON_2021-01-28.rds")
    diagnoses.last25.filename <- paste0(datadir, "CC_SMR01_ICD10_25_ANON_2021-01-28.rds")
    shielding.full.filename <- paste0(datadir,
                                      "CC_shielding_patients_anon_20210128.rds")
    smr00.filename <- paste0(datadir, "CC_SMR00__ANON_2021-01-28.rds")
    smr04.filename <- paste0(datadir, "CC_SMR04__ANON_2021-01-28.rds")
    smr06.filename <- paste0(datadir, "CC_SMR06_ICD10_ANON_2021-01-28.rds")
    scrips.filename <- paste0(datadir, "CC_PIS_x15_ANON_2021-01-28.rds")
    scrips.last15.filename <- paste0(datadir, "CC_PIS_15_ANON_2021-01-28.rds")
    procedures.filename <- paste0(datadir, "CC_SMR01_OPCS4_MAIN.x25_ANON_2021-01-28.rds")
    procedures.last25.filename <-
        paste0(datadir, "CC_SMR01_OPCS4_MAIN.25_ANON_2021-01-28.rds")
    ## this is a character vector containing the names of the objects, not the filenames
    raw.filenames <- grep("\\.filename", objects(), value=TRUE)
    raw.filenames <- raw.filenames[raw.filenames != "scrips.filename"]
    if(cleanid) {
        for(f in raw.filenames) {
            cleanID(f)
        }
    } else if(checkid) {
        for(f in raw.filenames) {
            checkID(f)
        }
    }
}

cc.all <- RDStodt(cc.filename, keyname="ANON_ID")

CoDvars <- grep("CAUSE_OF_DEATH_CODE_", names(cc.all), value=TRUE)
cc.all[, (CoDvars) := NULL]
print(summary(cc.all$SPECIMENDATE))

## replace years before 2020 with 2020
cc.all[lubridate::year(SPECIMENDATE) < 2020, SPECIMENDATE := year.to2020(SPECIMENDATE)]
## set year to 2021 for dates before Feb 2020
cc.all[lubridate::year(SPECIMENDATE) == 2020 & SPECIMENDATE < as.Date("2020-02-01"),
       SPECIMENDATE := year.to2021(SPECIMENDATE)]
## set year to 2020 for dates after system date
cc.all[SPECIMENDATE > Sys.Date(), SPECIMENDATE := year.to2020(SPECIMENDATE)]

lastdate.specimen <- max(cc.all$SPECIMENDATE)
cat("Last specimen date", format(lastdate.specimen,"%d %B %Y"), "\n")

setnames(cc.all, "CAREHOME", "care.home", skip_absent=TRUE)

cc.all[is.na(INSTITUTION_CODE==93) |
       (INSTITUTION_CODE!=93 & INSTITUTION_CODE!=98), care.home := 0]
cc.all[!is.na(INSTITUTION_CODE) & (INSTITUTION_CODE==93 | INSTITUTION_CODE==98),
       care.home := 1]

cc.all[, care.home := as.factor(car::recode(care.home, "0='Independent';
                                                        1='Care/nursing home'"))]
cc.all[, care.home := relevel(care.home, ref="Independent")]

setnames(cc.all, "CASE_NO", "stratum", skip_absent=TRUE)
setnames(cc.all, "SEX", "sex", skip_absent=TRUE)
setnames(cc.all, "imumune", "immune", skip_absent=TRUE)
setnames(cc.all, "simd", "SIMD.quintile", skip_absent=TRUE)
setnames(cc.all, "DATE_OF_DEATH", "Date.Death", skip_absent=TRUE)
setnames(cc.all, "age", "AGE", skip_absent=TRUE)
setnames(cc.all, "AgeYear", "AGE", skip_absent=TRUE)
setnames(cc.all, "simd2020_sc_quintile", "qSIMD.integer", skip_absent=TRUE)

if(!controls) { # cases only
    cc.all[, CASE := 1]
} else { # case-control
    cc.all[, CASE := as.integer(is.case)]
    cc.all[, stratum := as.integer(stratum)]
}
cc.all[CASE==1, testpositive.case := diag.case==0 & cod.case==0]

## include cases ascertained through discharge diagnosis
    ## case ascertainment uses three sources in order:
    ##   test positives, discharge diagnoses, death certs
    ## is.case is defined to assign case status
    ## diag.case is a 0-1 variable encoding cases ascertained through discharge diagnosis
    ## cod.case encodes ascertainment through death cert
with(cc.all, table(is.case, cod.case)) # 137 have is.case==false and cod.case==1
with(cc.all, table(is.case, diag.case)) # 74 have is.fase==false and diag.case==1 and is.case==0
with(cc.all[diag.case==1 & cod.case==0],
     table(icu + hdu > 0, !is.na(Date.Death))) # 4 cases ascertained through discharge diagnosis were not in critical care and have a death certificate with no mention of COVID.
## Sharon & Jen have assigned these a dummy specimen date 7 days before date of admission
    ## we reconstruct the date of admission as 7 days after the specimen date,
    ## then code as severe if the date of death was within 28 days of admission.

cc.all[AGE < 0, AGE := 0]

maxdate.specimen <- max(cc.all[CASE==1, SPECIMENDATE])

cc.all[is.na(Prisoner_Flag), Prisoner_Flag := FALSE]
## all cases not ascertained through diag or cod must be testpositive cases
## 3052 cases not testpositive cases
cc.all[, region := factor(car::recode(hb2019name,
   recode="c('NHS Borders', 'NHS Fife', 'NHS Lothian', 'NHS Tayside', 'NHS Forth Valley') ='Eastern';
   'NHS Greater Glasgow and Clyde'='Glasgow';
   'NHS Grampian' = 'North Eastern'; 
 c('NHS Ayrshire and Arran', 'NHS Dumfries and Galloway', 'NHS Lanarkshire') = 'South Western';
 c('NHS Western Isles', 'NHS Shetland', 'NHS Highland', 'NHS Orkney') = 'Highlands and Islands'"))]
cc.all[, region3 := car::recode(region,
                                 "c('North Eastern', 'Highlands and Islands')='Northern';
                                  'Glasgow'='South Western'",
                                levels=c("Eastern", "South Western", "Northern"))]
cc.all[, region4 := car::recode(region,
                                "c('North Eastern', 'Highlands and Islands')='Northern'",
                                levels=c("Eastern", "Glasgow", "South Western", "Northern"))]
maxdate.death <- max(cc.all$Date.Death, na.rm=TRUE)

cc.all[, SIMD.quintile := as.factor(car::recode(SIMD.quintile, "'Unknown'=NA"))]

cc.all[, sex := car::recode(as.factor(sex), "1='Male'; 2='Female'")]
cc.all[, sex := factor(sex, levels=c("Female", "Male"))]

cc.all[, agegr20 := as.factor(car::recode(as.integer(AGE),
                                          "0:39='0-39'; 40:59='40-59';
                               60:74='60-74'; 75:hi='75 or more'"))]
cc.all[, agegr20plus := as.factor(car::recode(as.integer(AGE),
                                              "0:19='0-19'; 20:39='20-39'; 40:59='40-59';
                               60:74='60-74'; 75:hi='75 or more'"))]
cc.all[, agegr3 :=
             as.factor(car::recode(AGE,
                                   "0:59='0-60 years'; 60:74='60-74 years'; 75:hi='75+ years'"))]
cc.all[, agegr2 :=
             as.factor(car::recode(AGE,
                                   "0:74='0-74 years'; 75:hi='75+ years'"))]

if(linkdate != "jun18") {
    cc.all[, occup := "Other / undetermined"]
    if(linkdate=="jan28") {
        cc.all[TEACHER_FLAG==1, occup := "Teacher"]
        cc.all[HCW_FLAG==1, occup := "Health care, not PF / undetermined"]
    } 
    if(occup) {
        ## use hcw table to code PF health care workers
        cc.all[ANON_ID %in% hcw[role=="pf_any", ANON_ID], occup := "Health care PF"]
        cc.all[, occup := base::factor(occup,
                                       levels=c("Other / undetermined",
                                                "Teacher",
                                                "Health care, not PF / undetermined",
                                                "Health care PF"))]
    } else { # only 3 levels
        cc.all[, occup := base::factor(occup,
                                       levels=c("Other / undetermined",
                                                "Teacher",
                                                "Health care, not PF / undetermined"))]
    }
 
#######################################################################################

    source("ct_ecoss.R", verbose=TRUE)
    rm(ecoss.withdate)
    
#######################################################################################

    ## left join cc.all with first positive test from ct
    ct <- ct[!is.na(Sgene.dropout)]  ## drop records where the S gene signal is NA
    setkey(ct, SpecimenDate)
    ct.first <- ct[!duplicated(ANON_ID)]
    setnames(ct.first, "SpecimenDate", "Ct.SpecimenDate") # don't want these fields imported twice
    setnames(ct.first, "flag_lighthouse_labs_testing", "Ct.Lighthouse.flag") 
    #ct.first[, SourceLab := NULL]
    setkey(ct.first, ANON_ID)
    cc.all <- ct.first[, ][cc.all]
    ## some controls later tested positive in the Lighthouse lab but were not ascertained as cases by the cutoff date for case ascertainment 
    print(head(cc.all[Ct.SpecimenDate==as.Date("2021-02-04"), .(SPECIMENDATE, CASE)]))
    rm(ct.first)
    
    ## use max and min sgtf to recode undetermined values where possible
    cc.all[Sgene.dropout=="Undetermined" & max.Sgene.dropout==3, Sgene.dropout := "Definite dropout"]
    cc.all[Sgene.dropout=="Undetermined" & min.Sgene.dropout==1, Sgene.dropout := "No dropout"]

    if("flag_covid_symptomatic" %in% names(cc.all)) {
        cc.all[flag_covid_symptomatic=="", flag_covid_symptomatic := NA]
        cc.all[, symptomatic := as.integer(flag_covid_symptomatic=="true")]
    }

    ## left join cc.all with first positive test from ECOSS
    ## sort within ID by test result (positive before negative) and date within test result
    ecoss.pos <- ecoss[ecoss.result=="Positive"]
    setkey(ecoss.pos, SpecimenDate) 
    ecoss.first <- ecoss.pos[!duplicated(ANON_ID)]
    setkey(ecoss.first, ANON_ID)
    setkey(cc.all, ANON_ID)
    ## left join cc.all with ecoss on ANON_ID
    cc.all <- ecoss.first[cc.all]
    rm(ecoss.pos)
    rm(ecoss.first)

    ## SMR00 and SMR04 
    smr00 <- RDStodt(smr00.filename, keyname="ANON_ID") # outpatients
    smr00[, CLINIC_DATE := as.Date(CLINIC_DATE)] # convert PosixCT
    ## Mode of Clinical Interaction identifies the setting where contact between a Health Care Professional and a patient/carer takes place. This is further defined as a two-way interaction where the patient has the option to have subsequent dialogue with HCP.
    ## Codes and Values: 1 Face to Face, 2 Telephone, 3 Videolink, 4 Written
    if("MODE_OF_CLINICAL_INTERACTION" %in% names(smr00)) {
        smr00 <- smr00[MODE_OF_CLINICAL_INTERACTION==1]
    }
    smr04 <- RDStodt(smr04.filename, keyname="ANON_ID")
    setnames(smr04, "ADMISSION_DATE", "Admission.Date", skip_absent=TRUE)
    setnames(smr04, "DISCHARGE_DATE", "Discharge.Date", skip_absent=TRUE)
    smr04[, Admission.Date:= as.Date(Admission.Date, format="%m/%d/%Y")]
    smr04[, Discharge.Date:= as.Date(Discharge.Date, format="%m/%d/%Y")]
   ## 1654 records of psychiatric admissions during 2020: last record 9 Oct
    ## smr06 <- RDStodt(smr06.filename, keyname="ANON_ID") # cancer incidence date and site

    ## for test-positive cases in households where an earlier member tested positive, add
    ## index (2 or more) and interval in days since last case in household
    households.multicase <- cc.all[!is.na(HID) & testpositive.case==TRUE, .N, by=HID]
    setnames(households.multicase, "N", "Ncases.household")
    setkey(households.multicase, HID)
    setkey(cc.all, HID)
    cc.all <- households.multicase[cc.all]
    
    ## assign secondary case status in households
    multicases <- cc.all[!is.na(HID) & CASE==1 & Ncases.household > 1,
                         .(ANON_ID, SPECIMENDATE, HID)]
    setkey(multicases, SPECIMENDATE)
    multicases <- multicases[, add.index(.SD), by=HID, .SDcols=c("ANON_ID", "SPECIMENDATE")]
    multicases <- multicases[!is.na(ANON_ID)] # why does this return NA values for ANON_ID? 
    multicases <- multicases[, .(ANON_ID, case.index, case.interval)]

    setkey(multicases, ANON_ID)
    setkey(cc.all, ANON_ID)

    cc.all <- multicases[cc.all]
    cc.all[, secondarycase := as.integer(case.interval > 4 & case.interval <= 14)]
    
} # Ct, ECOSS, SMR04

##############################################################

## create specimen date table to be imported into other tables
cc.specimendate <- cc.all[, .(ANON_ID, SPECIMENDATE)]
setkey(cc.specimendate, ANON_ID)

diagnoses <- RDStodt(diagnoses.filename, keyname="ANON_ID")

rapid <- RDStodt(rapid.filename, keyname="ANON_ID")
setnames(rapid, "DISCHARGE_DATE", "Discharge.Date", skip_absent=TRUE)
lastdate.rapid <- max(c(rapid$Admission.Date, rapid$Discharge.Date), na.rm=TRUE) 

procedures <- RDStodt(procedures.filename, keyname="ANON_ID") # should have dates

if(linkdate != "jun18") {
    ## import specimendate field into rapid
    setkey(rapid, ANON_ID) # why isn't it still keyed? 
    rapid <- cc.specimendate[rapid] ## left join of rapid with specimendate field in cc.all

    #################################
    ## generate table of admissions for time to hospitalisation analysis
    ##exclude discharges on or before SPECIMENDATE
    admissions <- rapid[Admission.Date >= SPECIMENDATE |
                        Admission.Date < SPECIMENDATE & Discharge.Date > SPECIMENDATE]
    admissions[, daystoadmission := as.integer(Admission.Date - SPECIMENDATE)]
    ## set daystoadmission to 0 for those already in hospital
    admissions[daystoadmission < 0, daystoadmission := 0]
    setkey(admissions, daystoadmission)
    admissions <- admissions[!duplicated(ANON_ID),
                             .(ANON_ID, Admission.Date, daystoadmission)]
    setkey(admissions, ANON_ID)
    cc.all <- admissions[cc.all]
    cc.all[, censoringdays.admission := as.integer(
                 pmin(lastdate.rapid, Admission.Date, na.rm=TRUE) - SPECIMENDATE)]
    ## some cases have specimendate after lastdate.rapid
    cc.all[censoringdays.admission < 0, censoringdays.admission := 0]
  
    cc.all[, adm.within14 := as.integer(!is.na(daystoadmission) & daystoadmission <= 14)] 
    cc.all[, adm.within28 := as.integer(!is.na(daystoadmission) & daystoadmission <= 28)]

    ################  hosp.recent #################################
    ## select for each individual the most recent smr00 attendance before specimendate
    setkey(smr00, ANON_ID)
    setkey(cc.all, ANON_ID)
    smr00 <- cc.all[, .(ANON_ID, SPECIMENDATE)][smr00]
    smr00 <- smr00[CLINIC_DATE < SPECIMENDATE]
    setorder(smr00, -CLINIC_DATE)
    smr00 <- smr00[!duplicated(ANON_ID)]
    setkey(smr00, ANON_ID)
    smr00[, outpat.recent := as.integer(SPECIMENDATE - CLINIC_DATE > 4 &
                                      SPECIMENDATE - CLINIC_DATE <= 14)]
    recent.smr00 <- smr00[!is.na(outpat.recent) & outpat.recent == 1]
  
    ## select for each individual the most recent smr04 discharge date before specimendate
    setkey(smr04, ANON_ID)
    ## import SPECIMENDATE into smr04
    smr04 <- cc.all[, .(ANON_ID, SPECIMENDATE)][smr04]
    smr04 <- smr04[Discharge.Date < SPECIMENDATE]
    setkey(smr04, Admission.Date)
    smr04 <- smr04[!duplicated(ANON_ID)]
    setkey(smr04, ANON_ID)
    smr04[, hosp.recent := as.integer(
    (SPECIMENDATE - Admission.Date > 4 &
     (is.na(Discharge.Date) | Discharge.Date - SPECIMENDATE >= 0)) | 
    (SPECIMENDATE - Discharge.Date > 0 & SPECIMENDATE - Discharge.Date <= 14) 
    )]
    recent.smr04 <- smr04[!is.na(hosp.recent) & hosp.recent == 1]

    diagnoses.last25 <- RDStodt(diagnoses.last25.filename, keyname="ANON_ID")
    diagnoses.all <- rbind(diagnoses, diagnoses.last25)
    ## import specimendate into diagnoses.all
    diagnoses.all <- cc.specimendate[diagnoses.all]

    ## select most recent discharge date for each individual from diagnoses.all
    diagnoses.all[, revdate := -as.integer(DISCHARGE_DATE)]
    setkey(diagnoses.all, revdate, ANON_ID)
    diagnoses.all <- diagnoses.all[!duplicated(ANON_ID), .(ANON_ID, SPECIMENDATE, DISCHARGE_DATE)]
    table(duplicated(diagnoses.all$ANON_ID)) 

    ## this step imports DISCHARGE_DATE from diagnoses.all into rapid
    diagnoses.dischargedates <- diagnoses.all[, .(ANON_ID, DISCHARGE_DATE)]
    setkey(diagnoses.dischargedates, ANON_ID)

    rapid <- diagnoses.dischargedates[rapid]
    ## discharge dates from diagnoses (SMR01) are used to fill in missing discharge dates
    rapid[is.na(Discharge.Date) & Admission.Date <= DISCHARGE_DATE,
          Discharge.Date := DISCHARGE_DATE]
    
    ## create hosp.recent variable for diagnoses.all (have kept SPECIMENDATE)
    diagnoses.all[, hosp.recent := as.integer(is.na(DISCHARGE_DATE) &
                                              SPECIMENDATE - DISCHARGE_DATE  > 0 & 
                                              SPECIMENDATE - DISCHARGE_DATE  < 15)]  
    recent.diagnoses <- diagnoses.all[!is.na(hosp.recent) & hosp.recent == 1]

    ## hosp.recent if in hospital at any time from 5 to 14 days before specimen date
    ## this set is defined as
    ## in hospital on specdate with admission >= 5 days before
    ## discharged before specdate with admission <= 14 days before
    
    ## (disch date missing or >= specdate) AND specdate - admission date > 4
    ## OR
    ## discharge date nonmissing AND specdate - discharge date >= 5 and <= 14
    
    rapid[, hosp.recent := as.integer(
    (SPECIMENDATE - Admission.Date > 4 &
     (is.na(Discharge.Date) | Discharge.Date - SPECIMENDATE >= 0))
    | 
    (!is.na(Discharge.Date) &
     (SPECIMENDATE - Discharge.Date >= 5 & 
      SPECIMENDATE - Discharge.Date <= 14)))]

    ## prob.hcai is defined as
    ## (disch date missing or >= specdate) AND specdate - admission date > 7

    rapid[, prob.hcai := as.integer(SPECIMENDATE - Admission.Date > 7 &
                                    (is.na(Discharge.Date) |
                                     Discharge.Date - SPECIMENDATE >= 0))]
    
    recent.rapid <- rapid[!is.na(hosp.recent) & hosp.recent == 1]
    probhcai.rapid <- rapid[!is.na(prob.hcai) & prob.hcai ==1]
    
    cc.all[, hosp.recent := as.factor(ANON_ID %in% recent.rapid$ANON_ID |
                                      ANON_ID %in% recent.diagnoses$ANON_ID |  
                                      ANON_ID %in% recent.smr04$ANON_ID |
                                      ANON_ID %in% recent.smr00$ANON_ID)]
    
    cc.all[, prob.hcai := as.factor(ANON_ID %in% probhcai.rapid$ANON_ID)] 
}

## import specimen date into diagnoses
diagnoses <- cc.specimendate[diagnoses]

source("icdchapters.R")

## FIXED: move subchapter assignment to here
## use overlap join to assign chapter and subchapter to ICD diagnoses
    
diagnoses[, icdnum := icdToInt(ICD10)]
diagnoses[, icdnum2 := icdnum]
setkey(diagnoses, icdnum, icdnum2)
setkey(icd.subchapters, startnum, endnum)
diagnoses <- foverlaps(diagnoses, icd.subchapters[, .(chnum, subchnum, startnum, endnum)])

if(linkdate != "jan28") {
    diagnoses <- diagnoses[, .(ANON_ID, SPECIMENDATE, DISCHARGE_DATE,
                           ICD10, chnum, subchnum)]
    ids.icd.neoplasm.lastyear <-
        unique(diagnoses[as.integer(SPECIMENDATE - DISCHARGE_DATE) + 25 > 365 &
                     grepl("^C[0-9]|^D[0-4]", ICD10), ANON_ID])
}

###############################################

if(linkdate != "jun18") {
    if(controls) { # household composition from UPRN vars
    cc.all[, hh.upto4 := UPRN_0_4] 
    cc.all[is.na(hh.upto4), hh.upto4 := 0] 
    cc.all[, hh.upto4 := car::recode(hh.upto4, "NA=0; 3:hi = 3")]
    cc.all[, hh.5to11 := UPRN_5_11] 
    cc.all[, hh.5to11 := car::recode(hh.5to11, "NA=0; 4:hi = 4")]
    cc.all[, hh.12to17 := UPRN_12_17] 
    cc.all[, hh.12to17 := car::recode(hh.12to17, "4:hi = 4")]
    cc.all[, hh.over18 := UPRN_adult] ## mostly nonmissing
    cc.all[, hh.over18 := car::recode(hh.over18, "0=1; 10:hi = 10")]

    setnafill(cc.all, cols=c("hh.upto4", "hh.5to11", "hh.12to17", "hh.over18"), fill=0)

    cc.all[, hh.schoolage := hh.5to11 + hh.12to17]
    cc.all[, hh.schoolage := car::recode(hh.schoolage, "5:hi = 5")]
    cc.all[, hh.schoolagegr := car::recode(hh.schoolage,
                                        "0='0 school age';
                                         1='1 school age';
                                        2:hi='2 or more'",
                                        as.factor=TRUE)]
    cc.all[, hh.schoolage.any := hh.schoolage > 0]

    cc.all[, hh.over18gr := car::recode(hh.over18,
                                        "1='1 adult';
                                         2='2 adults';
                                         3:4='3 to 4';
                                         5:9='5 to 9';
                                         10:hi='10 or more'",
                                        as.factor=TRUE,
                                        levels=c("1 adult", "2 adults", "3 to 4",
                                                 "5 to 9", "10 or more"))]
    cc.all[, hh.over18gr4 := car::recode(hh.over18gr, "'5 to 9'='5 or more';
                                                 '10 or more'='5 or more'")]
    cc.all[, hh.over18gr3 := car::recode(hh.over18gr, "'3 to 4'='3 or more';
                                                       '5 to 9'='3 or more';
                                                       '10 or more'='3 or more'")]
    cc.all[, adultsgt2 := as.factor(hh.over18 > 2)]
    cc.all[, adultsgt1 := as.factor(hh.over18 > 1)]
    cc.all[, preschool.any := as.factor(hh.upto4 > 0)]

    table.numadults.care <- with(cc.all[AGE > 70],
                                 table(hh.over18 >= 10, care.home))
    
    care.home.reclassified <- table.numadults.care[, 2]
    
    cc.all[hh.over18 >= 10 & AGE >=70, care.home=="Care/nursing home"]
    }


    cc.all[, inhosp := as.factor(inhosp)]
    cc.all[, exp.group := "No exposure"]
    cc.all[care.home=="Care/nursing home", exp.group :="Care home"]
    cc.all[care.home=="Independent" & hosp.recent==TRUE,
           exp.group := "Recent hospital exposure"]
    cc.all[, exp.group := factor(exp.group,
                                 levels=c("No exposure", "Care home",
                                        "Recent hospital exposure"))]
} # household composition

################# SICSAG ##########################

## merge SICSAG data and overwrite icu and hdu values incorrectly coded 0
if(sicsag) {
sicsag <- RDStodt(sicsag.filename, keyname="ANON_ID")
sicsag[, AdmitHosp := as.Date(AdmitHosp)]
sicsag[, AdmitUnit := as.Date(AdmitUnit)]
sicsag[, DiscDate := as.Date(DiscDate)]
if(linkdate != "jun18") {
    sicsag[DateDiscHosp=="", DateDiscHosp := NA]
}
sicsag[, DateDiscHosp := as.Date(DateDiscHosp)]
## get first and last date of entry in SICSAG table 
mindate.sicsag <- min(sicsag$AdmitUnit, na.rm=TRUE) 
maxdate.sicsag <- max(sicsag$AdmitUnit, na.rm=TRUE) 

## what are covidICUorHDU codes 4 and 5?
sicsag <- sicsag[covidICUorHDU==1 | covidICUorHDU==3]

## left join with cc.all to import SPECIMENDATE
sicsag <- cc.all[, .(ANON_ID, SPECIMENDATE)][sicsag]  ## 2880 records
                 
## all 2287 cases coded as COVID critical care are in SICSAG
print(table(cc.all[icu==1 | hdu==1, ANON_ID] %in% sicsag$ANON_ID))

## some anon ids in SICSAG are not coded in cc.all as COVID critical care -- possibly because did not have covid diagnosis or outside time window
print(table(sicsag$ANON_ID %in% cc.all[icu==1 | hdu==1, ANON_ID]))

## records with date of discharge from critical care before first positive test
with(sicsag, table(!is.na(DiscDate) & DiscDate < SPECIMENDATE))
## records with date of admission more than 21 days after first positive test
with(sicsag, table(SPECIMENDATE - AdmitUnit > 21)) # 1181records

######  generate Date.firstcritical to be used to compute failure time ########## 
## restrict to first SICSAG record for each ID that is no earlier than SPECIMENDATE
## and left join cc.all
sicsag.firstcritical <- sicsag[AdmitUnit >= SPECIMENDATE]
setkey(sicsag.firstcritical, ANON_ID)
sicsag.firstcritical <- sicsag.firstcritical[!duplicated(ANON_ID)]
setnames(sicsag.firstcritical, "AdmitUnit", "Date.firstcritical")
sicsag.firstcritical <- sicsag.firstcritical[, .(ANON_ID, Date.firstcritical)]
setkey(sicsag.firstcritical, ANON_ID)
cc.all <- sicsag.firstcritical[cc.all]

x <- with(cc.all,
         as.integer(
             pmin(max(sicsag$AdmitUnit, na.rm=TRUE), Date.firstcritical, na.rm=TRUE) -
             SPECIMENDATE))
print(summary(x))
cc.all[x < 0, .(SPECIMENDATE, Date.firstcritical)]

print(summary(as.integer(cc.all[CASE==1, Date.firstcritical - SPECIMENDATE])))
# 75% of admissions to critical care are within 10 days of SPECIMENDATE

#######  generate a separate field for date of admission within 21-day time window
## exclude from sicsag records that are outside 21-day time window
sicsag.retained <- sicsag[
(
    is.na(DiscDate) |
    (!is.na(DiscDate) & DiscDate > SPECIMENDATE)
) &
(
    is.na(AdmitUnit) |
    (!is.na(AdmitUnit) & SPECIMENDATE - AdmitUnit <= 21)
)]

## 1160 out of 1393 records retained

# drop duplicate ID variables after sorting by AdmitUnit
setkeyv(sicsag.retained, c("ANON_ID", "AdmitUnit"))
sicsag.unique <- sicsag.retained[!duplicated(ANON_ID),
                                 c("ANON_ID", "AdmitUnit", "covidICUorHDU")]
setkey(sicsag.unique, ANON_ID)
cc.all <- sicsag.unique[cc.all]

## overwrite icu or hdu field where coded in sicsag
cc.all[covidICUorHDU==1, icu := 1]
cc.all[covidICUorHDU==3, hdu := 1]
## overwrite icu and hdu field where not coded in sicsag
cc.all[is.na(covidICUorHDU), icu := 0]
cc.all[is.na(covidICUorHDU), hdu := 0]

cc.all[, daystocritical := as.integer(Date.firstcritical - SPECIMENDATE)]
# cc.all[daystocritical < 0, daystocritical := 0] # shouldn't need this
min(cc.all$AdmitUnit, na.rm=TRUE)
max(sicsag$AdmitUnit, na.rm=TRUE)
      
cc.all[, censoringdays.critical := as.integer(
             pmin(maxdate.sicsag, Date.firstcritical, na.rm=TRUE) - SPECIMENDATE)]
  cc.all[, critical.within14 := as.integer(!is.na(daystocritical) & daystocritical <= 14)] 
cc.all[, critical.within28 := as.integer(!is.na(daystocritical) & daystocritical <= 28)] 

}

cc.all[, dispensing.days := as.integer(SPECIMENDATE - as.Date("2019-06-01"))]

## exclude 3 cases already dead before their SPECIMENDATE
cases.deceased <- with(cc.all, CASE==1 &
                                  !is.na(Date.Death) &
                                  Date.Death < SPECIMENDATE)

## imputed SPECIMENDATE is 14 days before death
## but note that mode for days to death is 8 days
cc.all[, daystodeath := as.integer(Date.Death - SPECIMENDATE)]
## those with SPECIMENDATE later than maxdate.death should have censoringdays.death set to NA
cc.all[, censoringdays.death := as.integer(
             pmin(maxdate.death, Date.Death, na.rm=TRUE) - SPECIMENDATE)]
cc.all[censoringdays.death < 0, censoringdays.death := NA] 
cc.all[, death.within28 := as.integer(!is.na(daystodeath) & daystodeath <= 28)] 
cc.all[, death.within56 := as.integer(!is.na(daystodeath) & daystodeath <= 56)] 

## exclude cases classified as unobservable, for consistency with controls
## FIXME -- make this work with the field now in cc.all
if(!is.null(chistatus.filename)) {
    cases.status <-
        RDStodt(chistatus.filename, keyname="ANON_ID")
    setnames(cases.status, "CHI_EXTENDED_STATUS", "EXTENDED_STATUS", skip_absent=TRUE)
    cc.all <- cases.status[cc.all]
    with(cc.all[CASE==1],
     table(EXTENDED_STATUS, is.na(Date.Death)))
    cc.all <- cc.all[is.na(EXTENDED_STATUS) |
                     EXTENDED_STATUS == "C" |
                     (EXTENDED_STATUS == "D" & !is.na(Date.Death))]
}

if(controls) { # case-control study 

    ## exclude controls already dead on date of test of case they were matched to
    controls.deceased <- with(cc.all, CASE==0 &
                                      !is.na(Date.Death) &
                                      Date.Death <= SPECIMENDATE)
    cc.all <- cc.all[!(controls.deceased | cases.deceased)]                         
    
    
                                        #with(cc.all[CASE==1], table(ANON_ID %in%     
if(linkdate=="jun18") {
    with(cc.all[CASE==1],
         table(Explanation, is.na(Date.Death)))
} ## tabulate Explanation field among those not recorded as dead

## rows have been replicated within strata containing N cases so that each case and control appears N^2 times
## remove these replicated rows
    cc.all <- distinct(.data=cc.all, ANON_ID, stratum, .keep_all = TRUE)
    setkey(cc.all, ANON_ID)
    cc.all <- copy(cc.all) # force a deep copy
    
    num.casectrl.strata <- cc.all[, .N, by=c("stratum", "CASE")]
    ## what we want is a table with one row per stratum, separate columns for numcases and numctrls
    num.casectrl.strata <- dcast(num.casectrl.strata, stratum ~ CASE, value.var="N")
    colnames(num.casectrl.strata) <- c("stratum", "controls", "cases")
    with(num.casectrl.strata, print(table(controls, cases)))
} 

######################################################################
## classification of cases into 4 groups

## 8 individuals without a positive test result are included as cases but do not have
## covid_cod==1 or diag.case==1

## all cases have nonmissing SPECIMENDATE
## controls should be assigned same SPECIMENDATE as the case they were matched to

## deathwithin28 is a logical variable that literally means death within 28 days of a positive test
## should be coded as FALSE for those who did not test positive
if(linkdate == "jun18") { # death within 28 days of positive test
    cc.all[, deathwithin28 := testpositive.case & 
                 !is.na(Date.Death) & Date.Death - SPECIMENDATE <= 28]
} else { # death within 28 days of imputed diagnosis date
    cc.all[, deathwithin28 := testpositive.case &
                 !is.na(Date.Death) & Date.Death - SPECIMENDATE <= 28]
}

## Sharon's variable dead28 is assigned by this line
## Covid_CC_linkage_Part2_desktop.R:
## cc$dead28 <- ifelse(!is.na(cc$DATE_OF_DEATH) & cc$DATE_OF_DEATH >= cc$SPECIMENDATE & as.numeric(cc$DATE_OF_DEATH - cc$SPECIMENDATE) <=28, 1, 0)
## this evaluates to 1 for anyone classified as a case who dies within 28 days of specimen
## date even if this is a dummy specimen date.
## but fatal cases without a positive test could have been ascertained only through death cert or diagnosis. 
with(cc.all[CASE==1], table(dead28, deathwithin28, exclude=NULL)) 
with(cc.all[CASE==1], table(icu, hdu, exclude=NULL)) 
cc.all[, criticalcare := icu==1 | hdu==1]

## adm28 should be 0 for cases who did not test positive
cc.all[, adm28 := adm28 & testpositive.case]

## new linkage does not have the variable nrs_covid_case
with(cc.all[CASE==1 & covid_ucod==1], table(deathwithin28, exclude=NULL))
## 3701 cases with covid_cod have deathwithin28==1
## all those with covid_ucod==1 have covid_cod==1 
with(cc.all[CASE==1], table(dead28, deathwithin28, exclude=NULL)) 

with(cc.all[CASE==1 & testpositive.case], table(criticalcare, deathwithin28, exclude=NULL))
with(cc.all[CASE==1], table(criticalcare, testpositive.case, exclude=NULL))
with(cc.all[CASE==1], table(testpositive.case, deathwithin28, exclude=NULL))

## coding of case groups
cc.all[, group := "Unclassified"]

## define group A (severe cases)
## group A includes
## anyone entering critical care within 21 days of positive test or
## death within 28 days of positive test or
## certified with certified with covid as underlying cause
## in later releases, includes all cases in critical care with discharge diagnosis of covid
cc.all[CASE==1 & (criticalcare | deathwithin28 | covid_ucod==1),
       group := "A"]

## assign remaining test-positive cases hospitalized within 28 days of positive test as group B
cc.all[testpositive.case & group=="Unclassified" & (adm28 > 0), group := "B"]
## assign all remaining test-positive cases to group C
cc.all[testpositive.case & group=="Unclassified", group := "C"]
## assign remaining cases with mention on death cert to group D
cc.all[CASE==1 & group=="Unclassified" & covid_cod==1, group := "D"]

print(table(cc.all$CASE, cc.all$group, exclude=NULL))

## assign a logical variable for fatalcase, and a binary variable for fatal.casegroup
cc.all$fatalcase <- with(cc.all, CASE==1 & group=="A" & (deathwithin28 | covid_ucod==1))
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
casegroups <- cc.all[CASE==1, .(stratum, group, fatal.casegroup)]
colnames(casegroups)[2] <- "casegroup"
setkey(casegroups, casegroup) # orders by casegroup so dropping duplicates retains most severe casegroup
casegroups <- casegroups[!duplicated(casegroups$stratum), ]
setkey(casegroups, stratum)
setkey(cc.all, stratum)
cc.all <- casegroups[cc.all]
setkey(cc.all, ANON_ID)

## for cases, overwrite the casegroup field with the group assigned above
cc.all[CASE==1, casegroup := group]
table(cc.all$CASE, cc.all$casegroup, exclude=NULL)

## drop records with missing casegroup (controls in strata with no remaining classified case)
cc.all <- cc.all[!is.na(casegroup)]
cc.all[, casegr := ifelse(is.case & is.na(severe.casegr),
                          "Not severe", as.character(severe.casegr))]
cc.all[casegr == "Not severe" & RAPID == 1, casegr := "Not severe, hospitalized"]
cc.all[, casegr := car::recode(casegr, "'Not severe'='Not severe, not hospitalized'",
                               as.factor=TRUE,
                               levels=c("Not severe, not hospitalized",
                                        "Not severe, hospitalized",
                                        "Critical care, non-fatal", 
                                        "Critical care, fatal",
                                        "No critical care, test-positive, fatal", 
                                        "No critical care, not test-positive, fatal"))]
cc.all[, casegr2 := as.factor(car::recode(as.integer(casegr),
                                          "1:2='No critical care, non-fatal'; 3:6='Severe';"))]
cc.all[, casegr3 := car::recode(as.integer(casegr),
                                "1:2='No critical care, non-fatal';
                                 3='Critical care, non-fatal';
                                 4:6='Fatal';",
                                as.factor=TRUE,
                                levels=c("No critical care, non-fatal",
                                        "Critical care, non-fatal", 
                                        "Fatal"))]
####################################################################

## import most recent SMR01 discharge date into rapid
if(linkdate != "jun18") {
    cc.all[qSIMD.integer==9, qSIMD.integer := NA]
} # import case status into shielded.full, assign shielding interval, group num adults

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

if(linkdate == "jun18") {
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
} # import ONOMAP ethnicity

setkey(cc.all, ANON_ID)
######################################################

rm(smr00)
rm(rapid)
rm(diagnoses.all)
rm(ecoss)
rm(cc.specimendate)
rm(numctrls.strata)

gc()

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

#################### import shielding data ####

if(shielding) { ## shielding table ###
        
        shielded.full <- RDStodt(shielding.full.filename, keyname="ANON_ID")
        setnames(shielded.full, "group", "shield.group", skip_absent=TRUE)
        setnames(shielded.full, "batch", "shield.batch", skip_absent=TRUE)
        setnames(shielded.full, "age", "AGE", skip_absent=TRUE)
        setnames(shielded.full, "care_home", "care.home", skip_absent=TRUE)
                                        #shielded.full[nursing_home==1, care.home := 1]
                                        #shielded.full[, care.home := as.factor(care.home)]
        shielded.full[, sex := as.factor(car::recode(sex, "1='Male'; 2='Female'"))]
        setkey(shielded.full, shield.batch)
        shielded.full <- batches[shielded.full]
        print(with(shielded.full, table(is.na(ANON_ID)))) # 17314 with nonmissing ANON_ID
        ## no duplicates among nonmissing ANON_IDs
        
        shielded.full[, shield.group := car::recode(shield.group,
                                                "1='Solid organ transplant';
                                      2='Specific cancers';
                                      3='Severe respiratory';
                                      4='Rare diseases';
                                      5='On immunosuppressants';
                                      6='Pregnant with heart disease';
                                      7='Additional conditions'",
                                     as.factor=TRUE,
                                     levels=c(
                                         "Solid organ transplant",
                                         "Specific cancers",
                                         "Severe respiratory",
                                         "Rare diseases",
                                         "On immunosuppressants",
                                         "Pregnant with heart disease",
                                         "Additional conditions"
                                     ))]
        shielded.full[, agegr20 := as.factor(car::recode(as.integer(AGE),
                                                     "0:39='0-39'; 40:59='40-59';
                               60:74='60-74'; 75:hi='75 or more'"))]
        ## remove obvious wrong assignments
    shielded.full <- shielded.full[!(shield.group == "Pregnant with heart disease" & AGE > 60)]
 #   shielded$ANON_ID)) # 2305 shield-eligible  cases
    with(cc.all[CASE==1], table(ANON_ID %in% shielded.full$ANON_ID)) # 2319 shield-eligible cases

        ## left join of cc.all with subset of shielded.full in which ANON_ID is nonmissing
    shielded.full.cc <- shielded.full[!is.na(ANON_ID),
                                      .(ANON_ID, shield.batch, Date.Sent,
                                        shielding_id, shield.group)]
    setkey(shielded.full.cc, ANON_ID)
    names(shielded.full.cc)[names(shielded.full.cc) %in% names(cc.all)]

    cc.all <- shielded.full.cc[cc.all] ## make sure that cc.all is keyed on ANON_ID

    cc.all[, shield.any := as.factor(!is.na(shield.group))]
    cc.all[is.na(shield.batch), shield.batch := 0]
    cc.all[, shield.batch := as.factor(shield.batch)]

    ## import case status into shielded.full
    ## cannot just do a left join on ANON_ID as this is missing for most records
    ## in shielded.full 
    cases <- cc.all[CASE==1, .(ANON_ID, casegr, casegr2, casegr3)]
    setkey(cases, ANON_ID)

    shielded.full[!is.na(ANON_ID), casegr := cases$casegr[match(ANON_ID, cases$ANON_ID)]]
    shielded.full[!is.na(ANON_ID), casegr2 := cases$casegr2[match(ANON_ID, cases$ANON_ID)]]
    shielded.full[!is.na(ANON_ID), casegr3 := cases$casegr3[match(ANON_ID, cases$ANON_ID)]]

    shielded.full[, casegr := as.character(casegr)]
    shielded.full[is.na(casegr), casegr := "Not diagnosed as case"]
    shielded.full[, casegr := factor(casegr, levels=c("Not diagnosed as case",
                                                      levels(cc.all$casegr)))]
    shielded.full[, casegr2 := as.character(casegr2)]
    shielded.full[is.na(casegr2), casegr2 := "Not diagnosed as case"]
    shielded.full[, casegr2 := factor(casegr2, levels=c("Not diagnosed as case",
                                                      levels(cc.all$casegr2)))]
 
    shielded.full[, casegr3 := as.character(casegr3)]
    shielded.full[is.na(casegr3), casegr3 := "Not diagnosed as case"]
    shielded.full[, casegr3 := factor(casegr3, levels=c("Not diagnosed as case",
                                                      levels(cc.all$casegr3)))]
 
    ## assign after.letter as 1 for those shielded if letter sent at least numdays.letter before specimen date 
    cc.all[, after.letter := ifelse(Date.Sent + numdays.letter <= SPECIMENDATE, 1, 0)]
    ## assign after.letter as 0 for those not shielded
    cc.all[is.na(Date.Sent), after.letter := 0]
    cc.all[, after.letter := as.factor(after.letter)]

    cc.all[, interval2 := 0]
    cc.all[, interval3 := 0]
    cc.all[SPECIMENDATE >=  as.Date("2020-04-03") + 14 &
           SPECIMENDATE < as.Date("2020-08-01") + 14, interval2 := 1]
    cc.all[SPECIMENDATE >=  as.Date("2020-08-01") + 14, interval3 := 1]
    
    cc.all[, interval := ifelse(SPECIMENDATE <  as.Date("2020-04-03") + 14, 1,
                         ifelse(SPECIMENDATE <  as.Date("2020-08-01") + 14, 2, 3))]
    
    cc.all[, interval2.shielding := interval2 * (as.integer(shield.any) - 1)]
    cc.all[, interval3.shielding := interval3 * (as.integer(shield.any) - 1)]

    }

crm(controls.deceased)
rm(diagnoses.dischargedates)


if(pis) { # read scrips file and import BNF chapter variables ############# 

    scrips.firsttime <- FALSE
    scripsobject.filename <- paste0(datadir, "scrips.last240days.RData") # all further reads should be from this object
    ## could drop item code
    if(scrips.firsttime) { 
        cat("Restricting scrips to last 240 days before (specimendate - 15) ... \n")
        scrips <- RDStodt(scrips.filename)
        scrips[, ANON_ID := as.integer(gsub(".* ", "", ANON_ID))] # remove date prefix
        scrips[, daysbefore := as.integer(SPECIMENDATE - dispensed_date)]
        ## restrict to last 240 days before cutoff of (specimendate - 15) 
        scrips <- scrips[daysbefore - 15 <= 240]
        ## keep only IDs of cases or controls matched to severe cases
        ids.keep <- cc.all[CASE==1 | group=="A", ANON_ID]
        scrips <- scrips[ANON_ID %in% ids.keep]
        ## save this object for re-use 
        save(scrips, file=scripsobject.filename)
        cat("done\n")
    }
    
    cat("Loading saved scrips file ...")
    load(scripsobject.filename)
    cat("done\n")

    paste.colpercent(with(cc.all, table(ANON_ID %in% scrips$ANON_ID, CASE)))
    paste.colpercent(with(cc.all, table(ANON_ID %in% diagnoses$ANON_ID, CASE)))
    
    cc.all[, scrip.any := as.factor(as.integer(ANON_ID %in% scrips$ANON_ID))]
    cc.all[, diag.any := as.factor(as.integer(ANON_ID %in% diagnoses$ANON_ID))]
    cc.all[, scripordiag := as.factor(as.integer(diag.any=="1" | scrip.any=="1"))]

    ## colchicine
    ## BNF Subparagraph Code  1001040
    ## BNF Chemical Substance 1001040G0
    ## always prescribed as 500 mcg tablets, quantity usually 28, 56 or 100
    ## dose 2-4 tablets/day so at 2 tabs/day that corresponds to 14, 28 or 50 days supply
    
    if(FALSE) {
        scrips.colchicine <- 
            subset(scrips, approved_name=="COLCHICINE",
                   select=c("ANON_ID", "SPECIMENDATE",
                            "dispensed_date", "item_strength",
                            "quantity"))            
        scrips.last15.colchicine <- 
            subset(scrips.last15, approved_name=="COLCHICINE",
                   select=c("ANON_ID", "SPECIMENDATE",
                            "dispensed_date", "item_strength",
                            "quantity"))
        scrips.colchicine <- rbind(scrips.colchicine, scrips.last15.colchicine) %>%
            as.data.table(key="ANON_ID")
        scrips.colchicine[, end.date := dispensed_date + floor(quantity / 2)]
        scrips.colchicine[, on.colchicine := dispensed_date < SPECIMENDATE &
                                SPECIMENDATE <= end.date]
        scrips.colchicine <- scrips.colchicine[on.colchicine==TRUE]
        cc.all[, colchicine.current := as.factor(as.integer(ANON_ID %in%
                                                            scrips.colchicine$ANON_ID))]
    }
}

########## restrict to severe cases (and matched controls) before importing drug variables
cat("Restricting to severe cases and matched controls\n")
cc.severe <- cc.all[casegroup=="A"]
save(cc.all, file="cc.all.RData")

if(pis) {
    source("bnfcodes.R")
    source("drugs.R") # memory
    rm(scrips.nonopiates)
    rm(dose.protonpump)
    rm(scrips.wide)
    gc()
    
    ids.bnf.diabetes <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 601])
    ids.icd.diabetes <- unique(diagnoses$ANON_ID[grep("^E1[0-4]", diagnoses$ICD10)])
    ids.diabetes.extra <- unique(c(ids.icd.diabetes, ids.bnf.diabetes))

    objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
    print(tail(objmem))
    
###########################################################################

###  diabetes based on combining Sci-Diabetes records with ICD codes and drugs 
## add in BNF codes 6.1 for diabetes drugs and
## E10 to E14 for diabetes

########## coding listed conditions ####################

    ## missing recoded as zero
    cc.severe[is.na(dm.type), dm.type := 0]
    
    ## add in extra cases notified directly from SCI-Diabetes register, without assignment
## of diabetes type from SDRN database
    cc.severe[dm.type==0 & diab.reg==1, dm.type := 3]
    
                                        #cat("Extra diabetes cases from SCI-Diabetes by diabetes type\n")
#print(table(cc.severe$dm.type, cc.severe$ANON_ID %in% ids.diabetes.extra))

## code diagnoses detected from discharges or BNF codes as unknown type
## we could classify those not on insulin as definite Type 2 but Helen says no
## REVISION: for consistency with the diabetes paper, we will not include the extra cases identified through diagnostic codes or drug codes - they may be transient/resolved

## cc.severe$dm.type[cc.severe$dm.type==0 & cc.severe$ANON_ID %in% ids.diabetes.extra] <- 3

    ## recode diabetes type
    cc.severe[, dm.type := recode.dmtype(dm.type)]
    
    ## define indicator variable for any diabetes
    cc.severe[, diabetes.any := as.integer(dm.type != "Not diabetic")]
    cc.severe[, diabetes.any := as.factor(car::recode(diabetes.any,
                                               "0='Not diabetic'; 1='Diabetic'"))]
    cc.severe[, diabetes.any := relevel(diabetes.any, ref="Not diabetic")]
    
    objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
    print(tail(objmem))
    
    source("comorbidity.R")
    
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
    
    cc.severe[is.na(shield.group), shield.group := ifelse(listed.any==1,
                                                       "Moderate risk condition",
                                                       "No risk condition")]
    cc.severe[, shield.group := car::recode(shield.group,
                                         "'Pregnant with heartdisease'='Additional conditions'",
                                         as.factor=TRUE,
                                         levels=c(
                                             "No risk condition",
                                             "Moderate risk condition",
                                             "Solid organ transplant",
                                             "Specific cancers",
                                             "Severe respiratory",
                                             "Rare diseases",
                                             "On immunosuppressants",
                                             "Additional conditions"
                                         ))]
    cc.severe[, shieldedonly.group := car::recode(shield.group,
                                               "'No risk condition'=NA; 
                                                'Moderate risk condition'=NA",
                                     as.factor=TRUE,
                                     levels=c(
                                         "Severe respiratory",
                                         "Solid organ transplant",
                                         "Specific cancers",
                                         "Rare diseases",
                                         "On immunosuppressants",
                                         "Additional conditions"
                                     ))]
    
    cc.severe[, listedgr3 := 0]
    cc.severe[listed.any=="1", listedgr3 := 1]
    cc.severe[shield.any==TRUE, listedgr3 := 2]
    cc.severe[, listedgr3 := car::recode(listedgr3,
                                  recodes=
                                      "0='No risk condition';
                                       1='Moderate risk condition';
                                       2='Eligible for shielding'",
                                  as.factor=TRUE, 
                                  levels=c("No risk condition",
                                           "Moderate risk condition",
                                           "Eligible for shielding"))]
    
    cc.severe[, diag.other := as.integer(listed.any==0 & diag.any==1)]
    cc.severe[, diag.other := as.factor(car::recode(diag.other,
                                  "0='Listed condition or no admission';
                                   1='No listed condition, but other admission diagnosis'"))]

    rm(scrips)
    rm(subset.laporte.scrips)
}
    
## if(linkdate != "jun18")  .Internal(.invokeRestart(list(NULL, NULL), NULL))

#####################################################################

## cast ICD chapters in wide format and merge with cc.severe
icdchapters.wide <- data.table::dcast(diagnoses, ANON_ID ~ chnum, fun.aggregate=length,
                                      value.var="chnum")
## use Roman numerals for ICD chapters
chnums <- as.integer(colnames(icdchapters.wide)[-1])
colnames(icdchapters.wide)[-1] <-
    paste0("Ch.", as.roman(chnums), "_",
           icdchapters$shortname[match(chnums, icdchapters$chnum)])
icdchapters.wide <- as.data.table(icdchapters.wide, key="ANON_ID")
## drop rare chapters
#cols.keep <- colSums(icdchapters.wide) >= 20
cols.keep <- colSums(icdchapters.wide) >= 60 ## drops perinatal chapter
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
icdcolsnew <- grep("^Ch\\.", colnames(cc.severe), value=TRUE)
cc.severe[, num.icdchapters := rowSums(matrix(as.integer(as.matrix(cc.severe[, ..icdcolsnew])),
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

## recode subchapters as indicators
icdcols <- grep("^Ch_", colnames(cc.severe), value=TRUE)
cc.severe[, (icdcols) := lapply(.SD, recode.indicator), .SDcols = icdcols]

cc.severe[AGE < 75] %>% save(file="data/cc.severe.lt75.RData")

rm(scrips)
rm(diagnoses)
gc()

###########################################

if(linkdate == "jun18") {
   source("paper1tables.R")
   ## ##### stepwise regressions use saved version #####################
    #nfold <- 10
    #source("stepwise.R")
    #source("pharmaco.R")
    #rmarkdown::render("pharmaco.Rmd", output_file="pharmacoBMED-D-20-02413_rev4.pdf")
} else if(linkdate != "jan28") {

    #tabulate.freqs.regressions(varnames="colchicine.current",
    #                           data=cc.all[casegroup=="A" | casegroup=="B"])

    source("shielding.R")
    rmarkdown::render("shielding.Rmd")
    source("ct.R")
    save(table.coeffs.exposure, table.symptoms, table.withCt, file="ct.tables_jan6.RData")
} else {
    source("ct.R")
    source("shielding.R")
    rmarkdown::render("shielding.Rmd")

}

## remove large objects from memory
rm(diagnoses)

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

Rprof()
print(summaryRprof(tmp)$by.total[1:20, ])
