## analysis  script for case-control study
if(exists("cl")) {
    # parallel::stopCluster(cl)
    showConnections(all = TRUE)
    closeAllConnections()
}
rm(list=ls())
gc()

linkdate <- "jul12"
linkdate <- "jul28"
## all saved data files should be saved to datadir defined by linkdate

controls <- TRUE
shielding <- FALSE
sicsag <- TRUE
pis <- TRUE
scrips.firsttime <- FALSE

## stepwise <- TRUE   
stepwise <- FALSE ## uncomment to save time if old version still valid
fatal.predict <- TRUE
fatal.predict <- FALSE

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
library(ggplot2)
#library(doParallel)
library(readxl)
library(DescTools)
library(icd)
library(gam)
library(data.table)

Rprof(tmp <- tempfile())

source("helperfunctions.R")

##################################################################
batches <- as.data.table(read.table("data/BatchReference.csv", sep="\t", header=TRUE))
setnames(batches, "Batch", "shield.batch")
batches$shield.batch <- as.integer(gsub("Batch", "", batches$shield.batch))
setkey(batches, shield.batch)
batches$Date.Sent <- gsub("Friday ", "", batches$Date.Sent)
batches$Date.Sent <- gsub("^([[:digit:]]+)([stndrdth]{2} )", "\\1 ", batches$Date.Sent)
batches$Date.Sent <- as.Date(batches$Date.Sent, "%d %B %Y")

controls <- TRUE
shielding <- TRUE
occup <- TRUE
chistatus.filename <- NULL

if(linkdate == "jul12") {
        datadir <- "./data/2021-07-12/"
        ct.filename <-  paste0(datadir, "CC_CT_Test_Results_2021-07-12.rds")
        ecoss.tests.filename <- paste0(datadir, "CC_ecoss_tests_2021-07-12.rds")

        teach.filename <- paste0(datadir, "CC_TEACHER_2021-07-12.csv")
        teach <- fread(teach.filename, key="anon_id")
        
        hcw.filename <- paste0(datadir, "CC_HCW_HHD_SWISS_NOV20_2021-07-12.csv")
        hcw <- fread(hcw.filename)
        colstokeep <- c("nov_key", "hid", "unique_id", "job_family", "sub_job_family",
                        "specialty", "patient_facing_category",
                            "hid_members_anon", "anon_id")
        hcw <- hcw[, ..colstokeep]
        hcw <- hcw[!is.na(anon_id)]
        hcw <- unique(hcw, by="anon_id")
        setkey(hcw, anon_id)
        
        rapid.filename <- paste0(datadir, "CC_RAPID_ANON_2021-07-12.rds")
        sicsag.filename <- paste0(datadir, "CC_SICSAG_ANON_2021-07-12.rds")
        cc.filename <- paste0(datadir, "CC_linked_ANON_2021-07-12.rds")
        diagnoses.filename <- paste0(datadir, "CC_SMR01_ICD10_x25_2021-07-12.rds")
        diagnoses.last25.filename <- paste0(datadir, "CC_SMR01_ICD10_25_2021-07-12.rds")
        
        shielding.full.filename <- paste0(datadir,
                                          "CC_shielding_patients_anon_20210712.csv")
        smr00.filename <- paste0(datadir, "CC_SMR00_2021-07-12.rds")
        smr04.filename <- paste0(datadir, "CC_SMR04_ANON_2021-07-12.rds")
        smr06.filename <- paste0(datadir, "CC_SMR06_ICD10_ANON_2021-07-12.rds")
        scrips.filename <- paste0(datadir, "CC_PIS_x15_ANON_2021-07-12.rds")
        scrips.last15.filename <- paste0(datadir, "CC_PIS_15_ANON_2021-07-12.rds")
        procedures.filename <- paste0(datadir, "CC_SMR01_OPCS4_MAIN.x25_2021-07-12.rds")
        procedures.last25.filename <-
            paste0(datadir, "CC_SMR01_OPCS4_MAIN.25_2021-07-12.rds")
        vacc.filename <- paste0(datadir, "CC_VACC_2021-07-12.rds")
} else if(linkdate=="jul28") {
       datadir <- "./data/2021-07-28/"
       ct.filename <-  paste0(datadir, "CC_CT_Test_Results_2021-07-28.rds")
       ecoss.tests.filename <- paste0(datadir, "CC_ecoss_tests_2021-07-28.rds")
  
       teach.filename <- paste0(datadir, "CC_TEACHER_2021-07-28.csv")
       teach <- fread(teach.filename, key="anon_id")
       
       hcw.filename <- paste0(datadir, "CC_HCW_HHD_SWISS_NOV20_2021-07-28.csv")
       hcw <- fread(hcw.filename)
       colstokeep <- c("nov_key", "hid", "unique_id", "job_family", "sub_job_family",
                            "specialty", "patient_facing_category",
                       "hid_members_anon", "anon_id")
       hcw <- hcw[, ..colstokeep]
       hcw <- hcw[!is.na(anon_id)]
       hcw <- unique(hcw, by="anon_id")
       setkey(hcw, anon_id)

       rapid.filename <- paste0(datadir, "CC_RAPID_ANON_2021-07-28.rds")
       sicsag.filename <- paste0(datadir, "CC_SICSAG_ANON_2021-07-28.rds")
       cc.filename <- paste0(datadir, "CC_linked_ANON_2021-07-28.rds")
       diagnoses.filename <- paste0(datadir, "CC_SMR01_ICD10_2021-07-28.rds")
        
       ## shielding
       shielding.full.filename <- paste0(datadir,
                                         "CC_shielding_patients_anon_2021-07-28.csv")
       smr00.filename <- paste0(datadir, "CC_SMR00_2021-07-28.rds")
       smr04.filename <- paste0(datadir, "CC_SMR04_ANON_2021-07-28.rds")
       smr06.filename <- paste0(datadir, "CC_SMR06_ICD10_ANON_2021-07-28.rds")
       scrips.filename <- paste0(datadir, "CC_PIS_ANON_2021-07-28.rds")
       scrips.last15.filename <- paste0(datadir, "CC_PIS_15_ANON_2021-07-28.rds")
       procedures.filename <- paste0(datadir, "CC_SMR01_OPCS4_MAIN2021-07-28.rds")
       vacc.filename <- paste0(datadir, "CC_VACC_2021-07-28.rds")
}

cc.all <- RDStodt(cc.filename, keyname=c("anon_id", "specimen_date"))
cc.all <- na.omit(cc.all, cols=c("anon_id", "specimen_date"))
cc.all <- unique(cc.all, by=key(cc.all)) # 3153433 rows

## create specimen date table for rolling joins with other tables
cc.specimendate <- cc.all[, .(anon_id, specimen_date)]
setkey(cc.specimendate, anon_id, specimen_date)
cc.specimendate[, joindate:= specimen_date]
save(cc.specimendate, file=paste0(datadir, "cc.specimendate.RData"))

CoDvars <- grep("cause_of_death_code_", names(cc.all), value=TRUE)
cc.all[, (CoDvars) := NULL]

## replace years before 2020 with 2020
cc.all[lubridate::year(specimen_date) < 2020, specimen_date := year.to2020(specimen_date)]
## set year to 2021 for dates before Feb 2020
cc.all[lubridate::year(specimen_date) == 2020 & specimen_date < as.Date("2020-02-01"),
       specimen_date := year.to2021(specimen_date)]
## set year to 2020 for dates after system date
cc.all[specimen_date > Sys.Date(), specimen_date := year.to2020(specimen_date)]

cc.all[, wave := 1 + as.integer(specimen_date >= as.Date("2020-06-01"))]
cc.all[specimen_date >= as.Date("2020-06-01") &  specimen_date < as.Date("2020-09-01"),
      wave := NA]

lastdate.specimen <- max(cc.all$specimen_date)
cat("Last specimen date", format(lastdate.specimen,"%d %B %Y"), "\n")

setnames(cc.all, "CASE_NO", "stratum", skip_absent=TRUE)
setnames(cc.all, "SEX", "sex", skip_absent=TRUE)
setnames(cc.all, "imumune", "immune", skip_absent=TRUE)
setnames(cc.all, "DATE_OF_DEATH", "date_of_death", skip_absent=TRUE)
setnames(cc.all, "age", "age_years", skip_absent=TRUE)
setnames(cc.all, "AgeYear", "age_years", skip_absent=TRUE)
setnames(cc.all, "simd2020_sc_quintile", "qSIMD.integer", skip_absent=TRUE)
#setnames(cc.all, "CAREHOME", "care.home", skip_absent=TRUE)

cc.all[is.na(institution_code==93) |
       (institution_code!=93 & institution_code!=98), care.home := 0]
cc.all[!is.na(institution_code) & (institution_code==93 | institution_code==98),
       care.home := 1]
cc.all[, care.home := as.factor(car::recode(care.home, "0='Independent';
                                                        1='Care/nursing home'"))]
cc.all[, care.home := relevel(care.home, ref="Independent")]

if(!controls) { # cases only
    cc.all[, CASE := 1]
} else { # case-control
    cc.all[, CASE := as.integer(is_case)]
    cc.all[, stratum := as.integer(stratum)]
}
cc.all[CASE==1, testpositive.case := diag_case==0 & cod_case==0]

## include cases ascertained through discharge diagnosis
    ## case ascertainment uses three sources in order:
    ##   test positives, discharge diagnoses, death certs
    ## is.case is defined to assign case status
    ## diag.case is a 0-1 variable encoding cases ascertained through discharge diagnosis
    ## cod.case encodes ascertainment through death cert
with(cc.all, table(is_case, cod_case)) # 137 have is_case==false and cod_case==1
with(cc.all, table(is_case, diag_case)) # 74 have is.fase==false and diag_case==1 and is_case==0
#with(cc.all[diag_case==1 & cod_case==0],
#     table(icu + hdu > 0, !is.na(date_of_death))) # 4 cases ascertained through discharge diagnosis were not in critical care and have a death certificate with no mention of COVID.
## Sharon & Jen have assigned these a dummy specimen date 7 days before date of admission
## we reconstruct the date of admission as 7 days after the specimen date,
## then code as severe if the date of death was within 28 days of admission.

cc.all[age_years < 0, age_years := 0]

#cc.all[is.na(Prisoner_Flag), Prisoner_Flag := FALSE]
## all cases not ascertained through diag or cod must be testpositive cases
## 3052 cases not testpositive cases
cc.all[, region := factor(car::recode(hb2019name, # NUTS 2 classification, but West Central should include North Lanarkshire
                                      recode=
                             "c('NHS Fife', 'NHS Lothian', 'NHS Tayside', 'NHS Forth Valley')='Eastern';
                               'NHS Greater Glasgow and Clyde'='West Central';
                               'NHS Grampian'='North Eastern'; 
  c('NHS Borders', 'NHS Ayrshire and Arran', 'NHS Dumfries and Galloway', 'NHS Lanarkshire')='Southern';
   c('NHS Western Isles', 'NHS Shetland', 'NHS Highland', 'NHS Orkney')='Highlands and Islands'"))]

cc.all[, region4 := car::recode(region,
                                "c('North Eastern', 'Highlands and Islands')='Northern'",
                                levels=c("Eastern", "West Central", "Southern", "Northern"))]
maxdate.death <- max(cc.all$date_of_death, na.rm=TRUE)

cc.all[, SIMD.quintile := as.factor(qSIMD.integer)]
cc.all[, sex := car::recode(as.factor(sex), "1='Male'; 2='Female'")]
cc.all[, sex := factor(sex, levels=c("Female", "Male"))]

cc.all[, agegr20 := as.factor(car::recode(as.integer(age_years),
                                          "0:39='0-39'; 40:59='40-59';
                               60:74='60-74'; 75:hi='75 or more'"))]
cc.all[, agegr20plus5 := car::recode(age_years,
                                    "0:19='0 to 19'; 20:29='20 to 29'; 30:39='30 to 39'; 40:49='40 to 49'; 50:hi='50 years and over'",
                                    as.factor=TRUE)]
cc.all[, agegr20plus := as.factor(car::recode(as.integer(age_years),
                                              "0:19='0-19'; 20:39='20-39'; 40:59='40-59';
                               60:74='60-74'; 75:hi='75 or more'"))]
cc.all[, agegr3 :=
             as.factor(car::recode(age_years,
                                   "0:59='0-60 years'; 60:74='60-74 years'; 75:hi='75+ years'"))]
cc.all[, agegr2 :=
             as.factor(car::recode(age_years,
                                   "0:74='0-74 years'; 75:hi='75+ years'"))]

cc.all <- hcw[cc.all] # make sure duplicate IDS are removed from hcw
rm(hcw)

cc.all[, occup := "Other / undetermined"] # default assignment
## use patient_facing_category field to code PF health care workers
cc.all[!is.na(patient_facing_category), occup := patient_facing_category]
cc.all <- teach[cc.all]
cc.all[!is.na(sector), occup := "Teacher"]
cc.all[, occup := car::recode(occup,
         "NA='Other / undetermined';
          c('Non Patient Facing', 'Unknown')='Health care, not PF / undetermined';
          'Patient Facing'='Health care PF'",
            as.factor=TRUE,
            levels=c("Other / undetermined",
                     "Teacher",
                     "Health care, not PF / undetermined",
                     "Health care PF"))]

vacc <- RDStodt(vacc.filename, key="anon_id") # ## long format, one record per dose
vacc <- vacc[vacc_status == "completed"]
setnames(vacc, "vacc_occurence_time", "vaxdate")
vacc[, vaxdate := as.Date(vaxdate)]
vacc <- vacc[, .(anon_id, vacc_dose_number, vaxdate, vacc_product_name,
                     vacc_batch_number)]
setkey(vacc, anon_id, vacc_dose_number) 
vacc <- unique(vacc, by=key(vacc)) # https://riptutorial.com/data-table/example/26336/unique-and-duplicated-no-longer-works-on-keyed-data-table
vacc.wide <- dcast(vacc, anon_id ~ vacc_dose_number,
                         value.var=c("vaxdate", "vacc_product_name",
                                     "vacc_batch_number"))
                       # fun.aggregate=toString)
rm(vacc)

vacc.wide <- unique(vacc.wide)
vacc.wide[, weekdose1 := floor(as.integer(vaxdate_1 - as.Date("2020-12-01")) / 7)]
vacc.wide[, weekdose1 := as.Date("2010-12-01") + 7 * weekdose1]
setkey(vacc.wide, anon_id)

cc.all <- vacc.wide[cc.all]
rm(vacc.wide)
cc.all[, dayssincedose1 := as.integer(specimen_date - vaxdate_1)]
cc.all[is.na(dayssincedose1) | dayssincedose1 < 0, dayssincedose1 := -1]
cc.all[, dayssincedose2 := as.integer(specimen_date - vaxdate_2)]
cc.all[is.na(dayssincedose2) | dayssincedose2 < 0, dayssincedose2 := -1]
cc.all[, weekssincedose1 := floor(dayssincedose1 / 7)]
cc.all[, vax14.dose := 0] # code as 0 unless there is a completed vaccine record 
cc.all[dayssincedose1 >=14, vax14.dose := 1]
cc.all[dayssincedose2 >= 14 , vax14.dose := 2]
cc.all[, vax14.factor := as.factor(vax14.dose)]

cc.all[is.na(vaxdate_1), vacc_product_name_1 := "No vaccine"]
cc.all[, vacc_product_name_1 := car::recode(vacc_product_name_1,
                                       recodes="'Covid-19 mRNA Vaccine Pfizer'='Pfizer';
                                             'Covid-19 Vaccine AstraZeneca'='AstraZeneca';
                                             'Covid-19 mRNA Vaccine Moderna'='Moderna'",
                                     as.factor=TRUE,
                                     levels=c("No vaccine", "Pfizer",
                                              "AstraZeneca", "Moderna"))]
## combined variable for vax dose and product
cc.all[vax14.dose==0, vaxgr := 1]
cc.all[vax14.dose==1 & (vacc_product_name_1=="Pfizer" | vacc_product_name_1=="Moderna"),
       vaxgr := 2]
cc.all[vax14.dose==1 & vacc_product_name_1=="AstraZeneca", vaxgr := 3]
cc.all[vax14.dose==2 & (vacc_product_name_1=="Pfizer" | vacc_product_name_1=="Moderna"),
       vaxgr := 4]
cc.all[vax14.dose==2 & vacc_product_name_1=="AstraZeneca", vaxgr := 5]
cc.all[, vaxgr := car::recode(vaxgr,
                                    "1='Not vaccinated';
                                       2='1 dose mRNA vaccine';
                                       3='1 dose AZ vaccine';
                                       4='2 doses mRNA vaccine';
                                       5='2 doses AZ vaccine'",
                                    as.factor=TRUE,
                                    levels=c("Not vaccinated", 
                                             "1 dose mRNA vaccine", 
                                             "1 dose AZ vaccine", 
                                             "2 doses mRNA vaccine", 
                                             "2 doses AZ vaccine"))]

#######################################################################################
    
# source("ct_ecoss.R") # run first time only
load(paste0(datadir, "ecoss.RData"))
     
#######################################################################################

if(FALSE) { # import Ct
    ## left join cc.all with first positive test from ct
ct <- ct[!is.na(Sgene.dropout)]  ## drop records where the S gene signal is NA
setkey(ct, SpecimenDate)
ct.first <- ct[!duplicated(anon_id)]
setnames(ct.first, "SpecimenDate", "Ct.SpecimenDate") # don't want these fields imported twice
setnames(ct.first, "flag_lighthouse_labs_testing", "Ct.Lighthouse.flag") 
                                        #ct.first[, SourceLab := NULL]
setkey(ct.first, anon_id)

cc.all <- ct.first[, ][cc.all]
## some controls later tested positive in the Lighthouse lab but were not ascertained as cases by the cutoff date for case ascertainment 
rm(ct.first)
rm(ct)

## use max and min sgtf to recode undetermined values where possible
cc.all[Sgene.dropout=="Undetermined" & max.Sgene.dropout==3, Sgene.dropout := "Definite dropout"]
cc.all[Sgene.dropout=="Undetermined" & min.Sgene.dropout==1, Sgene.dropout := "No dropout"]

if("flag_covid_symptomatic" %in% names(cc.all)) {
        cc.all[flag_covid_symptomatic=="", flag_covid_symptomatic := NA]
        cc.all[, symptomatic := as.integer(flag_covid_symptomatic=="true")]
    }
}

## left join cc.all with first positive test from ECOSS
    ## sort within ID by test result (positive before negative) and date within test result
ecoss.pos <- ecoss[ecoss.result=="Positive"]
setkey(ecoss.pos, SpecimenDate) 
ecoss.first <- ecoss.pos[!duplicated(anon_id)]
setkey(ecoss.first, anon_id)
setkey(cc.all, anon_id)
## left join cc.all with ecoss on anon_id
cc.all <- ecoss.first[cc.all]
rm(ecoss.pos)
rm(ecoss.first)
rm(ecoss)
gc()

## check that specimen_date of test-positive cases is first specimen date in ecoss

    ## FIXME: this section should be moved
    ## for test-positive cases in households where an earlier member tested positive, add
    ## index (2 or more) and interval in days since last case in household

## households.multicase has one record per household
households.multicase <- cc.all[!is.na(hid) & testpositive.case==TRUE, .N, by=hid]
setnames(households.multicase, "N", "Ncases.household")
setkey(households.multicase, hid)
setkey(cc.all, hid)
cc.all <- households.multicase[cc.all]
    
## assign secondary case status in households
multicases <- cc.all[!is.na(hid) & CASE==1 & Ncases.household > 1,
                         .(anon_id, specimen_date, hid)]
setkey(multicases, hid, specimen_date)
multicases[, household.casenum := 1:.N, by=hid]
multicases[, case.interval := c(NA, diff(as.integer(specimen_date))), by=hid]
multicases <- multicases[!is.na(anon_id)] # why does this return NA values for anon_id? 
multicases <- multicases[, .(anon_id, household.casenum, case.interval)]
setkey(multicases, anon_id)
setkey(cc.all, anon_id)

cc.all <- multicases[cc.all]
cc.all[, secondarycase := as.integer(case.interval > 4 & case.interval <= 14)]
    
procedures <- RDStodt(procedures.filename, keyname="anon_id") # FIXME: should have dates

##############################################################

timewin.admissions <- function(dt) {
    ## calculate exposure in recent and less recent time windows
    setkey(dt, anon_id, admissiondate, dischargedate)
    ## overlap join recent time intervals
    setkey(cc.specimendate, anon_id, daysbefore14, daysbefore5)
    cc.recent.inpat <- foverlaps(cc.specimendate,
                             dt[, .(anon_id, admissiondate, dischargedate)],
                             type="any", nomatch=NULL)
    cc.recent.inpat <- cc.recent.inpat[, .(anon_id, specimen_date)]
    cc.recent.inpat[, inpat.recent := 1]
    cc.recent.inpat <- unique(cc.recent.inpat)

    ## overlap join less recent time intervals
    setkey(cc.specimendate, anon_id, daysbefore24, daysbefore15)
    cc.lessrecent.inpat <- foverlaps(cc.specimendate,
                                 dt[, .(anon_id, admissiondate, dischargedate)],
                                 type="any", nomatch=NULL)
    cc.lessrecent.inpat <- cc.lessrecent.inpat[, .(anon_id, specimen_date)]
    cc.lessrecent.inpat[, inpat.lessrecent := 1]

    ## outer join of recent and less recent
    setkey(cc.recent.inpat, anon_id, specimen_date)
    setkey(cc.lessrecent.inpat, anon_id, specimen_date)
    inpat.timewin <- merge(cc.recent.inpat, cc.lessrecent.inpat, all=TRUE)
    setnafill(inpat.timewin, cols=c("inpat.lessrecent", "inpat.recent"), fill=0)
    inpat.timewin <- unique(inpat.timewin)
    setkeyv(inpat.timewin, c("anon_id", "specimen_date"))
    #inpat.timewin[, fD := .N > 1, by=key(inpat.timewin)]
    #inpat.timewin[fD==TRUE]
    
    ## define categoric variable for exposure time windows: less recent only, recent only, both
    inpat.timewin[, timewingr := group.tw(lessrecent=inpat.lessrecent,
                                                recent=inpat.recent)]
    return(inpat.timewin[, .(anon_id, specimen_date, timewingr)])
}

## define recent and less recent time windows before specimen date for overlap join with inpatient spells
cc.specimendate[, daysbefore24 := specimen_date - 24]
cc.specimendate[, daysbefore15 := specimen_date - 15]
cc.specimendate[, daysbefore14 := specimen_date - 14]
cc.specimendate[, daysbefore5  := specimen_date - 5]

###########################
## use rapid to derive time from specimen date to hospitalisation 
rapid <- RDStodt(rapid.filename, keyname="anon_id")
setnames(rapid, "admission_date", "AdmissionDate.rapid", skip_absent=TRUE)
setnames(rapid, "discharge_date", "DischargeDate.rapid", skip_absent=TRUE)
lastdate.rapid <- max(c(rapid$AdmissionDate.rapid, rapid$DischargeDate.rapid), na.rm=TRUE) 
rapid[is.na(DischargeDate.rapid), DischargeDate := lastdate.rapid]
rapid[, joindate := AdmissionDate.rapid]
setkey(rapid, anon_id, joindate)
setkey(cc.specimendate, anon_id, joindate)

## left join of cc.specimendate with next admission date in rapid
## so roll backwards up to 28 days from admission date
cc.nexthosp <- rapid[cc.specimendate, roll=-28] # roll backwards from admission date
## retain only records for which an admission date after specimendate was found in rapid
cc.nexthosp <- cc.nexthosp[!is.na(AdmissionDate.rapid)]
cc.nexthosp[, daystoadmission := as.integer(AdmissionDate.rapid - specimen_date)]
setorder(cc.nexthosp, anon_id, specimen_date, daystoadmission)
setkey(cc.nexthosp, anon_id, specimen_date)
cc.nexthosp <- unique(cc.nexthosp, by=key(cc.nexthosp))
setkey(cc.all, anon_id, specimen_date)
cc.all <- cc.nexthosp[, .(anon_id, AdmissionDate.rapid, specimen_date,
                          daystoadmission)][cc.all]
dim(cc.all)

cc.all[, censoringdays.admission := as.integer(
             pmin(lastdate.rapid, AdmissionDate.rapid, na.rm=TRUE) - specimen_date)]
## some cases have specimendate after lastdate.rapid
cc.all[censoringdays.admission < 0, censoringdays.admission := 0]

cc.all[, adm.within14 := as.integer(!is.na(daystoadmission) & daystoadmission <= 14)] 
cc.all[, adm.within28 := as.integer(!is.na(daystoadmission) & daystoadmission <= 28)]

#################################################

diagnoses <- RDStodt(diagnoses.filename, keyname="anon_id")
#diagnoses.last25 <- RDStodt(diagnoses.last25.filename, keyname="anon_id")
#diagnoses <- rbind(diagnoses, diagnoses.last25)
#rm(diagnoses.last25)
diagnoses <- unique(diagnoses)
setnames(diagnoses, "admission_date", "admissiondate", skip_absent=TRUE)
setnames(diagnoses, "discharge_date", "dischargedate", skip_absent=TRUE)
diagnoses[, admissiondate := as.Date(admissiondate)]
diagnoses[, dischargedate := as.Date(dischargedate)]
## set missing discharge dates to last date in table
lastdate.diagnoses <-  max(c(diagnoses$admissiondate, diagnoses$dischargedate), na.rm=TRUE)
diagnoses[is.na(dischargedate), dischargedate := lastdate.diagnoses]

################  recent hospital exposure  #################################
## 2 time windows of exposure: 5 to 14 days, 15 to 24 days before specimen date
## 
## recent hospital exposure: any exposure with
## admission.daysbefore >= 5 AND admission.daysbefore <= 14
## OR
## discharge.daysbefore >= 5 AND discharge.daysbefore <= 14
## OR
## admission.daysbefore > 14 AND discharge.daysbefore < 5  

##### 1. recent SMR00 outpatient visits
smr00 <- RDStodt(smr00.filename, keyname="anon_id") # outpatients
smr00[, clinic_date := as.Date(clinic_date)] # convert PosixCT
## Mode of Clinical Interaction identifies setting where contact between Health Care Professional and patient/carer takes place. further defined as a two-way interaction where the patient has the option to have subsequent dialogue with HCP.
## Codes and Values: 1 Face to Face, 2 Telephone, 3 Videolink, 4 Written
smr00 <- smr00[mode_of_clinical_interaction==1]
smr00[, joindate := clinic_date]
setkey(smr00, anon_id, joindate)
setkey(cc.specimendate, anon_id, joindate)
## left join of smr00 clinic date, rolling backwards from specimendate,
## so we specify a roll backwards of cc.specimendate with roll=-25
cc.recentsmr00 <- cc.specimendate[smr00, roll=-25] # roll backwards from specimendate
## retain only records for which a specimen date was found in the 25-day interval
cc.recentsmr00 <- cc.recentsmr00[!is.na(specimen_date)]
cc.recentsmr00[, outpatient.daysbefore := as.integer(specimen_date - clinic_date)]
cc.recentsmr00[outpatient.daysbefore >14 & outpatient.daysbefore <=24, outpat.recent := "days15to24"]
cc.recentsmr00[outpatient.daysbefore >4 & outpatient.daysbefore <=14, outpat.recent := "days5to14"]
## cast from long to wide
outpat.timewin <- dcast(cc.recentsmr00, anon_id + specimen_date ~ outpat.recent)
outpat.timewin[, opd.timewingr := group.tw(lessrecent=days15to24, recent=days5to14)]
outpat.timewin <- outpat.timewin[, .(anon_id, specimen_date, opd.timewingr)]
setkey(outpat.timewin, anon_id, specimen_date)
setkey(cc.all, anon_id, specimen_date)
cc.all <- outpat.timewin[cc.all]
with(cc.all, table(CASE, opd.timewingr))
rm(smr00)

## 2. recent daycase exposure
daycases <- diagnoses[inpatient_daycase_identifier=="D"]
daycases[, joindate := admissiondate]
setkey(daycases, anon_id, joindate)
setkey(cc.specimendate, anon_id, joindate)
## we want a right join of daycases admission date with next cc.specimendate,
## so we specify a roll backwards of cc.specimendate with roll=-25
cc.recentdaycases <- cc.specimendate[daycases, roll=-25] # roll backwards from specimendate
## retain only records for which a specimen date was found in the 25-day interval
cc.recentdaycases <- cc.recentdaycases[!is.na(specimen_date)]
cc.recentdaycases[, daycase.daysbefore := as.integer(specimen_date - admissiondate)]
cc.recentdaycases[daycase.daysbefore >14 & daycase.daysbefore <=24, daycase.recent := "days15to24"]
cc.recentdaycases[daycase.daysbefore >4 & daycase.daysbefore <=14, daycase.recent := "days5to14"]
## cast from long to wide
daycase.timewin <- dcast(cc.recentdaycases, anon_id + specimen_date ~ daycase.recent)
daycase.timewin[, daycase.timewingr := group.tw(lessrecent=days15to24, recent=days5to14)]
daycase.timewin <- daycase.timewin[, .(anon_id, specimen_date, daycase.timewingr)]
setkey(daycase.timewin, anon_id, specimen_date)
cc.all <- daycase.timewin[cc.all]
rm(daycases)

## 3. recent SMR04 psychiatric inpatient spells
smr04 <- RDStodt(smr04.filename, keyname="anon_id")
setnames(smr04, "admission_date", "admissiondate", skip_absent=TRUE)
setnames(smr04, "discharge_date", "dischargedate", skip_absent=TRUE)
smr04[, admissiondate := as.Date(admissiondate)]
smr04[, dischargedate := as.Date(dischargedate)]
## set missing discharge dates to last date in table
lastdate.smr04 <-  max(c(smr04$admissiondate, smr04$dischargedate), na.rm=TRUE)
smr04[is.na(dischargedate), dischargedate := lastdate.smr04]
psych.timewin <- timewin.admissions(smr04)
setnames(psych.timewin, "timewingr", "psych.timewingr")
setkey(psych.timewin, anon_id, specimen_date)
cc.all <- psych.timewin[cc.all] # FIXME: extra records

## 4. recent SMR01 inpatient spells
smr01.timewin <- timewin.admissions(diagnoses[inpatient_daycase_identifier=="I"])
setnames(smr01.timewin, "timewingr", "disch.timewingr")
setkey(smr01.timewin, anon_id, specimen_date)
cc.all <- smr01.timewin[cc.all]

## 5. recent RAPID inpatient spells 
## FIXME: use diagnoses to fill missing discharge dates in rapid where possible
setnames(rapid, "AdmissionDate.rapid", "admissiondate")
setnames(rapid, "DischargeDate.rapid", "dischargedate")
## set missing discharge dates to last date in table
lastdate.rapid <-  max(c(rapid$admissiondate, rapid$dischargedate), na.rm=TRUE)
rapid[is.na(dischargedate), dischargedate := lastdate.rapid]
rapid.timewin <- timewin.admissions(rapid)
setnames(rapid.timewin, "timewingr", "rapid.timewingr")
setkey(rapid.timewin, anon_id, specimen_date)
cc.all <- rapid.timewin[cc.all]

#######################################################

cc.all[, daycase.recent := as.logical(daycase.timewingr=="Recent TW only" | daycase.timewingr=="Both TWs")]
cc.all[, disch.recent := as.logical(disch.timewingr=="Recent TW only" | disch.timewingr=="Both TWs")]
cc.all[, rapid.recent := as.logical(rapid.timewingr=="Recent TW only" | rapid.timewingr=="Both TWs")]   
cc.all[, psych.recent := as.logical(psych.timewingr=="Recent TW only" | psych.timewingr=="Both TWs")]
cc.all[, opd.recent := as.logical(opd.timewingr=="Recent TW only" | opd.timewingr=="Both TWs")]

cc.all[is.na(daycase.recent), daycase.recent := FALSE]
cc.all[is.na(disch.recent), disch.recent := FALSE]
cc.all[is.na(rapid.recent), rapid.recent := FALSE]
cc.all[is.na(psych.recent), psych.recent := FALSE]
cc.all[is.na(opd.recent), opd.recent := FALSE]

cc.all[, hosp.recent := daycase.recent | disch.recent | rapid.recent | psych.recent | opd.recent]
cc.all[, inpat.recent := disch.recent | rapid.recent | psych.recent]

#########################################################


## ECDPC definitions
## definite hcai: admission.daysbefore >= 15 AND discharge.daysbefore < 0
## probable hcai: admission.daysbefore >= 8 AND discharge.daysbefore <= 14
## indeterminate hcai: admission.daysbefore >=3 AND discharge.daysbefore < 0

## PHE / LSHTM definition
## Community-Onset Suspected Healthcare-Associated (COSHA)
## admission.daysbefore >= -2 AND discharge.daysbefore <= 14

## our definition
## these definitions are mutually exclusive so can be coded as a categoric variable with 0 = no hcai
## definite hcai: admission.daysbefore >= 15 AND discharge.daysbefore < 0
## probable hcai: admission.daysbefore >= 8 AND discharge.daysbefore < 0
## indeterminate hcai: admission.daysbefore >=3 AND discharge.daysbefore < 0

## we want specimen dates that are within interval from admission to discharge in rapid
## x table is cc.specimendate, y table is rapid, type="within"

setkey(cc.specimendate, anon_id, specimen_date, joindate)
setkey(rapid, anon_id, admissiondate, dischargedate)
cc.inhosp <- foverlaps(cc.specimendate[, .(anon_id, specimen_date, joindate)],
                          rapid[, .(anon_id, admissiondate, dischargedate)],
                          type="within", nomatch=NULL)
cc.inhosp[, admission.daysbefore := as.integer(specimen_date - admissiondate)]
## define variable inhosp: tested positive on or after admission date and before discharge date
cc.inhosp[, inhosp := as.integer(admission.daysbefore >= 0 &
                                 as.integer(dischargedate - specimen_date) > 1)]
cc.inhosp[, hcai := car::recode(admission.daysbefore,
                           "0:2='Non-hospital onset'; 3:7='Indeterminate'; 8:14='Probable'; 15:hi='Definite'",
                           as.factor=TRUE,
                           levels=c("Definite", "Probable", "Indeterminate", "Non-hospital onset"))] ## ordering
setorder(cc.inhosp, -hcai) ## order so that most definite code for hcai is first
cc.hcai <- unique(cc.inhosp, by=c("anon_id", "specimen_date"))
## reorder levels now that we have retained most definite code
cc.hcai[, hcai := factor(hcai, levels=levels(hcai)[4:1])]

setkey(cc.hcai, anon_id, specimen_date)
print(nrow(cc.hcai[cc.all]))
cc.all <- cc.hcai[cc.all]
rm(cc.hcai)

cc.all[is.na(hcai), hcai := "Community onset"]
cc.all[, hcai := factor(hcai,
                               levels=c("Community onset",
                                        "Non-hospital onset", 
                                        "Indeterminate",
                                        "Probable", 
                                        "Definite"))]
cc.all[, prob.hcai := hcai == "Probable" | hcai =="Definite"]
print(nrow(cc.all))


cc.inhosp <- unique(cc.inhosp[inhosp==1, .(anon_id, specimen_date, inhosp)])
setkey(cc.inhosp, anon_id, specimen_date)
print(nrow(cc.inhosp[cc.all]))
cc.all <- cc.inhosp[cc.all]
rm(cc.inhosp)
setnafill(cc.all, cols="inhosp", fill=0)

#####################################################

source("icdchapters.R")

## use overlap join to assign chapter and subchapter to ICD diagnoses
## convert 3-char ICD code to integer for overlap join
diagnoses[, icdnum := icdToInt(icd10)]
diagnoses[, icdnum2 := icdnum]
setkey(diagnoses, icdnum, icdnum2)
setkey(icd.subchapters, startnum, endnum)
diagnoses <- foverlaps(diagnoses, icd.subchapters[, .(chnum, subchnum, startnum, endnum)])
diagnoses <- diagnoses[, .(anon_id, admissiondate, dischargedate, icd10, chnum, subchnum,
                           inpatient_daycase_identifier)]

## create table cc.diagnoses with 5-year lookback from specimen date for each anon_id / specimendate
cc.specimendate[, lookback5yr := specimen_date - floor(5 * 365.25)]
setkey(cc.specimendate, anon_id, lookback5yr, specimen_date)
setkey(diagnoses, anon_id, admissiondate, dischargedate)
cc.diagnoses <- foverlaps(diagnoses,
                          cc.specimendate[, .(anon_id, specimen_date, lookback5yr)],
                          type="within", nomatch=NULL)
rm(diagnoses)
cc.diagnoses[, days.sincedischarge := as.integer(specimen_date - dischargedate)]
cc.diagnoses <- cc.diagnoses[days.sincedischarge >= 15]
setorder(cc.diagnoses, days.sincedischarge)
cc.diagnoses <- unique(cc.diagnoses, by=c("anon_id", "specimen_date", "icd10")) 
setkey(cc.diagnoses, anon_id, specimen_date)
gc()

######################################
## smr06 <- RDStodt(smr06.filename, keyname="anon_id") # cancer incidence date and site
## FIXME -- move to morbidity coding section
#ids.icd.neoplasm.lastyear <-
#        unique(diagnoses[as.integer(specimen_date - DISCHARGE_DATE) + 25 > 365 &
#                     grepl("^C[0-9]|^D[0-4]", icd10), anon_id])

###############################################

cc.all[, hh.upto4 := uprn_0_4] 
cc.all[is.na(hh.upto4), hh.upto4 := 0] 
cc.all[, hh.upto4 := car::recode(hh.upto4, "NA=0; 3:hi = 3")]
cc.all[, hh.5to11 := uprn_5_11] 
cc.all[, hh.5to11 := car::recode(hh.5to11, "NA=0; 4:hi = 4")]
cc.all[, hh.12to17 := uprn_12_17] 
cc.all[, hh.12to17 := car::recode(hh.12to17, "4:hi = 4")]
cc.all[, hh.over18 := uprn_adult] ## mostly nonmissing
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

table.numadults.care <- with(cc.all[age_years > 70],
                             table(hh.over18 >= 10, care.home))

care.home.reclassified <- table.numadults.care[, 2]
cc.all[hh.over18 >= 10 & age_years >=70, care.home :="Care/nursing home"]

## set number of adults in household to 1 for care home residents 
cc.all[care.home=="Care/nursing home", hh.over18 := 1]

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
     
#cc.all[, in_hosp := as.factor(in_hosp)]
cc.all[, exp.group := "No exposure"]
cc.all[care.home=="Care/nursing home", exp.group :="Care home"]
cc.all[care.home=="Independent" & hosp.recent==TRUE,
           exp.group := "Recent hospital exposure"]
cc.all[, exp.group := factor(exp.group,
                                 levels=c("No exposure", "Care home",
                                        "Recent hospital exposure"))]

################# SICSAG ##########################

## merge SICSAG data and overwrite icu and hdu values incorrectly coded 0
sicsag <- RDStodt(sicsag.filename, keyname="anon_id")
sicsag[, admit_hosp := as.Date(admit_hosp)]
sicsag[, admit_unit := as.Date(admit_unit)]
sicsag[, disc_date := as.Date(disc_date)]
sicsag[date_disc_hosp=="", date_disc_hosp := NA]
sicsag[, date_disc_hosp := as.Date(date_disc_hosp)]
## get first and last date of entry in SICSAG table 
mindate.sicsag <- min(sicsag$admit_unit, na.rm=TRUE) 
maxdate.sicsag <- max(sicsag$admit_unit, na.rm=TRUE) 
## what are covid_icu_or_hdu codes 4 and 5?
sicsag <- sicsag[covid_icu_or_hdu==1 | covid_icu_or_hdu==3]
sicsag[, joindate := admit_hosp]
setkey(sicsag, anon_id, joindate)
setkey(cc.specimendate, anon_id, joindate)
## left join cc.specimendate, rolling backward from hospital admission date
cc.sicsag <- sicsag[cc.specimendate[, .(anon_id, specimen_date, joindate)], roll=-28]  ## 2880 records
cc.sicsag <- cc.sicsag[!is.na(admit_unit)]
cc.sicsag[, daystocriticalcare := as.integer(admit_unit - specimen_date)]

## about 40% of cases coded as COVID critical care in cc.all are not matched in SICSAG
#print(table(cc.all[icu==1 | hdu==1, anon_id] %in% cc.sicsag$anon_id))

setorder(cc.sicsag, by="daystocriticalcare")
setkey(cc.sicsag, anon_id, specimen_date)
cc.sicsag <- unique(cc.sicsag, by=key(cc.sicsag))
cc.all <- cc.sicsag[cc.all]
print(nrow(cc.all))

## 75% of admissions to critical care are within 10 days of specimen_date

## overwrite icu or hdu field where coded in sicsag
cc.all[covid_icu_or_hdu==1, icu := 1]
cc.all[covid_icu_or_hdu==3, hdu := 1]
## overwrite icu and hdu field where not coded in sicsag
cc.all[is.na(covid_icu_or_hdu), icu := 0]
cc.all[is.na(covid_icu_or_hdu), hdu := 0]

mindate.sicsag <- min(cc.all$admit_unit, na.rm=TRUE)
maxdate.sicsag <- max(sicsag$admit_unit, na.rm=TRUE)
      
cc.all[, censoringdays.critical := as.integer(
             pmin(maxdate.sicsag, admit_unit, na.rm=TRUE) - specimen_date)]
cc.all[, critical.within14 := as.integer(!is.na(daystocriticalcare) & daystocriticalcare <= 14)] 
cc.all[, critical.within28 := as.integer(!is.na(daystocriticalcare) & daystocriticalcare <= 28)] 

####################################################################################

cc.all[, dispensing.days := as.integer(specimen_date - as.Date("2019-06-01"))]

## exclude cases already dead before their specimen_date, and controls already dead on date of test of case they were matched to
cc.all <- cc.all[is.na(date_of_death) | date_of_death >= specimen_date]

## imputed specimen_date is 14 days before death
## but note that mode for days to death is 8 days
cc.all[, daystodeath := as.integer(date_of_death - specimen_date)]
## those with specimen_date later than maxdate.death should have censoringdays.death set to NA
cc.all[, censoringdays.death := as.integer(
             pmin(maxdate.death, date_of_death, na.rm=TRUE) - specimen_date)]
cc.all[censoringdays.death < 0, censoringdays.death := NA] 
cc.all[, death.within28 := as.integer(!is.na(daystodeath) & daystodeath <= 28)] 
cc.all[, death.within56 := as.integer(!is.na(daystodeath) & daystodeath <= 56)] 

## exclude cases classified as unobservable, for consistency with controls
## FIXME -- make this work with the field now in cc.all
if(!is.null(chistatus.filename)) {
    cases.status <-
        RDStodt(chistatus.filename, keyname="anon_id")
    setnames(cases.status, "CHI_EXTENDED_STATUS", "EXTENDED_STATUS", skip_absent=TRUE)
    cc.all <- cases.status[cc.all]
    with(cc.all[CASE==1],
     table(EXTENDED_STATUS, is.na(date_of_death)))
    cc.all <- cc.all[is.na(EXTENDED_STATUS) |
                     EXTENDED_STATUS == "C" |
                     (EXTENDED_STATUS == "D" & !is.na(date_of_death))]
}
    
## FIXME: move this to beginning, ?delete
## rows have been replicated within strata containing N cases so that each case and control appears N^2 times
## remove these replicated rows
cc.all <- distinct(.data=cc.all, anon_id, stratum, .keep_all = TRUE)
setkey(cc.all, anon_id)
cc.all <- copy(cc.all) # force a deep copy
gc()
print(nrow(cc.all))

num.casectrl.strata <- cc.all[, .N, by=c("stratum", "CASE")]
## we want a table with one row per stratum, separate columns for numcases and numctrls
num.casectrl.strata <- dcast(num.casectrl.strata, stratum ~ CASE, value.var="N")
colnames(num.casectrl.strata) <- c("stratum", "controls", "cases")
with(num.casectrl.strata, print(table(controls, cases)))
 
######################################################################
## classification of cases into 4 groups
cat("Classifying cases into 4 groups ...")
## 8 individuals without a positive test result are included as cases but do not have
## covid_cod==1 or diag.case==1

## all cases have nonmissing SPECIMENDATE
## controls should be assigned same SPECIMENDATE as the case they were matched to

## deathwithin28 is a logical variable that literally means death within 28 days of a positive test
## should be coded as FALSE for those who did not test positive
cc.all[, deathwithin28 := testpositive.case &
             !is.na(date_of_death) & date_of_death - specimen_date <= 28]

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
## should include those admitted within 28 days (adm.within28=1) AND those already in hospital on date they tested positive

cc.all[, adm28 := as.integer((adm.within28==1 | inhosp==1) & testpositive.case)]

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
setkey(cc.all, anon_id, specimen_date)

## for cases, overwrite the casegroup field with the group assigned above
cc.all[CASE==1, casegroup := group]
table(cc.all$CASE, cc.all$casegroup, exclude=NULL)

## drop records with missing casegroup (controls in strata with no remaining classified case)
cc.all <- cc.all[!is.na(casegroup)]
cc.all[, casegr := ifelse(is_case & is.na(severe.casegr),
                          "Not severe", as.character(severe.casegr))]
cc.all[casegr == "Not severe" & !is.na(admit_unit), casegr := "Not severe, hospitalized"]
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
cat("done\n")

with(cc.all[casegroup=="A"], print(table(CASE, daycase.timewingr, exclude=NULL)))
with(cc.all[casegroup=="A"], print(table(CASE, disch.timewingr, exclude=NULL)))
with(cc.all[casegroup=="A"], print(table(CASE, rapid.timewingr, exclude=NULL)))
with(cc.all[casegroup=="A"], print(table(CASE, psych.timewingr, exclude=NULL)))
with(cc.all[casegroup=="A"], print(table(CASE, opd.timewingr, exclude=NULL)))

####################################################
## create hosp.recent variable from disch, rapid, smr04, smr00
    ## hosp.recent if any exposure in recent time window

####################################################################

## import most recent SMR01 discharge date into rapid
cc.all[qSIMD.integer==9, qSIMD.integer := NA]

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
narrow.case.freqs <- with(cc.all[narrow.case], table(age_years, sex, exclude=NULL))
broad.case.freqs <- with(cc.all[broad.case], table(age_years, sex, exclude=NULL))
narrow.death.freqs <- with(cc.all[narrow.fatalcase], table(age_years, sex, exclude=NULL))
broad.death.freqs <- with(cc.all[broad.fatalcase], table(age_years, sex, exclude=NULL))
save(narrow.case.freqs, broad.case.freqs,
     narrow.death.freqs, broad.death.freqs,
     file=paste0(datadir, "casefreqs.4cats.agesex.RData"))

#source("incidencemortality.R")
rm(narrow.case)
rm(narrow.fatalcase)
rm(broad.fatalcase)
rm(broad.case)

######## coding ethnicity ##############################

source("ethnic_assign.R")

cc.all[, ethnic5.smr := collapseto5.ethnicsmr(ethnic_smr_last)]
## recode SMR ethnicity to 4 categories: White, Black, South Asian, Other
cc.all[, ethnic4.smr := as.factor(car::recode(ethnic5.smr,
                                              "'Chinese'='Other'"))]
cc.all[, ethnic4.smr := factor(ethnic4.smr,
                               levels=levels(ethnic4.smr)[c(4, 3, 1, 2)])]   

setkey(cc.all, anon_id, specimen_date)
######################################################

rm(rapid)
gc()

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

#################### import shielding data ####
        
shielded.full <- fread(shielding.full.filename)
setkey(shielded.full, anon_id)
setnames(shielded.full, "group", "shield.group", skip_absent=TRUE)
setnames(shielded.full, "batch", "shield.batch", skip_absent=TRUE)
setnames(shielded.full, "age", "age_years", skip_absent=TRUE)
shielded.full[, sex := as.factor(car::recode(sex, "1='Male'; 2='Female'"))]
setkey(shielded.full, shield.batch)
shielded.full <- batches[shielded.full]
shielded.full[, removal := nchar(removal)]

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
## remove obvious wrong assignments
shielded.full <- shielded.full[!(shield.group == "Pregnant with heart disease" &
                                 (sex=="Male" | age_years > 55))]

## left join of cc.all with subset of shielded.full in which anon_id is nonmissing
setkey(shielded.full, anon_id)
cc.all <- shielded.full[, .(anon_id, shield.batch, Date.Sent,
                            shielding_id, shield.group)][cc.all]

cc.all[, shield.any := as.factor(!is.na(shield.group))]
cc.all[is.na(shield.group), shield.group := "Ineligible for shielding"]
cc.all[, shieldelig.group := car::recode(shield.group,
                                          "'Pregnant with heartdisease'='Additional conditions'",
                                          as.factor=TRUE,
                                          levels=c(
                                              "Ineligible for shielding",
                                              "Solid organ transplant",
                                              "Specific cancers",
                                              "Severe respiratory",
                                              "Rare diseases",
                                              "On immunosuppressants",
                                              "Additional conditions"
                                          ))]

cc.all[is.na(shield.batch), shield.batch := 0]
cc.all[, shield.batch := as.factor(shield.batch)]

## import case status into shielded.full

cases <- cc.all[CASE==1, .(anon_id, casegr, casegr2, casegr3)]
setkey(cases, anon_id)
shielded.full <- cases[shielded.full]

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
shielded.full[, agegr20 := as.factor(car::recode(as.integer(age_years),
                                                     "0:39='0-39'; 40:59='40-59';
                               60:74='60-74'; 75:hi='75 or more'"))]

## rm(controls.deceased)

## read raw scrips file in separate script
## read scrips file and import BNF chapter variables #############
rm(rapid)
rm(cc.recent.smr00)
colsToDelete <- grep("daysbefore", names(cc.specimendate), value=TRUE)
set(cc.specimendate, , colsToDelete, NULL)
gc()
objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

bnfcodes <- as.data.table(read_excel("./BNF_Code_Information.xlsx", sheet=4))
colnames(bnfcodes) <- c("chaptername", "chapternum", "sectionname", "sectioncode")
bnfcodes[, `:=`(chapternum=as.integer(chapternum), sectioncode=as.integer(sectioncode))]
bnfchapters <- bnfcodes[, 1:2]
bnfchapters <- base::unique(bnfchapters)
   
# source("cc.scrips.R")
load(paste0(datadir, "cc.scripsbnf.RData"))

## BNF chapter codes: up to 2 digits, leading 0 may be stripped
## section codes: 2 digits
## paragraph codes: 2 digits
## subpara codes: 1 digit
## chemical substance code: 2 characters, may include a letter

## number of distinct BNF subparas with at least one prescription in chapter 02 cardiovascular
cc.numdrugs.cv <- cc.scripsbnf[substr(bnf_item_code, 1, 2) == "02",
                                  .(numdrugs.cv = uniqueN(substr(bnf_item_code, 1, 7))),
                                  by=c("anon_id", "specimen_date")]
setorder(cc.numdrugs.cv, -numdrugs.cv)
setkey(cc.numdrugs.cv, anon_id, specimen_date)
cc.numdrugs.cv <- unique(cc.numdrugs.cv, by=key(cc.numdrugs.cv))
setkey(cc.all, anon_id, specimen_date)
cc.all <- cc.numdrugs.cv[cc.all]
setnafill(cc.all, cols="numdrugs.cv", fill=0) 

## number of distinct BNF subparas with at least one prescription, excluding chapter 02 cardiovascular
cc.numdrugs.notcv <- cc.scripsbnf[substr(bnf_item_code, 1, 2) != "02",
                                  .(numdrugs.notcv = uniqueN(substr(bnf_item_code, 1, 7))),
                                  by=c("anon_id", "specimen_date")]
setorder(cc.numdrugs.notcv, -numdrugs.notcv)
setkey(cc.numdrugs.notcv, anon_id, specimen_date)
cc.numdrugs.notcv <- unique(cc.numdrugs.notcv, by=key(cc.numdrugs.notcv))
setkey(cc.all, anon_id, specimen_date)
cc.all <- cc.numdrugs.notcv[cc.all]
rm(cc.numdrugs.notcv)
setnafill(cc.all, cols="numdrugs.notcv", fill=0) 

cc.all[, numdrugs := numdrugs.cv + numdrugs.notcv]

cc.numdrugs.cns <- cc.scripsbnf[substr(bnf_item_code, 1, 2) == "04",
                                  .(numdrugs.cns = uniqueN(substr(bnf_item_code, 1, 7))),
                                  by=c("anon_id", "specimen_date")]
setorder(cc.numdrugs.cns, -numdrugs.cns)
setkey(cc.numdrugs.cns, anon_id, specimen_date)
cc.numdrugs.cns <- unique(cc.numdrugs.cns, by=key(cc.numdrugs.cns))
setkey(cc.all, anon_id, specimen_date)

cc.all <- cc.numdrugs.cns[cc.all]
setnafill(cc.all, cols="numdrugs.cns", fill=0) 

## 103 is proton pump inhibitors 01 03
## 1001030 is rheumatic disease suppressants 10 01 03 0
## 1001030C0 is hydroxychloroquine
## 1001030U0 is methotrexate

## recode scrips$bnf.chapter values > 14 or NA to 14
#scrips[is.na(chapternum), chapternum := 14]
#scrips[chapternum > 14, chapternum := 14] 
## colchicine
## BNF Subparagraph Code  1001040
## BNF Chemical Substance 1001040G0
## always prescribed as 500 mcg tablets, quantity usually 28, 56 or 100
## dose 2-4 tablets/day so at 2 tabs/day that corresponds to 14, 28 or 50 days supply
## gabapentin 0408010G0 pregabalin 0408010AE - classified in subpara 0408010 Control of epilepsy

source("comorbidity.R")

setkey(cc.all, anon_id, specimen_date)
cc.all <- cc.comorbid[cc.all]
rm(cc.comorbid)
rm(cc.diabetes) 
gc()

########## coding diabetes ####################

###  diabetes based on combining Sci-Diabetes records with ICD codes and drugs 
## add in BNF codes 6.1 for diabetes drugs and
## E10 to E14 for diabetes
    
## missing recoded as zero
cc.all[is.na(dm_type), dm_type := 0]
## add in extra cases notified directly from SCI-Diabetes register, without assignment
## of diabetes type from SDRN database
cc.all[dm_type==0 & diab_reg==1, dm_type := 3]
## code diagnoses detected from discharges or BNF codes as unknown type
## we could classify those not on insulin as definite Type 2 
## REVISION: for consistency with the diabetes paper, we will not include the extra cases identified through diagnostic codes or drug codes - they may be transient/resolved
## recode diabetes type
cc.all[, dm_type := recode.dmtype(dm_type)]

## define indicator variable for any diabetes
cc.all[, diabetes.any := as.integer(dm_type != "Not diabetic")]
cc.all[, diabetes.any := as.factor(car::recode(diabetes.any,
                                                "0='Not diabetic'; 1='Diabetic'"))]
cc.all[, diabetes.any := relevel(diabetes.any, ref="Not diabetic")]
cc.all[, dm_type3 := car::recode(dm_type, "'Other/unknown type'='Type 2 diabetes'",
                                  levels=c("Not diabetic", "Type 1 diabetes",
                                           "Type 2 diabetes"))]

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))
gc()

cc.all[, listed.any := 
              as.factor(as.integer(diabetes.any=="Diabetic" | ihd==1 |
                                   heart.other==1 | ckd==1 | chronresp==1 |
                                   neuro==1 | liver==1 | immune==1))]

cc.all[shield.group=="Ineligible for shielding",
            shield.group := ifelse(listed.any==1,
                                   "Moderate risk condition",
                                   "No risk condition")]
cc.all[, shield.group := car::recode(shield.group,
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
cc.all[, shieldedonly.group := car::recode(shield.group,
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

cc.all[, listedgr3 := 0]
cc.all[listed.any=="1", listedgr3 := 1]
cc.all[shield.any==TRUE, listedgr3 := 2]
cc.all[, listedgr3 := car::recode(listedgr3,
                                       recodes=
                                           "0='No risk condition';
                                       1='Moderate risk condition';
                                       2='Eligible for shielding'",
                                       as.factor=TRUE, 
                                       levels=c("No risk condition",
                                                "Moderate risk condition",
                                                "Eligible for shielding"))]

save(cc.all, file=paste0(datadir, "cc.all.RData"))

###################################

Rprof()
print(summaryRprof(tmp)$by.total[1:20, ])

source("rrtimewin.R")

rmarkdown::render("vaxshieldupdate.Rmd")

## TODO: procedures, tabulate reinfections,  

if(FALSE) {

## generate and merge indicators for any diag in each  ICD chapter
icdchapters.anydiag <- diagnoses[, .N, by=c("anon_id", "chnum")] %>%
    dcast(anon_id ~ chnum, value.var="N")
## recode as 0, 1
setnafill(icdchapters.anydiag, cols=2:ncol(icdchapters.anydiag), fill=0)
for (j in 2:ncol(icdchapters.anydiag)) set(icdchapters.anydiag, j=j,
                                           value=as.integer(icdchapters.anydiag[[j]] > 0))
## use Roman numerals for ICD chapters
chnums <- as.integer(colnames(icdchapters.anydiag)[-1])
colnames(icdchapters.anydiag)[-1] <-
    paste0("Ch.", as.roman(chnums), "_",
           icdchapters$shortname[match(chnums, icdchapters$chnum)])
## drop rare ICD chapters (freq < 1%)
keep.cols <- colSums(icdchapters.anydiag) >= 0.01 * nrow(icdchapters.anydiag) ## drops perinatal chapter
icdchapters.anydiag <- icdchapters.anydiag[, ..keep.cols]
setkey(icdchapters.anydiag, anon_id)
cc.kept <- icdchapters.anydiag[cc.kept]
rm(icdchapters.anydiag)
icdcols <- grep("^Ch\\.", colnames(cc.kept)) # begins with Ch._
setnafill(cc.kept, cols=icdcols, fill=0)

## cc.kept[, (icdcols)] returns the column names as a vector
## cc.kept[, .(icdcols)] returns the column names as a data.table
## cc.kept[, ..icdcols] returns a data.table selected by column

## derive number of chapters with at least one diag 
x <- cc.kept[, ..icdcols]
cc.kept[, num.icdchapters := rowSums(x)]
rm(x)
cc.kept[, num.icdchapters.gr := as.factor(car::recode(num.icdchapters,
                                                      "0='No discharge records';
                         1:2='1-2 ICD-10 chapters';
                        3:hi='3 or more chapters'"))]
cc.kept[, num.icdchapters.gr := relevel(num.icdchapters.gr, ref="No discharge records")]

################################################################################

## generate and merge indicators for any diag in each  ICD subchapter
icdsubchapters.anydiag <- diagnoses[, .N, by=c("anon_id", "subchnum")] %>%
    dcast(anon_id ~ subchnum, value.var="N")
## recode as 0, 1
setnafill(icdsubchapters.anydiag, cols=2:ncol(icdsubchapters.anydiag), fill=0)
for (j in 2:ncol(icdsubchapters.anydiag)) set(icdsubchapters.anydiag, j=j,
                                              value=as.integer(icdsubchapters.anydiag[[j]] > 0))
subchnums <- as.integer(colnames(icdsubchapters.anydiag)[-1])
colnames(icdsubchapters.anydiag)[-1] <-
    paste0("Ch_", 
           as.roman(icd.subchapters$chnum[match(subchnums, icd.subchapters$subchnum)]), ": ",
           icdsubchapters$start[match(subchnums, icd.subchapters$subchnum)], "-",
           icdsubchapters$end[match(subchnums, icd.subchapters$subchnum)], " ", 
           icdsubchapters$name[match(subchnums, icd.subchapters$subchnum)])

## drop rare ICD subchapters (freq < 0.1%)
keep.cols <- colSums(icdsubchapters.anydiag) >= 0.001 * nrow(icdsubchapters.anydiag)  
icdsubchapters.anydiag <- icdsubchapters.anydiag[, ..keep.cols]
setkey(icdsubchapters.anydiag, anon_id)
cc.kept <- icdsubchapters.anydiag[cc.kept]
rm(icdsubchapters.anydiag)
icdcols <- grep("^Ch_", colnames(cc.kept)) ## no period after Ch 
setnafill(cc.kept, cols=icdcols, fill=0)

## derive number of subchapters with at least one diag 
x <- cc.kept[, ..icdcols]
cc.kept[, num.icdsubchapters := rowSums(x)]
rm(x)
rm(diagnoses)
rm(procedures)
gc()

## calculate risk score in controls and set intercept to equate observed and expected severe cases
if(TRUE) {
    #source("trainmodel.R")
    source("riskscore.R")

    ## thin beta.draws to have  < 1000 samples
    keep.draws <- seq(1, nrow(beta.draws), by=ceiling(nrow(beta.draws) / 1000))
    beta.draws <- beta.draws[keep.draws, ]
    
    ## calculate unnormalized risk score
    ## FIXME: impute missing values so that risk score can be calculated for everyone 
    covariate.names.agesex <- c("anon_id", "age_years", "sex", covariate.names) 
    cc.predict <- na.omit(cc.kept[, ..covariate.names.agesex])
    setnames(cc.predict, "age_years", "Age")
    setnames(cc.predict, "sex", "Sex")
    
    cat("calculating unnormalized logits with matrix multiplication ...")
    X.kept <- model.matrix(object=~ ., data=cc.predict[, ..covariateonly.names])[, -1]
    cc.predict <- cc.predict[, .(anon_id, Age, Sex)] # to save memory
    unnorm.lograte.draws <- X.kept %*% t(beta.draws) # this is a large object
    rm(X.kept)
    cat("done\n")
 
    objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
    print(tail(objmem))
    ## without thinning, this will use a lot of memory
    cat("calculating unnormalized logits of mean in arithmetic basis ...")
    unnorm.lograte <- matrixStats::rowLogSumExps(unnorm.lograte.draws) -
    log(ncol(unnorm.lograte.draws))
    rm(unnorm.lograte.draws)
    cc.predict[, logit.unnorm := ..unnorm.lograte]

    ## merge with dt.lookup to get gam.logit and calculate logit.norm 
    setkey(cc.predict, Age, Sex)
    setkey(dt.lookup, COVID.age, Sex)
    cc.predict <- dt.lookup[cc.predict]
    cc.predict[, logit.norm := ..intercept.gam.unnorm + gam.logit + logit.unnorm]
    setnames(cc.predict, "COVID.age", "Age")
    cc.predict <- cc.predict[, .(anon_id, Sex, logit.norm)]
    
    covidage <- getCOVIDage(x=cc.predict, y=dt.lookup)
    setkey(covidage, anon_id)
    cc.kept <- covidage[cc.kept]
}

#######################################################################################

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

print(table(cc.kept$CASE, cc.kept$casegroup))

cc.severe <- cc.kept[casegroup=="A"]
cc.lt75 <- cc.severe[age_years < 75]
save(cc.lt75, file=paste0(datadir, "cc.severe.lt75.RData"))
rm(cc.lt75)

###########################################

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))


#####################################################

## reinfections appear in wave==1 with SPECIMENDATE matching SpecimenDate1

#load(paste0(datadir, "cc.all.RData"))
with(cc.all, table(wave, specimen_date==SpecimenDate1))

cc.all[, reinfected := !is.na(SpecimenDate1)]

cc.all[, occup3 := car::recode(occup, "'Teacher'='Other / undetermined'", as.factor=TRUE, 
                               levels=levels(occup)[-2])]

## get reinfections into wave==2
wave2 <- cc.all[wave==2 | specimen_date==SpecimenDate1]
wave2[specimen_date==SpecimenDate1, specimen_date := SpecimenDate2]
wave2[, reinfected := !is.na(SpecimenDate2) & specimen_date == SpecimenDate2]
wave2[, death28 := !is.na(date_of_death) & date_of_death - specimen_date < 28]

table.wave2 <- with(wave2[CASE==1], table(reinfected))

## what factors in wave 1 predict reinfection in wave 2
varnames=c("agegr3", "sex", "care.home", "shield.any", "occup3")
table.reinfection <- tabulate.freqs.regressions(varnames=varnames, outcome="reinfected", 
                           data=cc.all[wave==1 & CASE==1 & fatalcase==0],
                           model="logistic")
colnames(table.reinfection) <- gsub("FALSE", "No reinfection", colnames(table.reinfection))
colnames(table.reinfection) <- gsub("TRUE", "Reinfection", colnames(table.reinfection))

## reinfection predicts fatal cases in wave 2
varnames.fatal <- c("age_years", "sex", "care.home", "shield.any", "occup3", "reinfected")
table.fatal <- tabulate.freqs.regressions(varnames=varnames.fatal, outcome="death28", 
                           data=wave2[age_years >= 60],
                           model="logistic")
colnames(table.fatal) <- gsub("FALSE", "Non-fatal", colnames(table.fatal))
colnames(table.fatal) <- gsub("TRUE", "Fatal", colnames(table.fatal))

cc.kept[CASE==1 & adm.within28==1 & specimen_date >= as.Date("2021-01-01"), mean(COVID.age, na.rm=TRUE), by=lubridate::month(specimen_date)]

paste.colpercent(with(cc.kept[CASE==1 & adm.within28==1 & specimen_date >= as.Date("2021-01-01")],
     table(criticalcare, lubridate::month(specimen_date))), digits=1)
}
