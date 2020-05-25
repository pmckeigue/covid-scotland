## analysis script for case-control study
rm(list=ls())

library(car)
library(survival)
library(MASS)
library(wevid)
library(rmarkdown)
library(pander)
library(ggplot2)
library(doParallel)
library(reshape2)
library(readxl)
library(DescTools)
library(icd.data)
library(gam)
# library(plyr) # ? needed to access .() function
library(dplyr)

registerDoParallel(cores=2)

source("helperfunctions.R")

scotpop <- read_excel("./Scotland_midyearpop_est2019.xlsx")

old <- TRUE
#old <- FALSE

if(old) {
    cc.all <- readRDS("./data/CC_linked_ANON_20200501 (2).rds")
    diagnoses <- readRDS("./data/CC_SMR01_ICD10_x25_ANON_20200501.rds") # 842 records ? excluded
    procedures <- readRDS("./data/CC_SMR01_OPCS4_MAIN.x25_ANON_20200501.rds")
    scrips <- readRDS("./data/CC_PIS_x15_ANON_20200501.rds")
} else {
    cc.all <- readRDS("./data/CC_linked_ANON_20200515.rds")
    diagnoses <- readRDS("./data/CC_SMR01_ICD10_x25_ANON_20200515.rds")
    procedures <- readRDS("./data/CC_SMR01_OPCS4_MAIN.x25_ANON_20200515.rds") 
    scrips <- readRDS("./data/CC_PIS_x15_ANON_20200515.rds")[, c("ANON_ID",
                                                           "bnf_paragraph_code",
                                                           "bnf_paragraph_description")] 
    ## scrips should be 2 tables to save space
    ## one record per scrip
    scripvars <- c("ANON_ID", "dispensed_date", "bnf_paragraph_code",
                   "formulation_code",
                   "item_strength", "item_strength_uom", "item_code",
                   "num_items", "quantity")
    ## lookup table for subpara code and item code
    drugvars <- c("bnf_paragraph_code", "bnf_paragraph_description",
                  "item_code", "approved_name")
    ## other files for diagnoses, procedures and scrips within the time limit
    ## data/CC_SMR01_ICD10_25_ANON_20200515.rds
    ## data/CC_SMR01_OPCS4_MAIN.25_ANON_20200501.rds
    ## data/CC_PIS_15_ANON_20200515.rds

    scrips.protonpump <-
        readRDS("./data/CC_PIS_x15_ANON_20200515.rds") %>%
        subset(bnf_paragraph_code=="0103050")
}

###################################################################

## for now, just keep id, paragraph code, paragraph_description
scrips$bnf_paragraph_description <- as.factor(scrips$bnf_paragraph_description)
cat("scrips object uses", object.size(scrips) * 1E-6, "MB\n")

## we have 7 digits on scrips, giving resolution to subpara level only

length(table(substr(scrips$bnf_paragraph_code, 1, 2))) # chapter
length(table(substr(scrips$bnf_paragraph_code, 1, 4))) # chapter, section
length(table(substr(scrips$bnf_paragraph_code, 1, 6)))  # chapter, section, paragraph
length(table(scrips$bnf_paragraph_code)) # 537 groups

## we need integer variables chapter, sectioncode, paracode for use with reshape2::dcast
scrips$chapternum <- as.integer(substr(scrips$bnf_paragraph_code, 1, 2))
scrips$sectioncode <- as.integer(substr(scrips$bnf_paragraph_code, 1, 4))
scrips$paracode <- as.integer(substr(scrips$bnf_paragraph_code, 1, 6))

## recode scrips$bnf.chapter values > 14 or NA to 14
scrips$chapternum[is.na(scrips$chapternum)] <- 14
scrips$chapternum[scrips$chapternum > 14] <- 14

################################################################################

icdchapters <- data.frame(names(icd10_chapters),
                             t(matrix(as.character(unlist(icd10_chapters)), nrow=2)))
colnames(icdchapters) <- c("name", "start", "end")
icdchapters$shortname <- gsub("Diseases of the ", "", icdchapters$name)
icdchapters$shortname <- gsub("Certain ", "", icdchapters$shortname)
icdchapters$shortname <- gsub("conditions originating in the ", "",
                              icdchapters$shortname)
icdchapters$shortname <- gsub("Factors influencing ", "", icdchapters$shortname)
truncate.at <- StrPos(icdchapters$shortname, " |,|-") - 1
truncate.at[is.na(truncate.at)] <- nchar(icdchapters$shortname[is.na(truncate.at)])
icdchapters$shortname <- substr(icdchapters$shortname, 1, truncate.at) 
icdchapters$start <- as.character(icdchapters$start)
icdchapters$end <- as.character(icdchapters$end)

icdsubchapters <- data.frame(names(icd10_sub_chapters),
                             t(matrix(as.character(unlist(icd10_sub_chapters)), nrow=2)))
colnames(icdsubchapters) <- c("name", "start", "end")
icdsubchapters$start <- as.character(icdsubchapters$start)
icdsubchapters$end <- as.character(icdsubchapters$end)

source("bnfcodes.R")

###############################################

names(cc.all) <- gsub("CASE_NO", "stratum", names(cc.all))
names(cc.all) <- gsub("^SEX$", "sex", names(cc.all))
names(cc.all) <- gsub("imumune", "immune", names(cc.all))
names(cc.all) <- gsub("^ethnic$", "ethnic.old", names(cc.all))
names(cc.all) <- gsub("^CAREHOME$", "care.home", names(cc.all))
names(cc.all) <- gsub("^simd$", "SIMD.quintile", names(cc.all))
names(cc.all) <- gsub("DATE_OF_DEATH", "Date.Death", names(cc.all))
names(cc.all) <- gsub("^age$", "AGE", names(cc.all))

## HAI is based on the ECDC definition of nosocomial infection

## exclude controls already dead on date of test of case they were matched to
controls.deceased <- with(cc.all, CASE==0 &
                                  !is.na(Date.Death) &
                                  Date.Death <= SPECDATE)
cc.all <- cc.all[!controls.deceased, ]                         

cc.all$stratum <- as.integer(cc.all$stratum)
## check that each stratum contains a single case
cat("checking that each stratum contains a single case ...")
table.strata <- tapply(cc.all$CASE, cc.all$stratum, sum) == 1
strata.onecase <- as.integer(names(table.strata)[as.integer(table.strata)==1])
keep <- cc.all$stratum %in% strata.onecase
cc.all <- cc.all[keep, ]
cat("done:", length(which(!keep)), "observations dropped\n")

cc.all$SIMD.quintile <- car::recode(cc.all$SIMD.quintile, "'Unknown'=NA")

cc.all$scrip.any <- as.factor(as.integer(cc.all$ANON_ID %in% scrips$ANON_ID))
cc.all$diag.any <- as.factor(as.integer(cc.all$ANON_ID %in% diagnoses$ANON_ID))

cc.all$scripordiag <- as.integer(with(cc.all, as.integer(diag.any)==2 |
                                                    as.integer(scrip.any)==2))
                  
cc.all$sex <- car::recode(as.factor(cc.all$sex), "1='Male'; 2='Female'")
cc.all <- within(cc.all, sex <- relevel(sex, ref="Female"))

cc.all$agegr20 <- as.factor(car::recode(as.integer(cc.all$AGE),
                              "0:39='0-39'; 40:59='40-59';
                                  60:74='60-74'; 75:hi='75 or more'"))

cc.all$agegr3 <-
    as.factor(car::recode(cc.all$AGE,
                          "0:59='0-60 years'; 60:74='60-74 years'; 75:hi='75+ years'"))

cc.all$care.home[cc.all$NURSINGHOME==1] <- 1
cc.all$care.home <- as.factor(car::recode(cc.all$care.home, "0='Independent'; 1='Care/nursing home'"))
cc.all <- within(cc.all, care.home <- relevel(care.home, ref="Independent"))

## all cases have nonmissing SPECDATE
## controls are recoded to have same SPECDATE as the case they were matched to
cc.all$deathwithin28 <- 0
cc.all$deathwithin28[with(cc.all,
                          CASE==1 &
                          !is.na(Date.Death) & 
                          Date.Death - SPECDATE >= 0 &
                          Date.Death - SPECDATE <= 28)] <- 1

## fatality rates by type of unit
## 34% in icu, 38% in icu.hdu.ccu,
## 24% in those in hdu but never in icu or icu.hdu.ccu (a small group)
cc.all$unitcategory <- numeric(nrow(cc.all))
cc.all$unitcategory[cc.all$icu==0 & cc.all$hdu==0 & cc.all$adm28==1] <- 0
cc.all$unitcategory[cc.all$icu==0 & cc.all$hdu==1] <- 1
cc.all$unitcategory[cc.all$icu==1 & cc.all$hdu==0] <- 2
cc.all$unitcategory[cc.all$icu==1 & cc.all$hdu==1] <- 3
cc.all$unitcategory[cc.all$CASE==0] <- NA
cc.all$unitcategory <- car::recode(cc.all$unitcategory,
                                   "0='No HDU or ICU'; 1='HDU only'; 2='ICU only'; 3='HDU and ICU'")

print(with(cc.all[cc.all$CASE==1, ], paste.colpercent(table(deathwithin28, unitcategory))))

## check this is correct: all those with icu==1 or icu.hdu.ccu==1 should be coded as severe

## integer values > 1 for icu and inhosp may represent days from test to entry
## values of 0 must be for those not admitted, as there are no missing values
with(cc.all[cc.all$CASE==1, ], table(adm28, exclude=NULL))
with(cc.all[cc.all$CASE==1, ], table(icu, exclude=NULL))
if(!old) {
    with(cc.all[cc.all$CASE==1, ], table(hdu, exclude=NULL))
    with(cc.all[cc.all$CASE==1, ], table(icu, hdu, exclude=NULL))
}

## with(cc.all[cc.all$CASE==1, ], table(icu.hdu.ccu, exclude=NULL))

## coding of case groups -- check this is correct
cc.all$group <- NA
cc.all$group[cc.all$CASE==1 & cc.all$adm28== 0 & cc.all$inhosp==0] <- "C"
cc.all$group[cc.all$CASE==1 & (cc.all$adm28==1 | cc.all$inhosp==1)] <- "B"
cc.all$group[cc.all$CASE==1 & (cc.all$icu==1 | cc.all$deathwithin28==1)]  <- "A"
table(cc.all$deathwithin28, cc.all$group, exclude=NULL)

## assign controls to same group as matched case i.e. A, B, C and create a new variable named casegroup
casegroups <- cc.all[cc.all$CASE==1, ][, c("stratum", "group")]
colnames(casegroups)[2] <- "casegroup"
cc.all <- merge(cc.all, casegroups, by=c("stratum"), all.x=T)
table(cc.all$CASE, cc.all$casegroup)
with(cc.all[cc.all$CASE==1, ], table(casegroup, deathwithin28, exclude=NULL))

cc.all$fatalcase <- as.integer(cc.all$CASE==1 & cc.all$deathwithin28==1)

cc.all$casegroup <- car::recode(cc.all$casegroup,
                           "'A'='Critical care or fatal'; 'B'='Hospitalised, not severe'; 'C'='Test-positive, not hospitalised'")
cc.all$casegroup <- as.factor(cc.all$casegroup)
cc.all <- within(cc.all, casegroup <- relevel(casegroup, ref="Critical care or fatal"))

######## coding ethnicity ##############################

OnolyticsType <- cc.all$OnolyticsType
GeographicalArea <- cc.all$GeographicalArea
ethnic.smr <- as.character(cc.all$ETHNIC_SMR_LAST)

source("ethnic_assign.R")

cc.all$ethnic5.smr <- ethnic5.smr

## recode SMR ethnicity to 4 categories: White, Black, South Asian, Other
cc.all$ethnic4.smr <- as.factor(car::recode(cc.all$ethnic5.smr, "'Chinese'='Other'"))
cc.all <- within(cc.all, ethnic4.smr <- factor(ethnic4.smr,
                                               levels=levels(ethnic4.smr)[c(4, 3, 1, 2)]))   
if(length(OnolyticsType) > 0) {
    cc.all$ethnic5 <- ethnic5
   
    ## tabulate ONOMAP ethnicity against SMR ethnicity
    table.ethnic <- table(cc.all$ethnic5, cc.all$ethnic5.smr, exclude=NULL)
    
    tn <- table(cc.all$ethnic5, cc.all$ethnic5.smr)
    SouthAsian.sensitivity <- 100 * tn[5, 2] / sum(tn[, 2])
    SouthAsian.specificity <- 100 * (sum(tn[, -2]) - sum(tn[5, ]) + tn[5, 2]) / sum(tn[, -2])
    sum.xtabulate <- sum(tn)
                                                       
    table.ethnic <- paste.colpercent(table.ethnic)

    ## recode ONOMAP ethnicity to 4 categories: White, South Asian, Chinese, Other
    cc.all$ethnic4 <- car::recode(cc.all$ethnic5, "'Black'='Other'")
    cc.all <- within(cc.all, ethnic4 <- relevel(as.factor(ethnic4), ref="White"))
    cc.all$ethnic4 <- factor(cc.all$ethnic4, levels=levels(cc.all$ethnic4)[c(1, 4, 2, 3)])
   
    ## recode ONOMAP ethnicity to 3 categories: White, South Asian, Other
    cc.all$ethnic3 <- car::recode(cc.all$ethnic4, "'Chinese'='Other'")
    cc.all <- within(cc.all, ethnic3 <- relevel(as.factor(ethnic3), ref="White"))
    cc.all$ethnic3 <- factor(cc.all$ethnic3, levels=levels(cc.all$ethnic3)[c(1, 3, 2)])
}

####################################################################

########################################################################################
###  diabetes based on combining Sci-Diabetes records with ICD codes and drugs 
## add in BNF codes 6.1 for diabetes drugs and
## E10 to E14 for diabetes

## immune.any includes primary immunodeficiency and secondary immunosuppression

ids.icd.diabetes <- unique(diagnoses$ANON_ID[grep("^E1[0-4]", diagnoses$ICD10)])

## 802 other immunomodulating drugs
## Methotrexate and chloroquine appear in musculoskeletal chapter 
ids.bnf.diabetes <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 601])

ids.diabetes.extra <- unique(c(ids.icd.diabetes, ids.bnf.diabetes))

# recode diabetes type
cc.all$dm.type <- as.integer(cc.all$dm.type)
## missing recoded as zero
cc.all$dm.type[is.na(cc.all$dm.type)] <- 0
## code diagnoses detected from discharges or BNF codes as unknown type
## we could classify those not on insulin as definite Type 2
cc.all$dm.type[cc.all$dm.type==0 & cc.all$ANON_ID %in% ids.diabetes.extra] <- 3

cc.all$dm.type <- 
  as.factor(car::recode(cc.all$dm.type, 
                        "c(0, 10)='Not diabetic';
                         c(1, 101, 102)='Type 1 diabetes';
                         c(2, 202, 203)='Type 2 diabetes'; 
                         3:9='Other/unknown type';
                         11:100='Other/unknown type'"))

cc.all <- within(cc.all, dm.type <- relevel(dm.type, ref="Not diabetic"))
cc.all$dm.type <- factor(cc.all$dm.type, levels=levels(cc.all$dm.type)[c(1, 3, 4, 2)])

## define an indicator variable for any diabetes
cc.all$diabetes.any <- as.integer(cc.all$dm.type != "Not diabetic")
cc.all$diabetes.any <- as.factor(car::recode(cc.all$diabetes.any,
                                        "0='Not diabetic'; 1='Diabetic'"))
cc.all <- within(cc.all, diabetes.any <- relevel(diabetes.any, ref="Not diabetic"))

####### drug exposure specifically relevant to proton pump ##################

ids.protonpump <- unique(scrips$ANON_ID[scrips$bnf_paragraph_code == "0103050"])
cc.all$protonpump <- as.factor(as.integer(cc.all$ANON_ID %in% ids.protonpump))
cc.all$y.protonpump <- as.integer(cc.all$protonpump =="1")

TW <- 120 # time window in days

if(!old) {
    ## get SPECDATE into scrips.protonpump
    scrips.protonpump <- merge(scrips.protonpump, cc.all[, c("ANON_ID", "SPECDATE")],
                               by="ANON_ID", all.x=TRUE) 
    scrips.protonpump$daysbefore <- as.integer(scrips.protonpump$SPECDATE -
                                               scrips.protonpump$dispensed_date)
    ## minimum value is 16 days -- the 15-day cutoff has been applied to the scrips table

    scrips.protonpump$dispensing.days <- as.integer(scrips.protonpump$SPECDATE - 15 - 
                                                    as.Date("2019-06-01"))
    scrips.protonpump$intervalTWday <- ceiling((scrips.protonpump$daysbefore - 15) / TW)
    scrips.protonpump$intervalTWday[scrips.protonpump$intervalTWday > 3] <- 3

## https://www.whocc.no/atc_ddd_index/?code=A02BC&showdescription=yes gives defined daily doses
##
##ATC code  	Name  	DDD 	 U 	Adm.R	 Note
##A02BC05 	esomeprazole 	30 	mg 	O
##A02BC03 	lansoprazole 	30 	mg 	O 
##A02BC01 	omeprazole 	20 	mg 	O 
##A02BC02 	pantoprazole 	40 	mg 	O 	
##A02BC04 	rabeprazole 	20 	mg 	O 	
    DDD5 <- c(30, 30, 20, 40, 20)

    scrips.protonpump$dose <- scrips.protonpump$item_strength * scrips.protonpump$quantity 

    ## calculate dose of each drug over entire period
    dose.protonpump <-  reshape2::dcast(data=scrips.protonpump,
                                        formula=ANON_ID + dispensing.days ~ approved_name,
                                        value.var="dose",
                                        fun.aggregate=sum)

    for(j in 1:5) { # loop over approved names to divide by DDD x dispensing.days 

        dose.protonpump[, j + 2] <- dose.protonpump[, j + 2] /
            (DDD5[j] * dose.protonpump$dispensing.days) 
    }
    dose.protonpump$DDDs.all <- rowSums(dose.protonpump[, -(1:2)])

    ###############  average dose in each TW-day interval ######################## 
    doseTWday.protonpump <-  reshape2::dcast(data=scrips.protonpump,
                                    formula=ANON_ID + intervalTWday ~ approved_name,
                                    value.var="dose",
                                    fun.aggregate=sum)
    ## strictly, each interval should be divided by the number of days in that interval for that individual
    ## for now, use TW days
    for(j in 1:5) { # loop over approved names to divide by DDD x TW 
        doseTWday.protonpump[, j + 2] <- doseTWday.protonpump[, j + 2] /
            (DDD5[j] * TW)
    }
    doseTWday.protonpump$DDDsTWday.all <- rowSums(doseTWday.protonpump[, -(1:2)])
    print(summary(doseTWday.protonpump))

    ## drop columns for individual drugs, and cast again to get one column for each interval
    doseTWday.protonpump <- doseTWday.protonpump[, c("ANON_ID", "intervalTWday", "DDDsTWday.all")]
    doseTWday.protonpump <- doseTWday.protonpump[!is.na(doseTWday.protonpump$intervalTWday), ]

    doseTWday.protonpump.wide <-  reshape2::dcast(data=doseTWday.protonpump,
                                    formula=ANON_ID ~ intervalTWday,
                                    value.var="DDDsTWday.all",
                                    fun.aggregate=sum)
    colnames(doseTWday.protonpump.wide)[-1] <- c("DDD.interval1", "DDD.interval2", "DDD.interval3") 

    cc.all <- merge(cc.all, dose.protonpump, by="ANON_ID", all.x=TRUE)
    cc.all$DDDs.all[is.na(cc.all$DDDs.all)] <- 0
    cc.all$ESOMEPRAZOLE[is.na(cc.all$ESOMEPRAZOLE)] <- 0
    cc.all$LANSOPRAZOLE[is.na(cc.all$LANSOPRAZOLE)] <- 0
    cc.all$OMEPRAZOLE[is.na(cc.all$OMEPRAZOLE)] <- 0
    cc.all$PANTOPRAZOLE[is.na(cc.all$PANTOPRAZOLE)] <- 0
    cc.all$RABEPRAZOLE[is.na(cc.all$RABEPRAZOLE)] <- 0
    cc.all$dispensing.days <- as.integer(cc.all$SPECDATE - as.Date("2019-06-01"))
                                        #cc.all$DDDs.average <- cc.all$DDDs.all / cc.all$dispensing.days
    cc.all$DDDsgr <- 0.5 * ceiling(2 * cc.all$DDDs.all)
    cc.all$DDDsgr <- as.factor(car::recode(cc.all$DDDsgr, "2:hi='2 or more'"))
    
    cc.all <- merge(cc.all, doseTWday.protonpump.wide, by="ANON_ID", all.x=TRUE)
    cc.all$DDD.interval1[is.na(cc.all$DDD.interval1)] <- 0
    cc.all$DDD.interval2[is.na(cc.all$DDD.interval2)] <- 0
    cc.all$DDD.interval3[is.na(cc.all$DDD.interval3)] <- 0
}

ids.nonopioid.analgesic <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 40701])
cc.all$nonopioid.analgesic <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                   ids.nonopioid.analgesic))
ids.antiplatelet <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 209])
cc.all$antiplatelet <- as.factor(as.integer(cc.all$ANON_ID %in% ids.antiplatelet))

ids.nsaid <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 100101])
cc.all$nsaid <- as.factor(as.integer(cc.all$ANON_ID %in% ids.nsaid))

ids.opioid.analgesic <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 40702])
cc.all$opioid.analgesic <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                      ids.opioid.analgesic))

ids.antipsychotic <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 40201])
cc.all$antipsychotic <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                      ids.antipsychotic))

ids.osmotic.laxative <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 10604])
cc.all$osmotic.laxative <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                   ids.osmotic.laxative))

###############################################################################
########################################################################
testpositives.ethnic.smr <- paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                                  table(ethnic5.smr, casegroup)), 1)

table.testpositives.demog.ethnicsmr <-
    tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home", "SIMD.quintile"),
                               data=cc.all[!is.na(cc.all$ethnic5.smr), ])


table.hospitalized.demog.ethnicsmr <-
    tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home", "SIMD.quintile"),
                               data=cc.all[(cc.all$casegroup=="Hospitalised, not severe" |
                                            cc.all$casegroup=="Critical care or fatal") &
                                           !is.na(cc.all$ethnic5.smr), ])

if(length(OnolyticsType) > 0) {
## tabulate ethnicity by case group
testpositives.ethnic <- paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                              table(ethnic4, casegroup)), 1)

testpositives.carehome <- paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                                table(ethnic4, care.home)), 0)

testpositives.healthboard <- t(paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                                table(ethnic4, HBRES_NAME)), 0))

table.testpositives.demog <-
    tabulate.freqs.regressions(varnames=c("ethnic4", "care.home", "SIMD.quintile"),
                               data=cc.all)

table.hospitalized.demog <-
    tabulate.freqs.regressions(varnames=c("ethnic4", "care.home", "SIMD.quintile"),
                               data=cc.all[cc.all$casegroup=="Hospitalised, not severe" |
                                           cc.all$casegroup=="Critical care or fatal", ])
}

#######################################################################################

hosp <- cc.all$casegroup=="Hospitalised, not severe"
## merge BNF chapters, one variable per subpara
chnums = 1:13
cc.hosp <- merge.bnfsubparas(chnums=chnums, data=cc.all[hosp, ])

########## restrict to severe cases and matched controls ###################### 
cat("Restricting to severe cases and matched controls\n")
cc.severe <- cc.all[cc.all$casegroup=="Critical care or fatal", ]
rm(cc.all)

## merge drugs, one variable per chapter
scrips.wide <- reshape2::dcast(scrips, ANON_ID ~ chapternum, fun.aggregate=length, 
                               value.var="chapternum")
shortnames.cols <-  bnfchapters$shortname[match(as.integer(colnames(scrips.wide)[-1]),
                                                as.integer(bnfchapters$chapternum))]
colnames(scrips.wide)[-1] <- paste("BNF", colnames(scrips.wide)[-1], shortnames.cols,
                                   sep="_")
cc.severe <- merge(cc.severe, scrips.wide, by="ANON_ID", all.x=TRUE)

bnfcols <- grep("^BNF", colnames(cc.severe))
for(j in bnfcols) {
    cc.severe[, j][is.na(cc.severe[, j])] <- 0
    cc.severe[, j][cc.severe[, j] > 1] <- 1
    cc.severe[, j] <- as.factor(cc.severe[, j])
}

## merge BNF chapters, one variable per subpara
chnums = 1:13
cc.severe <- merge.bnfsubparas(chnums=chnums, data=cc.severe)

## merge ICD diagnoses

## chapter is assigned as the position of the first element in icdchapters$start that x is greater than or equal to
unique.diagnoses <- as.character(unique(diagnoses$ICD10))
chapter <- integer(length(unique.diagnoses))
subchapter <- integer(length(unique.diagnoses))

for(i in 1:length(unique.diagnoses)) {
    chapter[i] <- min(which(substr(unique.diagnoses[i], 1, 3) <= icdchapters$end))
    ## subchapter is row in icdsubchapters table 
    subchapter[i] <- min(which(substr(unique.diagnoses[i], 1, 3) <= icdsubchapters$end))
}
unique.diagnoses <- data.frame(ICD10=unique.diagnoses, chapter=chapter, subchapter=subchapter)
diagnoses <- merge(diagnoses, unique.diagnoses, by="ICD10", all.x=TRUE)

diagnoses.wide <- reshape2::dcast(diagnoses, ANON_ID ~ chapter, fun.aggregate=length,
                                  value.var="chapter")
colnames(diagnoses.wide)[-1] <-
    paste0("Ch.", as.integer(colnames(diagnoses.wide)[-1]), "_", 
           icdchapters$shortname[as.integer(colnames(diagnoses.wide)[-1])])
## drop rare chapters
diagnoses.wide <- diagnoses.wide[, colSums(diagnoses.wide) > 20]

cc.severe <- merge(cc.severe, diagnoses.wide, by="ANON_ID", all.x=TRUE)
icdcols <- grep("^Ch.", colnames(cc.severe))
for(j in icdcols) {
    cc.severe[, j][is.na(cc.severe[, j])] <- 0
    cc.severe[, j][cc.severe[, j] > 1] <- 1
    cc.severe[, j] <- as.factor(cc.severe[, j])
}

cc.severe$num.icdchapters <- rowSums(matrix(as.integer(as.matrix(cc.severe[, icdcols])),
                                            nrow=nrow(cc.severe)))
cc.severe$num.icdchapters <-
    as.factor(car::recode(cc.severe$num.icdchapters,
                          "0='No discharge records'; 1:2='1-2 ICD-10 chapters'; 3:4='3-4 chapters'; 5:hi='5 or more chapters'")
                         )
cc.severe <- within(cc.severe, num.icdchapters <- relevel(num.icdchapters, ref="No discharge records"))

###########################################

source("comorbidity.R")

## 8 listed conditions designated by NHS
listed.conditions <- c("dm.type", "IHD.any", "heart.other.any", "oad.any",
                       "ckd.any", "neuro.any", "liver.any", "immune.any")

############ extract predefined disease categories #################
## as these are coded as factors, lowest level will be 1

cc.severe$listed.any <-
    as.factor(as.integer(with(cc.severe,
                              diabetes.any==1 | IHD.any==1 |
                              heart.other.any==1 |
                              ckd.any==1 | oad.any==1 |
                              neuro.any==1 | liver.any==1 | immune.any==1)))

########### variable lists for tabulating

demog <- c("ethnic3", "SIMD.quintile", "care.home")

if(length(OnolyticsType)==0) {
      demog <- c("SIMD.quintile", "care.home")
}

demog.smr <- c("ethnic4.smr", "SIMD.quintile", "care.home")

bnf.chapternames <- colnames(cc.severe)[bnfcols]
drugs <- bnf.chapternames

## chapters 1 and 2 subpara codes begin with 1 or 2 and have 6 digits
subparanames <- colnames(cc.severe)[grep("^subpara\\.[12][0-9]{5}\\.", colnames(cc.severe))]
icd.chapternames <- colnames(cc.severe)[icdcols]
conditions <- icd.chapternames

### incidence and mortality using national population estimates #####

source("incidencemortality.R")

###############################################################

if(FALSE) { # coding antihypertensives
antihypertensive.classes <- c("vasodilator_ah6.bnf", 
                              "centrally_acting_ah6.bnf",  
                              "adrenergic_neurone_block6.bnf",  
                              "alpha_adrenoceptor6.bnf", 
                              "ace6.bnf",  
                              "angio6.bnf",  
                              "renin_angiotensin6.bnf",  
                              "thiazides6.bnf",  
                              "calcium_channel6.bnf")
antihypertensives <- c(antihypertensive.classes[4:9], "antihypertensive.other")
}

########################################################

table.severe.demog <-
    tabulate.freqs.regressions(varnames=demog,
                               data=cc.severe)

cc.severe$diag.other <- as.integer(cc.severe$listed.any==0 & cc.severe$diag.any==1)
cc.severe$diag.other <- as.factor(car::recode(cc.severe$diag.other,
                                  "0='Listed condition or no diagnosis';
                                   1='No listed condition, other diagnosis'"))

table.scripordiag <- NULL
for(agegr in levels(cc.severe$agegr20)) {
    x <- with(cc.severe[cc.severe$agegr20==agegr, ],
              table(scripordiag, CASE))
    colnames(x) <- c("Controls", "Cases")
    x <- paste.colpercent(x)
    y <- with(cc.severe[cc.severe$agegr20==agegr, ],
              table(listed.any, CASE))
    y <- paste.colpercent(y)
    y <- y[2, , drop=FALSE]
    rownames(y) <- "Any listed condition"
    z <- with(cc.severe[cc.severe$agegr20==agegr, ],
              table(listed.any==0 & diag.any==1, CASE))
    z <- paste.colpercent(z)
    z <- z[2, , drop=FALSE]
    rownames(z) <- "No listed condition, other diagnosis"
    
    x <- rbind(x[1, , drop=FALSE], y, z, x[2, , drop=FALSE])
    table.scripordiag <- cbind(table.scripordiag, x)
}


table.scripordiag.fatal <- NULL
for(agegr in levels(cc.severe$agegr3)) {
    x <- with(cc.severe[cc.severe$agegr3==agegr, ],
              table(scripordiag, fatalcase))
    colnames(x) <- c("Controls", "Fatal cases")
    x <- paste.colpercent(x)
    x <- x[1, , drop=FALSE]
    rownames(x) <- agegr
    table.scripordiag.fatal <- rbind(table.scripordiag.fatal, x)
}

varnames.listed <- c("care.home", "scrip.any", "diag.any", "listed.any", "diag.other",
                     "scripordiag", listed.conditions)

table.agegr <- NULL
for(agegr in levels(cc.severe$agegr20)) {
    x <- univariate.tabulate(varnames=c("deathwithin28", varnames.listed), 
                             outcome="CASE",
                             data=cc.severe[cc.severe$agegr20==agegr, ],
                             drop.reflevel=FALSE)
    table.agegr <- cbind(table.agegr, x)
}
#table.agegr <- rbind(table.agegr, table.scripordiag)

freqs.all <- univariate.tabulate(varnames=c("deathwithin28", varnames.listed, "listed.any"),
                             outcome="CASE",
                             data=cc.severe,
                             drop.reflevel=FALSE)

keep.varnames <- logical(length(varnames.listed))
for(i in 1:length(varnames.listed)) {
    x <- cc.severe[, match(varnames.listed[i], colnames(cc.severe))]
    exposed <- as.integer(x) > 1
    a <- with(cc.severe[exposed, ], table(agegr3, CASE))
    keep.varnames[i] <- !any(as.integer(a)==0)
}

tables.agegr <- vector("list", length(levels(cc.severe$agegr3)))
for(i in 1:length(levels(cc.severe$agegr3))) {
    agegr <- levels(cc.severe$agegr3)[i]
    tables.agegr[[i]] <-
        tabulate.freqs.regressions(varnames=c("care.home", "scrip.any",
                                              "diag.any", listed.conditions),
                                   data=cc.severe[cc.severe$agegr3==agegr, ])
}

table.agegr.all <- tabulate.freqs.regressions(varnames=c("care.home", "scrip.any",
                                              "diag.any", listed.conditions),
                                               data=cc.severe)

## demographic vars
table.demog.aug <- tabulate.freqs.regressions(varnames=demog, data=cc.severe)

## separate analysis using SMR ethnicity 
table.ethnicsmr <- univariate.tabulate(varnames="ethnic4.smr", outcome="CASE",
                                       data=cc.severe[!is.na(cc.severe$ethnic4.smr), ],
                                       drop.reflevel=FALSE)
univariate.ethnicsmr <-
    univariate.clogit(varnames="ethnic4.smr",
                      data=cc.severe[!is.na(cc.severe$ethnic4.smr), ],
                      add.reflevel=TRUE)
table.ethnicsmr.aug <- combine.tables2(table.ethnicsmr, univariate.ethnicsmr)
rownames(table.ethnicsmr.aug) <- replace.names(rownames(table.ethnicsmr.aug))

#### 5 ethnic groups for HPS report
table.ethnic5smr <-
    tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home",
                                          "SIMD.quintile"),
                               outcome="CASE",
                               data=cc.severe[!is.na(cc.severe$ethnic5.smr), ])
rownames(table.ethnic5smr) <- replace.names(rownames(table.ethnic5smr))

## listed conditions
table.listed.conditions.lt60 <-
    tabulate.freqs.regressions(varnames=listed.conditions,
                               data=cc.severe[cc.severe$AGE < 60, ])
table.listed.conditions.ge60 <-
    tabulate.freqs.regressions(varnames=listed.conditions,
                               data=cc.severe[cc.severe$AGE >= 60, ])

## full multivariate model variables -- for this use dm.type rather than diabetes.any 
multivariate.all <-
    multivariate.clogit(varnames=c(demog, "dm.type", listed.conditions[-1],
                                   "diag.any", conditions, "scrip.any", drugs,
                                   "protonpump"),
                        data=cc.severe, add.reflevel=TRUE)

################# restrict to those without listed conditions #############

nocare <- cc.severe$care.home=="Independent"
notlisted <- cc.severe$listed.any == 0
nocare.notlisted <- nocare & notlisted

## conditions
table.conditions.aug <- tabulate.freqs.regressions(varnames=conditions, 
                                                   data=cc.severe[notlisted, ])
cat("Tabulating ICD subchapter diagnoses ...")
## tabulate subchapters in ICD chapters of interest
table.icdchapter2 <- tabulate.icdchapter(chnum=2, data=cc.severe[notlisted, ])
table.icdchapter7 <-  tabulate.icdchapter(chnum=7, data=cc.severe[notlisted, ])
table.icdchapter11 <- tabulate.icdchapter(chnum=11, data=cc.severe[notlisted, ])

table.icdsubchapters <- NULL
for(i in 1:20) {
    table.icdsubchapters <-
        rbind(table.icdsubchapters,
              tabulate.icdchapter(chnum=i, data=cc.severe[notlisted, ], minrowsum=50))
}
table.icdsubchapters <- table.icdsubchapters[grep("ensuremath",
                                                  table.icdsubchapters$u.pvalue), ]
cat("done\n")

#########################################################################

## drugs 
table.drugs.aug <- tabulate.freqs.regressions(varnames=drugs, 
                                              data=cc.severe[notlisted, ])

#############################################################################

## tabulate scrip.any effect by carehome and listed.any

## code care home residents as 1, other as 0
cc.severe$cats3 <-  as.integer(cc.severe$care.home == "Care/nursing home")
cc.severe$cats3[cc.severe$cats3==0 & cc.severe$listed.any==1] <- 2 
cc.severe$cats3[cc.severe$cats3==0 & cc.severe$listed.any==0] <- 3
cc.severe$cats3 <- as.factor(car::recode(cc.severe$cats3, "1='Care/nursing home';
                               2='Independent living, listed condition';
                               3='Independent living, no listed condition'"))     

table.scrip.cats3 <- NULL
for(cat in levels(cc.severe$cats3)) {
    x <- tabulate.freqs.regressions(varnames="scrip.any",
                                    data=cc.severe[cc.severe$cats3==cat, ])[, 1:4]
    rownames(x) <- cat
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.scrip.cats3 <- rbind(table.scrip.cats3, x)
}

################################################################

## backwards selection of smallest subset of BNF chapters that explains most of the scrip.any effect
    
x <- cc.severe[nocare.notlisted, ][, drugs]
x <- matrix(as.integer(as.matrix(x)), nrow=nrow(x))
colnames(x) <- drugs
y <- cc.severe[nocare.notlisted, ]$CASE
stratum <- cc.severe[nocare.notlisted, ]$stratum

stepwise.drop <- stepwise.union.dropcols(x=x, y=y, stratum=stratum)

x <- cc.severe[nocare.notlisted, ][, subparanames]
x <- matrix(as.integer(as.matrix(x)), nrow=nrow(x))
colnames(x) <- subparanames

cat("Stepwise drop procedure over subparas in BNF chapters 1 and 2 ...")
stepwise.drop.subparas <- stepwise.union.dropcols(x=x, y=y, stratum=stratum)
cat("done\n")

## tabulate associations with drug chapters in those not in care homes and without listed conditions 
table.drugs.nocare.notlisted <- tabulate.freqs.regressions(varnames=drugs, 
                                                           data=cc.severe[nocare.notlisted, ])

## tabulate proportion of effect of scrip.any that is explained by each chapter
## use cc.severe[nocare.notlisted, ]
table.anyscrip.chapter <-
    summary(clogit(formula=as.formula(paste("CASE ~ scrip.any + strata(stratum)")),
                                      data=cc.severe[nocare.notlisted, ]))$coefficients[1, 1:2, drop=FALSE]
for(bnfchapter in drugs) {
    ch.formula=as.formula(paste("CASE ~ scrip.any +", bnfchapter, "+ strata(stratum)"))
    table.anyscrip.chapter <-
        rbind(table.anyscrip.chapter,
              summary(clogit(ch.formula, data=cc.severe[nocare.notlisted, ]))$coefficients[1, 1:2])
}
rownames(table.anyscrip.chapter) <- c("Unadjusted", drugs)
table.anyscrip.chapter <- as.data.frame(table.anyscrip.chapter)
table.anyscrip.chapter$prop.explained <- round(with(table.anyscrip.chapter,
                                                    c(0, 1 - coef[-1] / coef[1])), 2)

## tabulate para or subpara codes in BNF chapters of interest
table.bnfchapter1 <- tabulate.bnfparas(chnum=1, data=cc.severe[nocare.notlisted, ])
table.bnfchapter2 <- tabulate.bnfsubparas(chnum=2, data=cc.severe[nocare.notlisted, ])
table.bnfchapter4 <- tabulate.bnfsubparas(chnum=4, data=cc.severe[nocare.notlisted, ])
table.bnfchapter9 <- tabulate.bnfsubparas(chnum=9, data=cc.severe[nocare.notlisted, ])
table.bnfchapter10 <- tabulate.bnfsubparas(chnum=10, data=cc.severe[nocare.notlisted, ])

## fix to tabulate BNF chemical substance

if(!old) {
############# proton pump #########################

## calculate propensity score trained on cc.hosp

source("propensity.R")

cc.severe$propensity <- propensity

####### effects of scrip and protonpump by care home / listed condition status ######################

table.protonpump.cats3 <- NULL
for(condition in levels(cc.severe$cats3)) {
    x <- tabulate.freqs.regressions(varnames="protonpump",
                                    data=cc.severe[cc.severe$cats3==condition, ])[, 1:4]
    rownames(x) <- condition
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.protonpump.cats3 <- rbind(table.protonpump.cats3, x)
}

####### tabulate dose-response effect  ##########################################

## tabulate.freqs.regressions() can only be used with factor variables

table.dose.protonpump <- NULL
withcovariates.formula <- as.formula("CASE ~ DDDs.all + propensity + strata(stratum)")
coeff.row <- 1

for(agegr in levels(cc.severe$agegr3)) {
    x <- summary(clogit(formula=CASE ~ DDDs.all + strata(stratum), 
                        data=cc.severe[cc.severe$agegr3==agegr, ]))$coefficients[coeff.row, , drop=FALSE]
    x <- as.data.frame(x)
    x$u.ci <- or.ci(x[, 1], x[, 3])
    x$u.pvalue <- pvalue.latex(x[, 5])
    x <- x[, c("u.ci", "u.pvalue")]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.severe[cc.severe$agegr3==agegr, ]))$coefficients[coeff.row, , drop=FALSE]
    y <- as.data.frame(y)
    y$m.ci <- or.ci(y[, 1], y[, 3])
    y$m.pvalue <- pvalue.latex(y[, 5])
    y <- y[, c("m.ci", "m.pvalue")]
    x <- data.frame(x, y)               
    table.dose.protonpump <- rbind(table.dose.protonpump, x)
}

rownames(table.dose.protonpump) <- NULL
freqs <- paste.colpercent(with(cc.severe, table(agegr3, CASE)))
colnames(freqs)[1:2] <- c("Controls", "Cases")
freqs <- as.data.frame(freqs)
table.dose.protonpump  <- cbind(freqs, table.dose.protonpump)

print(table.dose.protonpump)

####### tabulate dose-response effect  ##########################################

## tabulate.freqs.regressions() can only be used with factor variables

table.dosegr.protonpump <- NULL
withcovariates.formula <- as.formula("CASE ~ DDDsgr + propensity + strata(stratum)")
coeff.rows <- 1:4

for(agegr in levels(cc.severe$agegr3)) {
    x <- summary(clogit(formula=CASE ~ DDDsgr + strata(stratum), 
                        data=cc.severe[cc.severe$agegr3==agegr, ]))$coefficients
    x <- as.data.frame(x)
    x$u.ci <- or.ci(x[, 1], x[, 3])
    x$u.pvalue <- pvalue.latex(x[, 5])
    x <- x[, c("u.ci", "u.pvalue")]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.severe[cc.severe$agegr3==agegr, ]))$coefficients[coeff.rows, ] # drop=FALSE]
    y <- as.data.frame(y)
    y$m.ci <- or.ci(y[, 1], y[, 3])
    y$m.pvalue <- pvalue.latex(y[, 5])
    y <- y[, c("m.ci", "m.pvalue")]
    x <- data.frame(x, y)               
   
    x <- rbind(rep(NA, ncol(x)), x)
    freqs <- paste.colpercent(with(cc.severe[cc.severe$agegr3==agegr, ], table(DDDsgr, CASE)))
    colnames(freqs)[1:2] <- c("Controls", "Cases")
    x <- data.frame(DDD.average=rownames(freqs), freqs, x)
    rownames(x) <- paste0(agegr, ": ", rownames(x))
    table.dosegr.protonpump <- rbind(table.dosegr.protonpump, x)
}
print(table.dosegr.protonpump)

################# tabulate ever-use effect

withcovariates.formula <- as.formula("CASE ~ protonpump + propensity + strata(stratum)")

coeff.row <- 1
table.everuse.protonpump <- NULL
for(agegr in levels(cc.severe$agegr20)) {
    x <- tabulate.freqs.regressions(varnames="protonpump",
                                    data=cc.severe[cc.severe$agegr20==agegr, ])[, 1:4]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.severe[cc.severe$agegr20==agegr, ]))$coefficients[coeff.row, , drop=FALSE]
    x$m.ci <- or.ci(y[, 1], y[, 3])
    x$m.pvalue <- pvalue.latex(y[, 5])
    rownames(x) <- agegr
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.everuse.protonpump <- rbind(table.everuse.protonpump, x)
}
x <- tabulate.freqs.regressions(varnames="protonpump", data=cc.severe)[, 1:4]
y <- summary(clogit(formula=withcovariates.formula,
                    data=cc.severe))$coefficients[coeff.row, , drop=FALSE]
x$m.ci <- or.ci(y[, 1], y[, 3])
x$m.pvalue <- pvalue.latex(y[, 5])
rownames(x) <- "All"
colnames(table.everuse.protonpump)[1:2] <- colnames(x)[1:2]
table.everuse.protonpump <- rbind(table.everuse.protonpump, x)

#### time window #######################

cc.severe$exposure.nonrecent <- as.integer(cc.severe$DDD.interval2 > 0 | cc.severe$DDD.interval3 > 0)
cc.severe$exposure.recent <- as.integer(cc.severe$DDD.interval1 > 0)

cc.severe$exposurecat <- integer(nrow(cc.severe))
cc.severe$exposurecat[cc.severe$exposure.nonrecent==0 & cc.severe$exposure.recent==0] <- 0
cc.severe$exposurecat[cc.severe$exposure.nonrecent==1 & cc.severe$exposure.recent==0] <- 1
cc.severe$exposurecat[cc.severe$exposure.nonrecent==0 & cc.severe$exposure.recent==1] <- 2
cc.severe$exposurecat[cc.severe$exposure.nonrecent==1 & cc.severe$exposure.recent==1] <- 3

cc.severe$exposurecat <- as.factor(car::recode(cc.severe$exposurecat,
                                "0='No prescriptions'; 1='Non-recent only'; 2='Recent only';
                                 3='Prescriptions in both time windows'"))
cc.severe$exposurecat <- factor(cc.severe$exposurecat,
                                levels=levels(cc.severe$exposurecat)[c(4, 2, 3, 1)])

summary(clogit(formula=CASE ~ exposurecat + strata(stratum), 
               data=cc.severe[cc.severe$AGE < 75, ]))$coefficients

## FIXME: this function reorders levels of the factor exposurecat
table.timewindow.protonpump <- tabulate.freqs.regressions(varnames="exposurecat",
                                                          data=cc.severe)[, 1:4]
 
###########################################

## tabulate effect of each drug
ppinames <- grep("PRAZOLE$", names(cc.severe), value=TRUE)
ppinames.formula <- as.formula(paste("CASE ~", paste(ppinames, collapse="+"), "+ strata(stratum)"))

ppi.means <- cbind(
    round(apply(cc.severe[cc.severe$CASE==0, ppinames], 2, mean), 4),
    round(apply(cc.severe[cc.severe$CASE==1, ppinames], 2, mean), 4)
)

table.ppinames <- NULL
for(ppi in ppinames) {
    x  <- summary(clogit(formula=as.formula(paste("CASE ~", ppi, "+ strata(stratum)")),
                         data=cc.severe))$coefficients
    table.ppinames <- rbind(table.ppinames, x)
}
y <- summary(clogit(formula=ppinames.formula, 
                    data=cc.severe))$coefficients
x <- as.data.frame(table.ppinames)
x$u.ci <- or.ci(x[, 1], x[, 3])
x$u.pvalue <- pvalue.latex(x[, 5])
x$m.ci <- or.ci(y[, 1], y[, 3])
x$m.pvalue <- pvalue.latex(y[, 5])
x <- x[, c("u.ci", "u.pvalue", "m.ci", "m.pvalue")]
table.ppinames <- data.frame(ppi.means, x)

################# tabulate effect of any proton pump exposure by age group #############################

withcovariates.formula <- as.formula("CASE ~ protonpump + propensity + strata(stratum)")
coeff.row <- 1

table.protonpump <- NULL
for(agegr in levels(cc.severe$agegr20)) {
    x <- tabulate.freqs.regressions(varnames="protonpump",
                                    data=cc.severe[cc.severe$agegr20==agegr, ])[, 1:4]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.severe[cc.severe$agegr20==agegr, ]))$coefficients[coeff.row, , drop=FALSE]
    x$m.ci <- or.ci(y[, 1], y[, 3])
    x$m.pvalue <- pvalue.latex(y[, 5])
    rownames(x) <- agegr
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.protonpump <- rbind(table.protonpump, x)
}
x <- tabulate.freqs.regressions(varnames="protonpump", data=cc.severe)[, 1:4]
y <- summary(clogit(formula=withcovariates.formula,
                    data=cc.severe))$coefficients[coeff.row, , drop=FALSE]
x$m.ci <- or.ci(y[, 1], y[, 3])
x$m.pvalue <- pvalue.latex(y[, 5])
rownames(x) <- "All"
colnames(table.protonpump)[1:2] <- colnames(x)[1:2]
table.protonpump <- rbind(table.protonpump, x)
################################################################

## tabulate fatal cases  by age group
table.fatal.protonpump <- NULL
for(agegr in levels(cc.severe$agegr20)) {
    x <- tabulate.freqs.regressions(varnames="protonpump", outcome="fatalcase",
                                    data=cc.severe[cc.severe$agegr20==agegr, ])[, 1:4]
    rownames(x) <- agegr
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.fatal.protonpump <- rbind(table.fatal.protonpump, x)
}
x <- tabulate.freqs.regressions(varnames="protonpump", data=cc.severe)[, 1:4]
rownames(x) <- "All"
colnames(x)[1:2] <- c("Controls", "Cases")
table.fatal.protonpump <- rbind(table.fatal.protonpump, x)
}
######## stepwise regressions use saved version #####################
nfold <- 4
#stepwise <- TRUE
stepwise <- FALSE

source("stepwise.R")

#####################################################################

if(old) {
    rmarkdown::render("casecontrol.Rmd", output_file="casecontrol.pdf")
} else  {
    rmarkdown::render("pharmaco.Rmd", output_file="pharmaco.pdf")
} 

#rmarkdown::render("Covid_ethnicity_Scotland.Rmd")

## remove large objects from memory
objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(objmem)

