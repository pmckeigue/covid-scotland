## analysis script for case-control study
if(exists("cl")) {
    parallel::stopCluster(cl)
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
                     "reshape2",  # deprecated
                     "readxl", 
                     "DescTools", 
                     "icd.data", 
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
library(reshape2)
library(readxl)
library(DescTools)
library(icd.data)
library(gam)
library(data.table)

Rprof(tmp <- tempfile())

source("helperfunctions.R")

old <- TRUE
old <- FALSE # uncomment to use 8 June linkage

include.nrs.only <- TRUE
#include.nrs.only <- FALSE # uncomment to exclude NRS-only deaths

stepwise <- TRUE   
stepwise <- FALSE ## uncomment to save time if old version still valid

if(old) {
    cc.all <- as.data.table(readRDS("./data/CC_linked_ANON_20200515.rds"))
    diagnoses <- readRDS("./data/CC_SMR01_ICD10_x25_ANON_20200515.rds")
    procedures <- readRDS("./data/CC_SMR01_OPCS4_MAIN.x25_ANON_20200515.rds")
    scrips.filename <- "./data/CC_PIS_x15_ANON_20200515.rds" 
    controls.status <- readRDS("./data/CC_CHI_CHECK_2020-05-15.rds")
} else {
    cc.all <- as.data.table(readRDS("./data/CC_linked_ANON_2020-06-18.rds"), key="ANON_ID")
    diagnoses <- readRDS("./data/CC_SMR01_ICD10_x25_ANON_2020-06-18.rds")
    procedures <- readRDS("./data/CC_SMR01_OPCS4_MAIN.x25_ANON_2020-06-18.rds")
    scrips.filename <- "./data/CC_PIS_x15_ANON_2020-06-18.rds"
    onomap <- as.data.table(readRDS("./data/ONOMAP_ANON_2020-06-18.rds"), key="ANON_ID")
    rapid <- readRDS("./data/CC_RAPID_ANON_2020-06-18.rds")
    cases.status <-
        as.data.table(readRDS("./data/CC_EXTENDED_STATUS_ANON_2020-06-18.rds"), key="ANON_ID")
}

scrips <- as.data.table(readRDS(scrips.filename)[, c("ANON_ID", "bnf_paragraph_code",
                                                     "bnf_paragraph_description")])

################### short names for ICD chapters ########################

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
icdchapters$shortname <- gsub("^health$", "Health_factors", icdchapters$shortname)
icdchapters$start <- as.character(icdchapters$start)
icdchapters$end <- as.character(icdchapters$end)
icdsubchapters <- data.frame(names(icd10_sub_chapters),
                             t(matrix(as.character(unlist(icd10_sub_chapters)), nrow=2)))
colnames(icdsubchapters) <- c("name", "start", "end")
icdsubchapters$start <- as.character(icdsubchapters$start)
icdsubchapters$end <- as.character(icdsubchapters$end)

###############################################

names(cc.all) <- gsub("CASE_NO", "stratum", names(cc.all))
names(cc.all) <- gsub("^SEX$", "sex", names(cc.all))
names(cc.all) <- gsub("imumune", "immune", names(cc.all))
#names(cc.all) <- gsub("^ethnic$", "ethnic.old", names(cc.all))
names(cc.all) <- gsub("^CAREHOME$", "care.home", names(cc.all))
names(cc.all) <- gsub("^simd$", "SIMD.quintile", names(cc.all))
names(cc.all) <- gsub("DATE_OF_DEATH", "Date.Death", names(cc.all))
names(cc.all) <- gsub("SPECIMENDATE", "SPECDATE", names(cc.all))
names(cc.all) <- gsub("^age$", "AGE", names(cc.all))
names(cc.all) <- gsub("^AgeYear$", "AGE", names(cc.all))
if(!old) {
    cc.all$CASE <- as.integer(cc.all$is.case)
}
cc.all$stratum <- as.integer(cc.all$stratum)

cc.all <- mutate(cc.all, dispensing.days = as.integer(SPECDATE - as.Date("2019-06-01")))

## HAI is based on the ECDC definition of nosocomial infection

paste.colpercent(with(cc.all, table(ANON_ID %in% scrips$ANON_ID, CASE)))
paste.colpercent(with(cc.all, table(ANON_ID %in% diagnoses$ANON_ID, CASE)))

cc.all$scrip.any <- as.factor(as.integer(cc.all$ANON_ID %in% scrips$ANON_ID))
cc.all$diag.any <- as.factor(as.integer(cc.all$ANON_ID %in% diagnoses$ANON_ID))

ids.noscrip.agegt75 <- subset(cc.all, CASE==0 & scrip.any==0 & AGE > 75)[["ANON_ID"]][1:20]
ids.nodiag.agegt75 <- subset(cc.all, CASE==0 & diag.any==0 & AGE > 75)[["ANON_ID"]][1:20]

save(ids.noscrip.agegt75, ids.nodiag.agegt75, file="anon_ids_notmatched.RData")

## exclude controls already dead on date of test of case they were matched to
controls.deceased <- with(cc.all, CASE==0 &
                                  !is.na(Date.Death) &
                                  Date.Death <= SPECDATE)
cc.all <- subset(cc.all, !controls.deceased)                         

if(old) {
    ## tabulate controls not observable
    controls.status <- controls.status[, c("ANON_ID", "CHI.EXTENDED.STATUS",
                                           "CHI.DATE.OF.DEATH",
                                           "CHI.Explanation")]
    cc.all <- merge(cc.all, controls.status, by="ANON_ID", all.x=TRUE)
    table.exclusions.scrip <- with(cc.all, table(CHI.Explanation, scrip.any))
    table.exclusions.diag <- with(cc.all, table(CHI.Explanation, diag.any))
    ## exclude controls classified on CHI database as no longer current
    cc.all <- subset(cc.all, subset=is.na(CHI.Explanation) |
                                 CHI.Explanation=="Current and no history" |
                                 CHI.Explanation=="Current with history")
} else { # exclude cases classified as unobservable, for consistency with controls 
    cc.all <- cases.status[cc.all]
    with(subset(cc.all, CASE==1),
         table(Explanation, is.na(Date.Death)))
    with(subset(cc.all, CASE==1),
         table(EXTENDED_STATUS, is.na(Date.Death)))
    cc.all <- subset(cc.all,
                     is.na(EXTENDED_STATUS) |
                     EXTENDED_STATUS == "C" |
                     (EXTENDED_STATUS == "D" & !is.na(Date.Death)))
}

## rows have been replicated within strata containing N cases so that each case and control appears N^2 times
## remove these replicated rows
cc.all <- distinct(.data=cc.all, ANON_ID, stratum, .keep_all = TRUE)
setkey(cc.all, ANON_ID)

numcases.strata <- tapply(cc.all$CASE, cc.all$stratum, sum) 
numctrls.strata <- tapply(1 - cc.all$CASE, cc.all$stratum, sum)
print(table(numctrls.strata, numcases.strata))

cc.all$SIMD.quintile <- car::recode(cc.all$SIMD.quintile, "'Unknown'=NA")

cc.all$scripordiag <- as.integer(with(cc.all, as.integer(diag.any)==2 |
                                              as.integer(scrip.any)==2))

cc.all$sex <- car::recode(as.factor(cc.all$sex), "1='Male'; 2='Female'")
cc.all[, sex := relevel(sex, ref="Female")]

cc.all$agegr20 <- as.factor(car::recode(as.integer(cc.all$AGE),
                              "0:39='0-39'; 40:59='40-59';
                                  60:74='60-74'; 75:hi='75 or more'"))
cc.all$agegr3 <-
    as.factor(car::recode(cc.all$AGE,
                          "0:59='0-60 years'; 60:74='60-74 years'; 75:hi='75+ years'"))

cc.all$agegr2 <-
    as.factor(car::recode(cc.all$AGE,
                          "0:74='0-74 years'; 75:hi='75+ years'"))

cc.all$care.home[cc.all$NURSINGHOME==1] <- 1
cc.all$care.home <- as.factor(car::recode(cc.all$care.home, "0='Independent'; 1='Care/nursing home'"))
cc.all[, care.home := relevel(care.home, ref="Independent")]

if(old) {
    cc.all <- mutate(cc.all, testpositive.case = CASE==1 & nrs_covid_case==0)
} else {
    cc.all$ECOSS_POSITIVE[is.na(cc.all$ECOSS_POSITIVE)] <- "No test"
    cc.all$ECOSS_POSITIVE[cc.all$ECOSS_POSITIVE==""] <- "No result"
    ## 91 individuals are included as cases but they do not have covid on death cert and their test result is negative or empty string (labelled "No result")
    ## probably have a later test result that was positive
    ## assigned for now as test-positive cases
    cc.all$testpositive.case <- cc.all$CASE==1 &
        (cc.all$ECOSS_POSITIVE=="Positive" |
         cc.all$ECOSS_POSITIVE=="Negative" | 
         cc.all$ECOSS_POSITIVE=="No result")
    with(subset(cc.all, CASE==1), table(ECOSS_POSITIVE, covid_cod))
    with(cc.all[cc.all$CASE==1, ], table(ECOSS_POSITIVE, adm28, exclude=NULL))
}


## all cases have nonmissing SPECDATE
## controls should be assigned same SPECDATE as the case they were matched to
## but deathwithin28 should be 0 for those who did not test positive
cc.all <- mutate(cc.all, deathwithin28 = testpositive.case &
                             !is.na(Date.Death) & 
                             Date.Death - SPECDATE >= 0 & Date.Death - SPECDATE <= 28)

## Sharon's variable dead28 is assigned by this line
## Covid_CC_linkage_Part2_desktop.R:
## cc$dead28 <- ifelse(!is.na(cc$DATE_OF_DEATH) & cc$DATE_OF_DEATH >= cc$SPECIMENDATE & as.numeric(cc$DATE_OF_DEATH - cc$SPECIMENDATE) <=28, 1, 0)
## evaluates to 1 for cases without positive test result 

with(subset(cc.all, CASE==1), table(dead28, deathwithin28, exclude=NULL)) 

with(subset(cc.all, CASE==1), table(icu, hdu, exclude=NULL)) 
## only 854 in critical care -- must be wrong
cc.all$criticalcare <- cc.all$icu==1 | cc.all$hdu==1

## what does a value of 2 represent? 
## Covid_CC_linkage_Part2_desktop.R:
## icu$icu <- ifelse(icu$covidICUorHDU ==1, 1, 0)
## icu$hdu <- ifelse(icu$covidICUorHDU ==3, 1, 0)

## adm28 should be 0 for cases who did not test positive
cc.all <- mutate(cc.all, adm28 = adm28 & testpositive.case)

with(cc.all[cc.all$CASE==1, ], table(inhosp, adm28, exclude=NULL))
## Covid_CC_linkage_Part2_desktop.R:rapid$days <- as.numeric(rapid$Admission.Date - rapid$SPECIMENDATE)
## Covid_CC_linkage_Part2_desktop.R:rapid$adm28 <- ifelse(rapid$days %in% c(0:28), 1, 0)


## new linkage does not have the variable nrs_covid_case
with(subset(cc.all, CASE==1 & covid_ucod==1), table(deathwithin28, exclude=NULL))
## 3701 cases with covid_cod have deathwithin28==1
## all those with covid_ucod==1 have covid_cod==1 

with(subset(cc.all, CASE==1 & testpositive.case), table(criticalcare, deathwithin28, exclude=NULL))
with(subset(cc.all, CASE==1), table(criticalcare, testpositive.case, exclude=NULL))
with(subset(cc.all, CASE==1), table(testpositive.case, deathwithin28, exclude=NULL))

## coding of case groups
cc.all <- mutate(cc.all, group = "Unclassified")

## define group A (severe cases) 
if(include.nrs.only) { # group A includes anyone certified with COVID as underlying cause
    cc.all[CASE==1 & (criticalcare | deathwithin28 | covid_ucod==1),
           group := "A"]
} else { # group A includes test-negative cases ascertained through NRS only if they were in critical care
    cc.all[CASE==1 & (testpositive.case & (criticalcare | deathwithin28==1)) |
                     (criticalcare & covid_ucod==1), group := "A"]
}

## assign remaining test-positive cases hospitalized within 28 days of positive test as group B
cc.all[testpositive.case & group=="Unclassified" & (adm28 > 0 | inhosp > 0), group := "B"]
## assign all remaining test-positive cases to group C
cc.all[testpositive.case & group=="Unclassified", group := "C"]
## assign remaining cases with mention on death cert to group D
cc.all[CASE==1 & group=="Unclassified" & covid_cod==1, group := "D"]

print(table(cc.all$CASE, cc.all$group, exclude=NULL))

#print(subset(cc.all, CASE==1 & group=="Unclassified", 
#             select=c(covid_cod, deathwithin28, ECOSS_POSITIVE, criticalcare)))

#print(with(subset(cc.all, CASE==1 & group=="Unclassified"), table(ECOSS_POSITIVE, covid_cod)))

## remove unclassified cases 
#cc.all <- subset(cc.all, !(CASE==1 & group == "Unclassified"))

## create a new variable named casegroup, assign controls in each stratum to same group as case
## for strata with 2 or more cases, assign all controls to highest group of case
casegroups <- subset(cc.all, subset=CASE==1, select=c("stratum", "group"))
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
cc.all <- subset(cc.all, !is.na(casegroup))

with(cc.all[cc.all$CASE==1, ], table(casegroup, deathwithin28, exclude=NULL))

#cc.all$casegroup <- car::recode(cc.all$casegroup,
#                                "'A'='Critical care or fatal'; 'B'='Hospitalised, not severe'; 'C'='Test-positive, not hospitalised'; 'D'='Contributing cause on death cert, no test'")

print(paste.colpercent(with(cc.all[cc.all$CASE==1, ], table(care.home, casegroup))))

### incidence and mortality using national population estimates #####
narrow.case <- cc.all$CASE==1 & cc.all$casegroup=="A"
narrow.fatalcase <- narrow.case & (cc.all$deathwithin28 | cc.all$covid_ucod==1)
table(narrow.case, narrow.fatalcase)
cc.all$fatalcase <- narrow.fatalcase # for later use

## broad.case adds in all those with mention of covid on death cert
broad.case <- cc.all$CASE==1 & (cc.all$casegroup=="A" | cc.all$casegroup=="D")
broad.fatalcase <- broad.case & (cc.all$deathwithin28 | cc.all$covid_cod==1)
table(broad.case, broad.fatalcase)

## for graphing, we want three categories
## severe cases as defined
## broad cases as used for diabetes report: severe, hospitalized or NRS
## fatal cases including NRS deaths
narrow.case.freqs <- with(cc.all[narrow.case, ], table(AGE, sex, exclude=NULL))
broad.case.freqs <- with(cc.all[broad.case, ], table(AGE, sex, exclude=NULL))
narrow.death.freqs <- with(cc.all[narrow.fatalcase, ], table(AGE, sex, exclude=NULL))
broad.death.freqs <- with(cc.all[broad.fatalcase, ], table(AGE, sex, exclude=NULL))

save(narrow.case.freqs, broad.case.freqs,
     narrow.death.freqs, broad.death.freqs, file="casefreqs.4cats.agesex.RData")

source("incidencemortality.R")

######## coding ethnicity ##############################

source("ethnic_assign.R")

cc.all <- mutate(cc.all, ethnic5.smr = collapseto5.ethnicsmr(ETHNIC_SMR_LAST))

## ANON_ID may be duplicated in cc.all but not in onomap

if(!old) {
    ## some IDs in onomap are not in cc.all
    onomap <- subset(onomap, ANON_ID %in% cc.all$ANON_ID)
    setkey(onomap, ANON_ID)
    setkey(cc.all, ANON_ID)
    cc.all <- onomap[cc.all]
}

cc.all <- mutate(cc.all,
                 group.onomap = group.onomap(OnolyticsType, GeographicalArea))
cc.all <- mutate(cc.all,
                 ethnic5.onomap = collapseto5.onomap.group(group.onomap))

## tabulate ONOMAP ethnicity against SMR ethnicity
table.ethnic <- table(cc.all$ethnic5.onomap, cc.all$ethnic5.smr, exclude=NULL)
tn <- as.data.frame.matrix(table(cc.all$ethnic5.onomap, cc.all$ethnic5.smr))
SouthAsian.sensitivity <- 100 * tn["South Asian", "South Asian"] / sum(tn[, "South Asian"])
SA.trueneg <- sum(tn) - sum(tn[, "South Asian"])
SA.trueneg.correct <- SA.trueneg - (sum(tn[, "South Asian"]) - tn["South Asian", "South Asian"]) 
SouthAsian.specificity <- 100 * SA.trueneg.correct / SA.trueneg
sum.xtabulate <- sum(tn)

## recode SMR ethnicity to 4 categories: White, Black, South Asian, Other
cc.all <- mutate(cc.all, ethnic4.smr = as.factor(car::recode(ethnic5.smr,
                                                             "'Chinese'='Other'")))
cc.all[, ethnic4.smr := factor(ethnic4.smr,
                               levels=levels(ethnic4.smr)[c(4, 3, 1, 2)])]   

## recode ONOMAP ethnicity to 4 categories: White, South Asian, Chinese, Other
cc.all <- mutate(cc.all, ethnic4.onomap = car::recode(ethnic5.onomap, "'Black'='Other'"))

cc.all[, ethnic4.onomap := factor(ethnic4.onomap,
                                  levels=levels(ethnic4.onomap)[c(4, 3, 1, 2)])]


## recode to 3 categories: White, South Asian, Other
cc.all <- mutate(cc.all, ethnic3.onomap = car::recode(ethnic4.onomap, "'Chinese'='Other'"))
cc.all[, ethnic3.onomap := factor(ethnic3.onomap,
                                  levels=levels(ethnic3.onomap)[c(3, 2, 1)])]
setkey(cc.all, ANON_ID)
######################################################
                 
source("bnfcodes.R")

source("drugs.R")

###########################################################################

###  diabetes based on combining Sci-Diabetes records with ICD codes and drugs 
## add in BNF codes 6.1 for diabetes drugs and
## E10 to E14 for diabetes

## immune.any includes primary immunodeficiency and secondary immunosuppression

########## coding listed conditions ####################

ids.icd.diabetes <- unique(diagnoses$ANON_ID[grep("^E1[0-4]", diagnoses$ICD10)])
ids.bnf.diabetes <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 601])
ids.diabetes.extra <- unique(c(ids.icd.diabetes, ids.bnf.diabetes))

# recode diabetes type
cc.all$dm.type <- as.integer(cc.all$dm.type)
## missing recoded as zero
cc.all$dm.type[is.na(cc.all$dm.type)] <- 0

## add in extra cases notified directly from SCI-Diabetes register, without assignment
## of diabetes type from SDRN database
cc.all$dm.type[cc.all$dm.type==0 & cc.all$diab.reg==1] <- 3

cat("Extra diabetes cases from SCI-Diabetes by diabetes type\n")
print(table(cc.all$dm.type, cc.all$ANON_ID %in% ids.diabetes.extra))

## code diagnoses detected from discharges or BNF codes as unknown type
## we could classify those not on insulin as definite Type 2 but Helen says no

## REVISION: for consistency with the diabetes paper, we will not include the extra cases identified through diagnostic codes or drug codes - they may be transient/resolved

## cc.all$dm.type[cc.all$dm.type==0 & cc.all$ANON_ID %in% ids.diabetes.extra] <- 3

cc.all$dm.type <- 
  as.factor(car::recode(cc.all$dm.type, 
                        "c(0, 10, 17, 97)='Not diabetic';
                         c(1, 101, 102)='Type 1 diabetes';
                         c(2, 202, 203)='Type 2 diabetes'; 
                         3:9='Other/unknown type';
                         11:16='Other/unknown type';
                         18:96='Other/unknown type';
                         98:100='Other/unknown type'"))

cc.all <- within(cc.all, dm.type <- relevel(dm.type, ref="Not diabetic"))
cc.all$dm.type <- factor(cc.all$dm.type, levels=levels(cc.all$dm.type)[c(1, 3, 4, 2)])

## define an indicator variable for any diabetes
cc.all$diabetes.any <- as.integer(cc.all$dm.type != "Not diabetic")
cc.all$diabetes.any <- as.factor(car::recode(cc.all$diabetes.any,
                                        "0='Not diabetic'; 1='Diabetic'"))
cc.all <- within(cc.all, diabetes.any <- relevel(diabetes.any, ref="Not diabetic"))

####################################################################

## code care home residents as 1, other as 0

cc.all$cats3 <-  as.integer(cc.all$care.home == "Care/nursing home")
cc.all$cats3[cc.all$cats3==0 & cc.all$listed.any==1] <- 2 
cc.all$cats3[cc.all$cats3==0 & cc.all$listed.any==0] <- 3
cc.all$cats3 <- as.factor(car::recode(cc.all$cats3, "1='Care/nursing home';
                               2='Independent living, listed condition';
                               3='Independent living, no listed condition'"))     

############ save hospitalized non-severe for later use as training dataset #############

hosp <- cc.all$casegroup=="B"
## merge BNF chapters, one variable per subpara
chnums = 1:13
cc.hosp <- merge.bnfsubparas(chnums=chnums, data=cc.all[hosp, ])

saveRDS(cc.hosp, file="cchosp.rds")
rm(cc.hosp)

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

########## restrict to severe cases and matched controls ###################### 

cat("Restricting to severe cases and matched controls\n")
cc.severe <- subset(cc.all, subset=casegroup=="A")

############ FIXME: with enough memory, this can be done on cc.all ###

## merge drugs, one variable per chapter
scrips.wide <- reshape2::dcast(scrips, ANON_ID ~ chapternum, fun.aggregate=length, 
                               value.var="chapternum")
shortnames.cols <-  bnfchapters$shortname[match(as.integer(colnames(scrips.wide)[-1]),
                                                as.integer(bnfchapters$chapternum))]
colnames(scrips.wide)[-1] <- paste("BNF", colnames(scrips.wide)[-1], shortnames.cols,
                                   sep="_")
scrips.wide <- as.data.table(scrips.wide, key="ANON_ID")
cc.severe <- scrips.wide[cc.severe]

bnfcols <- grep("^BNF", colnames(cc.severe))
for(j in bnfcols) {
    cc.severe[[j]][is.na(cc.severe[[j]])] <- 0
    cc.severe[[j]][cc.severe[[j]] > 1] <- 1
    cc.severe[[j]] <- as.factor(cc.severe[[j]])
}

## merge BNF chapters, one variable per subpara
cat("Merging BNF subparagraph codes ...")
chnums = 1:13
cc.severe <- merge.bnfsubparas(chnums=chnums, data=cc.severe)
cat("done\n")

subparacols <- grep("^subpara\\.", names(cc.severe))
x <- apply(subset(cc.severe, select=subparacols), 2, as.character)
x <- apply(x, 2, as.integer)
cc.severe$numdrugs.subpara <- rowSums(x)
cc.severe$numdrugsgr <- 3 * ceiling(cc.severe$numdrugs.subpara / 3)
cc.severe$numdrugsgr <- as.factor(car::recode(cc.severe$numdrugsgr,
                                              "3='1 to 3'; 6='4 to 6'; 9='7 to 9';
                                              12='10 to 12'; 15:hi='>12'"))
cc.severe[, numdrugsgr := factor(numdrugsgr,
                                 levels=levels(numdrugsgr)[c(2, 3, 5, 6, 4, 1)])]
table(cc.severe$numdrugsgr)

cc.severe$numdrugs.notppi <- cc.severe$numdrugs.subpara - cc.severe$y.protonpump
cc.severe$numdrugs.notppi.gr <- 3 * ceiling(cc.severe$numdrugs.notppi / 3)
cc.severe$numdrugs.notppi.gr <- as.factor(car::recode(cc.severe$numdrugs.notppi.gr, "15:hi='>12'"))
cc.severe[, numdrugs.notppi.gr := factor(numdrugs.notppi.gr,
                                        levels=levels(numdrugs.notppi.gr)[c(2, 4:6, 3, 1)])]
table(cc.severe$numdrugs.notppi.gr)

cardiovasc.subparacols <- grep("^subpara\\.2", names(cc.severe))
x <- apply(subset(cc.severe, select=cardiovasc.subparacols), 2, as.character)
x <- apply(x, 2, as.integer)
cc.severe$numdrugs.cardiovasc <- rowSums(x)
cc.severe$numdrugs.notcardiovasc <- cc.severe$numdrugs.subpara - cc.severe$numdrugs.cardiovasc

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

diagnoses.wide <- as.data.table(diagnoses.wide, key="ANON_ID")
cc.severe <- diagnoses.wide[cc.severe]

icdcols <- grep("^Ch.", colnames(cc.severe))
for(j in icdcols) {
    cc.severe[[j]][is.na(cc.severe[[j]])] <- 0
    cc.severe[[j]][cc.severe[[j]] > 1] <- 1
    cc.severe[[j]] <- as.factor(cc.severe[[j]])
}

cc.severe$num.icdchapters <-
    rowSums(matrix(as.integer(as.matrix(subset(cc.severe, select=icdcols))),
                   nrow=nrow(cc.severe)))
cc.severe$num.icdchapters.gr <-
    as.factor(car::recode(cc.severe$num.icdchapters,
                          "0='No discharge records'; 1:2='1-2 ICD-10 chapters';
                           3:hi='3 or more chapters'"))
cc.severe[, num.icdchapters.gr := relevel(num.icdchapters.gr, ref="No discharge records")]

###########################################

source("comorbidity.R")

## 8 listed conditions designated by NHS
listed.conditions <- c("dm.type", "IHD.any", "heart.other.any", "oad.any",
                       "ckd.any", "neuro.any", "liver.any", "immune.any")

############ extract predefined disease categories #################
## as these are coded as factors, lowest level will be 1

cc.severe$listed.any <-
    as.factor(as.integer(with(cc.severe,
                              diabetes.any=="Diabetic" | IHD.any==1 |
                              heart.other.any==1 |
                              ckd.any==1 | oad.any==1 |
                              neuro.any==1 | liver.any==1 | immune.any==1)))

cc.severe$diag.other <- as.integer(cc.severe$listed.any==0 & cc.severe$diag.any==1)
cc.severe$diag.other <- as.factor(car::recode(cc.severe$diag.other,
                                  "0='Listed condition or no admission';
                                   1='No listed condition, but other admission diagnosis'"))


####################################################################

## code care home residents as 1, other as 0
cc.severe$cats3 <-  as.integer(cc.severe$care.home == "Care/nursing home")
cc.severe$cats3[cc.severe$cats3==0 & cc.severe$listed.any==1] <- 2 
cc.severe$cats3[cc.severe$cats3==0 & cc.severe$listed.any==0] <- 3
cc.severe$cats3 <- as.factor(car::recode(cc.severe$cats3, "1='Care/nursing home';
                               2='Independent living, listed condition';
                               3='Independent living, no listed condition'"))     

########### variable lists for tabulating

demog <- c("ethnic3.onomap", "SIMD.quintile", "care.home")

demog.smr <- c("ethnic4.smr", "SIMD.quintile", "care.home")

varnames.extended <- c("care.home", "scrip.any", "diag.any", "listed.any", "diag.other",
                     "scripordiag", listed.conditions)

varnames.listedpluscovs <- c("care.home", "scrip.any", "diag.any", listed.conditions)

bnfcols <- grep("^BNF", colnames(cc.severe))
bnf.chapternames <- colnames(cc.severe)[bnfcols]

icdcols <- grep("^Ch.", colnames(cc.severe))
icd.chapternames <- colnames(cc.severe)[icdcols]

full.varnames <- c(demog, listed.conditions,
                   "diag.any", icd.chapternames, "scrip.any", bnf.chapternames)

## logical vectors for subsetting
nocare <- cc.severe$care.home=="Independent"
notlisted <- cc.severe$listed.any == 0
nocare.notlisted <- nocare & notlisted

##### tables for first paper  ######################

table.severe.demog <- tabulate.freqs.regressions(varnames=demog,
                                                 data=cc.severe)

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
    rownames(z) <- "No listed condition, but other diagnosis"
    
    x <- rbind(x[1, , drop=FALSE], y, z, x[2, , drop=FALSE])
    table.scripordiag <- cbind(table.scripordiag, x)
}

table.scripordiag.fatal <- NULL
for(agegr in levels(cc.severe$agegr3)) {
    x <- with(cc.severe[cc.severe$agegr3==agegr, ],
              table(scripordiag, fatalcase))
    colnames(x) <- c("Controls", "Fatal cases")
    x <- paste.colpercent(x)
    # x <- x[1, , drop=FALSE]
    rownames(x) <- c("No scrip or diagnosis", "Scrip or diagnosis")
    table.scripordiag.fatal <- rbind(table.scripordiag.fatal, x)
}

cc.severe$scripordiag <- as.factor(cc.severe$scripordiag)

## table.agegr shown in paper 1: frequencies by age group but no rate ratios
## varnames.listed
table.agegr <- NULL
for(agegr in levels(cc.severe$agegr20)) {
    x <- univariate.tabulate(varnames=varnames.extended, 
                             outcome="CASE",
                             data=cc.severe[cc.severe$agegr20==agegr, ],
                             drop.reflevel=FALSE)
    table.agegr <- cbind(table.agegr, x)
}
## all ages including 28-day fatality
freqs.all <- univariate.tabulate(varnames=c("deathwithin28", varnames.extended),
                             outcome="CASE",
                             data=cc.severe,
                             drop.reflevel=FALSE)

## list of tables by age group shown in paper 1
## varnames.listedpluscovs
tables.agegr <- vector("list", length(levels(cc.severe$agegr3)))
for(i in 1:length(levels(cc.severe$agegr3))) {
    agegr <- levels(cc.severe$agegr3)[i]
    tables.agegr[[i]] <-
        tabulate.freqs.regressions(varnames=varnames.listedpluscovs,
                                   data=cc.severe[cc.severe$agegr3==agegr, ])
}

## table.agegr.all shown in paper 1: rate ratios for all age groups combined
## varnames.listedpluscovs
table.agegr.all <- tabulate.freqs.regressions(varnames=varnames.listedpluscovs,
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

## listed conditions
table.listed.conditions.lt60 <-
    tabulate.freqs.regressions(varnames=listed.conditions,
                               data=cc.severe[cc.severe$AGE < 60, ])
table.listed.conditions.ge60 <-
    tabulate.freqs.regressions(varnames=listed.conditions,
                               data=cc.severe[cc.severe$AGE >= 60 ])

## full multivariate model variables -- for this use dm.type rather than diabetes.any 
multivariate.all <-
    multivariate.clogit(varnames=full.varnames,
                        data=cc.severe, add.reflevel=TRUE)

################# restrict to those without listed conditions #############

## conditions
table.conditions.aug <- tabulate.freqs.regressions(varnames=icd.chapternames, 
                                                   data=cc.severe[notlisted, ])

cat("Tabulating ICD subchapter diagnoses ...")
## tabulate subchapters in ICD chapters of interest
table.icdchapter2 <- tabulate.icdchapter(chnum=2, data=cc.severe[notlisted, ])
table.icdchapter7 <-  tabulate.icdchapter(chnum=7, data=cc.severe[notlisted, ])
table.icdchapter11 <- tabulate.icdchapter(chnum=11, data=cc.severe[notlisted, ])

## table icdsubchapters shown in paper 1
table.icdsubchapters <- NULL
for(i in 1:20) {
    table.icdsubchapters <-
        rbind(table.icdsubchapters,
              tabulate.icdchapter(chnum=i, data=cc.severe[notlisted, ], minrowsum=50))
}

#table.icdsubchapters <- table.icdsubchapters[grep("ensuremath",
#                                                  table.icdsubchapters$u.pvalue), ]
cat("done\n")

#########################################################################

## drugs 
table.drugs.aug <- tabulate.freqs.regressions(varnames=bnf.chapternames, 
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

################################################################

## tabulate associations with drug chapters in those not in care homes and without listed conditions 
table.drugs.nocare.notlisted <- tabulate.freqs.regressions(varnames=bnf.chapternames, 
                                                           data=cc.severe[nocare.notlisted, ])

######## stepwise regressions use saved version #####################

nfold <- 4
source("stepwise.R")

#####################################################
if(old) {
    rmarkdown::render("casecontrol.Rmd", output_file="casecontrol15May.pdf")
} else {
    rmarkdown::render("casecontrol.Rmd", output_file="casecontrol8June.pdf")
    source("pharmaco.R")
    rmarkdown::render("pharmaco.Rmd", output_file="pharmaco.pdf")
} 

## remove large objects from memory

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

if(old & include.nrs.only) {
    saveRDS(table.agegr, file="table.agegr.withNRS.rds")
}

rmarkdown::render("casecontrol_script.Rmd", output_file="casecontrol_script.pdf")

Rprof()
print(summaryRprof(tmp)$by.total[1:20, ])
