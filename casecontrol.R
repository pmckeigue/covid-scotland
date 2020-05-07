## analysis script for case-control study

library(car)
library(survival)
library(MASS)
library(wevid)
library(rmarkdown)
library(pander)
library(ggplot2)
library(doParallel)
registerDoParallel(cores=2)
library(reshape2)
library(readxl)
library(DescTools)
library(icd.data)

source("helperfunctions.R")
        
cc.all <- readRDS("./data/CC_linked_ANON_20200501.rds")

diagnoses <- readRDS("./data/CC_SMR01_ICD10_x25_ANON_20200501.rds") # 842 records ? excluded
procedures <- readRDS("./data/CC_SMR01_OPCS4_MAIN.x25_ANON_20200501.rds")
scrips <- readRDS("./data/CC_PIS_x15_ANON_20200501.rds")
bnfcodes <- read_excel("./data/BNF_Code_Information.xlsx", sheet=4)
colnames(bnfcodes) <- c("chaptername", "chapternum", "sectionname", "sectioncode")
    
bnfchapters <- bnfcodes[, 1:2]
bnfchapters <- base::unique(bnfchapters)
#colnames(bnfchapters) <- c("fullname", "chapternum")
truncate.at <- StrPos(bnfchapters$chaptername, " |,|-") -1
truncate.at[is.na(truncate.at)] <- nchar(bnfchapters$chaptername[is.na(truncate.at)])
bnfchapters$shortname <- substr(bnfchapters$chaptername, 1, truncate.at) 

names(cc.all) <- gsub("CASE_NO", "stratum", names(cc.all))
names(cc.all) <- gsub("^SEX$", "sex", names(cc.all))
names(cc.all) <- gsub("imumune", "immune", names(cc.all))
names(cc.all) <- gsub("^ethnic$", "ethnic.old", names(cc.all))
names(cc.all) <- gsub("^CAREHOME$", "care.home", names(cc.all))
names(cc.all) <- gsub("^simd$", "SIMD.quintile", names(cc.all))

## exclude controls already dead on date of testing case they were matched to
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

cc.all$SIMD.quintile <- recode(cc.all$SIMD.quintile, "'Unknown'=NA")

cc.all$scrip.any <- as.factor(as.integer(cc.all$ANON_ID %in% scrips$ANON_ID))
cc.all$diag.any <- as.factor(as.integer(cc.all$ANON_ID %in% diagnoses$ANON_ID))

cc.all$emerg <- as.factor(cc.all$emerg.x25)
cc.all$icu.hdu.ccu <- as.factor(cc.all$icu.hdu.ccu.x25)
cc.all$inpat <- as.factor(cc.all$inpat.x25)

######## coding ethnicity ##############################

OnolyticsType <- cc.all$OnolyticsType
GeographicalArea <- cc.all$GeographicalArea
ethnic.smr <- cc.all$ETHNIC_SMR_LAST

source("ethnic_assign.R")

cc.all$ethnic5 <- ethnic5
cc.all$ethnic5.smr <- ethnic5.smr

## tabulate ONOMAP ethnicity against SMR ethnicity
table.ethnic <- table(cc.all$ethnic5, cc.all$ethnic5.smr, exclude=NULL)
table.ethnic <- paste.colpercent(table.ethnic)

## recode to 4 categories to pool groups not distinguished by name
cc.all$ethnic4 <- car::recode(cc.all$ethnic5, "'Black'='Other'; 'Chinese'='Other'")
cc.all <- within(cc.all, ethnic4 <- relevel(as.factor(ethnic4), ref="White"))

####################################################################

cc.all$sex <- car::recode(as.factor(cc.all$sex), "1='Male'; 2='Female'")
cc.all <- within(cc.all, sex <- relevel(sex, ref="Female"))

cc.all$agegr20 <- 20 * floor(0.05 * (cc.all$AGE - 15)) + 15
cc.all$agegr20 <- as.factor(car::recode(cc.all$agegr20,
                              "-5:15='0-34'; 35='35-54';
                                  55='55-74'; 75:hi='75 or more'"))

cc.all$care.home <- car::recode(as.factor(cc.all$care.home), "0='Independent'; 1='Care home resident'")
cc.all <- within(cc.all, care.home <- relevel(care.home, ref="Independent"))

## all cases have nonmissing SPECDATE
cc.all$deathwithin28 <- with(cc.all,
                             CASE==1 &
                             Date.Death - SPECDATE >= 0 &
                             Date.Death - SPECDATE <= 28)
cc.all$deathwithin28[is.na(cc.all$deathwithin28)] <- 0

## integer values > 1 for icu and inhosp may represent days from test to admission
## values of 0 must be for those not admitted, as there are no missing values
with(cc.all[cc.all$CASE==1, ], table(inhosp, exclude=NULL))
with(cc.all[cc.all$CASE==1, ], table(icu, exclude=NULL))

## temporary coding of case groups -- check this is correct
## integer values 0 for icu and inhosp may represent days from test to entry
cc.all$group <- NA
cc.all$group[cc.all$CASE==1 & cc.all$inhosp== 0] <- "C"
cc.all$group[cc.all$CASE==1 & cc.all$inhosp > 0] <- "B"
cc.all$group[(cc.all$CASE==1 & cc.all$icu > 0) | cc.all$deathwithin28==1]  <- "A"
table(cc.all$deathwithin28, cc.all$group, exclude=NULL)

## assign controls to same group as cases i.e. A, B, C and create a new variable named casegroup
casegroups <- cc.all[cc.all$CASE==1, ][, c("stratum", "group")]
colnames(casegroups)[2] <- "casegroup"
cc.all <- merge(cc.all, casegroups, by=c("stratum"), all.x=T)
table(cc.all$CASE, cc.all$casegroup)

with(cc.all[cc.all$CASE==1, ], table(casegroup, deathwithin28, exclude=NULL))

cc.all$casegroup <- car::recode(cc.all$casegroup,
                           "'A'='Critical care or fatal'; 'B'='Hospitalised, not severe'; 'C'='Test-positive, not hospitalised'")
cc.all$casegroup <- as.factor(cc.all$casegroup)
cc.all <- within(cc.all, casegroup <- relevel(casegroup, ref="Critical care or fatal"))

## recode diabetes type
cc.all$dm.type <- as.integer(cc.all$dm.type)
## missing recoded as zero
cc.all$dm.type[is.na(cc.all$dm.type)] <- 0
cc.all$diabetes.any <- as.factor(car::recode(cc.all$dm.type, 
                                             "0='Not diabetic'; 1:hi='Diabetic'"))
cc.all <- within(cc.all, diabetes.any <- relevel(diabetes.any, ref="Not diabetic"))
cc.all$dm.type <- 
  as.factor(car::recode(cc.all$dm.type, 
                        "0='Not diabetic'; 1='Type 1'; 2='Type 2';
                         3:hi='Other diabetes type'"))
cc.all <- within(cc.all, dm.type <- relevel(dm.type, ref="Not diabetic"))


########################################################################

## tabulate ethnicity by case group
testpositives.ethnic <- paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                              table(ethnic5, casegroup)), 1)
testpositives.ethnic.smr <- paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                                  table(ethnic5.smr, casegroup)), 1)

########### restrict to severe cases and matched controls ###################### 
 
cc.severe <- cc.all[cc.all$casegroup=="Critical care or fatal", ]

## merge drugs
length(table(substr(scrips$bnf_paragraph_code, 1, 2))) # chapter
length(table(substr(scrips$bnf_paragraph_code, 1, 4))) # chapter, section
length(table(substr(scrips$bnf_paragraph_code, 1, 6)))  # chapter, section, paragraph
length(table(scrips$bnf_paragraph_code)) # 537 groups

scrips$chapternum <- as.integer(substr(scrips$bnf_paragraph_code, 1, 2))
scrips$sectioncode <- as.integer(substr(scrips$bnf_paragraph_code, 1, 4))

## recode scrips$bnf.chapter values > 15 to 15
scrips$chapternum[is.na(scrips$chapternum)] <- 15
scrips$chapternum[scrips$chapternum > 15] <- 15

## drop records in bnfchapters with chapternum > 15 and label this category "Other"
bnfchapters <- bnfchapters[bnfchapters$chapternum <=15, ]
bnfchapters$shortname[bnfchapters$chapternum==15] <- "Other" 

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
cc.all <- within(cc.severe, num.icdchapters <- relevel(num.icdchapters, ref="No discharge records"))

if(FALSE) {

cc.all$kidney.any <- as.integer(cc.all$kidney.advanced==1 | cc.all$kidney.other==1)

cc.all$asthma.any <- as.integer(cc.all$asthma12.bnf==1 | cc.all$asthma==1)

cc.all$cysticfibrosis.any <- as.integer(cc.all$cystic_fibrosis12.bnf==1 |
                                        cc.all$cystic.fibrosis==1)

cc.all$tuberculosis.any <- as.integer(cc.all$tuberculosis==1 |
                                      cc.all$tuberculosis12.bnf==1)

cc.all$otherresp.any <- as.integer(cc.all$resp.other==1 |
                                  cc.all$other_chronic_lrd12.bnf)

cc.all$respinf.orTB <- as.integer(cc.all$resp.inf==1 |
                                  cc.all$tuberculosis.any)

cc.all$HIV.any <- as.integer(cc.all$hiv12.bnf==1 |
                             cc.all$HIV==1)

cc.all$immune.any <- as.integer(cc.all$immune==1 |
                                cc.all$transplant.not.kidney==1 |
                                cc.all$cytotoxic_immune6.bnf==1 |
                                cc.all$HIV.any)

cc.all$IHD.any <- as.integer(cc.all$ihd12.bnf==1 | cc.all$IHD==1)

cc.all$connectivetissue.any <- as.integer(cc.all$connective_tissue12.bnf==1 | cc.all$connective==1)

cc.all$epilepsy.any <- as.integer(cc.all$epilepsy12.bnf==1 | cc.all$epilepsy==1)

cc.all$otherneuro.any <- as.integer(cc.all$ms12.bnf==1 | cc.all$neuro.other==1)

cc.all$antihypertensive.any <- with(cc.all,
                                    as.integer(vasodilator_ah6.bnf==1 |
                                               centrally_acting_ah6.bnf==1 |
                                               adrenergic_neurone_block6.bnf==1 |
                                               alpha_adrenoceptor6.bnf==1 |
                                               ace6.bnf==1 |
                                               angio6.bnf==1 |
                                               renin_angiotensin6.bnf==1 |
                                               thiazides6.bnf==1 |
                                               calcium_channel6.bnf==1))
cc.all$antihypertensive.other <- with(cc.all,
                                      as.integer(vasodilator_ah6.bnf==1 |
                                                 centrally_acting_ah6.bnf==1))

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

demog <- c("ethnic4", "SIMD.quintile", "care.home")

drugs <- c("scrip.any", colnames(cc.severe)[bnfcols])
conditions <- c("diag.any", "dm.type", colnames(cc.severe)[icdcols])

if(FALSE) {
demog.smr <- c("ETHNIC_smr", "simd2020_sc_decile", "care.home")

conditions <- c("dm.type", "antihypertensive.any", "IHD.any", "heart.other",
                "CVD",  "circulatory.other",
                "asthma.any", "respinf.orTB", "resp.other",
                "cysticfibrosis.any",
                "kidney.any", "connectivetissue.any", 
                "epilepsy.any", "mono.poly.neuro", "otherneuro.any",
                "blood.cancer", "lung.cancer", "other.cancer",
                "immune.any") # , "HIV.any")

drugs <- eval(c("antihypertensive.any",
                antihypertensives, # "antihypertensive.other",
                "anticoagulants6.bnf", "nsaids6.bnf",
                "lipid_regulating6.bnf", "statins6.bnf",
                "hydroxychloroquine6.bnf"))

admissions <- "SMR01"
}

lookup.names <- data.frame(varname=c("scrip.any", "diag.any", "care.home",
                                     "emerg", "icu.hdu.ccu", "inpat",
                                     "diabetes.any",
                                     "antihypertensive.any",
                                     "IHD.any",
                                     "CVD",
                                     "heart.other",
                                     "circulatory.other",
                                     "asthma.any",
                                     "respinf.orTB",
                                     "resp.other",
                                     "cysticfibrosis.any",
                                     "kidney.any",
                                     "connectivetissue.any",
                                     "mono.poly.neuro",
                                     "epilepsy.any",
                                     "otherneuro.any",
                                     "blood.cancer",
                                     "lung.cancer",
                                     "other.cancer",
                                     "immune.any",
                                     "alpha_adrenoreceptor6.bnf",
                                     "ace6.bnf", "angio6.bnf",
                                     "renin_angiotensin6.bnf",
                                     "thiazides6.bnf",
                                     "calcium_channel6.bnf",
                                     "antihypertensive.other",
                                     "anticoagulants6.bnf",
                                     "nsaids6.bnf",
                                     "lipid_regulating6.bnf",
                                     "statins6.bnf",
                                     "hydroxychloroquine6.bnf"
                                     ),
                           longname=c("Any prescription", "Any admission", "Care home",
                                      "Emergency admission last year",
                                      "Critical care admission last year",
                                      "Any admission last year",
                                      "Diabetes (any type)",
                                      "Any antihypertensive",
                                      "Ischaemic heart disease",
                                      "Cerebrovascular disease",
                                      "Other heart disease",
                                      "Other circulatory disease",
                                      "Asthma",
                                      "Respiratory infections",
                                      "Other respiratory disease",
                                      "Cystic fibrosis",
                                      "Kidney disease",
                                      "Connective tissue disease",
                                      "Neuropathy (mono- or poly-)",
                                      "Epilepsy",
                                      "Other neurological conditions",
                                      "Cancer of blood-forming organs",
                                      "Lung cancer",
                                      "Other cancer",
                                      "Any immune deficiency / suppression", 
                                      "alpha-adrenoreceptor blocker",
                                      "ACE inhibitor", "Angiotensin-II receptor blocker",
                                      "ACE or A-IIR inhibitor",
                                      "Thiazides",
                                      "Calcium channel blocker",
                                      "Other antihypertensive",
                                      "Anticoagulants",
                                      "Non-steroidal anti-inflammatory drugs",
                                      "Lipid-regulating agents",
                                      "Statins",
                                      "Hydroxychloroquine"
                                      ))



############### ethnicity for report ########################################

## case-control analysis for ONOMAP ethnicity 

table.ethnic <- paste.colpercent(table(cc.severe$ethnic5, cc.severe$CASE), 1)
univariate.ethnic <- summary(clogit(formula=CASE ~ ethnic5 + strata(stratum),
                                    data=cc.severe))$coefficients
univariate.ethnic <- rbind(rep(NA, 2), univariate.ethnic[, c(2, 5)])
table.ethnic <- data.frame(table.ethnic, univariate.ethnic)
colnames(table.ethnic)[1:2] <- paste(c("Controls", "Cases"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(cc.severe$CASE[!is.na(cc.severe$ethnic)]))))
colnames(table.ethnic)[3:4] <- c("Rate ratio", "p-value")

## case-control analysis for SMR ethnicity
table.ethnic.smr <- paste.colpercent(table(cc.severe$ethnic5.smr, cc.severe$CASE), 1)
univariate.ethnic.smr <- summary(clogit(formula=CASE ~ ethnic5.smr + strata(stratum),
                                        data=cc.severe))$coefficients
univariate.ethnic.smr <- rbind(rep(NA, 2), univariate.ethnic.smr[, c(2, 5)])
table.ethnic.smr <- data.frame(table.ethnic.smr, univariate.ethnic.smr)
colnames(table.ethnic.smr)[1:2] <- paste(c("Controls", "Cases"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(cc.severe$CASE[!is.na(cc.severe$ethnic5.smr)]))))
colnames(table.ethnic.smr)[3:4] <- c("Rate ratio", "p-value")

###############################################################
                                        #
table.agegr <- NULL
for(agegr in levels(cc.severe$agegr20)) {
    x <- univariate.tabulate(varnames=c("care.home", "scrip.any",
                                        "num.icdchapters", "diabetes.any"),
                             outcomevar="CASE",
                             data=cc.severe[cc.severe$agegr20==agegr, ],
                             drop.referencelevel=FALSE)
    table.agegr <- cbind(table.agegr, x)
}

###################### restrict by age
cc.severe <- cc.severe[cc.severe$AGE < 75, ] ## restrict by age 

############################# demographic vars ##########

table.demog <- univariate.tabulate(varnames=demog, outcomevar="CASE", data=cc.severe)
colnames(table.demog) <- c("Controls", "Cases")
univariate.demog <- NULL
for(i in 1:length(demog)) {
    univariate.formula <- as.formula(paste("CASE ~ ", demog[i], "+ strata(stratum)"))
    x <- summary(clogit(formula=univariate.formula, data=cc.severe))$coefficients
    univariate.demog <- rbind(univariate.demog, x)
}

################# drugs ######################################

table.drugs <- univariate.tabulate(varnames=drugs, outcomevar="CASE", data=cc.severe)
colnames(table.drugs) <- c("Controls", "Cases")
univariate.drugs <- NULL
for(i in 1:length(drugs)) {
    univariate.formula <- as.formula(paste("CASE ~ ", drugs[i], "+ strata(stratum)"))
    x <- summary(clogit(formula=univariate.formula, data=cc.severe))$coefficients
    univariate.drugs <- rbind(univariate.drugs, x)
}

############################ conditions #################

table.conditions <- univariate.tabulate(varnames=conditions, outcomevar="CASE", data=cc.severe)
colnames(table.conditions) <- c("Controls", "Cases")
univariate.conditions <- NULL
for(i in 1:length(conditions)) {
    univariate.formula <- as.formula(paste("CASE ~ ", conditions[i], "+ strata(stratum)"))
    x <- summary(clogit(formula=univariate.formula, data=cc.severe))$coefficients
    univariate.conditions <- rbind(univariate.conditions, x)
}
    
########## multivariate regression models #####################

## demog
demog.formula <- as.formula(paste("CASE ~",
                                  paste(demog, collapse=" + "),
                                "+ strata(stratum)"))
cc.severe.demog.nonmissing <- nonmissing.obs(cc.severe, (demog))
demog.model <- clogit(formula=demog.formula, data=cc.severe.demog.nonmissing)
multivariate.demog <- summary(demog.model)$coefficients

## drugs
drugs.formula <- as.formula(paste("CASE ~",
                                  paste(drugs, collapse=" + "),
                                "+ strata(stratum)"))
cc.severe.drugs.nonmissing <- nonmissing.obs(cc.severe, (drugs))
drugs.model <- clogit(formula=drugs.formula, data=cc.severe.drugs.nonmissing)
multivariate.drugs <- summary(drugs.model)$coefficients

## conditions
conditions.formula <- as.formula(paste("CASE ~",
                                  paste(conditions, collapse=" + "),
                                "+ strata(stratum)"))
cc.severe.conditions.nonmissing <- nonmissing.obs(cc.severe, (conditions))
conditions.model <- clogit(formula=conditions.formula, data=cc.severe.conditions.nonmissing)
multivariate.conditions <- summary(conditions.model)$coefficients

## all
all.formula <- as.formula(paste("CASE ~",
                                paste(demog, collapse=" + "), "+",
                                paste(conditions, collapse=" + "), "+", 
                                paste(drugs, collapse=" + "), 
                                "+ strata(stratum)"))
cc.severe.all.nonmissing <- nonmissing.obs(cc.severe, c(demog, conditions, drugs))
all.model <- clogit(formula=all.formula, data=cc.severe.all.nonmissing)
multivariate.all <- summary(all.model)$coefficients

#####################################################################

lower.formula <- "CASE ~ care.home + strata(stratum)"

nfold <- 2
run.stepwise <- FALSE
if(run.stepwise) { # stepwise variable selection ######################

    ## FIXME: redefine demog.model etc to include a baseline variable
## demog
stepwise.demog <- step(demog.model,
                       scope=list(lower=lower.formula, upper=demog.formula),
                       direction="both", method="approximate", trace=-1)
stepwise.demog <- summary(stepwise.demog)$coefficients
rownames(stepwise.demog) <- replace.names(rownames(stepwise.demog))
print(stepwise.demog)

## drugs
stepwise.drugs <- step(drugs.model,
                       scope=list(lower=lower.formula, upper=drugs.formula),
                       direction="both", method="approximate", trace=-1)
stepwise.drugs <- summary(stepwise.drugs)$coefficients
rownames(stepwise.drugs) <- replace.names(rownames(stepwise.drugs))
print(stepwise.drugs)

## conditions
stepwise.conditions <- step(conditions.model,
                       scope=list(lower=lower.formula, upper=conditions.formula),
                       direction="both", method="approximate", trace=-1)
stepwise.conditions <- summary(stepwise.conditions)$coefficients
rownames(stepwise.conditions) <- replace.names(rownames(stepwise.conditions))
print(stepwise.conditions)

## all
stepwise.all <- step(all.model,
                       scope=list(lower=lower.formula, upper=all.formula),
                       direction="both", method="approximate", trace=-1)
stepwise.all <- summary(stepwise.all)$coefficients
rownames(stepwise.all) <- replace.names(rownames(stepwise.all))
print(stepwise.all)

############# cross-validation of stepwise regression #########################

## demog
test.folds <- testfolds.bystratum(stratum=cc.severe.demog.nonmissing$stratum,
                                  y=cc.severe.demog.nonmissing$CASE, nfold=nfold)
cv.data <- merge(cc.severe.demog.nonmissing, test.folds, by="stratum")

##demog.predicted <- foreach(i=1:nfold, .combine=rbind) %do% {
demog.predicted <- NULL
for(i in 1:nfold) {
    test.data  <- cv.data[cv.data$test.fold == i, ]
    train.data <- cv.data[cv.data$test.fold != i, ]
    start.model <- clogit(formula=demog.formula, data=train.data)
    stepwise.model <- step(start.model,
                           scope=list(lower=lower.formula, upper=demog.formula),
                           direction="both", trace=-1, method="approximate")
    ## normalize within each stratum
    unnorm.p <- predict(object=stepwise.model, newdata=test.data,
                        na.action="na.pass", 
                        type="risk", reference="sample")
    norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                            stratum=test.data$stratum,
                                            y=test.data$CASE)
    demog.predicted <- rbind(demog.predicted, norm.predicted)
}
demog.densities <- with(demog.predicted, Wdensities(y, posterior.p, prior.p))
pander(summary(demog.densities), table.style="multiline",
       split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
       caption="Prediction of severe COVID-19 from demographic variables only")

## drugs
test.folds <- testfolds.bystratum(stratum=cc.severe.drugs.nonmissing$stratum,
                                  y=cc.severe.drugs.nonmissing$CASE, nfold=nfold)
cv.data <- merge(cc.severe.drugs.nonmissing, test.folds, by="stratum")

##drugs.predicted <- foreach(i=1:nfold, .combine=rbind) %do% {
drugs.predicted <- NULL
for(i in 1:nfold) {
    test.data  <- cv.data[cv.data$test.fold == i, ]
    train.data <- cv.data[cv.data$test.fold != i, ]
    start.model <- clogit(formula=drugs.formula, data=train.data)
    stepwise.model <- step(start.model,
                           scope=list(lower=lower.formula, upper=drugs.formula),
                           direction="both", trace=-1, method="approximate")
    ## normalize within each stratum
    unnorm.p <- predict(object=stepwise.model, newdata=test.data,
                        na.action="na.pass", 
                        type="risk", reference="sample")
    norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                            stratum=test.data$stratum,
                                            y=test.data$CASE)
    drugs.predicted <- rbind(drugs.predicted, norm.predicted)
}
drugs.densities <- with(drugs.predicted, Wdensities(y, posterior.p, prior.p))
   pander(summary(drugs.densities), table.style="multiline",
           split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
           caption="Prediction of severe COVID-19 from drugs only")

## conditions
test.folds <- testfolds.bystratum(stratum=cc.severe.conditions.nonmissing$stratum,
                                  y=cc.severe.conditions.nonmissing$CASE, nfold=nfold)
cv.data <- merge(cc.severe.conditions.nonmissing, test.folds, by="stratum")

conditions.predicted <- NULL
for(i in 1:nfold) {
    test.data  <- cv.data[cv.data$test.fold == i, ]
    train.data <- cv.data[cv.data$test.fold != i, ]
    start.model <- clogit(formula=conditions.formula, data=train.data)
    stepwise.model <- step(start.model,
                           scope=list(lower=lower.formula, upper=conditions.formula),
                           direction="both", trace=-1, method="approximate")
    ## normalize within each stratum
    unnorm.p <- predict(object=stepwise.model, newdata=test.data,
                        na.action="na.pass", 
                        type="risk", reference="sample")
    norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                            stratum=test.data$stratum,
                                            y=test.data$CASE)
    conditions.predicted <- rbind(conditions.predicted, norm.predicted)
}
conditions.densities <- with(conditions.predicted, Wdensities(y, posterior.p, prior.p))
   pander(summary(conditions.densities), table.style="multiline",
           split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
           caption="Prediction of severe COVID-19 from conditions only")

## all
test.folds <- testfolds.bystratum(stratum=cc.severe.all.nonmissing$stratum,
                                  y=cc.severe.all.nonmissing$CASE, nfold=nfold)
cv.data <- merge(cc.severe.all.nonmissing, test.folds, by="stratum")

all.predicted <- NULL
for(i in 1:nfold) {
    test.data  <- cv.data[cv.data$test.fold == i, ]
    train.data <- cv.data[cv.data$test.fold != i, ]
    start.model <- clogit(formula=all.formula, data=train.data)
    stepwise.model <- step(start.model,
                           scope=list(lower=lower.formula, upper=all.formula),
                           direction="both", trace=-1, method="approximate")
    ## normalize within each stratum
    unnorm.p <- predict(object=stepwise.model, newdata=test.data,
                        na.action="na.pass", 
                        type="risk", reference="sample")
    norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                            stratum=test.data$stratum,
                                            y=test.data$CASE)
    all.predicted <- rbind(all.predicted, norm.predicted)
}
    all.densities <- with(all.predicted, Wdensities(y, posterior.p, prior.p))
    pander(summary(all.densities), table.style="multiline",
           split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
           caption="Prediction of severe COVID-19 from all variables")

save(stepwise.demog, stepwise.drugs, stepwise.conditions, stepwise.all,
     demog.densities, drugs.densities, conditions.densities, all.densities,
     file="./data/stepwise.RData")
} else {
    load("./data/stepwise.RData")
}

#################################################

tables.icd.list <- vector("list", 14)
for(i in 1:14) {
    tables.icd.list[[i]] <- tabulate.icdsubchapter(i)
}
    
tables.bnf.list <- vector("list", 14)
for(i in 1:14) {
    tables.bnf.list[[i]] <- tabulate.bnfsection(i)
}

    
