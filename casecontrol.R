## analysis script for case-control study

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

registerDoParallel(cores=2)

source("helperfunctions.R")
        
cc.all <- readRDS("./data/CC_linked_ANON_20200501.rds")

diagnoses <- readRDS("./data/CC_SMR01_ICD10_x25_ANON_20200501.rds") # 842 records ? excluded
procedures <- readRDS("./data/CC_SMR01_OPCS4_MAIN.x25_ANON_20200501.rds")
scotpop <- read_excel("./Scotland_midyearpop_est2019.xlsx")

## we have only 7 digits on scrips, giving resolution to subpara level only
scrips <- readRDS("./data/CC_PIS_x15_ANON_20200501.rds")
scrips$paracode <- as.integer(substr(scrips$bnf_paragraph_code, 1, 6))

source("bnfcodes.R")

names(cc.all) <- gsub("CASE_NO", "stratum", names(cc.all))
names(cc.all) <- gsub("^SEX$", "sex", names(cc.all))
names(cc.all) <- gsub("imumune", "immune", names(cc.all))
names(cc.all) <- gsub("^ethnic$", "ethnic.old", names(cc.all))
names(cc.all) <- gsub("^CAREHOME$", "care.home", names(cc.all))
names(cc.all) <- gsub("^simd$", "SIMD.quintile", names(cc.all))

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
cc.all <- within(cc.all, ethnic5.smr <- relevel(as.factor(ethnic5.smr), ref="White"))

## tabulate ONOMAP ethnicity against SMR ethnicity
table.ethnic <- table(cc.all$ethnic5, cc.all$ethnic5.smr, exclude=NULL)
table.ethnic <- paste.colpercent(table.ethnic)

## recode SMR ethnicity to 4 categories: White, Black, South Asian, Other
cc.all$ethnic4.smr <- car::recode(cc.all$ethnic5.smr, "'Chinese'='Other'")
cc.all <- within(cc.all, ethnic4.smr <- relevel(as.factor(ethnic4.smr), ref="White"))

## recode ONOMAP ethnicity to 4 categories: White, South Asian, Chinese, Other
cc.all$ethnic4 <- car::recode(cc.all$ethnic5, "'Black'='Other'")
cc.all <- within(cc.all, ethnic4 <- relevel(as.factor(ethnic4), ref="White"))

## recode ONOMAP ethnicity to 3 categories: White, South Asian, Other
cc.all$ethnic3 <- car::recode(cc.all$ethnic4, "'Chinese'='Other'")
cc.all <- within(cc.all, ethnic3 <- relevel(as.factor(ethnic3), ref="White"))

####################################################################

cc.all$sex <- car::recode(as.factor(cc.all$sex), "1='Male'; 2='Female'")
cc.all <- within(cc.all, sex <- relevel(sex, ref="Female"))

#cc.all$agegr20 <- 20 * floor(0.05 * (cc.all$AGE - 15)) + 15
cc.all$agegr20 <- as.factor(car::recode(as.integer(cc.all$AGE),
                              "0:39='0-39'; 40:59='40-59';
                                  60:74='60-74'; 75:hi='75 or more'"))

cc.all$care.home <- car::recode(as.factor(cc.all$care.home), "0='Independent'; 1='Care home resident'")
cc.all <- within(cc.all, care.home <- relevel(care.home, ref="Independent"))

## all cases have nonmissing SPECDATE
## FIXME: recode this for controls also
cc.all$deathwithin28 <- with(cc.all,
                             CASE==1 &
                             Date.Death - SPECDATE >= 0 &
                             Date.Death - SPECDATE <= 28)
cc.all$deathwithin28[is.na(cc.all$deathwithin28)] <- 0
cc.all$deathwithin28 <- as.factor(recode(cc.all$deathwithin28, "0='Alive'; 1='Dead'"))

## integer values > 1 for icu and inhosp may represent days from test to entry
## values of 0 must be for those not admitted, as there are no missing values
with(cc.all[cc.all$CASE==1, ], table(inhosp, exclude=NULL))
with(cc.all[cc.all$CASE==1, ], table(icu, exclude=NULL))

## coding of case groups -- check this is correct
cc.all$group <- NA
cc.all$group[cc.all$CASE==1 & cc.all$inhosp== 0] <- "C"
cc.all$group[cc.all$CASE==1 & cc.all$inhosp > 0] <- "B"
cc.all$group[(cc.all$CASE==1 & cc.all$icu > 0) | cc.all$deathwithin28=="Dead"]  <- "A"
table(cc.all$deathwithin28, cc.all$group, exclude=NULL)

## assign controls to same group as matched case i.e. A, B, C and create a new variable named casegroup
casegroups <- cc.all[cc.all$CASE==1, ][, c("stratum", "group")]
colnames(casegroups)[2] <- "casegroup"
cc.all <- merge(cc.all, casegroups, by=c("stratum"), all.x=T)
table(cc.all$CASE, cc.all$casegroup)

with(cc.all[cc.all$CASE==1, ], table(casegroup, deathwithin28, exclude=NULL))

cc.all$fatalcase <- as.integer(cc.all$CASE==1 & cc.all$deathwithin28=="Dead")

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
                        "c(0, 10)='Not diabetic';
                         c(1, 101, 102)='Type 1';
                         c(2, 202, 203)='Type 2'; 
                         3:9='Other diabetes type';
                         11:100='Other diabetes type'"))

cc.all <- within(cc.all, dm.type <- relevel(dm.type, ref="Not diabetic"))

########################################################################

## tabulate ethnicity by case group
testpositives.ethnic <- paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                              table(ethnic4, casegroup)), 1)
testpositives.ethnic.smr <- paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                                  table(ethnic5.smr, casegroup)), 1)

testpositives.carehome <- paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                                table(ethnic4, care.home)), 0)

testpositives.healthboard <- t(paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                                table(ethnic4, HBRES_NAME)), 0))

table.severe.demog <-
    tabulate.freqs.regressions(varnames=c("ethnic4", "care.home", "SIMD.quintile"),
                               data=cc.severe)


table.testpositives.demog <-
    tabulate.freqs.regressions(varnames=c("ethnic4", "care.home", "SIMD.quintile"),
                               data=cc.all)
table.testpositives.demog.ethnicsmr <-
    tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home", "SIMD.quintile"),
                               data=cc.all[!is.na(cc.all$ethnic5.smr), ])

table.hospitalized.demog <-
    tabulate.freqs.regressions(varnames=c("ethnic4", "care.home", "SIMD.quintile"),
                               data=cc.all[cc.all$casegroup=="Hospitalized, not severe" |
                                           cc.all$casegroup=="Critical care or fatal", ])

table.hospitalized.demog.ethnicsmr <-
    tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home", "SIMD.quintile"),
                               data=cc.all[(cc.all$casegroup=="Hospitalized, not severe" |
                                            cc.all$casegroup=="Critical care or fatal") &
                                           !is.na(cc.all$ethnic5.smr), ])
 
########### restrict to severe cases and matched controls ###################### 
 
cc.severe <- cc.all[cc.all$casegroup=="Critical care or fatal", ]

## merge drugs
length(table(substr(scrips$bnf_paragraph_code, 1, 2))) # chapter
length(table(substr(scrips$bnf_paragraph_code, 1, 4))) # chapter, section
length(table(substr(scrips$bnf_paragraph_code, 1, 6)))  # chapter, section, paragraph
length(table(scrips$bnf_paragraph_code)) # 537 groups

scrips$chapternum <- as.integer(substr(scrips$bnf_paragraph_code, 1, 2))
scrips$sectioncode <- as.integer(substr(scrips$bnf_paragraph_code, 1, 4))

## recode scrips$bnf.chapter values > 14 or NA to 14
scrips$chapternum[is.na(scrips$chapternum)] <- 14
scrips$chapternum[scrips$chapternum > 14] <- 14

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
cc.severe <- within(cc.severe, num.icdchapters <- relevel(num.icdchapters, ref="No discharge records"))

############ extract predefined disease categories #################

## nitrates are BNF code 020601
ids.icd.IHD <- unique(diagnoses$ANON_ID[grep("^I2[0-5]", diagnoses$ICD10)])
ids.bnf.IHD <- unique(scrips$ANON_ID[substr(as.character(scrips$bnf_paragraph_code), 1, 6) == "020601"])
table(ids.bnf.IHD %in% ids.icd.IHD)
## should add procedure codes for angioplasty etc
ids.IHD <- unique(c(ids.icd.IHD, ids.bnf.IHD))
cc.severe$IHD.any <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.IHD))

ids.icd.heart.other <- unique(diagnoses$ANON_ID[grep("^I0[01256789]|^I1[0-5]|^I2[6-8]|^I3[0-9]|^I4[0-9]|^I5[0-2]",
                                  diagnoses$ICD10)])
ids.bnf.heart.other <- unique(scrips$ANON_ID[substr(as.character(scrips$bnf_paragraph_code), 1, 4) == "0203"])  # anti-arrhythmics
table(ids.bnf.heart.other %in% ids.icd.heart.other)
ids.heart.other <- unique(c(ids.icd.heart.other, ids.bnf.heart.other))
cc.severe$heart.other.any <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.heart.other))

ids.icd.ckd <- unique(diagnoses$ANON_ID[grep("^N18[3-5]|^Z49[0-2]|^Z94[02]",
                                             diagnoses$ICD10)])
ids.kidneytransplant <- unique(procedures$ANON_ID[grep("^M01[1234589]",
                                                       procedures$MAIN_OPERATION)])
table(ids.kidneytransplant %in% ids.icd.ckd)
ids.ckd.any <- unique(c(ids.icd.ckd, ids.kidneytransplant))
cc.severe$ckd.any <-  as.factor(as.integer(cc.severe$ANON_ID %in% ids.ckd.any))

ids.icd.asthma <- unique(diagnoses$ANON_ID[grep("^J4[56]", diagnoses$ICD10)])
ids.icd.chronresp <- unique(diagnoses$ANON_ID[grep("^J4[012347]|^J6[0-9]|^J70|^J8[0-6]|^J9[0-9]|^G47\\.?3",
                                                   diagnoses$ICD10)])
ids.bnf.broncho <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) >= 301 &
                                         as.integer(scrips$sectioncode) <= 303])
table(ids.icd.asthma %in% ids.bnf.broncho)
table(ids.icd.chronresp %in% ids.bnf.broncho)
ids.oad.any <- unique(c(ids.icd.asthma, ids.icd.chronresp, ids.bnf.broncho))
cc.severe$oad.any <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.oad.any))

## include all Nervous chapter except G40 "Episodic and Paroxysmal Disorders"
## also include F03 dementia NOS
ids.icd.neuro <- unique(diagnoses$ANON_ID[grep("^F03|^G[012356789]", diagnoses$ICD10)])
ids.bnf.neuro <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 409 |
                                         as.integer(scrips$sectioncode) == 411])
table(ids.bnf.neuro %in% ids.icd.neuro)
## interferon beta 080204M, Glatiramer acetate 0802040U0, Natalizumab 0802040W0
## Dimethyl fumar 0802040AK, Teriflunomide 0802040AL, Alemtuzumab 0802030
## no records in scrips for these drugs
## 526 records in scrips[substr(scrips$bnf_paragraph_code, 1, 5) == "08020", ]
ids.neuro.any <- unique(c(ids.icd.neuro, ids.bnf.neuro))
cc.severe$neuro.any <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.neuro.any))

liver.grep.string <- "^C22\\.?0|^I85\\.?0|^I98\\.?3|^K70\\.?[234|^K71\\.?7|^K72\\.?[019]|^K72\\.?[019|^K73|^K74\\.?[023456]|^K76\\.?7|^R18"
table(grep(liver.grep.string, diagnoses$ICD10, value=TRUE))
ids.icd.liver <- unique(diagnoses$ANON_ID[grep(liver.grep.string, diagnoses$ICD10)])
cc.severe$liver.any <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.icd.liver))

ids.esoph.stomach.duod <-  unique(diagnoses$ANON_ID[grep("^K2[0-9]|^K3[01]", diagnoses$ICD10)])
cc.severe$esoph.stomach.duod <-
    as.factor(as.integer(cc.severe$ANON_ID %in% ids.esoph.stomach.duod))

## immune.any includes primary immunodeficiency and secondary immunosuppression
ids.icd.immune <- unique(diagnoses$ANON_ID[grep("^B2[0-3|^D8[0-9]", diagnoses$ICD10)])

## 802 other immunomodulating drugs
## Methotrexate and chloroquine appear in musculoskeletal chapter 
ids.bnf.immune <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 802 |
                                        as.integer(scrips$sectioncode) == 411])

ids.immune.any <- unique(c(ids.icd.immune, ids.bnf.immune))
cc.severe$immune.any <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.immune.any))

listed.conditions <- c("IHD.any", "heart.other.any", "oad.any",
                       "ckd.any", "neuro.any", "liver.any", "immune.any")
cc.severe$listed.any <- as.factor(as.integer(with(cc.severe,
                                        ! dm.type == "Not diabetic" | IHD.any==1 | heart.other.any==1 |
                                        ckd.any==1 | oad.any==1 |
                                        neuro.any==1 | liver.any==1 | immune.any==1)))

ids.protonpump <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 103])
cc.severe$protonpump <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.protonpump))

ids.nonopioid.analgesic <- unique(scrips$ANON_ID[as.integer(scrips$paracode) == 40701])
cc.severe$nonopioid.analgesic <- as.factor(as.integer(cc.severe$ANON_ID %in%
                                                      ids.nonopioid.analgesic))

ids.icd.neoplasm <- unique(diagnoses$ANON_ID[grep("^C[0-9]|^D[0-4]", diagnoses$ICD10)])
ids.bnf.neoplasm <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 801])
ids.neoplasm.any <- unique(c(ids.icd.neoplasm, ids.bnf.neoplasm))
cc.severe$neoplasm.any <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.neoplasm.any))


ids.antiplatelet <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 209])
cc.severe$antiplatelet <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.antiplatelet))

ids.nsaid <- unique(scrips$ANON_ID[as.integer(scrips$paracode) == 100101])
cc.severe$nsaid <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.nsaid))

cc.severe$y.protonpump <- as.integer(cc.severe$protonpump =="1")


if(FALSE) { # coding antihypertensives
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

demog <- c("ethnic3", "SIMD.quintile", "care.home")
demog.smr <- c("ethnic4.smr", "SIMD.quintile", "care.home")

bnf.chapternames <- colnames(cc.severe)[bnfcols]
drugs <- bnf.chapternames
icd.chapternames <- colnames(cc.severe)[icdcols]
conditions <- icd.chapternames

lookup.names <- data.frame(varname=c("deathwithin28", "scrip.any", "diag.any", "care.home",
                                     "emerg", "icu.hdu.ccu", "inpat",
                                     "protonpump",
                                     "diabetes.any",
                                     "antihypertensive.any",
                                     "IHD.any",
                                     "CVD",
                                     "heart.other.any",
                                     "circulatory.other",
                                     "oad.any",
                                     "respinf.orTB",
                                     "resp.other",
                                     "cysticfibrosis.any",
                                     "ckd.any",
                                     "neuro.any",
                                     "liver.any",
                                     "listed.any",
                                     "connectivetissue.any",
                                     "mono.poly.neuro",
                                     "epilepsy.any",
                                     "otherneuro.any",
                                     "blood.cancer",
                                     "lung.cancer",
                                     "other.cancer",
                                     "immune.any",
                                     "alpha_adrenoreceptor6.bnf",
                                     "ace6.bnf",
                                     "angio6.bnf",
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
                           longname=c("Death within 28 days of test",
                                      "Any prescription", "Any admission", "Care home",
                                      "Emergency admission last year",
                                      "Critical care admission last year",
                                      "Any admission last year",
                                      "Proton pump inhibitor",
                                      "Diabetes (any type)",
                                      "Any antihypertensive",
                                      "Ischaemic heart disease",
                                      "Cerebrovascular disease",
                                      "Other heart disease",
                                      "Other circulatory disease",
                                      "Asthma or chronic airway disease",
                                      "Respiratory infections",
                                      "Other respiratory disease",
                                      "Cystic fibrosis",
                                      "Chronic kidney disease or transplant recipient",
                                      "Neurological (except epilepsy) or dementia",
                                      "Liver disease",
                                      "Any listed condition", 
                                      "Connective tissue disease",
                                      "Neuropathy (mono- or poly-)",
                                      "Epilepsy",
                                      "Other neurological conditions",
                                      "Cancer of blood-forming organs",
                                      "Lung cancer",
                                      "Other cancer",
                                      "Immune deficiency or suppression", 
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

####### incidence and mortality using national population estimates ######################

case.freqs <- with(cc.severe[cc.severe$CASE==1, ], table(AGE, sex, exclude=NULL))
case.freqs <- data.frame(Age=as.integer(rownames(case.freqs)),
                         Females=as.integer(case.freqs[, 1]),
                         Males=as.integer(case.freqs[, 2]))
case.long <- reshape2::melt(case.freqs, id="Age")
colnames(case.long) <- c("Age", "Sex", "Cases")

death.freqs <- with(cc.severe[cc.severe$fatalcase==1, ], table(AGE, sex, exclude=NULL))
death.freqs <- data.frame(Age=as.integer(rownames(death.freqs)),
                         Females=as.integer(death.freqs[, 1]),
                         Males=as.integer(death.freqs[, 2]))
death.long <- reshape2::melt(death.freqs, id="Age")
colnames(death.long) <- c("Age", "Sex", "Deaths")

scotpop.long <- reshape2::melt(scotpop[, -2], id="Age")
colnames(scotpop.long) <- c("Age", "Sex", "Population")

discrim <- merge(scotpop.long, case.long, by=c("Age", "Sex"), all.x=TRUE)
discrim <- merge(discrim, death.long, by=c("Age", "Sex"), all.x=TRUE)
discrim$Cases[is.na(discrim$Cases)] <- 0
discrim$Deaths[is.na(discrim$Deaths)] <- 0

discrim$Sex <- as.factor(discrim$Sex)

discrim$Noncases <- discrim$Population - discrim$Cases
y.cases <- cbind(as.integer(discrim$Cases), as.integer(discrim$Noncases))

discrim$Survivors <- discrim$Population - discrim$Deaths
y.deaths <- cbind(as.integer(discrim$Deaths), as.integer(discrim$Survivors))


cases.model <- glm(formula=y.cases ~ Sex + Age, family="binomial", data=discrim)
deaths.model <- glm(formula=y.deaths ~ Sex + Age, family="binomial", data=discrim)

cases.model.coeffs <- summary(cases.model)$coefficients
deaths.model.coeffs <- summary(deaths.model)$coefficients
logistic.coeffs <- data.frame(severecase=cases.model.coeffs[, 1],
                              death=deaths.model.coeffs[, 1])

male <- discrim$Sex=="Males"
female <- discrim$Sex=="Females"
gam.model.MaleDeaths <- gam::gam(formula=y.deaths[male, ] ~ s(Age), family=binomial("logit"),
                                 data=discrim[male, ])
gam.model.FemaleDeaths <- gam::gam(formula=y.deaths[female, ] ~ s(Age), family=binomial("logit"),
                                   data=discrim[female, ])
gam.model.MaleCases<- gam::gam(formula=y.cases[male, ] ~ s(Age), family=binomial("logit"),
                               data=discrim[male, ])
gam.model.FemaleCases <- gam::gam(formula=y.cases[female, ] ~ s(Age), family=binomial("logit"),
                                  data=discrim[female, ])

gam.male <- data.frame(Cases=car::logit(gam.model.MaleCases$fitted.values),
                       Deaths=car::logit(gam.model.MaleDeaths$fitted.values),
                       Age=discrim$Age[male])
gam.male.long <- reshape2::melt(data=gam.male, id="Age")
colnames(gam.male.long)[2] <- "Status"
gam.male.long$Sex <- "Males"

gam.female <- data.frame(Cases=car::logit(gam.model.FemaleCases$fitted.values),
                       Deaths=car::logit(gam.model.FemaleDeaths$fitted.values),
                       Age=discrim$Age[female])
gam.female.long <- reshape2::melt(data=gam.female, id="Age")
colnames(gam.female.long)[2] <- "Status"
gam.female.long$Sex <- "Females"
gam <- rbind(gam.male.long, gam.female.long)
     
###############################################################

logodds.posterior <- predict(object=cases.model, newdata=discrim, type="link")
logodds.prior <- log(sum(discrim$Cases) / sum(discrim$Noncases))
log.likratio <- logodds.posterior - logodds.prior
discrim$W <- log.likratio / log(2)
lambda1 <- sum(discrim$W * discrim$Cases) / sum(discrim$Cases)
lambda0 <- sum(-discrim$W * discrim$Noncases) / sum(discrim$Noncases)
cases.Lambda.agesex <- 0.5 * (lambda0 +  lambda1)


logodds.posterior <- predict(object=deaths.model, newdata=discrim, type="link")
logodds.prior <- log(sum(discrim$Deaths) / sum(discrim$Survivors))
log.likratio <- logodds.posterior - logodds.prior
discrim$W <- log.likratio / log(2)
lambda1 <- sum(discrim$W * discrim$Deaths) / sum(discrim$Deaths)
lambda0 <- sum(-discrim$W * discrim$Survivors) / sum(discrim$Survivors)
deaths.Lambda.agesex <- 0.5 * (lambda0 +  lambda1)

########################################################
varnames.listed <- c("care.home", "scrip.any", "diag.any",
                  "diabetes.any", listed.conditions)

table.agegr <- NULL
for(agegr in levels(cc.severe$agegr20)) {
    x <- univariate.tabulate(varnames=c("deathwithin28", varnames.listed, "listed.any"), 
                             outcome="CASE",
                             data=cc.severe[cc.severe$agegr20==agegr, ],
                             drop.reflevel=FALSE)
    table.agegr <- cbind(table.agegr, x)
}

freqs.all <- univariate.tabulate(varnames=c("deathwithin28", varnames.listed, "listed.any"), 
                             outcome="CASE",
                             data=cc.severe,
                             drop.reflevel=FALSE)


  
cc.severe$agegr3 <-
    as.factor(car::recode(cc.severe$AGE,
                          "0:59='0-60 years'; 60:74='60-74 years'; 75:hi='75+ years'"))

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
        tabulate.freqs.regressions(varnames=varnames.listed,
                                   data=cc.severe[cc.severe$agegr3==agegr, ])
}

table.agegr.all <- tabulate.freqs.regressions(varnames=varnames.listed,
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
                               data=cc.severe[cc.severe$AGE >= 60, ])
## all variables
multivariate.all <-
    multivariate.clogit(varnames=c(demog, "dm.type", listed.conditions,
                                   "diag.any", conditions, "scrip.any", drugs,
                                   "protonpump"),
                        data=cc.severe, add.reflevel=TRUE)
                                                          
################# restrict to those without listed conditions #############
cc.nocare <- cc.severe[cc.severe$care.home=="Independent", ]
cc.notlisted <- cc.severe[cc.severe$listed.any == 0, ]

## conditions
table.conditions.aug <- tabulate.freqs.regressions(varnames=conditions, 
                                                   data=cc.notlisted)

## tabulate subchapters in ICD chapters of interest
table.icdchapter2 <- tabulate.icdchapter(chnum=2, data=cc.notlisted)

table.icdchapter11 <- tabulate.icdchapter(chnum=11, data=cc.notlisted)

table.icdchapter7 <-  tabulate.icdchapter(chnum=7, data=cc.notlisted)

tabulate.icdchapter(chnum=8, data=cc.notlisted)

table.icdsubchapters <- NULL
for(i in 1:20) {
    table.icdsubchapters <-
        rbind(table.icdsubchapters,
              tabulate.icdchapter(chnum=i, data=cc.notlisted, minrowsum=50))
}
table.icdsubchapters <- table.icdsubchapters[grep("ensuremath",
                                                  table.icdsubchapters$u.pvalue), ]

#########################################################################

## drugs 
table.drugs.aug <- tabulate.freqs.regressions(varnames=drugs, 
                                              data=cc.notlisted)


## tabulate scrip.any effect by carehome
table.scrip.any.carehome <- NULL
for(residence in levels(cc.severe$care.home)) {
    x <- tabulate.freqs.regressions(varnames="scrip.any",
                                    data=cc.severe[cc.severe$care.home==residence, ])[, 1:4]
    rownames(x) <- residence
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.scrip.any.carehome <- rbind(table.scrip.any.carehome, x)
}


## tabulate proportion of effect of scrip.any that is explained by each chapter
## use cc.nocare
table.anyscrip.chapter <-
    summary(clogit(formula=as.formula(paste("CASE ~ scrip.any + strata(stratum)")),
                                      data=cc.nocare))$coefficients[1, 1:2, drop=FALSE]
for(bnfchapter in drugs) {
    ch.formula=as.formula(paste("CASE ~ scrip.any +", bnfchapter, "+ strata(stratum)"))
    table.anyscrip.chapter <-
        rbind(table.anyscrip.chapter,
              summary(clogit(ch.formula, data=cc.nocare))$coefficients[1, 1:2])
}
rownames(table.anyscrip.chapter) <- c("Unadjusted", drugs)
table.anyscrip.chapter <- as.data.frame(table.anyscrip.chapter)
table.anyscrip.chapter$prop.explained <- round(with(table.anyscrip.chapter,
                                                    c(0, 1 - coef[-1] / coef[1])), 2)

## tabulate para or subpara codes in BNF chapters of interest
table.bnfchapter1 <- tabulate.bnfparas(chnum=1, data=cc.notlisted)

############# proton pump #########################

## tabulate associations with and without covariate adjustment, excluding care home

table.protonpump <- NULL
withcovariates.formula <- as.formula("CASE ~ nsaid + antiplatelet + esoph.stomach.duod + protonpump + strata(stratum)")

for(agegr in levels(cc.nocare$agegr20)) {
    x <- tabulate.freqs.regressions(varnames="protonpump",
                                    data=cc.nocare[cc.nocare$agegr20==agegr, ])[, 1:4]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.nocare[cc.nocare$agegr20==agegr, ]))$coefficients[4, , drop=FALSE]
    x$m.ci <- or.ci(y[, 1], y[, 3])
    x$m.pvalue <- pvalue.latex(y[, 5])
    rownames(x) <- agegr
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.protonpump <- rbind(table.protonpump, x)
}
x <- tabulate.freqs.regressions(varnames="protonpump", data=cc.nocare)[, 1:4]
y <- summary(clogit(formula=withcovariates.formula,
                    data=cc.nocare))$coefficients[4, , drop=FALSE]
x$m.ci <- or.ci(y[, 1], y[, 3])
x$m.pvalue <- pvalue.latex(y[, 5])
rownames(x) <- "All"
colnames(x)[1:2] <- c("Controls", "Cases")
table.protonpump <- rbind(table.protonpump, x)


####### effects of scrip and protonpump by care home status ############################


table.protonpump.carehome <- NULL
for(residence in levels(cc.severe$care.home)) {
    x <- tabulate.freqs.regressions(varnames="protonpump",
                                    data=cc.severe[cc.severe$care.home==residence, ])[, 1:4]
    rownames(x) <- residence
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.protonpump.carehome <- rbind(table.protonpump.carehome, x)
}

##################################################################################
## tabulate proton pump by age group, excluding care home residents

table.nocare.protonpump <- NULL
withcovariates.formula <- as.formula("CASE ~ nsaid + antiplatelet + esoph.stomach.duod + protonpump + strata(stratum)")

for(agegr in levels(cc.severe$agegr20)) {
    x <- tabulate.freqs.regressions(varnames="protonpump",
                                    data=cc.nocare[cc.nocare$agegr20==agegr, ])[, 1:4]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.nocare[cc.nocare$agegr20==agegr, ]))$coefficients[4, , drop=FALSE]
    x$m.ci <- or.ci(y[, 1], y[, 3])
    x$m.pvalue <- pvalue.latex(y[, 5])
    rownames(x) <- agegr
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.nocare.protonpump <- rbind(table.nocare.protonpump, x)
}
x <- tabulate.freqs.regressions(varnames="protonpump", data=cc.nocare)[, 1:4]
y <- summary(clogit(formula=withcovariates.formula,
                    data=cc.nocare))$coefficients[4, , drop=FALSE]
x$m.ci <- or.ci(y[, 1], y[, 3])
x$m.pvalue <- pvalue.latex(y[, 5])
rownames(x) <- "All"
colnames(x)[1:2] <- c("Controls", "Cases")
table.nocare.protonpump <- rbind(table.nocare.protonpump, x)

########################################################################

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




## other BNF chapters
table.bnfchapter4 <- tabulate.bnfsubparas(chnum=4, data=cc.notlisted)
table.bnfchapter2 <- tabulate.bnfsubparas(chnum=2, data=cc.notlisted)
table.bnfchapter9 <- tabulate.bnfsubparas(chnum=9, data=cc.notlisted)

## fix to tabulate BNF chemical substance

#########################################################
nfold <- 4
#stepwise <- TRUE
stepwise <- FALSE

source("stepwise.R")

#####################################################################

rmarkdown::render("casecontrol.Rmd", output_file="casecontrol.pdf")
rmarkdown::render("pharmaco.Rmd", output_file="pharmaco.pdf")
rmarkdown::render("Covid_ethnicity_Scotland.Rmd")
  

library(infotheo)

X <- with(cc.severe, data.frame(protonpump, scrip.any))

mi <- mean(unlist(lapply(X=split(X, cc.severe$stratum),
       FUN=function(X)
           infotheo::mutinformation(X)[2, 1])))

y.analgesic <- as.integer(cc.severe$nonopioid.analgesic)

analgesic.formula=as.formula(paste("CASE ~ care.home +",
                                   paste(conditions, collapse="+"),
                                   "+ nonopioid.analgesic + protonpump + strata(stratum)"))

summary(clogit(formula=analgesic.formula,
               data=cc.severe[cc.severe$listed.any=="0", ]))
