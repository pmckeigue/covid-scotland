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
        
#cc.all <- readRDS("./data/CC_linked_ANON_20200501.rds")
cc.all <- readRDS("./data/CC_linked_ANON_20200501 (2).rds")

diagnoses <- readRDS("./data/CC_SMR01_ICD10_x25_ANON_20200501.rds") # 842 records ? excluded
procedures <- readRDS("./data/CC_SMR01_OPCS4_MAIN.x25_ANON_20200501.rds")
scotpop <- read_excel("./Scotland_midyearpop_est2019.xlsx")

## we have only 7 digits on scrips, giving resolution to subpara level only
scrips <- readRDS("./data/CC_PIS_x15_ANON_20200501.rds")
scrips$paracode <- as.integer(substr(scrips$bnf_paragraph_code, 1, 6))
length(table(substr(scrips$bnf_paragraph_code, 1, 2))) # chapter
length(table(substr(scrips$bnf_paragraph_code, 1, 4))) # chapter, section
length(table(substr(scrips$bnf_paragraph_code, 1, 6)))  # chapter, section, paragraph
length(table(scrips$bnf_paragraph_code)) # 537 groups
scrips$chapternum <- as.integer(substr(scrips$bnf_paragraph_code, 1, 2))
scrips$sectioncode <- as.integer(substr(scrips$bnf_paragraph_code, 1, 4))
## recode scrips$bnf.chapter values > 14 or NA to 14
scrips$chapternum[is.na(scrips$chapternum)] <- 14
scrips$chapternum[scrips$chapternum > 14] <- 14


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
cc.all$ethnic4 <- factor(cc.all$ethnic4, levels=levels(cc.all$ethnic4)[c(1, 4, 2, 3)])

## recode ONOMAP ethnicity to 3 categories: White, South Asian, Other
cc.all$ethnic3 <- car::recode(cc.all$ethnic4, "'Chinese'='Other'")
cc.all <- within(cc.all, ethnic3 <- relevel(as.factor(ethnic3), ref="White"))
cc.all$ethnic3 <- factor(cc.all$ethnic3, levels=levels(cc.all$ethnic3)[c(1, 3, 2)])

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
cc.all$dm.type <- factor(cc.all$dm.type, levels=levels(cc.all$dm.type)[c(1, 3, 4, 2)])

###############################

ids.protonpump <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 103])
cc.all$protonpump <- as.factor(as.integer(cc.all$ANON_ID %in% ids.protonpump))
cc.all$y.protonpump <- as.integer(cc.all$protonpump =="1")

scrips.protonpump <- scrips$ANON_ID[as.integer(scrips$sectioncode) == 103]
scrips.protonpump <- table(scrips.protonpump)
scrips.protonpump <- data.frame(ANON_ID=as.integer(names(scrips.protonpump)),
                                scrips.protonpump=as.integer(scrips.protonpump))
cc.all <- merge(cc.all, scrips.protonpump, by="ANON_ID", all.x=TRUE)
cc.all$scrips.protonpump[is.na(cc.all$scrips.protonpump)] <- 0
cc.all$scrips.protonpump <- as.factor(cc.all$scrips.protonpump)
                                               
cc.all$protonpump <- as.factor(as.integer(cc.all$ANON_ID %in% ids.protonpump))
cc.all$y.protonpump <- as.integer(cc.all$protonpump =="1")


ids.nonopioid.analgesic <- unique(scrips$ANON_ID[as.integer(scrips$paracode) == 40701])
cc.all$nonopioid.analgesic <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                   ids.nonopioid.analgesic))

ids.antiplatelet <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 209])
cc.all$antiplatelet <- as.factor(as.integer(cc.all$ANON_ID %in% ids.antiplatelet))

ids.nsaid <- unique(scrips$ANON_ID[as.integer(scrips$paracode) == 100101])
cc.all$nsaid <- as.factor(as.integer(cc.all$ANON_ID %in% ids.nsaid))

ids.opioid.analgesic <- unique(scrips$ANON_ID[as.integer(scrips$paracode) == 40702])
cc.all$opioid.analgesic <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                      ids.opioid.analgesic))

ids.antipsychotic <- unique(scrips$ANON_ID[as.integer(scrips$paracode) == 40201])
cc.all$antipsychotic <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                      ids.antipsychotic))

ids.osmotic.laxative <- unique(scrips$ANON_ID[as.integer(scrips$paracode) == 10604])
cc.all$osmotic.laxative <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                   ids.osmotic.laxative))



###############################################################################

table.protonpump.testpos <- NULL
withcovariates.formula <- as.formula("CASE ~ protonpump + antiplatelet + nsaid + antipsychotic +  nonopioid.analgesic + strata(stratum)")
coeff.row <- 1

for(agegr in levels(cc.all$agegr20)) {
    x <- tabulate.freqs.regressions(varnames="protonpump",
                                    data=cc.all[cc.all$agegr20==agegr, ])[, 1:4]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.all[cc.all$agegr20==agegr, ]))$coefficients[coeff.row, , drop=FALSE]
    x$m.ci <- or.ci(y[, 1], y[, 3])
    x$m.pvalue <- pvalue.latex(y[, 5])
    rownames(x) <- agegr
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.protonpump.testpos <- rbind(table.protonpump.testpos, x)
}
x <- tabulate.freqs.regressions(varnames="protonpump", data=cc.all)[, 1:4]
y <- summary(clogit(formula=withcovariates.formula,
                    data=cc.all))$coefficients[coeff.row, , drop=FALSE]
x$m.ci <- or.ci(y[, 1], y[, 3])
x$m.pvalue <- pvalue.latex(y[, 5])
rownames(x) <- "All"
colnames(table.protonpump.testpos)[1:2] <- colnames(x)[1:2]
table.protonpump.testpos <- rbind(table.protonpump.testpos, x)

tabulate.freqs.regressions(varnames=c("care.home", "SIMD.quintile",
                                      "protonpump", "antiplatelet", "nsaid",
                                      "opioid.analgesic", "nonopioid.analgesic"), 
                           data=cc.all[cc.all$AGE < 40, ])


###############################################################################

cc.hosp <- cc.all[cc.all$casegroup=="Critical care or fatal" |
                  cc.all$casegroup=="Hospitalised, not severe", ]

table.protonpump.hosp <- NULL
withcovariates.formula <- as.formula("CASE ~ protonpump + nonopioid.analgesic + strata(stratum)")
coeff.row <- 1

for(agegr in levels(cc.hosp$agegr20)) {
    x <- tabulate.freqs.regressions(varnames="protonpump",
                                    data=cc.hosp[cc.hosp$agegr20==agegr, ])[, 1:4]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.hosp[cc.hosp$agegr20==agegr, ]))$coefficients[coeff.row, , drop=FALSE]
    x$m.ci <- or.ci(y[, 1], y[, 3])
    x$m.pvalue <- pvalue.latex(y[, 5])
    rownames(x) <- agegr
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.protonpump.hosp <- rbind(table.protonpump.hosp, x)
}
x <- tabulate.freqs.regressions(varnames="protonpump", data=cc.hosp)[, 1:4]
y <- summary(clogit(formula=withcovariates.formula,
                    data=cc.hosp))$coefficients[coeff.row, , drop=FALSE]
x$m.ci <- or.ci(y[, 1], y[, 3])
x$m.pvalue <- pvalue.latex(y[, 5])
rownames(x) <- "All"
colnames(table.protonpump.hosp)[1:2] <- colnames(x)[1:2]
table.protonpump.hosp <- rbind(table.protonpump.hosp, x)

tabulate.freqs.regressions(varnames=c("SIMD.quintile",
                                      "protonpump", "antiplatelet", "nsaid",
                                      "opioid.analgesic", "nonopioid.analgesic"), 
                           data=cc.hosp[cc.hosp$AGE < 40 & cc.hosp$HAI==0, ])

########## restrict to severe cases and matched controls ###################### 
 
cc.severe <- cc.all[cc.all$casegroup=="Critical care or fatal", ]

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

## merge BNF chapters 1 and 2, one variable per subpara

chnums = c(1, 2)
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


############ extract predefined disease categories #################
cc.severe$listed.any <-
    as.factor(as.integer(with(cc.severe,
                              !(dm.type == "Not diabetic") | IHD.any==1 |
                              heart.other.any==1 |
                              ckd.any==1 | oad.any==1 |
                              neuro.any==1 | liver.any==1 | immune.any==1)))



########### variable lists for tabulating

demog <- c("ethnic3", "SIMD.quintile", "care.home")
demog.smr <- c("ethnic4.smr", "SIMD.quintile", "care.home")

bnf.chapternames <- colnames(cc.severe)[bnfcols]
drugs <- bnf.chapternames

subparanames <- colnames(cc.severe)[grep("subpara", colnames(cc.severe))]

icd.chapternames <- colnames(cc.severe)[icdcols]
conditions <- icd.chapternames

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
                               data=cc.all[cc.all$casegroup=="Hospitalised, not severe" |
                                           cc.all$casegroup=="Critical care or fatal", ])

table.hospitalized.demog.ethnicsmr <-
    tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home", "SIMD.quintile"),
                               data=cc.all[(cc.all$casegroup=="Hospitalised, not severe" |
                                            cc.all$casegroup=="Critical care or fatal") &
                                           !is.na(cc.all$ethnic5.smr), ])

### incidence and mortality using national population estimates #####

source("incidencemortality.R")

###############################################################

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

#### 5 ethnic groups for report
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
## all variables
multivariate.all <-
    multivariate.clogit(varnames=c(demog, "dm.type", listed.conditions,
                                   "diag.any", conditions, "scrip.any", drugs,
                                   "protonpump"),
                        data=cc.severe, add.reflevel=TRUE)

                                                          
################# restrict to those without listed conditions #############
cc.nocare <- cc.severe[cc.severe$care.home=="Independent", ]
cc.notlisted <- cc.severe[cc.severe$listed.any == 0, ]
cc.nocare.notlisted <- cc.nocare[cc.nocare$listed.any == 0, ]

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

#############################################################################
## tabulate scrip.any effect by carehome
table.scrip.any.carehome <- NULL
for(residence in levels(cc.severe$care.home)) {
    x <- tabulate.freqs.regressions(varnames="scrip.any",
                                    data=cc.severe[cc.severe$care.home==residence, ])[, 1:4]
    rownames(x) <- residence
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.scrip.any.carehome <- rbind(table.scrip.any.carehome, x)
}

## tabulate scrip.any effect by listed.any in those not resident in care homes
table.scrip.any.nocare.listed <- NULL
for(listed in levels(cc.nocare$listed.any)) {
    x <- tabulate.freqs.regressions(varnames="scrip.any",
                                    data=cc.nocare[cc.nocare$listed.any==listed, ])[, 1:4]
    rownames(x) <- listed
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.scrip.any.nocare.listed <- rbind(table.scrip.any.nocare.listed, x)
}

################################################################

## backwards selection of smallest subset of BNF chapters that explains most of the scrip.any effect
    
x <- cc.nocare.notlisted[, drugs]
x <- matrix(as.integer(as.matrix(x)), nrow=nrow(x))
colnames(x) <- drugs
y <- cc.nocare.notlisted$CASE
stratum <- cc.nocare.notlisted$stratum

stepwise.drop <- stepwise.union.dropcols(x=x, y=y, stratum=stratum)

x <- cc.nocare.notlisted[, subparanames]
x <- matrix(as.integer(as.matrix(x)), nrow=nrow(x))
colnames(x) <- subparanames

cat("Stepwise drop procedure over subparas in BNF chapters 1 and 2 ...")
stepwise.drop.subparas <- stepwise.union.dropcols(x=x, y=y, stratum=stratum)
cat("done\n")

## tabulate associations with drug chapters in those not in care homes and without listed conditions 
table.drugs.nocare.notlisted <- tabulate.freqs.regressions(varnames=drugs, 
                                                           data=cc.nocare.notlisted)

## tabulate proportion of effect of scrip.any that is explained by each chapter
## use cc.nocare.notlisted
table.anyscrip.chapter <-
    summary(clogit(formula=as.formula(paste("CASE ~ scrip.any + strata(stratum)")),
                                      data=cc.nocare.notlisted))$coefficients[1, 1:2, drop=FALSE]
for(bnfchapter in drugs) {
    ch.formula=as.formula(paste("CASE ~ scrip.any +", bnfchapter, "+ strata(stratum)"))
    table.anyscrip.chapter <-
        rbind(table.anyscrip.chapter,
              summary(clogit(ch.formula, data=cc.nocare.notlisted))$coefficients[1, 1:2])
}
rownames(table.anyscrip.chapter) <- c("Unadjusted", drugs)
table.anyscrip.chapter <- as.data.frame(table.anyscrip.chapter)
table.anyscrip.chapter$prop.explained <- round(with(table.anyscrip.chapter,
                                                    c(0, 1 - coef[-1] / coef[1])), 2)


## tabulate para or subpara codes in BNF chapters of interest

table.bnfchapter1 <- tabulate.bnfparas(chnum=1, data=cc.nocare.notlisted)
table.bnfchapter2 <- tabulate.bnfsubparas(chnum=2, data=cc.nocare.notlisted)
table.bnfchapter4 <- tabulate.bnfsubparas(chnum=4, data=cc.nocare.notlisted)
table.bnfchapter9 <- tabulate.bnfsubparas(chnum=9, data=cc.nocare.notlisted)
table.bnfchapter10 <- tabulate.bnfsubparas(chnum=10, data=cc.nocare.notlisted)

## fix to tabulate BNF chemical substance


####################################################

tabulate.freqs.regressions(varnames=c("care.home",
                                      "neoplasm.any", "neuro.any", 
                                      "esoph.stomach.duod",
                                      "antiplatelet",
                                      "antipsychotic", 
                                      "nsaid", "opioid.analgesic", "nonopioid.analgesic",
                                      "osmotic.laxative", "protonpump"),
                           data=cc.severe[cc.severe$AGE < 60, ])

############# proton pump #########################

####### effects of scrip and protonpump by care home status ######################

table.protonpump.carehome <- NULL
for(residence in levels(cc.severe$care.home)) {
    x <- tabulate.freqs.regressions(varnames="protonpump",
                                    data=cc.severe[cc.severe$care.home==residence, ])[, 1:4]
    rownames(x) <- residence
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.protonpump.carehome <- rbind(table.protonpump.carehome, x)
}

## tabulate associations with and without covariate adjustment, excluding care home residents
table.scrips.protonpump <- tabulate.freqs.regressions(varnames="scrips.protonpump",
                                                      data=cc.severe)


table.protonpump <- NULL
withcovariates.formula <- as.formula("CASE ~ care.home + esoph.stomach.duod + nsaid + antiplatelet +  protonpump + strata(stratum)")
coeff.row <- 5

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

withcovariates.formula <- as.formula("CASE ~ neoplasm.any + nsaid + antiplatelet + protonpump + nonopioid.analgesic + strata(stratum)")

table.nonopioid.analgesic <- NULL
for(agegr in levels(cc.nocare$agegr20)) {
    x <- tabulate.freqs.regressions(varnames="nonopioid.analgesic",
                                    data=cc.nocare[cc.nocare$agegr20==agegr, ])[, 1:4]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.nocare[cc.nocare$agegr20==agegr, ]))$coefficients[coeff.row, , drop=FALSE]
    x$m.ci <- or.ci(y[, 1], y[, 3])
    x$m.pvalue <- pvalue.latex(y[, 5])
    rownames(x) <- agegr
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.nonopioid.analgesic <- rbind(table.nonopioid.analgesic, x)
}
x <- tabulate.freqs.regressions(varnames="nonopioid.analgesic", data=cc.nocare)[, 1:4]
y <- summary(clogit(formula=withcovariates.formula,
                    data=cc.nocare))$coefficients[coeff.row, , drop=FALSE]
x$m.ci <- or.ci(y[, 1], y[, 3])
x$m.pvalue <- pvalue.latex(y[, 5])
rownames(x) <- "All"
colnames(x)[1:2] <- colnames(table.nonopioid.analgesic)[1:2]
table.nonopioid.analgesic <- rbind(table.nonopioid.analgesic, x)

##################################################################################

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


######## stepwise regressions use saved version #####################
nfold <- 4
#stepwise <- TRUE
stepwise <- FALSE

source("stepwise.R")

#####################################################################

rmarkdown::render("casecontrol.Rmd", output_file="casecontrol.pdf")
rmarkdown::render("pharmaco.Rmd", output_file="pharmaco.pdf")
#rmarkdown::render("Covid_ethnicity_Scotland.Rmd")
 
library(infotheo)

X <- with(cc.severe, data.frame(protonpump, scrip.any))

mi <- mean(unlist(lapply(X=split(X, cc.severe$stratum),
       FUN=function(X)
           infotheo::mutinformation(X)[2, 1])))

y.analgesic <- as.integer(cc.severe$nonopioid.analgesic)

