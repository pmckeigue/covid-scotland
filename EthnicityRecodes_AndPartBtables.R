#PART B##
## The changes in ethnicity coding need to run before assigning cc.severe so it doesn't need to run against cc.all and cc.severe seperately
##--------------------------------------------------------------------------------------------------------

#RELEVEL ETHNIC5.SMR TO MATCH ONOMAP
cc.all$ethnic5.smr <- factor(cc.all$ethnic5.smr, levels=levels(cc.all$ethnic5.smr)[c(1,4,3,2,5)])
cc.severe$ethnic5.smr <- factor(cc.severe$ethnic5.smr, levels=levels(cc.severe$ethnic5.smr)[c(1,4,3,2,5)])

cc.all$ethnic5 <- factor(cc.all$ethnic5, levels=levels(cc.all$ethnic5)[c(1,2,3,5,4)])
cc.severe$ethnic5 <- factor(cc.severe$ethnic5, levels=levels(cc.severe$ethnic5)[c(1,2,3,5,4)])

######## coding ethnicity PART B dissaggregate South Asian and Black ethnicities ##############################

## dataset is named cc.all

## this script should be edited as a function that will work with any dataset, with user-specified names for source variables and output variables

## Source variables: 

## ethnic.smr - raw SMR categories
## OnolyticsType - raw Onomap types based on name classification
## Geographical Area - regional groupings of OnolyticsType

## Output variables: ethnic5.smr, eth, eth5

## this script

## (1) recodes ETHNIC_smr to specify Chinese as separate category, then collapses to 5 categories: White, South Asian, Chinese, Black, Other

## (2) derives a new variable "eth" from Onomap types

## (3) collapses eth to eth5, which has five categories as above

## We cannot combine SMR coding with Onomap coding because this will introduce non-differential missclassification of ethnicity between cases and controls, especially for those categorised as Black.


### SMR -----------------------------------------------------------------------------------

## label the codes used in the raw SMR variable ethnic.smr
##  https://www.ndc.scot.nhs.uk/Dictionary-A-Z/Definitions/index.asp?Search=E&ID=243&Title=Ethnic%20Group
#'1A'='Scottish';
#'1B'='Other British';
#'1C'='Irish';
#'1K'='Gypsy/ Traveller';
#'1L'='Polish';
#'1Z'='Any other white ethnic group';
#'2A'='Any mixed or multiple ethnic groups';
#'3F'='Pakistani, Pakistani Scottish or Pakistani British';
#'3G'='Indian, Indian Scottish or Indian British';
#'3H'='Bangladeshi, Bangladeshi Scottish or Bangladeshi British';
#'3J'='Chinese, Chinese Scottish or Chinese British';
#'3Z'='Other Asian, Asian Scottish or Asian British';
#'4D'='African, African Scottish or African British';
#'4Y'='Other African';
#'5C'='Caribbean, Caribbean Scottish or Caribbean British';
#'5D'='Black, Black Scottish or Black British';
#'5Y'='Other Caribbean or Black';
#'6A'='Arab, Arab Scottish or Arab British';
#'6Z'='Other ethnic group';
#'98'='Refused/Not provided by patient';
#'99'='Not Known'")

## extra codes: 1[DEFGHJ] White, 3[ABC] South Asian, 3[E] Chinese, 3D Other Asian, 4[ABCE]  African

##### Part B SMR ----------------------------------------------------------------------------

ethnic.smr <- as.character(cc.all$ETHNIC_SMR_LAST)

ethnic9.smr <- character(length(ethnic.smr))
ethnic9.smr[grep("^1[A-Z]$", ethnic.smr)] <- "White"
ethnic9.smr[grep("^1[L]$", ethnic.smr)] <- "White Polish"
ethnic9.smr[grep("^2[A-Z]$", ethnic.smr)] <- "Other"
ethnic9.smr[grep("^3[BCFH]$", ethnic.smr)] <- "Pakistani/Bangladeshi"
ethnic9.smr[grep("^3[AG]$", ethnic.smr)] <- "Indian"
ethnic9.smr[grep("^3[EJ]$", ethnic.smr)] <- "Chinese"
ethnic9.smr[grep("^3[DZ]$", ethnic.smr)] <- "Other"
ethnic9.smr[grep("^4[AE]$", ethnic.smr)] <- "Caribbean"
ethnic9.smr[grep("^4[BDY]$", ethnic.smr)] <- "African"
ethnic9.smr[grep("^4[C]$", ethnic.smr)] <- "Black"
#ethnic9.smr[grep("^4[F]$", ethnic.smr)] <- "Black" #4F not included in coding above
ethnic9.smr[grep("^5[C]$", ethnic.smr)] <- "Caribbean"
ethnic9.smr[grep("^5[ABDY]$", ethnic.smr)] <- "Black"
#ethnic9.smr[grep("^4[Z]$", ethnic.smr)] <- "Other" #not included in coding above
ethnic9.smr[grep("^6[AZ]$", ethnic.smr)] <- "Other"
ethnic9.smr[grep("^9", ethnic.smr)] <- NA
ethnic9.smr[ethnic9.smr==""] <- NA
ethnic9.smr <- as.factor(ethnic9.smr)
ethnic9.smr <- factor(ethnic9.smr, levels=levels(ethnic9.smr)[c(8,9,3,1,2,4,7,5,6)])

### ONOLYTICS -----------------------------------------------------------------------------------

if(length(OnolyticsType) > 0) {
  OnolyticsType <- car::recode(OnolyticsType,
                               "'NOT FOUND'=NA; 'INTERNATIONAL'=NA; 'UNCLASSIFIED'=NA; 'VOID'=NA; 'VOID - FORENAME'=NA; 'VOID INITIAL'=NA")
  
  table(OnolyticsType[GeographicalArea=="SOUTH ASIA"])
  table(OnolyticsType[GeographicalArea=="AFRICA"])
  table(OnolyticsType[GeographicalArea=="BRITISH ISLES"])
  table(OnolyticsType[GeographicalArea=="EAST ASIA"])
  table(OnolyticsType[GeographicalArea=="MIDDLE EAST"])
  
##### Part B ONOLYTICS ------------------------------------------------------------------------------------
  
  eth8 <- rep("Other", nrow(cc.all))
  eth8[is.na(OnolyticsType)] <- NA
  eth8[OnolyticsType=="CHINESE" |
        OnolyticsType=="HONG KONGESE" |
        OnolyticsType=="SINGAPORESE"] <- "Chinese"
  eth8[GeographicalArea=="EAST ASIA" &
        eth != "Chinese"] <- "Other Asia & Pacific"
  
  eth8[GeographicalArea=="SOUTH ASIA"] <- "Other South Asian"
  
  eth8[GeographicalArea=="BRITISH ISLES"] <- "Britain&Ireland"
  eth8[GeographicalArea=="CENTRAL EUROPE" |
       GeographicalArea=="EASTERN EUROPE" |
        GeographicalArea=="NORTHERN EUROPE" |
        GeographicalArea=="SOUTHERN EUROPE" |
        OnolyticsType=="AFRIKAANS"] <- "Other Europe"
  
  eth8[GeographicalArea=="AFRICA" &
        OnolyticsType != "AFRIKAANS" &
        OnolyticsType != "LIBYAN"] <- "Black African"
  eth8[GeographicalArea=="MIDDLE EAST" &
        OnolyticsType != "MUSLIM"] <- "East Med"
  eth8[OnolyticsType == "MUSLIM"] <- "South Asian" # "Muslim, not localized"
  
  eth8[OnolyticsType == "BLACK CARIBBEAN"] <- "Black Caribbean"
  
  eth8[OnolyticsType == "BANGLADESHI"] <- "Muslim South Asian"
  eth8[OnolyticsType == "MUSLIM INDIAN"] <- "Muslim South Asian"
  eth8[OnolyticsType == "PAKISTANI"] <- "Muslim South Asian"
  eth8[OnolyticsType == "PAKISTANI KASHMIR"] <- "Muslim South Asian"
  eth8[OnolyticsType == "MUSLIM"] <- "Muslim South Asian"

  eth8[OnolyticsType == "POLISH"] <- "White Polish"
  
  ## increase categories to 8 categories: White, White Polish, Black Caribbean, Black African, Chinese, Muslim South Asian,
  ## Oter South Asian, Other
  ethnic8 <- car::recode(eth8, "'Britain&Ireland'='White';  'Other Europe'='White';  'East Med'='Other'; 'Other Asia & Pacific'='Other'") 
  ethnic8 <- as.factor(ethnic8)
  ethnic8 <- factor(ethnic8, levels=levels(ethnic8)[c(7,8,2,1,3,4,6,5)])
  ethnic8 <- relevel(as.factor(ethnic8), ref="White")  
  
  
}


## PART B Tables##
## CASE CONTROL TABLES USING ethnic9.smr and ethnic8--------------------------------------------------------------------

cc.all$ethnic9.smr <- ethnic9.smr
cc.all$ethnic8 <- ethnic8

## tabulate ONOMAP ethnicity against SMR ethnicity
table.ethnicB <- table(cc.all$ethnic8, cc.all$ethnic9.smr, exclude=NULL)
table.ethnicB <- paste.colpercent(table.ethnicB)


testpositives.ethnic.smrB <- paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                                  table(ethnic9.smr, casegroup)), 1)

table.testpositives.demog.ethnicsmrB <-
  tabulate.freqs.regressions(varnames=c("ethnic9.smr", "care.home", "SIMD.quintile"),
                             data=cc.all[!is.na(cc.all$ethnic9.smr), ])


table.hospitalized.demog.ethnicsmrB <-
  tabulate.freqs.regressions(varnames=c("ethnic9.smr", "care.home", "SIMD.quintile"),
                             data=cc.all[(cc.all$casegroup=="Hospitalised, not severe" |
                                            cc.all$casegroup=="Critical care or fatal") &
                                           !is.na(cc.all$ethnic9.smr), ])

if(length(OnolyticsType) > 0) {
  ## tabulate ethnicity by case group
  testpositives.ethnicB <- paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                                table(ethnic8, casegroup)), 1)
  
  testpositives.carehomeB <- paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                                  table(ethnic8, care.home)), 0)
  
  testpositives.healthboardB <- t(paste.colpercent(with(cc.all[cc.all$CASE==1, ],
                                                       table(ethnic8, HBRES_NAME)), 0))
  
  table.testpositives.demogB <-
    tabulate.freqs.regressions(varnames=c("ethnic8", "care.home", "SIMD.quintile"),
                               data=cc.all)
  
  table.hospitalized.demogB <-
    tabulate.freqs.regressions(varnames=c("ethnic8", "care.home", "SIMD.quintile"),
                               data=cc.all[cc.all$casegroup=="Hospitalised, not severe" |
                                             cc.all$casegroup=="Critical care or fatal", ])
}

table.ethnic9smrB <-
  tabulate.freqs.regressions(varnames=c("ethnic9.smr", "care.home",
                                        "SIMD.quintile"),
                             outcome="CASE",
                             data=cc.severe[!is.na(cc.severe$ethnic9.smr), ])
rownames(table.ethnic9smrB) <- replace.names(rownames(table.ethnic9smrB))


table.severe.demogB <-
  tabulate.freqs.regressions(varnames=c("ethnic8", "care.home",
                             "SIMD.quintile"),
                             data=cc.severe)


###--------------------------------------------------------------------------------------------------
### PART A - remove NRS deaths from severe tables and re-run
###--------------------------------------------------------------------------------------------------

## tabulate ONOMAP ethnicity against SMR ethnicity
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

testpositives.ethnic.smr <- paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                  table(ethnic5.smr, casegroup)), 1)

table.testpositives.demog.ethnicsmr <-
  tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home", "SIMD.quintile"),
                             data=cc.all[cc.all$nrs_covid_case!=1&!is.na(cc.all$ethnic5.smr), ])


table.hospitalized.demog.ethnicsmr <-
  tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home", "SIMD.quintile"),
                             data=cc.all[(cc.all$casegroup=="Hospitalised, not severe" |
                                            cc.all$casegroup=="Critical care or fatal") & cc.all$nrs_covid_case!=1 &
                                           !is.na(cc.all$ethnic5.smr), ])

if(length(OnolyticsType) > 0) {
  ## tabulate ethnicity by case group
  testpositives.ethnic <- paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                table(ethnic4, casegroup)), 1)
  
  testpositives.carehome <- paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                  table(ethnic4, care.home)), 0)
  
  testpositives.healthboard <- t(paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                       table(ethnic4, HBRES_NAME)), 0))
  
  table.testpositives.demog <-
    tabulate.freqs.regressions(varnames=c("ethnic4", "care.home", "SIMD.quintile"),
                               data=cc.all[cc.all$nrs_covid_case!=1,])
  
  table.hospitalized.demog <-
    tabulate.freqs.regressions(varnames=c("ethnic4", "care.home", "SIMD.quintile"),
                               data=cc.all[(cc.all$casegroup=="Hospitalised, not severe" |
                                             cc.all$casegroup=="Critical care or fatal") &cc.all$nrs_covid_case!=1, ])
}


table.ethnic5smr <-
  tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home",
                                        "SIMD.quintile"),
                             outcome="CASE",
                             data=cc.severe[cc.severe$nrs_covid_case!=1&!is.na(cc.severe$ethnic5.smr), ])
rownames(table.ethnic5smr) <- replace.names(rownames(table.ethnic5smr))


table.severe.demog <-
  tabulate.freqs.regressions(varnames=c("ethnic3", "care.home",
                                        "SIMD.quintile"),
                             data=cc.severe[cc.severe$nrs_covid_case!=1,])



table.ethnic5smrNRS <-
  tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home",
                                        "SIMD.quintile"),
                             outcome="CASE",
                             data=cc.severe[!is.na(cc.severe$ethnic5.smr), ])
rownames(table.ethnic5smrNRS) <- replace.names(rownames(table.ethnic5smrNRS))


table.severe.demogNRS <-
  tabulate.freqs.regressions(varnames=c("ethnic3", "care.home",
                                        "SIMD.quintile"),
                             data=cc.severe)


###--------------------------------------------------------------------------------------------------
###--------------------------------------------------------------------------------------------------
### PART B - remove NRS deaths from severe tables and re-run using disaggregated ethnicity
###--------------------------------------------------------------------------------------------------
###--------------------------------------------------------------------------------------------------

## tabulate ONOMAP ethnicity against SMR ethnicity
if(length(OnolyticsType) > 0) {
  #cc.all$ethnic5 <- ethnic5
  
  ## tabulate ONOMAP ethnicity against SMR ethnicity
  table.ethnicB <- table(cc.all$ethnic8, cc.all$ethnic9.smr, exclude=NULL)
  
  #tn <- table(cc.all$ethnic5, cc.all$ethnic9.smr)
  #SouthAsian.sensitivity <- 100 * tn[5, 2] / sum(tn[, 2])
  #SouthAsian.specificity <- 100 * (sum(tn[, -2]) - sum(tn[5, ]) + tn[5, 2]) / sum(tn[, -2])
  #sum.xtabulate <- sum(tn)
  
  table.ethnicB <- paste.colpercent(table.ethnicB)
  
  ## recode ONOMAP ethnicity to 4 categories: White, South Asian, Chinese, Other
  #cc.all$ethnic4 <- car::recode(cc.all$ethnic5, "'Black'='Other'")
  #cc.all <- within(cc.all, ethnic4 <- relevel(as.factor(ethnic4), ref="White"))
  #cc.all$ethnic4 <- factor(cc.all$ethnic4, levels=levels(cc.all$ethnic4)[c(1, 4, 2, 3)])
  
  ## recode ONOMAP ethnicity to 3 categories: White, South Asian, Other
  #cc.all$ethnic3 <- car::recode(cc.all$ethnic4, "'Chinese'='Other'")
  #cc.all <- within(cc.all, ethnic3 <- relevel(as.factor(ethnic3), ref="White"))
  #cc.all$ethnic3 <- factor(cc.all$ethnic3, levels=levels(cc.all$ethnic3)[c(1, 3, 2)])
}

testpositives.ethnic.smrB <- paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                  table(ethnic9.smr, casegroup)), 1)

table.testpositives.demog.ethnicsmrB <-
  tabulate.freqs.regressions(varnames=c("ethnic9.smr", "care.home", "SIMD.quintile"),
                             data=cc.all[cc.all$nrs_covid_case!=1&!is.na(cc.all$ethnic9.smr), ])


table.hospitalized.demog.ethnicsmrB <-
  tabulate.freqs.regressions(varnames=c("ethnic9.smr", "care.home", "SIMD.quintile"),
                             data=cc.all[(cc.all$casegroup=="Hospitalised, not severe" |
                                            cc.all$casegroup=="Critical care or fatal") & cc.all$nrs_covid_case!=1 &
                                           !is.na(cc.all$ethnic9.smr), ])

if(length(OnolyticsType) > 0) {
  ## tabulate ethnicity by case group
  testpositives.ethnicB <- paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                table(ethnic8, casegroup)), 1)
  
  testpositives.carehomeB <- paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                  table(ethnic8, care.home)), 0)
  
  testpositives.healthboardB <- t(paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                       table(ethnic8, HBRES_NAME)), 0))
  
  table.testpositives.demogB <-
    tabulate.freqs.regressions(varnames=c("ethnic8", "care.home", "SIMD.quintile"),
                               data=cc.all[cc.all$nrs_covid_case!=1,])
  
  table.hospitalized.demogB <-
    tabulate.freqs.regressions(varnames=c("ethnic8", "care.home", "SIMD.quintile"),
                               data=cc.all[(cc.all$casegroup=="Hospitalised, not severe" |
                                              cc.all$casegroup=="Critical care or fatal") &cc.all$nrs_covid_case!=1, ])
}


table.ethnic9smrB <-
  tabulate.freqs.regressions(varnames=c("ethnic9.smr", "care.home",
                                        "SIMD.quintile"),
                             outcome="CASE",
                             data=cc.severe[cc.severe$nrs_covid_case!=1&!is.na(cc.severe$ethnic9.smr), ])
rownames(table.ethnic9smrB) <- replace.names(rownames(table.ethnic9smrB))


table.severe.demogB <-
  tabulate.freqs.regressions(varnames=c("ethnic8", "care.home",
                                        "SIMD.quintile"),
                             data=cc.severe[cc.severe$nrs_covid_case!=1,])



table.ethnic9smrNRSB <-
  tabulate.freqs.regressions(varnames=c("ethnic9.smr", "care.home",
                                        "SIMD.quintile"),
                             outcome="CASE",
                             data=cc.severe[!is.na(cc.severe$ethnic9.smr), ])
rownames(table.ethnic9smrNRSB) <- replace.names(rownames(table.ethnic9smrNRSB))


table.severe.demogNRSB <-
  tabulate.freqs.regressions(varnames=c("ethnic8", "care.home",
                                        "SIMD.quintile"),
                             data=cc.severe)


