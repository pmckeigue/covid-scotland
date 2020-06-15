# Covid-19 case control linkage
# HPS stats support team
# r desktop

# ### set up --------------------------------------------------------------

library(readr)
library(haven)
library(dplyr)
library(lubridate)
library(reshape2)

# rm(list=ls())

stats_path <- "//stats/cl-out/HPS/Covid19/case_control"
westdata_path <- "//westdata01/dept/HPS/private/Statistics/StatsSupport/Coronavirus"
shared_path <- "//stats/HPS/Covid19/Case_control/data"

filedate <- "2020-06-08" # date input files used for matching

# ### working file --------------------------------------------------------

cc <- readRDS(paste0(stats_path, "/CC_working_", filedate, ".rds"))

# ### RAPID ---------------------------------------------------------------
# admissions in 28 days of test and hospitalised when tested
# NB this data is dynamic, updates overwrite the previous raw file but always go back to feb - use latest submission

rapid <- read_csv(paste0(westdata_path, "/Hospitalisations/RAPID_admissions.csv"), trim_ws = T)
rapid <- as.data.frame(rapid)

rapid <- rename(rapid, CHI=`CHI Number [C]`, Admission.Date=`Admission Date`, Discharge.Date=`Discharge Date`)
table(is.na(rapid$CHI))
rapid <- subset(rapid, !is.na(CHI))

min(nchar(rapid$CHI))
rapid$CHI <- ifelse((nchar(rapid$CHI) == 9) & !is.na(rapid$CHI), paste0("0", rapid$CHI), rapid$CHI) # lpad with 0

tp <- cc[c("upi", "ANON_ID", "SPECIMENDATE")]
rapid <- merge(x=rapid, y=tp, by.x="CHI", by.y="upi", all.x=T)
rapid <- subset(rapid, !is.na(SPECIMENDATE))
rm(tp)

rapid$Admission.Date <- as.Date(substr(rapid$Admission.Date, 1, 10), "%Y/%m/%d")
rapid$Discharge.Date <- as.Date(substr(rapid$Discharge.Date, 1, 10), "%Y/%m/%d")
rapid$days <- as.numeric(rapid$Admission.Date - rapid$SPECIMENDATE)
rapid$adm28 <- ifelse(rapid$days %in% c(0:28), 1, 0)
rapid$inhosp <- ifelse(rapid$SPECIMENDATE >= rapid$Admission.Date &
    (rapid$SPECIMENDATE <= rapid$Discharge.Date | is.na(rapid$Discharge.Date)), 1, 0)

rapid$no.disdate <- ifelse(is.na(rapid$Discharge.Date), 1, 0)
rapid <- rapid[order(rapid$CHI, rapid$Admission.Date, rapid$no.disdate, rapid$Discharge.Date),]

rapid.agg <- rapid %>%
  group_by(CHI, SPECIMENDATE) %>%
  summarise(RAPID=1, inhosp=max(inhosp), adm28=max(adm28)) %>%
  ungroup
cc <- merge(x=cc, y=rapid.agg, by.x=c("upi", "SPECIMENDATE"), by.y=c("CHI", "SPECIMENDATE"), all.x=T)
cc$RAPID[is.na(cc$RAPID)] <- 0
cc$inhosp[is.na(cc$inhosp)] <- 0
cc$adm28[is.na(cc$adm28)] <- 0

## admission and discharge date long format files

# match anon_id to rapid

rapid.long <- rapid[c("ANON_ID", "Admission.Date", "Discharge.Date")]
rapid.long <- rapid.long[order(rapid.long$ANON_ID, rapid.long$Admission.Date),]

saveRDS(rapid.long, paste0(stats_path, "/CC_RAPID_ANON_", filedate, ".rds"))

rm(rapid, rapid.agg, rapid.long)

# ### nosocomial infection -----------------------------------------------------------------
# this is derived from RAPID and provided by Chris R

nos <- read_csv(paste0(westdata_path, "/Hospitalisations/COVID19_nosocomial/analysis/outputs/CHIs_HAIs_ECDC.csv"), trim_ws = T)
nos$HAI <- 1

nos <- nos[c("chi_number", "specimen_date", "HAI", "HAI_category")]

cc <- merge(x=cc, y=nos, by.x=c("upi", "SPECIMENDATE"), by.y=c("chi_number", "specimen_date"), all.x=T)
cc$HAI[is.na(cc$HAI)] <- 0

nos.uk <- read_csv(paste0(westdata_path, "/Hospitalisations/COVID19_nosocomial/analysis/outputs/CHIs_HAIs_UK.csv"), trim_ws = T)
nos.uk$HAI_UK <- 1

nos.uk <- nos.uk[c("chi_number", "specimen_date", "HAI_UK", "HAI_category_UK")]

cc <- merge(x=cc, y=nos.uk, by.x=c("upi", "SPECIMENDATE"), by.y=c("chi_number", "specimen_date"), all.x=T)
cc$HAI_UK[is.na(cc$HAI_UK)] <- 0

rm(nos, nos.uk)

# ### SICSAG ICU ----------------------------------------------------------

icu <- read_spss("//freddy/dept/PHIBCS/PHI/HPS-OA/Collaborative Working/2019-nCoV/01 Surveillance and investigation/24 ICU/CHI Linkage/2020-06-08 Full episode level extract - Input.sav")
icu <- as.data.frame(icu)

icu <- rename(icu, CHI=ChiNo)
table(is.na(icu$CHI) | icu$CHI=="")
icu <- subset(icu, CHI != "")

min(nchar(icu$CHI))
icu$CHI <- as.character(icu$CHI)
icu$CHI <- ifelse((nchar(icu$CHI) == 9) & !is.na(icu$CHI), paste0("0", icu$CHI), icu$CHI) # lpad with 0

tp <- cc[c("upi", "SPECIMENDATE")]
icu <- merge(x=icu, y=tp, by.x="CHI", by.y="upi", all.x=T)
icu <- subset(icu, !is.na(SPECIMENDATE) & SPECIMENDATE <= AdmitUnit & !is.na(covidICUorHDU))

summary(icu$AdmitUnit)
icu$icu <- ifelse(icu$covidICUorHDU ==1, 1, 0)
icu$hdu <- ifelse(icu$covidICUorHDU ==3, 1, 0)

icu.agg <- icu %>%
  group_by(CHI, SPECIMENDATE) %>%
  summarise(icu=max(icu), hdu=max(hdu)) %>%
  ungroup

cc <- merge(x=cc, y=icu.agg, by.x=c("upi", "SPECIMENDATE"), by.y=c("CHI", "SPECIMENDATE"), all.x=T)
cc$icu[is.na(cc$icu)] <- 0
cc$hdu[is.na(cc$hdu)] <- 0

rm(tp, icu, icu.agg)

# ### diabetes register ---------------------------------------------------
# this data is static - each person represented once

diab <- read_csv(paste0(westdata_path, "/Diabetes/Diabetes_type.csv"), trim_ws = T)
diab <- as.data.frame(diab)
diab <- subset(diab, !is.na(CHI) & CHI !="")

min(nchar(diab$CHI))
diab$CHI <- as.character(diab$CHI)
diab$CHI <- ifelse((nchar(diab$CHI) == 9) & !is.na(diab$CHI), paste0("0", diab$CHI), diab$CHI) # lpad with 0

diab <- diab %>% group_by(CHI, dm.type) %>% summarise(diab.reg=1) %>% ungroup

cc <- merge(cc, diab, by.x="upi", by.y="CHI", all.x=T)
cc$diab.reg[is.na(cc$diab.reg)] <- 0


# ### save working file ---------------------------------------------------

saveRDS(cc, paste0(stats_path, "/CC_working_", filedate, ".rds"))
# cc <- readRDS(paste0(stats_path, "/CC_working_", filedate, ".rds"))

# ### save out final anonymised files ----------------------------------------------

cc.anon <- cc[c("ANON_ID", "PATID", "stratum", "is.case", "SPECIMENDATE", "SEX", "AgeYear",
  "INSTITUTION_CODE", "Prisoner_Flag", "DATE_OF_DEATH", "DATE_OF_REGISTRATION", "UNDERLYING_CAUSE_OF_DEATH",
  "CAUSE_OF_DEATH_CODE_0", "CAUSE_OF_DEATH_CODE_1", "CAUSE_OF_DEATH_CODE_2", "CAUSE_OF_DEATH_CODE_3",
  "CAUSE_OF_DEATH_CODE_4", "CAUSE_OF_DEATH_CODE_5", "CAUSE_OF_DEATH_CODE_6", "CAUSE_OF_DEATH_CODE_7",
  "CAUSE_OF_DEATH_CODE_8", "CAUSE_OF_DEATH_CODE_9", "covid_ucod", "covid_cod", "dead28",
  "ETHNIC_SMR_LAST", "can.reg", "CAREHOME", "NURSINGHOME", "simd2020_sc_quintile", "simd", "hb2019name",
  "RAPID", "inhosp", "adm28", "HAI", "HAI_category", "HAI_UK", "HAI_category_UK",
  "icu", "hdu", "dm.type", "diab.reg")]
cc.anon <- cc.anon[order(cc.anon$stratum, -cc.anon$is.case, cc.anon$ANON_ID),]

saveRDS(cc, paste0(stats_path, "/CC_linked_final_", filedate, ".rds"))
saveRDS(cc.anon, paste0(shared_path, "/CC_linked_ANON_", filedate, ".rds"))


# long format files:
diag.25 <- readRDS(paste0(stats_path, "/CC_SMR01_ICD10_25_", filedate, ".rds"))
diag.25 <- diag.25[c("ANON_ID", "ICD10")]
saveRDS(diag.25, paste0(shared_path, "/CC_SMR01_ICD10_25_ANON_", filedate, ".rds"))

diag.x25 <- readRDS(paste0(stats_path, "/CC_SMR01_ICD10_x25_", filedate, ".rds"))
diag.x25 <- diag.x25[c("ANON_ID", "ICD10")]
saveRDS(diag.x25, paste0(shared_path, "/CC_SMR01_ICD10_x25_ANON_", filedate, ".rds"))

op.x25 <- readRDS(paste0(stats_path, "/CC_SMR01_OPCS4_MAIN.x25_", filedate, ".rds"))
op.x25 <- op.x25[c("ANON_ID", "MAIN_OPERATION")]
saveRDS(op.x25, paste0(shared_path, "/CC_SMR01_OPCS4_MAIN.x25_ANON_", filedate, ".rds"))

op.25 <- readRDS(paste0(stats_path, "/CC_SMR01_OPCS4_MAIN.25_", filedate, ".rds"))
op.25 <- op.25[c("ANON_ID", "MAIN_OPERATION")]
saveRDS(op.25, paste0(shared_path, "/CC_SMR01_OPCS4_MAIN.25_ANON_", filedate, ".rds"))

rapid.long <- readRDS(paste0(stats_path, "/CC_RAPID_ANON_", filedate, ".rds"))
rapid.long <- rapid.long[c("ANON_ID", "Admission.Date", "Discharge.Date")]
saveRDS(rapid.long, paste0(shared_path, "/CC_RAPID_ANON_", filedate, ".rds"))


#PIS.x15 <- readRDS(paste0(stats_path, "/CC_PIS_x15_ANON_", filedate, ".rds"))
#PIS.x15 <- select(PIS.x15, -c("CHI))
#saveRDS(PIS.x15, paste0(shared_path, "/CC_PIS_x15_ANON_", filedate, ".rds"))

#PIS.15 <- readRDS(paste0(stats_path, "/CC_PIS_15_ANON_", filedate, ".rds"))
#PIS.15 <- select(PIS.x15, -c("CHI))
#saveRDS(PIS.15, paste0(shared_path, "/CC_PIS_15_ANON_", filedate, ".rds"))




