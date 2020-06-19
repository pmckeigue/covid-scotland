# Covid-19 case control linkage
# HPS stats support team
# r server

# ### set up --------------------------------------------------------------

library(readr)
library(dplyr)
library(odbc)
library(lubridate)
library(reshape2)

# rm(list=ls())

chili_path <- "/conf/linkage/output/y2k_cat_check/conf"
stats_path <- "//conf/linkage/output/HPS/Covid19/case_control"
lookups_path <- "//conf/linkage/output/lookups/Unicode"

filedate <- "2020-06-08" # date input files used for matching 

# ### returns from case control matching -----------------------------------

cc <- read_csv(paste0(chili_path, "/case_and_controls_June_2.csv"), trim_ws = TRUE)
cc.raw <- cc #copy

nrow(cc)
table(cc$is.case)
length(unique(cc$stratum))
with(subset(cc, is.case==TRUE), length(unique(stratum)))
with(subset(cc, is.case==FALSE), length(unique(stratum)))
length(unique(cc$upi))
with(subset(cc, is.case==TRUE), length(unique(upi)))
with(subset(cc, is.case==FALSE), length(unique(upi)))

# remove duplicate controls within same stratum and controls that are also a case and have matched to the same stratum

cc <- cc %>% 
  arrange(stratum, -is.case, upi) %>%
  group_by(stratum, upi) %>% 
  summarise(is.case=first(is.case),
            SPECIMENDATE=first(SPECIMENDATE),
            SEX=first(SEX),
            AgeYear=first(AgeYear),
            CURRENT_POSTCODE=first(CURRENT_POSTCODE),
            GP_PRAC_NO=first(GP_PRAC_NO),
            INSTITUTION_CODE=first(INSTITUTION_CODE),
            Prisoner_Flag=first(Prisoner_Flag)) %>%
  ungroup %>%
  arrange(stratum, -is.case, upi) 

nrow(cc)
table(cc$is.case)
length(unique(cc$stratum))
with(subset(cc, is.case==TRUE), length(unique(stratum)))
with(subset(cc, is.case==FALSE), length(unique(stratum)))
length(unique(cc$upi))
with(subset(cc, is.case==TRUE), length(unique(upi)))
with(subset(cc, is.case==FALSE), length(unique(upi)))

# remove stratum where case has matched only theirself (n=2)

cc <- cc %>% group_by(stratum) %>% mutate(num_upi=n_distinct(upi)) %>%  ungroup
cc <- subset(cc, num_upi >1)

# assign control specimen date based on stratum case specimen date

cc$SPECIMENDATE[cc$is.case==FALSE] <- NA
cc <- cc %>% 
  arrange(stratum, -is.case, upi) %>% 
  group_by(stratum) %>% 
  tidyr::fill(SPECIMENDATE, .direction = "down") %>%
  ungroup %>% 
  arrange(stratum, -is.case, upi)

# baseline case control file

nrow(cc)
table(cc$is.case)
with(subset(cc, is.case==TRUE), length(unique(upi)))
length(unique(cc$stratum))
with(subset(cc, is.case==TRUE), length(unique(stratum)))
with(subset(cc, is.case==FALSE), length(unique(stratum)))
table(is.na(cc$stratum))
table(is.na(cc$upi) | cc$upi=="")
table(is.na(cc$SPECIMENDATE))

# add IDs

cc$ANON_ID <- seq(from = 1, to = nrow(cc), by=1)
cc <- cc %>% mutate(PATID = match(upi, unique(upi)))


# save base file

cc <- cc[c("ANON_ID", "PATID", "stratum", "upi", "is.case", "SPECIMENDATE", "SEX", "AgeYear", "CURRENT_POSTCODE",
           "GP_PRAC_NO", "INSTITUTION_CODE", "Prisoner_Flag")]

saveRDS(cc, paste0(stats_path, "/CC_base_file_", filedate, ".rds"))


# ### connect to SMRA -----------------------------------------------------

SMRAconnection <- suppressWarnings(dbConnect(odbc(),
                                             dsn="SMRA",
                                             uid=.rs.askForPassword("Enter Username:"),
                                             pwd=.rs.askForPassword("Enter Password:") ))


# ### NRS deaths ----------------------------------------------------------
# weekly deaths added to monthly deaths before matching

# odbcListObjects(SMRAconnection, schema="ANALYSIS")
# odbcPreviewObject(SMRAconnection, table="ANALYSIS.GRO_DEATHS_WEEKLY_C", rowLimit=0)

NRS.wk <- tbl_df(dbGetQuery(SMRAconnection, statement=
                              'SELECT
                            T1.CHI AS NRS_CHI, T1.DERIVED_CHI, T1.DATE_OF_DEATH, T1.DATE_OF_REGISTRATION, 
                            T1.UNDERLYING_CAUSE_OF_DEATH, T1.CAUSE_OF_DEATH_CODE_0, T1.CAUSE_OF_DEATH_CODE_1,
                            T1.CAUSE_OF_DEATH_CODE_2, T1.CAUSE_OF_DEATH_CODE_3, T1.CAUSE_OF_DEATH_CODE_4, T1.CAUSE_OF_DEATH_CODE_5,
                            T1.CAUSE_OF_DEATH_CODE_6, T1.CAUSE_OF_DEATH_CODE_7, T1.CAUSE_OF_DEATH_CODE_8, T1.CAUSE_OF_DEATH_CODE_9
                            FROM ANALYSIS.GRO_DEATHS_WEEKLY_C T1' ))

# POSIXct to date
cols <- c("DATE_OF_DEATH", "DATE_OF_REGISTRATION")

for (i in cols) {
  NRS.wk[i] <- lapply(NRS.wk[i], as.Date, "%Y-%m-%d", tz="GMT")
}

summary(NRS.wk$DATE_OF_REGISTRATION)

# fix derived_CHI for matching

NRS.wk$DERIVED_CHI_FIX <- gsub("\\D", "", NRS.wk$DERIVED_CHI)
NRS.wk$DERIVED_CHI_FIX <- ifelse(nchar(NRS.wk$DERIVED_CHI_FIX) ==9, paste0("0", NRS.wk$DERIVED_CHI_FIX), NRS.wk$DERIVED_CHI_FIX)
NRS.wk$DERIVED_CHI_FIX[NRS.wk$DERIVED_CHI_FIX==""] <- NA
summary(nchar(NRS.wk$DERIVED_CHI_FIX))

NRS.wk$DERIVED_CHI <- NULL

NRS.wk <- melt(NRS.wk,
               measure.vars=c("NRS.wk_CHI", "DERIVED_CHI_FIX"),
               value.name="CHI",
               na.rm = FALSE )

NRS.wk <- subset(NRS.wk, !is.na(CHI) & CHI !="")
NRS.wk <- NRS.wk %>% filter(!is.na(CHI) & CHI !="") %>% group_by(CHI, DATE_OF_REGISTRATION, DATE_OF_DEATH, 
                                                                 UNDERLYING_CAUSE_OF_DEATH, CAUSE_OF_DEATH_CODE_0, 
                                                                 CAUSE_OF_DEATH_CODE_1, CAUSE_OF_DEATH_CODE_2, 
                                                                 CAUSE_OF_DEATH_CODE_3, CAUSE_OF_DEATH_CODE_4,
                                                                 CAUSE_OF_DEATH_CODE_5, CAUSE_OF_DEATH_CODE_6, 
                                                                 CAUSE_OF_DEATH_CODE_7, CAUSE_OF_DEATH_CODE_8, 
                                                                 CAUSE_OF_DEATH_CODE_9) %>% summarise() %>% ungroup

length(unique(NRS.wk$CHI)) 
nrow(NRS.wk)
# will be same unless people have >1 death record


## previous coded deaths

# odbcPreviewObject(SMRAconnection, table="ANALYSIS.GRO_DEATHS_C", rowLimit=0)

NRS <- tbl_df(dbGetQuery(SMRAconnection, statement=
                           "SELECT
                         T1.UPI_NUMBER, T1.DATE_OF_DEATH, T1.DATE_OF_REGISTRATION,
                         T1.UNDERLYING_CAUSE_OF_DEATH, T1.CAUSE_OF_DEATH_CODE_0, T1.CAUSE_OF_DEATH_CODE_1,
                         T1.CAUSE_OF_DEATH_CODE_2, T1.CAUSE_OF_DEATH_CODE_3, T1.CAUSE_OF_DEATH_CODE_4, T1.CAUSE_OF_DEATH_CODE_5,
                         T1.CAUSE_OF_DEATH_CODE_6, T1.CAUSE_OF_DEATH_CODE_7, T1.CAUSE_OF_DEATH_CODE_8,
                         T1.CAUSE_OF_DEATH_CODE_9
                         FROM ANALYSIS.GRO_DEATHS_C T1
                         WHERE DATE_OF_DEATH >= '01 January 2020' "  ))

NRS <- subset(NRS, !is.na(UPI_NUMBER) & UPI_NUMBER !="")

# POSIXct to date
cols <- c("DATE_OF_DEATH", "DATE_OF_REGISTRATION")

for (i in cols) {
  NRS[i] <- lapply(NRS[i], as.Date, "%Y-%m-%d", tz="GMT")
}

summary(NRS$DATE_OF_REGISTRATION)

names(NRS)
NRS <- rename(NRS, CHI=UPI_NUMBER)

length(unique(NRS$CHI)) 
nrow(NRS)
# will be same unless people have >1 death record

## combine and merge to cc 
# NRS.all <- NRS
NRS.all <- rbind(NRS.wk, NRS)
NRS.all <- NRS.all %>% 
  arrange(CHI, DATE_OF_DEATH) %>%
  group_by(CHI) %>% 
  summarise(
    DATE_OF_DEATH = first(DATE_OF_DEATH),
    DATE_OF_REGISTRATION = first(DATE_OF_REGISTRATION),
    UNDERLYING_CAUSE_OF_DEATH = first(UNDERLYING_CAUSE_OF_DEATH),
    CAUSE_OF_DEATH_CODE_0 = first(CAUSE_OF_DEATH_CODE_0),
    CAUSE_OF_DEATH_CODE_1 = first(CAUSE_OF_DEATH_CODE_1),
    CAUSE_OF_DEATH_CODE_2 = first(CAUSE_OF_DEATH_CODE_2),
    CAUSE_OF_DEATH_CODE_3 = first(CAUSE_OF_DEATH_CODE_3),
    CAUSE_OF_DEATH_CODE_4 = first(CAUSE_OF_DEATH_CODE_4),
    CAUSE_OF_DEATH_CODE_5 = first(CAUSE_OF_DEATH_CODE_5),
    CAUSE_OF_DEATH_CODE_6 = first(CAUSE_OF_DEATH_CODE_6),
    CAUSE_OF_DEATH_CODE_7 = first(CAUSE_OF_DEATH_CODE_7),
    CAUSE_OF_DEATH_CODE_8 = first(CAUSE_OF_DEATH_CODE_8),
    CAUSE_OF_DEATH_CODE_9 = first(CAUSE_OF_DEATH_CODE_9)) %>% 
  ungroup

cc <- merge(cc, NRS.all, by.x="upi", by.y="CHI", all.x=T)

## derived death fields - covid deaths U071 U072

# underlying cause
cc$covid_ucod <- ifelse(cc$UNDERLYING_CAUSE_OF_DEATH %in% c("U071", "U072"), 1, 0)

# any cause
vars <- cc[c("UNDERLYING_CAUSE_OF_DEATH", "CAUSE_OF_DEATH_CODE_0", "CAUSE_OF_DEATH_CODE_1", "CAUSE_OF_DEATH_CODE_2",
             "CAUSE_OF_DEATH_CODE_3", "CAUSE_OF_DEATH_CODE_4", "CAUSE_OF_DEATH_CODE_5", "CAUSE_OF_DEATH_CODE_6",
             "CAUSE_OF_DEATH_CODE_7", "CAUSE_OF_DEATH_CODE_8", "CAUSE_OF_DEATH_CODE_9")]

cc$covid_cod <- 0
for (i in vars) {
  cc$covid_cod <- ifelse(!is.na(i) & (substr(i, 1, 3) %in% c("U07")), 1, cc$covid_cod)
}


# dead in 28 days
cc$dead28 <- ifelse(!is.na(cc$DATE_OF_DEATH) & cc$DATE_OF_DEATH >= cc$SPECIMENDATE & as.numeric(cc$DATE_OF_DEATH - cc$SPECIMENDATE) <=28, 1, 0)


# check for and remove stratum where case death is before specimen date - likely recording error or bad link
cc$death.b4.test <- ifelse(cc$DATE_OF_DEATH < cc$SPECIMENDATE, 1, 0)
table(cc$death.b4.test, cc$is.case)
ck <- subset(cc, death.b4.test==1 & is.case==TRUE) # (n=14)
summary(ck$DATE_OF_DEATH)
cc <- subset(cc, !(stratum %in% c(ck$stratum)))

# check for and remove remaining controls with death before matched case specimen date
with(subset(cc, death.b4.test==1 & is.case==FALSE), length(upi)) # (n=1002)
with(subset(cc, death.b4.test==1 & is.case==FALSE), summary(DATE_OF_DEATH)) 

cc <- subset(cc, death.b4.test ==0 | is.na(death.b4.test))


# updated baseline case control file

nrow(cc)
table(cc$is.case)
with(subset(cc, is.case==TRUE), length(unique(upi)))
length(unique(cc$stratum))
with(subset(cc, is.case==TRUE), length(unique(stratum)))
with(subset(cc, is.case==FALSE), length(unique(stratum)))
ck <- cc %>% filter(is.case==FALSE) %>% group_by(stratum) %>% summarise(num_controls=n()) %>% ungroup
table(ck$num_controls)

# saveRDS(cc, paste0(stats_path, "/CC_base_file_", filedate, ".rds"))

rm(NRS.wk, NRS, NRS.all, ck, vars)


# ### SMR ethnicity -------------------------------------------------------

### list of CHIs and period for uploading to SMRA

tp <- cc

# lookback period start date - must be numeric format
tp$RFSTART <- as.numeric(format(tp$SPECIMENDATE - years(10), "%Y%m%d"))
tp$SPECDATE <- as.numeric(format(tp$SPECIMENDATE, "%Y%m%d"))
tp$CHI <- tp$upi

# temp file for uploading - character or numeric formats only
tp <- tp[c("CHI", "ANON_ID", "SPECDATE", "RFSTART")]
str(tp)

# tp <- tp[1:10, ]

### upload CHIs and extract all episodes in period

dbWriteTable(SMRAconnection, "TEMP", tp, overwrite=T)
# TEMP <- dbGetQuery(SMRAconnection, statement='SELECT * FROM SHAROK01."TEMP"')
# View(TEMP)
# rm(TEMP)

SMR01 <- tbl_df(dbGetQuery(SMRAconnection, statement=
                             'SELECT
                           T0.CHI, T0.ANON_ID, T0.SPECDATE, T0.RFSTART,
                           T1.UPI_NUMBER, T1.ETHNIC_GROUP, T1.DISCHARGE_DATE
                           FROM SHAROK01."TEMP" T0, ANALYSIS.SMR01_PI T1
                           WHERE T0.CHI = T1.UPI_NUMBER (+)
                           AND to_number(to_char(T1.ADMISSION_DATE,\'YYYYMMDD\')) >= T0.RFSTART 
                           AND to_number(to_char(T1.ADMISSION_DATE,\'YYYYMMDD\')) <= T0.SPECDATE  
                           ORDER BY T1.UPI_NUMBER, T1.CIS_MARKER, T1.ADMISSION_DATE, T1.DISCHARGE_DATE, T1.ADMISSION, T1.DISCHARGE, T1.URI'
                           ))

SMR00 <- tbl_df(dbGetQuery(SMRAconnection, statement=
                             'SELECT
                           T0.CHI, T0.ANON_ID, T0.SPECDATE, T0.RFSTART,
                           T1.UPI_NUMBER, T1.ETHNIC_GROUP, T1.CLINIC_DATE
                           FROM SHAROK01."TEMP" T0, ANALYSIS.SMR00_PI T1
                           WHERE T0.CHI = T1.UPI_NUMBER (+)
                           AND to_number(to_char(T1.CLINIC_DATE,\'YYYYMMDD\')) >= T0.RFSTART
                           AND to_number(to_char(T1.CLINIC_DATE,\'YYYYMMDD\')) <= T0.SPECDATE
                           ORDER BY T0.CHI, T0.SPECDATE, T1.CLINIC_DATE'  
                           ))

dbRemoveTable(SMRAconnection, "TEMP")


### format dates

# POSIXct to date
SMR01$DISCHARGE_DATE <- as.Date(SMR01$DISCHARGE_DATE, "%Y-%m-%d", tz="GMT")
SMR00$CLINIC_DATE <- as.Date(SMR00$CLINIC_DATE, "%Y-%m-%d", tz="GMT")

SMR01 <- rename(SMR01, DATE=DISCHARGE_DATE)
SMR00 <- rename(SMR00, DATE=CLINIC_DATE)

SMR <- rbind(SMR01, SMR00)

# ethnic group - last available 
SMR$ETHNIC_GROUP[SMR$ETHNIC_GROUP %in% c("00", "97", "98", "99")] <- NA
SMR <- subset(SMR, !is.na(ETHNIC_GROUP))
SMR <- SMR[order(SMR$CHI, SMR$DATE),]

SMR <- SMR %>% 
  group_by(CHI) %>% 
  summarise(ETHNIC_SMR_LAST=last(ETHNIC_GROUP)) %>% 
  ungroup

# merge to linked file
cc <- merge(cc, SMR, by.x="upi", by.y="CHI", all.x=T)
rm(tp, SMR01, SMR00, SMR)


# ### SMR01 ---------------------------------------------------------------
# long format SMR01 diagnosis files

# list of CHIs and period for uploading to SMRA

tp <- cc

# lookback period start date - must be numeric format
tp$RFSTART <- as.numeric(format(tp$SPECIMENDATE - years(5), "%Y%m%d"))
tp$SPECDATE <- as.numeric(format(tp$SPECIMENDATE, "%Y%m%d"))
tp$CHI <- tp$upi

# temp file for uploading - character or numeric formats only
tp <- tp[c("CHI", "ANON_ID", "SPECDATE", "RFSTART")]
str(tp)

# upload CHIs and extract all episodes in period

dbWriteTable(SMRAconnection, "TEMP", tp, overwrite=T)

SMR01 <- tbl_df(dbGetQuery(SMRAconnection, statement=
  'SELECT
  T0.CHI, T0.ANON_ID, T0.SPECDATE, T0.RFSTART,
  T1.UPI_NUMBER, T1.ADMISSION_DATE, T1.DISCHARGE_DATE, 
  T1.MAIN_CONDITION, T1.OTHER_CONDITION_1, T1.OTHER_CONDITION_2, T1.OTHER_CONDITION_3, T1.OTHER_CONDITION_4, 
  T1.OTHER_CONDITION_5, T1.MAIN_OPERATION
  FROM SHAROK01."TEMP" T0, ANALYSIS.SMR01_PI T1
  WHERE T0.CHI = T1.UPI_NUMBER (+)
  AND to_number(to_char(T1.DISCHARGE_DATE,\'YYYYMMDD\')) >= T0.RFSTART
  AND to_number(to_char(T1.DISCHARGE_DATE,\'YYYYMMDD\')) < T0.SPECDATE
  ORDER BY T1.UPI_NUMBER, T1.CIS_MARKER, T1.ADMISSION_DATE, T1.DISCHARGE_DATE, T1.ADMISSION, T1.DISCHARGE, T1.URI' 
  ))

# format dates
SMR01$SPECDATE <- as.Date(as.character(SMR01$SPECDATE), "%Y%m%d")
SMR01$DISCHARGE_DATE <- as.Date(SMR01$DISCHARGE_DATE, "%Y-%m-%d", tz="GMT")

# lookback period 
SMR01$pre25 <- ifelse(SMR01$DISCHARGE_DATE < SMR01$SPECDATE - days(25), 1, 0)

### diagnoses

tp <- subset(SMR01, pre25==1)
tp <- tp[c("ANON_ID", "MAIN_CONDITION", "OTHER_CONDITION_1", "OTHER_CONDITION_2", "OTHER_CONDITION_3", "OTHER_CONDITION_4", "OTHER_CONDITION_5")]
tp <- melt(tp,
           id.vars=c("ANON_ID"),
           variable.name="diag",
           value.name="ICD10",
           na.rm = TRUE )
tp <- tp %>% group_by(ANON_ID, ICD10) %>% summarise() %>% ungroup

saveRDS(tp, paste0(stats_path, "/CC_SMR01_ICD10_x25_", filedate, ".rds"))
rm(tp)

tp <- subset(SMR01, pre25!=1)
tp <- tp[c("ANON_ID", "MAIN_CONDITION", "OTHER_CONDITION_1", "OTHER_CONDITION_2", "OTHER_CONDITION_3", "OTHER_CONDITION_4", "OTHER_CONDITION_5")]
tp <- melt(tp,
           id.vars=c("ANON_ID"),
           variable.name="diag",
           value.name="ICD10",
           na.rm = TRUE )
tp <- tp %>% group_by(ANON_ID, ICD10) %>% summarise() %>% ungroup

saveRDS(tp, paste0(stats_path, "/CC_SMR01_ICD10_25_", filedate, ".rds"))
rm(tp)


### surgeries 

# main op only - keep pairing
tp <- subset(SMR01, pre25==1)
tp <- tp %>% group_by(ANON_ID, MAIN_OPERATION) %>% summarise() %>% ungroup
tp <- subset(tp, !is.na(MAIN_OPERATION))

saveRDS(tp, paste0(stats_path, "/CC_SMR01_OPCS4_MAIN.x25_", filedate, ".rds"))
rm(tp)

tp <- subset(SMR01, pre25!=1)
tp <- tp %>% group_by(ANON_ID, MAIN_OPERATION) %>% summarise() %>% ungroup
tp <- subset(tp, !is.na(MAIN_OPERATION))

saveRDS(tp, paste0(stats_path, "/CC_SMR01_OPCS4_MAIN.25_", filedate, ".rds"))
rm(tp, SMR01)


# ### SMR06 ---------------------------------------------------------------

SMR06 <- tbl_df(dbGetQuery(SMRAconnection, statement=
                             'SELECT
                           T0.CHI, T0.ANON_ID, T0.SPECDATE, T0.RFSTART,
                           T1.UPI_NUMBER, T1.INCIDENCE_DATE, to_number(to_char(T1.INCIDENCE_DATE,\'YYYYMMDD\')) AS INCIDENCE_DATE_NUM,
                           T1.ICD10S_CANCER_SITE
                           FROM SHAROK01."TEMP" T0, ANALYSIS.SMR06_PI T1
                           WHERE T0.CHI = T1.UPI_NUMBER (+)
                           AND to_number(to_char(T1.INCIDENCE_DATE,\'YYYYMMDD\')) >= T0.RFSTART
                           AND to_number(to_char(T1.INCIDENCE_DATE,\'YYYYMMDD\')) <= T0.SPECDATE
                           ORDER BY T0.CHI, T0.SPECDATE, T1.INCIDENCE_DATE'  ))

dbRemoveTable(SMRAconnection, "TEMP")

# format dates

SMR06$SPECDATE <- as.Date(as.character(SMR06$SPECDATE), "%Y%m%d")
SMR06$INCIDENCE_DATE <- as.Date(SMR06$INCIDENCE_DATE, "%Y-%m-%d", tz="GMT")

# long format file

tp <- SMR06[c("ANON_ID", "INCIDENCE_DATE", "ICD10S_CANCER_SITE")]
saveRDS(tp, paste0(stats_path, "/CC_SMR06_ICD10_", filedate, ".rds"))


# any diag
SMR06 <- SMR06 %>% group_by(ANON_ID) %>% summarise(can.reg=1) %>% ungroup

cc <- merge(cc, SMR06, by="ANON_ID", all.x=T)
cc$can.reg[is.na(cc$can.reg)] <- 0

rm(SMR06, tp)

# ### care home/nursing home  ---------------------------------------------------

cc$CAREHOME <- ifelse(cc$INSTITUTION_CODE %in% c(93), 1, 0)
cc$NURSINGHOME <- ifelse(cc$INSTITUTION_CODE %in% c(98), 1, 0)


# ### SIMD ----------------------------------------------------------------

pc <- readRDS(paste0(lookups_path, "/Deprivation/postcode_2020_1_simd2020.rds"))
pc <- pc[c("pc7", "simd2020_sc_quintile")]

cc <- merge(x=cc, y=pc, by.x="CURRENT_POSTCODE", by.y="pc7", all.x=T)

cc$simd2020_sc_quintile[is.na(cc$simd2020_sc_quintile)] <- 9
cc$simd <- factor(cc$simd2020_sc_quintile,
                  levels = c(1, 2, 3, 4, 5, 9),
                  labels = c("1 - most deprived", "2", "3", "4", "5 - least deprived", "Unknown"))
table (cc$simd)
table (is.na(cc$simd))
rm(pc)


# ### HB residence --------------------------------------------------------

pc <- readRDS(paste0(lookups_path, "/Geography/Scottish Postcode Directory/Scottish_Postcode_Directory_2020_1.rds"))
pc <- pc[c("pc7", "hb2019name")]

cc <- merge(x=cc, y=pc, by.x="CURRENT_POSTCODE", by.y="pc7", all.x=T)
table (cc$hb2019name)
table (is.na(cc$hb2019name))
rm(pc)


# ### save temp working file ---------------------------------------------------

saveRDS(cc, paste0(stats_path, "/CC_working_", filedate, ".rds"))
# cc <- readRDS(paste0(stats_path, "/CC_working_", filedate, ".rds"))


