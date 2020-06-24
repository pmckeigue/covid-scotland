### Covid19 - additional covid cases for CHILI matching
# NRS weekly+monthly deaths due to covid where no ECOSS test positive
# R studio Server Pro

library(odbc)
library(dplyr)
library(lubridate)
library(reshape2)

rm(list=ls())

clout_path <- "/conf/linkage/output/HPS/Covid19"

# ### weekly coded deaths -----------------------------------------

SMRAconnection <- suppressWarnings(dbConnect(odbc(),
                                             dsn="SMRA",
                                             uid=.rs.askForPassword("Enter Username:"),
                                             pwd=.rs.askForPassword("Enter Password:") ))

# odbcListObjects(SMRAconnection, schema="ANALYSIS")
# odbcPreviewObject(SMRAconnection, table="ANALYSIS.GRO_DEATHS_WEEKLY_C", rowLimit=0)

NRS <- tbl_df(dbGetQuery(SMRAconnection, statement=
                           'SELECT
                         T1.CHI AS NRS_CHI, T1.DERIVED_CHI, T1.DATE_OF_BIRTH, T1.SEX, T1.DATE_OF_DEATH, T1.DATE_OF_REGISTRATION, 
                         T1.UNDERLYING_CAUSE_OF_DEATH, T1.CAUSE_OF_DEATH_CODE_0, T1.CAUSE_OF_DEATH_CODE_1,
                         T1.CAUSE_OF_DEATH_CODE_2, T1.CAUSE_OF_DEATH_CODE_3, T1.CAUSE_OF_DEATH_CODE_4, T1.CAUSE_OF_DEATH_CODE_5,
                         T1.CAUSE_OF_DEATH_CODE_6, T1.CAUSE_OF_DEATH_CODE_7, T1.CAUSE_OF_DEATH_CODE_8,
                         T1.CAUSE_OF_DEATH_CODE_9
                         FROM ANALYSIS.GRO_DEATHS_WEEKLY_C T1'
))


# POSIXct to date
cols <- c("DATE_OF_BIRTH", "DATE_OF_DEATH", "DATE_OF_REGISTRATION")

for (i in cols) {
  NRS[i] <- lapply(NRS[i], as.Date, "%Y-%m-%d", tz="GMT")
}

summary(NRS$DATE_OF_REGISTRATION)

# fix derived_CHI for matching

NRS$DERIVED_CHI_FIX <- gsub("\\D", "", NRS$DERIVED_CHI)
NRS$DERIVED_CHI_FIX <- ifelse(nchar(NRS$DERIVED_CHI_FIX) ==9, paste0("0", NRS$DERIVED_CHI_FIX), NRS$DERIVED_CHI_FIX)
NRS$DERIVED_CHI_FIX[NRS$DERIVED_CHI_FIX==""] <- NA
summary(nchar(NRS$DERIVED_CHI_FIX))

NRS$DERIVED_CHI <- NULL

NRS <- melt(NRS,
            measure.vars=c("NRS_CHI", "DERIVED_CHI_FIX"),
            value.name="CHI",
            na.rm = FALSE )

NRS <- subset(NRS, !is.na(CHI) & CHI !="")
NRS <- NRS %>% filter(!is.na(CHI) & CHI !="") %>% group_by(CHI, DATE_OF_BIRTH, SEX, DATE_OF_REGISTRATION, DATE_OF_DEATH, 
                                                           UNDERLYING_CAUSE_OF_DEATH, CAUSE_OF_DEATH_CODE_0, 
                                                           CAUSE_OF_DEATH_CODE_1, CAUSE_OF_DEATH_CODE_2, 
                                                           CAUSE_OF_DEATH_CODE_3, CAUSE_OF_DEATH_CODE_4,
                                                           CAUSE_OF_DEATH_CODE_5, CAUSE_OF_DEATH_CODE_6, 
                                                           CAUSE_OF_DEATH_CODE_7, CAUSE_OF_DEATH_CODE_8, 
                                                           CAUSE_OF_DEATH_CODE_9) %>% summarise() %>% ungroup

length(unique(NRS$CHI)) 
nrow(NRS)
# will be same unless people have >1 death record

NRS.wk <- NRS
rm(NRS)

# ### previous coded deaths -----------------------------------------------

# odbcPreviewObject(SMRAconnection, table="ANALYSIS.GRO_DEATHS_C", rowLimit=0)

NRS <- tbl_df(dbGetQuery(SMRAconnection, statement=
                           "SELECT
                         T1.UPI_NUMBER, T1.DATE_OF_BIRTH, T1.SEX, T1.DATE_OF_DEATH, T1.DATE_OF_REGISTRATION,
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

# ### combine and merge to cc ---------------------------------------------

NRS.all <- rbind(NRS.wk, NRS)
NRS.all <- NRS.all[order(NRS.all$CHI, NRS.all$DATE_OF_DEATH),]
NRS.all <- NRS.all %>% group_by(CHI) %>% summarise(
  DATE_OF_BIRTH = first(DATE_OF_BIRTH),
  SEX = first(SEX),
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
  CAUSE_OF_DEATH_CODE_9 = first(CAUSE_OF_DEATH_CODE_9)) %>% ungroup

# ### all covid deaths U071 U072 ----------------------------------------------

# underlying cause
NRS.all$covid_ucod <- ifelse(NRS.all$UNDERLYING_CAUSE_OF_DEATH %in% c("U071", "U072"), 1, 0)

# any cause
vars <- NRS.all[c("UNDERLYING_CAUSE_OF_DEATH", "CAUSE_OF_DEATH_CODE_0", "CAUSE_OF_DEATH_CODE_1", "CAUSE_OF_DEATH_CODE_2",
             "CAUSE_OF_DEATH_CODE_3", "CAUSE_OF_DEATH_CODE_4", "CAUSE_OF_DEATH_CODE_5", "CAUSE_OF_DEATH_CODE_6",
             "CAUSE_OF_DEATH_CODE_7", "CAUSE_OF_DEATH_CODE_8", "CAUSE_OF_DEATH_CODE_9")]

NRS.all$covid_cod <- 0
for (i in vars) {
  NRS.all$covid_cod <- ifelse(!is.na(i) & (substr(i, 1, 3) %in% c("U07")), 1, NRS.all$covid_cod)
}

NRS.sub <- subset(NRS.all, covid_cod==1)
#saveRDS(NRS.sub, paste0(clout_path, "/NRS_Covid_Coded_Deaths_all.rds"))
# df <- readRDS(paste0(clout_path, "/NRS_Covid_coded_Deaths.rds"))


# ### additional covid coded cases for CHILI team -------------------------

ecoss <- readRDS(paste0(clout_path, "/ECOSS_deduped.rds"))

df.all <- merge(x=ecoss, y=NRS.sub, by="CHI", all=T)
table(!is.na(df.all$DATE_OF_DEATH))
with(subset(df.all, !is.na(DATE_OF_DEATH)), table(NCOV_RESULT))
with(subset(df.all, !is.na(DATE_OF_DEATH)), table(is.na(NCOV_RESULT)))

df.sub <- subset(df.all, !is.na(DATE_OF_DEATH) & (NCOV_RESULT=="Negative" | is.na(NCOV_RESULT)))

# dummy test date

df.sub$SpecimenDate.dummy <- as.Date(df.sub$DATE_OF_DEATH - days(14), "%Y-%m-%d")
df.sub$age <- trunc(as.numeric(df.sub$SpecimenDate.dummy - df.sub$DATE_OF_BIRTH)/365.25)

df.sub <- df.sub[c("CHI", "age", "SEX", "SpecimenDate.dummy", "DATE_OF_DEATH", "DATE_OF_REGISTRATION", "covid_cod", "covid_ucod")]

saveRDS(df.sub, paste0(clout_path, "/Additional_NRS_Covid_Cases_", Sys.Date(), ".rds"))
# df <- readRDS(paste0(clout_path, "/Additional_NRS_Covid_Cases.rds"))
