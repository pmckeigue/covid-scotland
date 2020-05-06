### Extract SMR01 and match SIMD
### Covid19 risk factors, comorbidities and ICD10 + ethnicity
# R studio Server Pro

library(odbc)
library(dplyr)
library(lubridate)
library(haven)

rm(list=ls())

dfe <- readRDS(paste0(stats_path, "/HPS/Covid19/case_control/CC_NEW.rds"))
dfe <- rename(dfe, CHI=UPI_NUMBER)


# ### list of CHIs for uploading -------------------------------------------------------

 # lookback period start date - must be numeric format to extract with
dfe$SPECDATE <- dfe$SpecimenDate
dfe$RFSTART <- as.numeric(format(dfe$SPECDATE - years(5), "%Y%m%d"))

 # temp file for uploading - character or numeric formats only
tp <- dfe[c("CHI", "SPECDATE", "RFSTART")]
tp$SPECDATE <- as.character(tp$SPECDATE)
str(tp)

# tp <- tp[1:10, ]

### upload CHIs and extract all episodes in period

SMRAconnection <- suppressWarnings(dbConnect(odbc(),
  dsn="SMRA",
  uid=.rs.askForPassword("Enter Username:"),
  pwd=.rs.askForPassword("Enter Password:") ))

dbWriteTable(SMRAconnection, "TEMP", tp, overwrite=T)
# TEMP <- dbGetQuery(SMRAconnection, statement='SELECT * FROM SHAROK01."TEMP"')
# View(TEMP)
# rm(TEMP)

extract <- tbl_df(dbGetQuery(SMRAconnection, statement=
  'SELECT
  T0.CHI, T0.SPECDATE, T0.RFSTART,
  T1.UPI_NUMBER, T1.DOB, T1.SEX, T1.ETHNIC_GROUP, T1.DR_POSTCODE, T1.LOCATION,
  T1.CIS_MARKER, T1.ADMISSION_DATE, to_number(to_char(T1.ADMISSION_DATE,\'YYYYMMDD\')) AS DOA,
  T1.ADMISSION_TYPE, T1.ADMISSION_TRANSFER_FROM, T1.ADMISSION_TRANSFER_FROM_LOC,
  T1.READY_FOR_DISCHARGE_DATE, T1.DISCHARGE_DATE, T1.DISCHARGE_TYPE, T1.DISCHARGE_TRANSFER_TO, T1.DISCHARGE_TRANSFER_TO_LOCATION,
  T1.SIGNIFICANT_FACILITY, T1.SPECIALTY, T1.LENGTH_OF_STAY, T1.CONSULTANT_HCP_RESPONSIBLE,
  T1.MAIN_CONDITION, T1.OTHER_CONDITION_1, T1.OTHER_CONDITION_2, T1.OTHER_CONDITION_3, T1.OTHER_CONDITION_4, T1.OTHER_CONDITION_5,
  T1.MAIN_OPERATION, T1.DATE_OF_MAIN_OPERATION, T1.OTHER_OPERATION_1, T1.DATE_OF_OTHER_OPERATION_1, T1.OTHER_OPERATION_2,
  T1.DATE_OF_OTHER_OPERATION_2, T1.OTHER_OPERATION_3, T1.DATE_OF_OTHER_OPERATION_3
  FROM SHAROK01."TEMP" T0, ANALYSIS.SMR01_PI T1
  WHERE T0.CHI = T1.UPI_NUMBER (+)
  AND to_number(to_char(T1.ADMISSION_DATE,\'YYYYMMDD\')) >= T0.RFSTART
  ORDER BY T1.UPI_NUMBER, T1.CIS_MARKER, T1.ADMISSION_DATE, T1.DISCHARGE_DATE, T1.ADMISSION, T1.DISCHARGE, T1.URI'
  )
)

dbRemoveTable(SMRAconnection, "TEMP")

### format dates
# character to date
extract$SPECDATE <- as.Date(extract$SPECDATE, "%Y-%m-%d") 

# POSIXct to date
cols <- c("ADMISSION_DATE", "DISCHARGE_DATE", "DOB", "DATE_OF_MAIN_OPERATION",
  "DATE_OF_OTHER_OPERATION_1", "DATE_OF_OTHER_OPERATION_2", "DATE_OF_OTHER_OPERATION_3")

for (i in cols) {
    extract[i] <- lapply(extract[i], as.Date, "%Y-%m-%d", tz="GMT")
    }


# CIS admission and discharge dates and episode order
extract <- mutate(extract %>% group_by(UPI_NUMBER, CIS_MARKER),
  cisdoa = min(ADMISSION_DATE),
  cisdodis = max(DISCHARGE_DATE),
  episode = 1:n())

# ethnic group - last available if missing - get last date too for comparison with SMR00
extract$ETHNIC_GROUP[extract$ETHNIC_GROUP %in% c("98", "99")] <- NA
extract <- extract[order(extract$CHI, extract$ADMISSION_DATE, extract$DISCHARGE_DATE),]
extract$last_ethnic_date <- ifelse(!is.na(extract$ETHNIC_GROUP), extract$cisdoa, NA)
extract <- extract %>% 
  group_by(CHI) %>% 
  tidyr::fill(ETHNIC_GROUP, last_ethnic_date,.direction = "down")
extract$last_ethnic_date <- as.Date(extract$last_ethnic_date, origin="1970-01-01", tz="GMT")

# working file
df <- extract




# ### Risk factors in stay containing test --------------------------------

# in hospital at time of test
df$hosp.test <- ifelse(df$SPECDATE >= df$cisdoa & df$SPECDATE <= df$cisdodis, 1, 0)
df$test.cis <- ifelse(df$hosp.test==1, df$CIS_MARKER, 0)

# in ICU/HDU/CCU during test stay
df$icuhduccu.test.cis <- ifelse(df$SIGNIFICANT_FACILITY %in% c("13", "1H", "14") & df$hosp.test==1, 1, 0)


# person level 
df <- df %>% group_by(UPI_NUMBER) %>% mutate(
  hosp.test=max(hosp.test), test.cis=first(test.cis), icuhduccu.test.cis=max(icuhduccu.test.cis)) %>% ungroup



# ### risk factors in lookback period ----------------------------------

df <- subset(df, ADMISSION_DATE < SPECDATE)


# ### Charlson Comorbidity Index score ------------------------------------

  # copy dataset
temp <- df

  # list ICD10 codes for required conditions

conheart3 <- c("I43", "I50")
conheart4 <- c("I099","I110","I130","I132","I255","I420","I425","I426","I427","I428","I429","P290")
demen3 <- c("F00","F01","F02","F03","G30")
demen4 <- c("F051","G311")
pulmon3 <- c("J40","J41","J42","J43","J44","J45","J46","J47","J60","J61","J62","J63","J64","J65","J66","J67")
pulmon4 <- c("I278","I279","J684","J701","J703")
conrheum3 <- c("M05","M32","M33","M34","M06")
conrheum4 <- c("M315","M351","M353","M360")
diabcomp4 <- c("E102","E103","E104","E105","E107","E112","E113","E114","E115","E117","E122","E123","E124","E125","E127","E132","E133","E134","E135","E137","E142","E143","E144","E145","E147")
parahemi3 <- c("G81","G82")
parahemi4 <- c("G041","G114","G801","G802","G830","G831","G832","G833","G834","G839")
renal3 <- c("N18","N19")
renal4 <- c("N052","N053","N054","N055","N056","N057","N250","I120","I131","N032","N033","N034","N035","N036","N037","Z490","Z491","Z492","Z940","Z992")
liver3 <- c("B18","K73","K74")
liver4 <- c("K700","K701","K702","K703","K709","K717","K713","K714","K715","K760","K762","K763","K764","K768","K769","Z944")
liverser4 <- c("K704","K711","K721","K729","K765","K766","K767","I850","I859","I864","I982")
hivaids3 <- c("B20","B21","B22","B24")
can3 <- c("C00","C01","C02","C03","C04","C05","C06","C07","C08","C09","C10","C11","C12","C13","C14","C15","C16","C17","C18","C19","C20","C21","C22","C23","C24","C25","C26","C30","C31","C32","C33","C34","C37","C38","C39","C40","C41","C43","C44","C45","C46","C47","C48","C49","C50","C51","C52","C53","C54","C55","C56","C57","C58","C60","C61","C62","C63","C64","C65","C66","C67","C68","C69","C70","C71","C72","C73","C74","C75","C76","C81","C82","C83","C84","C85","C88","C90","C91","C92","C93","C94","C95","C96","C97")
metcan3 <- c("C77","C78","C79","C80")

  # diagnosis variables for looping through

diagdata <- temp[c("MAIN_CONDITION", "OTHER_CONDITION_1", "OTHER_CONDITION_2", "OTHER_CONDITION_3", "OTHER_CONDITION_4", "OTHER_CONDITION_5")]

  # score for each diagnosis field

temp$conheart <- 0
for (i in diagdata){
  temp$conheart <- ifelse(
    substr(i, 1, 3) %in% conheart3 | substr(i, 1, 4) %in% conheart4
    , 2, temp$conheart
  )
}

temp$demen <- 0
for (i in diagdata){
  temp$demen <- ifelse(
    substr(i, 1, 3) %in% demen3 | substr(i, 1, 4) %in% demen4
    , 2, temp$demen
  )
}

temp$pulmon <- 0
for (i in diagdata){
  temp$pulmon <- ifelse(
    substr(i, 1, 3) %in% pulmon3 | substr(i, 1, 4) %in% pulmon4
    , 1, temp$pulmon
  )
}

temp$conrheum <- 0
for (i in diagdata){
  temp$conrheum <- ifelse(
    substr(i, 1, 3) %in% conrheum3 | substr(i, 1, 4) %in% conrheum4
    , 1, temp$conrheum
  )
}

temp$diabcomp <- 0
for (i in diagdata){
  temp$diabcomp <- ifelse(substr(i, 1, 4) %in% diabcomp4, 1, temp$diabcomp)
}

temp$parahemi <- 0
for (i in diagdata){
  temp$parahemi <- ifelse(
    substr(i, 1, 3) %in% parahemi3 | substr(i, 1, 4) %in% parahemi4
    , 2, temp$parahemi
  )
}

temp$renal <- 0
for (i in diagdata){
  temp$renal <- ifelse(
    substr(i, 1, 3) %in% renal3 | substr(i, 1, 4) %in% renal4
    , 1, temp$renal
  )
}

temp$hivaids <- 0
for (i in diagdata){
  temp$hivaids <- ifelse(substr(i, 1, 3) %in% hivaids3, 4, temp$hivaids)
}

temp$liver <- 0
for (i in diagdata){
  temp$liver <- ifelse(
    substr(i, 1, 3) %in% liver3 | substr(i, 1, 4) %in% liver4
    , 2,
      ifelse(substr(i, 1, 4) %in% liverser4, 4, temp$liver)
  )
}

temp$can <- 0
for (i in diagdata){
  temp$can <- ifelse(substr(i, 1, 3) %in% can3, 2,
    ifelse(substr(i, 1, 3) %in% metcan3, 6, temp$can)
  )
}

  # max scores for individual conditions
  # this should be at patient admission level

temp <- mutate(temp %>% group_by(CHI),
  conheart = max(conheart),
  demen = max(demen),
  pulmon = max(pulmon),
  conrheum = max(conrheum),
  diabcomp= max(diabcomp),
  parahemi= max(parahemi),
  renal= max(renal),
  hivaids= max(hivaids),
  liver= max(liver),
  can= max(can)
)

  # add scores for each condition to give overall score

temp$CCIscore <- rowSums(temp[,c("conheart", "demen", "pulmon", "conrheum", "diabcomp",
  "parahemi", "renal", "hivaids", "liver", "can")])

  # aggregate to give one row per person

temp <- temp %>% group_by(CHI) %>% summarise(
  CCIscore = max(CCIscore, na.rm = TRUE)
  ) %>% ungroup

  # match final score back to main file

df <- merge(df, temp , by=c("CHI"), all=TRUE)
df$CCIscore[is.na(df$CCIscore)] <- 0



# ### comorbidites --------------------------------------------------------


### emergency admissions - exclude 2 weeks prior to inf
df$SPECDATE <- as.Date(df$SPECDATE, "%Y-%m-%d")
df$emerg <- ifelse(df$ADMISSION_DATE < df$SPECDATE - days(14) & (df$ADMISSION_TYPE >=20 & df$ADMISSION_TYPE <= 22) | (df$ADMISSION_TYPE >= 30 & df$ADMISSION_TYPE <= 39), 1, 0)


### ICU/HDU/CCU (yes/no) - exclude 2 weeks prior to inf
df$SPECDATE <- as.Date(df$SPECDATE, "%Y-%m-%d")
df$icu.hdu.ccu <- ifelse(df$ADMISSION_DATE < df$SPECDATE - days(14) & df$SIGNIFICANT_FACILITY %in% c("13", "14", "1H"), 1, 0)


### Previous inpatient stay (yes/no)  - exclude 2 weeks prior to inf
df$SPECDATE <- as.Date(df$SPECDATE, "%Y-%m-%d")
df$inpat <- ifelse(df$ADMISSION_DATE < df$SPECDATE - days(14) & df$LENGTH_OF_STAY >=1, 1, 0)


### Covid19 case control comorbs (HC 23/4/2020)

vars <- df[c("MAIN_CONDITION", "OTHER_CONDITION_1", "OTHER_CONDITION_2", "OTHER_CONDITION_3", "OTHER_CONDITION_4", "OTHER_CONDITION_5")]

## CIRCULATORY DISORDERS	
# Ischaemic heart disease 	
df$IHD <- 0
for (i in vars){
  df$IHD <- ifelse(!is.na(i) &
  (substr(i, 1, 1) == "I" & as.numeric(substr(i, 2, 3)) %in% c(20:25)) 
    , 1, df$IHD)
  }

# Other heart disease 	
df$heart.other <- 0
for (i in vars){
  df$heart.other <- ifelse(!is.na(i) &
  (substr(i, 1, 1) == "I" & as.numeric(substr(i, 2, 3)) %in% c(00:02, 05:15, 26:28, 30:52)) 
    , 1, df$heart.other)
  }

# Cerebrovascular Disease	
df$CVD <- 0
for (i in vars){
  df$CVD <- ifelse(!is.na(i) &
  (substr(i, 1, 2) == "I6" | substr(i, 1, 3) %in% c("G45", "G46")) 
    , 1, df$CVD)
  }

# Other circulatory system diseases	
df$circulatory.other <- 0
for (i in vars){
  df$circulatory.other <- ifelse(!is.na(i) &
      (substr(i, 1, 2) %in% c("I7", "I8") | 
          (substr(i, 1, 1) == "I" & as.numeric(substr(i, 2, 3)) %in% c(95:99)) |
          substr(i, 1, 4) %in% c("Z958", "Z959")) 
    , 1, df$circulatory.other)
  }


## Neurological diseases except inflammatory	
# Epilepsy	
df$epilepsy <- 0
for (i in vars) {
  df$epilepsy <- ifelse(!is.na(i) &
  (substr(i, 1, 1) == "G" & as.numeric(substr(i, 2, 3)) %in% c(40:47))
    , 1, df$epilepsy)
  }

# Mono and polyneuropathies 
df$mono.poly.neuro <- 0
for (i in vars) {
  df$mono.poly.neuro <- ifelse(!is.na(i) & (
    substr(i, 1, 2) == "G5" | 
      (substr(i, 1, 1) == "G" & as.numeric(substr(i, 2, 3)) %in% c(60:64))   
    ), 1, df$mono.poly.neuro)
  }

# Other neurological conditions	
df$neuro.other <- 0
for (i in vars) {
  df$neuro.other <- ifelse(!is.na(i) & (
      substr(i, 1, 1) == "G" & as.numeric(substr(i, 2, 3)) %in% c(10:14, 20:26, 30:32, 35:37, 70:73, 80:83, 90:99)   
    ), 1, df$neuro.other)
  }


## Respiratory diseases	
# Acute respiratory infections
df$resp.inf <- 0
for (i in vars) {
  df$resp.inf <- ifelse(!is.na(i) & (
      substr(i, 1, 1) == "J" & as.numeric(substr(i, 2, 3)) %in% c(00:06, 09:18, 20:22)   
    ), 1, df$resp.inf)
  }

# Asthma	
df$asthma <- 0
for (i in vars) {
  df$asthma <- ifelse(!is.na(i) & (
      substr(i, 1, 1) == "J" & as.numeric(substr(i, 2, 3)) %in% c(45, 46)   
    ), 1, df$asthma)
  }

# Other Chronic lower respiratory disease
df$resp.other <- 0
for (i in vars) {
  df$resp.other <- ifelse(!is.na(i) & (
    substr(i, 1, 4) == "G473" | 
      (substr(i, 1, 1) == "J" & as.numeric(substr(i, 2, 3)) %in% c(40:44, 47, 60:70, 80:86, 90:99))   
    ), 1, df$resp.other)
  }

# Tuberculosis	
df$tuberculosis <- 0
for (i in vars) {
  df$tuberculosis <- ifelse(!is.na(i) & (
    substr(i, 1, 1) == "A" & as.numeric(substr(i, 2, 3)) %in% c(15:19)   
    ), 1, df$tuberculosis)
  }

## Connective tissue diseases	
# Connective tissue disorder	
df$connective <- 0
for (i in vars) {
  df$connective <- ifelse(!is.na(i) & (
    (substr(i, 1, 1) == "M" & as.numeric(substr(i, 2, 4)) %in% c(332, 353)) |
    (substr(i, 1, 1) == "M" & substr(i, 2, 2) == "0" & as.numeric(substr(i, 3, 4)) %in% c(50:53, 58:60, 63, 69)) |
    (substr(i, 1, 1) == "M" & as.numeric(substr(i, 2, 3)) %in% c(32, 34)) 
    ), 1, df$connective)
  }

##liver disease (SH 28/4/20)
df$liver <- 0
for (i in vars) {
  df$liver <- ifelse(!is.na(i) & (
    (substr(i, 1, 3) %in% c("K73", "R18")) |
    (substr(i, 1, 1) == "K" & as.numeric(substr(i, 2, 4)) %in% c(702:704, 717, 740, 742:746, 720:721, 729, 767 )) |
    (substr(i, 1, 4) %in% c("I850", "I983")) 
    ), 1, df$liver)
  }

## Kidney disease	
# Advanced Chronic kidney disease or RRT
df$kidney.advanced <- 0
for (i in vars) {
  df$kidney.advanced <- ifelse(!is.na(i) & (
    (substr(i, 1, 1) == "N" & as.numeric(substr(i, 2, 4)) %in% c(183:185))  |
    (substr(i, 1, 1) == "Z" & as.numeric(substr(i, 2, 4)) %in% c(490:492, 940, 992)) 
    ), 1, df$kidney.advanced)
  }


# Other chronic kidney disease	
df$kidney.other <- 0
for (i in vars) {
  df$kidney.other <- ifelse(!is.na(i) & (
    (substr(i, 1, 1) == "N" & as.numeric(substr(i, 2, 3)) %in% c(00:08, 10:16)) |
    (substr(i, 1, 1) == "N" & as.numeric(substr(i, 2, 4)) %in% c(181:182, 189))  
    ), 1, df$kidney.other)
  }

# Diabetes
df$diabetes <- 0
for (i in vars) {
  df$diabetes <- ifelse(!is.na(i) & (
    (substr(i, 1, 1) == "E" & as.numeric(substr(i, 2, 3)) %in% c(10:14)) 
    ), 1, df$diabetes)
  }

# Malignant neoplasms	
df$cancer <- 0
for (i in vars) {
  df$cancer <- ifelse(!is.na(i) & (
    (substr(i, 1, 1) == "C" & as.numeric(substr(i, 2, 3)) %in% c(00:97)) 
    ), 1, df$cancer)
  }

# Lung cancers	
df$lung.cancer <- 0
for (i in vars) {
  df$lung.cancer <- ifelse(!is.na(i) & (
    substr(i, 1, 3) == "C34"  
    ), 1, df$lung.cancer)
  }

# Blood cancers 	
df$blood.cancer <- 0
for (i in vars) {
  df$blood.cancer <- ifelse(!is.na(i) & (
    (substr(i, 1, 1) == "C" & substr(i, 2, 3) %in% c(81:96))  
    ), 1, df$blood.cancer)
  }


# Other cancers 
df$other.cancer <- 0
for (i in vars) {
  df$other.cancer <- ifelse(!is.na(i) & (
    (substr(i, 1, 1) == "C" & substr(i, 2, 3) %in% c(00:33, 35:80, 87))  
    ), 1, df$other.cancer)
}

## Immunological disease	
# HIV
df$HIV <- 0
for (i in vars) {
  df$HIV <- ifelse(!is.na(i) & (
    (substr(i, 1, 1) == "B" & substr(i, 2, 3) %in% c(20:23))  
    ), 1, df$HIV)
}


# Certain disorders involving the immune system	
df$immune <- 0
for (i in vars) {
  df$immune <- ifelse(!is.na(i) & (
    (substr(i, 1, 1) == "D" & substr(i, 2, 3) %in% c(80:89))  
    ), 1, df$immune)
  }

# Sickle cell disease	
df$SCD <- 0
for (i in vars) {
  df$SCD <- ifelse(!is.na(i) & (
    substr(i, 1, 3) == "D57"  
    ), 1, df$SCD)
  }

# Cystic fibrosis 	
df$cystic.fibrosis <- 0
for (i in vars) {
  df$cystic.fibrosis <- ifelse(!is.na(i) & (
    substr(i, 1, 3) == "E84"  
    ), 1, df$cystic.fibrosis)
  }

# Organ transplantation other than kidney 	
df$transplant.not.kidney <- 0
for (i in vars) {
  df$transplant.not.kidney <- ifelse(!is.na(i) & (
    substr(i, 1, 1) == "Z" & substr(i, 2, 4) %in% c(941:949) 
    ), 1, df$transplant.not.kidney)
  }


# ### ICD10 chapter flags -------------------------------------------------

vars_new <- as.character(01:22)
vars_new <- ifelse(nchar(vars_new)==1, paste0("0", vars_new), vars_new)
vars_new <- paste0("ICD10_", vars_new)

for (i in vars_new) { 
  df[i] <- 0 
  }

vars <- df[c("MAIN_CONDITION", "OTHER_CONDITION_1", "OTHER_CONDITION_2", "OTHER_CONDITION_3", "OTHER_CONDITION_4", "OTHER_CONDITION_5")]

for (i in vars) {
  df$ICD10_01 <- ifelse(!is.na(i) & (substr(i, 1, 1) %in% c("A", "B")), 1, df$ICD10_01)
  df$ICD10_02 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "C" | (substr(i, 1, 1) == "D" & as.numeric(substr(i, 2, 3)) <=48) ), 1, df$ICD10_02)  
  df$ICD10_03 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "D" & as.numeric(substr(i, 2, 3)) >=50 & as.numeric(substr(i, 2, 3)) <=89), 1, df$ICD10_03) 
  df$ICD10_04 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "E" & as.numeric(substr(i, 2, 3)) <=90), 1, df$ICD10_04) 
  df$ICD10_05 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "F"), 1, df$ICD10_05)
  df$ICD10_06 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "G"), 1, df$ICD10_06) 
  df$ICD10_07 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "H" & as.numeric(substr(i, 2, 3)) <=59), 1, df$ICD10_07)
  df$ICD10_08 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "H" & as.numeric(substr(i, 2, 3)) >=60 & as.numeric(substr(i, 2, 3)) <=95), 1, df$ICD10_08) 
  df$ICD10_09 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "I"), 1, df$ICD10_09) 
  df$ICD10_10 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "J"), 1, df$ICD10_10) 
  df$ICD10_11 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "K" & as.numeric(substr(i, 2, 3)) <=93), 1, df$ICD10_11)
  df$ICD10_12 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "L"), 1, df$ICD10_12) 
  df$ICD10_13 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "M"), 1, df$ICD10_13) 
  df$ICD10_14 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "N"), 1, df$ICD10_14)
  df$ICD10_15 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "O"), 1, df$ICD10_15)
  df$ICD10_16 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "P" & as.numeric(substr(i, 2, 3)) <=96), 1, df$ICD10_16)
  df$ICD10_17 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "Q"), 1, df$ICD10_17)
  df$ICD10_18 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "R"), 1, df$ICD10_18) 
  df$ICD10_19 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "S" | (substr(i, 1, 1) == "T" & as.numeric(substr(i, 2, 3)) <=98)), 1, df$ICD10_19)
  df$ICD10_20 <- ifelse(!is.na(i) & (substr(i, 1, 1) %in% c("V", "W", "X") | (substr(i, 1, 1) == "Y" & as.numeric(substr(i, 2, 3)) <=98)), 1, df$ICD10_20)
  df$ICD10_21 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "Z"), 1, df$ICD10_21)
  df$ICD10_22 <- ifelse(!is.na(i) & (substr(i, 1, 1) == "U"), 1, df$ICD10_22)
}

### aggregate to person level and match back to CC file
df <- df[order(df$CHI, df$ADMISSION_DATE, df$DISCHARGE_DATE),]
dff <- df %>% group_by(CHI) %>% summarise(
  LAST_CISDOA = max(cisdoa),
  LAST_CISDODIS = max(cisdodis),
  ETHNIC_GROUP = last(ETHNIC_GROUP),
  last_ethnic_date = max(last_ethnic_date),
  pc.smr = last(DR_POSTCODE),
  hosp.test = max(hosp.test),
  test.cis = max(test.cis), 
  icuhduccu.test.cis = max(icuhduccu.test.cis),
  CCIscore = max(CCIscore),
  emerg = max(emerg),
  icu.hdu.ccu = max(icu.hdu.ccu),
  inpat = max(inpat),
  IHD = max(IHD),
  heart.other = max(heart.other),
  CVD = max(CVD),
  circulatory.other = max(circulatory.other),
  epilepsy = max(epilepsy),
  mono.poly.neuro = max(mono.poly.neuro),
  neuro.other = max(neuro.other),
  resp.inf = max(resp.inf),
  asthma = max(asthma),
  resp.other = max(resp.other),
  tuberculosis = max(tuberculosis),
  connective = max(connective),
  liver = max(liver),
  kidney.advanced = max(kidney.advanced),
  kidney.other = max(kidney.other),
  diabetes = max(diabetes),
  cancer = max(cancer),
  lung.cancer = max(lung.cancer),
  blood.cancer = max(blood.cancer),
  other.cancer = max(other.cancer),
  HIV = max(HIV),
  immune = max(immune),
  SCD = max(SCD),
  cystic.fibrosis = max(cystic.fibrosis),
  transplant.not.kidney = max(transplant.not.kidney),
  SMR01 = 1
)

dff <- merge(x=dfe, y=dff, by="CHI", all.x=T)

# recode to 0 for those with no SMR in period

dff[, c(21:23, 25:ncol(dff))][is.na(dff[, c(21:23, 25:ncol(dff))])] <- 0


# ### SIMD -------------------------------------------------------

table(dff$pc7=="")
table(is.na(dff$pc.smr))

pc <- readRDS(paste0(lookups_path, "/Deprivation/postcode_2019_2_simd2020.rds"))
pc <- pc[c("pc7", "simd2020_sc_quintile")]

dff <- merge(x=dff, y=pc, by="pc7", all.x=T)

dff$simd2020_sc_quintile[is.na(dff$simd2020_sc_quintile)] <- 9

dff$simd <- factor(dff$simd2020_sc_quintile,
  levels = c(1, 2, 3, 4, 5, 9),
  labels = c("1 - most deprived", "2", "3", "4", "5 - least deprived", "Unknown"))
table (dff$simd)
table (is.na(dff$simd))
rm(pc)



# ### format and save ------------------------------------------


 #  CCI group

dff$CCIscore <- ifelse(dff$SMR01==0, NA, dff$CCIscore)
dff <- within(dff,{
  CCIgrp <- NA
  CCIgrp [is.na(CCIscore)] <- "[0]No_stay"
  CCIgrp [CCIscore ==0] <- "[0]No_score"
  CCIgrp [CCIscore >=1 & CCIscore <=2] <- "[1-2]Mild"
  CCIgrp [CCIscore >=3 & CCIscore <=5] <- "[3-5]Moderate"
  CCIgrp [CCIscore >=6] <- "[6+]Severe"  })
dff$CCIgrp <- as.factor(dff$CCIgrp)
table(dff$CCIgrp)


# save

saveRDS(dff, paste0(stats_path, "/HPS/Covid19/case_control/CC_NEW_linked_working.rds")) # server


### HC to send new codes for transplant 2020-05-01 ###



