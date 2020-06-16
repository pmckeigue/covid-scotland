#****************************************
#matching case and control records
#Martin 12-06-20
#run on RStudio Server Cluster, v3.6.1
#****************************************

#****************************************
#load packages
#****************************************


library(dplyr)
library(odbc)
library(lubridate)
library(tidylog)
library(data.table)
#library(tictoc)
library(stringdist)
library(glue)
library(phsmethods)
library(stringr)
# install.packages("remotes")
# remotes::install_github("Health-SocialCare-Scotland/phsmethods")


#########################################################
#New case matching function 11-06-20 from Luke Blackford#
#########################################################
SelectMatchedControls <- function(cases, popextract, first.stratum.number) {
  
  cases <- data.table::data.table(cases)
  popextract <- data.table::data.table(popextract)
  
  # Make a copy of cases so do not alter original by reference
  cases.mkey <- data.table::copy(cases)
  
  # For each mkey value in cases, count number n[k] of cases with this mkey value
  cases.mkey[, num.cases := .N, by=mkey]
  # For each mkey in population how many cases are there
  pop.m <- merge(popextract, unique(cases.mkey[, .(mkey, num.cases)]), how="left", on="mkey")
  # For each mkey in population how many potential controls are there
  pop.m[, num.pop.mkey := .N, by=mkey]
  
  # All cases are returned - then controls are appended
  ret.table <- cases
  ret.table$is.case <- TRUE
  
  # If <10*n[k] matched in population then take all available (no replacement)
  ret.pop.lt.10 <- pop.m[(num.cases > 0) & (num.pop.mkey < 10 * num.cases), .(upi, mkey)]
  ret.pop.lt.10$is.case <- FALSE
  ret.table <- rbind(ret.table, ret.pop.lt.10)
  
  # For >=10*n[k] in each mkey in popextract then randomly sample 10 * n[k] without replacement
  pop.m.gt.10 <- pop.m[(num.cases > 0) & (num.pop.mkey >= 10 * num.cases)]
  
  # Randomly order the data.table for sampling
  rand.order <- sample(nrow(pop.m.gt.10))
  ret.pop.gt.10 <- pop.m.gt.10[rand.order]
  
  # Select first 10*num.cases from each subgroup (.SD)
  ret.pop.gt.10 <- ret.pop.gt.10[, .SD[1:(10*unique(.SD$num.cases))], by=mkey][,.(upi, mkey)]

  
  ret.pop.gt.10$is.case <- FALSE
  ret.table <- rbind(ret.table, ret.pop.gt.10)
  
  # Add stratum - incrementing integer for each mkey starting at first.stratum.number
  mkey.vals <- unique(ret.table$mkey)
  stratum.vals <- seq(first.stratum.number,first.stratum.number+length(mkey.vals) - 1)
  ret.table[, stratum := as.integer(as.character(factor(mkey, mkey.vals, stratum.vals)))][]
  
  return(as.data.frame(ret.table))
}
##################
##################


#****************************************
#set up connection to database
#****************************************
con <- dbConnect(odbc(), dsn = "SMRA",
                 uid = rstudioapi::askForPassword("Database user"),
                 pwd = rstudioapi::askForPassword("Database password"),
                 port = "1527",
                 host = "nssstats01.csa.scot.nhs.uk",
                 SVC = "SMRA.nss.scot.nhs.uk")

#****************************************
#read/process incoming cases
#****************************************
Additional_20200608 <- readRDS("/conf/linkage/output/HPS/Covid19/Additional_NRS_Covid_Cases_2020-06-08.rds")
Additional_20200608 <-Additional_20200608 %>% 
  mutate(valid_extra=phsmethods::chi_check(Additional_20200608$CHI)) %>% 
  filter(valid_extra=="Valid CHI") %>% 
  select(CHI,SpecimenDate.dummy) %>% 
  rename(SpecimenDate=SpecimenDate.dummy,UPI_NUMBER=CHI) %>% 
  mutate(LabSpecimenNo=row_number()) %>% 
  mutate(LabSpecimenNo=as.character(LabSpecimenNo)) %>% 
  mutate(SpecimenDate=as.character(SpecimenDate))

ECOSS_deduped <- readRDS("/conf/linkage/output/HPS/Covid19/case_control/ECOSS_deduped_2020-06-08.rds")
#fwrite(ECOSS_deduped,"/chi/(1) Project Folders/Case Control/CCseed.csv")
ECOSS_deduped1 <-ECOSS_deduped %>% 
  #mutate(valid_extra=phsmethods::chi_check(CHI)) %>% 
  filter(result==1) %>% #ecoss filtered for positive results note some records without chi maybe dups and so could get 2 speciment dates - originally told only 1 specimen date per person in this file.
  mutate(PATID=as.character(PATID))

ECOSS_deduped2 <-ECOSS_deduped1 %>% 
  mutate(DateOB=as.character(DateOB)) %>%
  mutate(DateOB=str_replace_all(DateOB, "-", "")) %>% 
    mutate(PATID=as.character(PATID)) %>% 
  select(CHI,PATID,Forename,Surname,Sex,DateOB,PostCode)


fwrite(ECOSS_deduped2,"/chi/(1) Project Folders/Case Control/CCseed080620_2.csv")

#file from indexer after seeding
new_indexer_result<-fread("/chi/(1) Project Folders/Case Control/File_3118_UPI.csv")
#new_indexer_result<-fread("/chi/(1) Project Folders/Case Control/File_3093_UPI.csv")
new_indexer_result <- new_indexer_result %>% 
  select(SERIAL_NO,UPI) %>% 
  mutate(UPI=as.character(UPI)) %>% 
  mutate(UPI=phsmethods::chi_pad(UPI)) %>% 
  filter(!is.na(UPI)) %>% 
 rename(PATID=SERIAL_NO) %>% 
  mutate(PATID=as.character(PATID))

#seeded file with original payload
SEEDED<-left_join(new_indexer_result,ECOSS_deduped1,by="PATID")

SEEDED<-SEEDED %>% 
  select(UPI,LabSpecimenNo, SpecimenDate) %>% 
  rename(UPI_NUMBER=UPI) %>% 
  mutate(LabSpecimenNo=gsub(",","",LabSpecimenNo))  # commas in serial numbers!
  
#join the 2 incoming files
cases<-rbind(SEEDED,Additional_20200608)
# 12-06-20 teams conversation btween SK and MR
#if it is a duplicate across both then take the ecoss specimendate, 
#if it is a duplicate in ecoss then take the first specimendate. 
casesdups<-cases %>% 
  mutate(lab_len=nchar(LabSpecimenNo)) %>% 
  mutate(ecoss=if_else(lab_len<5,"NRS","ECOSS"))%>% #split file into origin of records
  group_by(UPI_NUMBER) %>% 
  mutate(n=n()) %>% 
ungroup() %>% 
filter(n=="2")#get dups

casesdups1<-casesdups %>% #duplicats internal to eccoss only
  group_by(UPI_NUMBER) %>% 
filter(ecoss=="ECOSS") %>% 
    mutate(min_date=min(SpecimenDate)) %>% 
  filter(min_date==SpecimenDate) %>% 
  ungroup() %>% 
  select(UPI_NUMBER,LabSpecimenNo,SpecimenDate)

casesdups2<-casesdups %>% #NRS duplicated in Ecoss file
  group_by(UPI_NUMBER) %>% 
  filter(ecoss=="NRS") %>% 
  ungroup() %>% 
  left_join(casesdups,by="UPI_NUMBER") %>% 
  filter(ecoss.y=="ECOSS") %>% #collect ECoss data rather than NRS.
  select(UPI_NUMBER,LabSpecimenNo.y,SpecimenDate.y) %>% 
  rename(LabSpecimenNo=LabSpecimenNo.y,SpecimenDate=SpecimenDate.y)

dups_sorted<-rbind(casesdups1,casesdups2)#join dups identified back together into 1 frame
dups_sorted<-dups_sorted %>% 
  distinct() 

cases_reduce<-cases %>% filter(!UPI_NUMBER %in% dups_sorted$UPI_NUMBER)#remove all dups from original df

dedupedcases<-rbind(cases_reduce,dups_sorted)#add corrected records from dups.
dedupedcases<-dedupedcases %>% 
  distinct()
cases<-dedupedcases 
  
#fwrite(cases,"/conf/linkage/output/HPS/Covid19/casecontrolupi_seeded.csv")

cases <- cases %>% rename(SPECIMENDATE = SpecimenDate)  

cases2 <-cases %>% 
  select (UPI_NUMBER,SPECIMENDATE) %>% 
  mutate(SPECIMENDATE=as_date(SPECIMENDATE))

#write to server for joining- added specimen date variable.
dbWriteTable(con,
             "GG1", 
             cases2, #data
             overwrite = TRUE,
             field.types = c(UPI_NUMBER ="VARCHAR2(10)"
             ))

#match to CHI to get (GP_practice)
#select only some columns, filter to include current CHI and arrange

cases_ex <- (dbGetQuery(con, statement = "SELECT DISTINCT L_UPI_DATA.GP_PRAC_NO,
                        L_UPI_DATA.CHI_STATUS,
                        L_UPI_DATA.SURNAME,
                        L_UPI_DATA.FIRST_FORENAME,
                        L_UPI_DATA.SEX,
                        L_UPI_DATA.DATE_OF_BIRTH,
                        L_UPI_DATA.CURRENT_POSTCODE,
                        L_UPI_DATA.INSTITUTION_CODE,
                        GG1_1.UPI_NUMBER
                        FROM MARTIR03.GG1 GG1_1 LEFT OUTER JOIN UPIP.L_UPI_DATA L_UPI_DATA
                        ON (GG1_1.UPI_NUMBER = L_UPI_DATA.UPI_NUMBER)
                        WHERE ( (L_UPI_DATA.CHI_STATUS = 'C') OR (L_UPI_DATA.CHI_STATUS IS NULL))"))

cases_ex1<-left_join(cases,cases_ex,by="UPI_NUMBER")
t=as_date(18353)

cases_ex2<-cases_ex1 %>% 
  mutate(AgeYear=floor((t-as_date(DATE_OF_BIRTH))/365.25)) %>% #calc age based on CHI DOB
  mutate(AgeYear=as.numeric(AgeYear)) %>% 
  select(UPI_NUMBER,AgeYear,SPECIMENDATE,GP_PRAC_NO,CURRENT_POSTCODE,SEX,INSTITUTION_CODE) 

cases_ex<-cases_ex2

#remove records with no GP practices

cases_ex <- cases_ex %>% 
    filter(!is.na(GP_PRAC_NO))#21 cases lost here due to no gp registration and therefore complete mkey isn't possible.

#Year age bandings - increased from last time not necessary was experimenting with need for age bandingin upper group sonly but not required this could be removed(12-06-20)
age_breaks <- c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,
                15,16,17,18,19,20,21,22,23, 24,25,26,27,
                28,29,30,31,32,33,34,35,36,37,38,39,
                40,41,42,43,44,45,46,47,48,49,50,51,52,
                53,54,55,56,57,58,59,60,61,62,63,64,65,
                66,67,68,69,70,71,72,73,74,75,76,77,78,
                79, 80,81,82,83,84,85,86,87,88,89,90,91,
                92,93,94,95,96,97,98,99,100,101,102,103,104,
                105,106,107,108,109,110,111,112,113,114,115,
                116,117,118,119)

age_labels <- c("0-1","2","3","4","5","6","7","8","9","10","11","12","13","14",
                "15","16","17","18","19","20","21","22","23","24","25","26","27",
                "28","29","30","31","32","33","34","35","36","37","38","39",
                "40","41","42","43","44","45","46","47","48","49","50","51","52",
                "53","54","55","56","57","58","59","60","61","62","63","64","65",
                "66","67","68","69","70","71","72","73", "74","75","76","77","78",
                "79", "80","81","82","83","84","85","86","87","88","89","90","91",
                "92","93","94","95","96","97","98","99","100","101","102","102","104",
                "105","106","107","108","109","110","111","112","113","114","115",
                "116","117","118","119")

#group ages into those groups above
#make key for matching later
cases_ex <- cases_ex  %>%
  mutate(age_group = as.character(cut(AgeYear, 
                                      breaks = age_breaks, 
                                      labels = age_labels, 
                                      include.lowest = TRUE))) %>% 
  mutate(mkey = paste(age_group, SEX, GP_PRAC_NO, sep = "_")) %>% 
  mutate(SPECIMENDATE=as.Date(SPECIMENDATE))

#how many keys are unique? 14768 of 16939 cases - 2227 not unique
n_distinct(cases_ex$mkey)

#**********************************
#extract matched info from UPIP####
#**********************************

#get list of practices to use, remove blanks and duplicates
#n = 883
gp <- cases_ex %>% 
  select(GP_PRAC_NO) %>% 
  filter(!is.na(GP_PRAC_NO)) %>%
  distinct()

#upload gp practices to use to connection
dbWriteTable(con, 
             "GP", 
             gp, #data
             overwrite = TRUE)

#select columns
#remove those with a date of death less than specimen date and keep only current CHI
#keep only unique rows
#join to gp practices from cases
#extracting the postcode changes the order so need to arrange
#get data for these practices only


upip_ex <- (dbGetQuery(con, statement ="SELECT DISTINCT L_UPI_DATA.UPI_NUMBER,
                       L_UPI_DATA.CHI_STATUS,
                       L_UPI_DATA.DATE_OF_BIRTH,
                       L_UPI_DATA.DATE_OF_DEATH,
                       L_UPI_DATA.SEX,
                       L_UPI_DATA.SURNAME,
                       L_UPI_DATA.FIRST_FORENAME,
                       L_UPI_DATA.DATE_TRANSFER_IN,
                       L_UPI_DATA.DATE_TRANSFER_OUT,
                       L_UPI_DATA.TRANSFER_OUT_CODE,
                       L_UPI_DATA.TRANSFER_IN_CODE,
                       L_UPI_DATA.CURRENT_POSTCODE,
                       L_UPI_DATA.GP_PRAC_NO,
                        L_UPI_DATA.INSTITUTION_CODE,
                       L_UPI_DATA.EXTENDED_STATUS
                       FROM MARTIR03.GP GP LEFT OUTER JOIN UPIP.L_UPI_DATA L_UPI_DATA
                       ON (GP.GP_PRAC_NO = L_UPI_DATA.GP_PRAC_NO)
                       WHERE ( (L_UPI_DATA.CHI_STATUS = 'C') OR (L_UPI_DATA.CHI_STATUS IS NULL))"))
#work out who is still in scotland.
EFFDATE=as_date(17986)
#9809706 in upip_ex 05-06-20
#05-06-20 only scottish residents added date transfer in / out and codes to allow filter in of those moving out of scotland or unknown.
upip_ex1<-upip_ex %>% 
  filter(EXTENDED_STATUS=="C"|is.na(EXTENDED_STATUS)|EXTENDED_STATUS=="d"|EXTENDED_STATUS=="D")
  #9827166records
  
  
  #remove people who are known to be outside scotland but don't remove deaths.
  #filter ((DATE_TRANSFER_IN <= EFFDATE|is.na(DATE_TRANSFER_IN)) & (EFFDATE <=DATE_TRANSFER_OUT|is.na(DATE_TRANSFER_OUT)))
  
  
  #work out age in years for matching
#maybe want to fix today to a specific date?
upip_ex1 <- upip_ex1 %>% 
  mutate(age = interval(start = DATE_OF_BIRTH, end = "2020-04-01"),
         age = floor(time_length(age, unit = "year")))

#add age groups
upip_ex1 <- upip_ex1 %>% 
  mutate(age_group = as.character(cut(age, 
                                      breaks = age_breaks, 
                                      labels = age_labels, 
                                      include.lowest = TRUE)))

#add key for matching - 43915 unique keys
upip_ex1 <- upip_ex1 %>%
  mutate(mkey = paste(age_group, SEX, GP_PRAC_NO, sep = "_"))
cases_temp<-cases_ex %>% 
  select(mkey,UPI_NUMBER,SPECIMENDATE)
upip_ex_a<- left_join(upip_ex1,cases_temp,by="mkey")
upip_ex_a<-upip_ex_a %>% 
  distinct()
#NRS death info to be filtered to ensure all deaths are after specimen date
NRS_Deaths <- readRDS("/conf/linkage/output/HPS/Covid19/deaths/NRS_Deaths.rds")
NRS_Deaths <-NRS_Deaths %>% 
  select(CHI,Date.Death)
upip_ex_b<-left_join(upip_ex_a,NRS_Deaths,by=c("UPI_NUMBER.x"="CHI"))
upip_ex_b<-upip_ex_b %>% 
  mutate(DATE_OF_DEATH=as.Date(DATE_OF_DEATH)) %>% 
  mutate(DATE_OF_DEATH=fifelse(is.na(DATE_OF_DEATH),Date.Death,DATE_OF_DEATH)) %>% 
  mutate(SPECIMENDATE=as.Date(SPECIMENDATE)) %>% 
  filter((DATE_TRANSFER_IN<=SPECIMENDATE |is.na(DATE_TRANSFER_IN)) & (SPECIMENDATE <=DATE_TRANSFER_OUT |is.na(DATE_TRANSFER_OUT)))

#filter extract to include only those that match the key: removed 4,825,082 ROWS
#no point searching the entire thing for controls
#NRS DATA UPDATES 'DATE_OF_DEATH' (249 fewer NA in file now)

#upip_ex_b <- upip_ex_b %>% filter(!UPI_NUMBER %in% cases_ex$UPI_NUMBER)Spec change cases can be controls too 09-06-20
upip_ex_b <- upip_ex_b %>% filter(mkey %in% cases_ex$mkey)

#remove duplicated rows: no rows removed
upip_ex_b <- upip_ex_b %>% distinct()

upip_ex_d <- upip_ex_b %>% 
  dplyr::rename(upi=UPI_NUMBER.x) %>% 
select(upi,mkey) 

cases_ex1<-cases_ex %>% 
    dplyr::rename(upi=UPI_NUMBER) %>%
  select(mkey,upi,SPECIMENDATE) %>% 
  filter(!is.na(SPECIMENDATE))  
#check duplicates
# cases_ex1_check<-cases_ex1 %>% 
#   group_by(upi) %>% 
#   mutate(n=n()) %>% 
#   filter(n>"1")


#loop over all specimen dates checking mkey matches 
cc.table<-NULL#assign a nul object to work with


startnum=1 
for( i in unique(cases_ex1$SPECIMENDATE)) {#take unique incoming specimen dates
 
  s.cases <- cases_ex1%>% filter(SPECIMENDATE == i) %>% select(mkey,upi)
  
  s.popextract <- upip_ex_d %>% filter(mkey %in% s.cases$mkey)#filter extract to include key for current case
  #cat("Selecting matched controls starting at stratum number", startnum, "\n")#debug line not needed
  s.cc <- SelectMatchedControls(cases=s.cases, popextract=s.popextract, first.stratum.number=startnum)#call function wich selects 10 controls and assigns stratum by specimen date
  
  cc.table<-rbind(cc.table, s.cc)#write result to table
  
  
  startnum=1 + max(cc.table$stratum)#go round loop until max specimen date reached
  #startnum=1 + cc.table[, .N]
}

#debug reporting
#x<-table(cc.table$stratum,cc.table$is.case)

cc.table_j<-cc.table %>%
  select(-mkey) %>%
  inner_join(upip_ex, by = c("upi"="UPI_NUMBER"))#add payload
cc.table_j<-cc.table_j %>%
  distinct() %>% 
  mutate(AgeYear=floor((t-as_date(DATE_OF_BIRTH))/365.25)) %>%
  mutate(AgeYear=as.numeric(AgeYear))
cc.table_e<-left_join(cc.table_j,cases_ex1,by="upi")#add original incomg variables

#****************************************
#clean up
#****************************************

#input - remove columns and rename
output <- cc.table_e %>% 
  dplyr::select(upi, is.case, stratum,SPECIMENDATE, SEX,AgeYear, CURRENT_POSTCODE,GP_PRAC_NO, INSTITUTION_CODE) 
 

#######get prisoner records identified and extracted from CHI database
PRISONS<-fread("/chi/(1) Project Folders/Case Control/Prisons.csv")
PRISONS<-PRISONS %>% 
  select(-X) %>% 
  mutate(PRISONPC=gsub(" ","",PRISONPC)) %>% 
  mutate(HC_PC=gsub(" ","",HC_PC)) %>% 
  filter(!is.na(GPPRISON))

con <- dbConnect(odbc(), dsn = "SMRA",
                 uid = rstudioapi::askForPassword("Database user"),
                 pwd = rstudioapi::askForPassword("Database password"),
                 port = "1527",
                 host = "nssstats01.csa.scot.nhs.uk",
                 SVC = "SMRA.nss.scot.nhs.uk")


dbWriteTable(con, "PRISONS",
             PRISONS,
             overwrite = TRUE,
             field.types = c(PRISON="VARCHAR2(25)",
                             LINE1="VARCHAR2(50)",
                             LINE2="VARCHAR2(50)",
                             PRISONPC="VARCHAR2(7)",
                             HC_LINE1="VARCHAR2(50)",
                             HC_LINE2="VARCHAR2(50)",
                             HC_PC="VARCHAR2(7)",
                             GPPRISON="VARCHAR(5)"
             ))

PRISONERS1<- (dbGetQuery(con, statement =glue::glue("SELECT DISTINCT PRISONS.PRISON,
                                                    PRISONS.LINE1,
                                                    PRISONS.LINE2,
                                                    PRISONS.PRISONPC,
                                                    PRISONS.HC_LINE1,
                                                    PRISONS.HC_LINE2,
                                                    PRISONS.HC_PC,
                                                    PRISONS.GPPRISON,
                                                    L_UPI_DATA.UPI_NUMBER,
                                                    L_UPI_DATA.CURRENT_LINE1,
                                                    L_UPI_DATA.CURRENT_LINE2,
                                                    L_UPI_DATA.CURRENT_LINE3,
                                                    L_UPI_DATA.CURRENT_POSTCODE,
                                                    L_UPI_DATA.CHI_STATUS,
                                                    L_UPI_DATA.DATE_OF_BIRTH,
                                                    L_UPI_DATA.SEX,
                                                    L_UPI_DATA.DATE_OF_DEATH,
                                                    L_UPI_DATA.GP_PRAC_NO,
                                                    L_UPI_DATA.INSTITUTION_CODE,
                                                    L_UPI_DATA.FIRST_FORENAME,
                                                    L_UPI_DATA.SURNAME
                                                    FROM {toupper(Sys.info()['user'])}.PRISONS PRISONS INNER JOIN UPIP.L_UPI_DATA L_UPI_DATA
                                                    ON (PRISONS.HC_PC = L_UPI_DATA.CURRENT_POSTCODE)
                                                    WHERE ( (L_UPI_DATA.CHI_STATUS = 'C') OR (L_UPI_DATA.CHI_STATUS IS NULL))
                                                    ")))

PRISONERS1<-PRISONERS1 %>% 
  mutate(HC_LINE1=toupper(HC_LINE1)) %>% 
  mutate(HC_LINE2=toupper(HC_LINE2)) %>% 
  mutate(prison_add=paste0(LINE1,LINE2,HC_PC)) %>% 
  mutate(prison_add=gsub("THE HEALTH CENTRE","",prison_add)) %>% 
  
  mutate(CHI_prison_add=paste0(CURRENT_LINE1,CURRENT_LINE2,CURRENT_LINE3,CURRENT_POSTCODE)) %>% 
  mutate(CHI_prison_add=gsub("THE HEALTH CENTRE","",prison_add)) %>% 
  mutate(ADD_SIM = (stringdist(prison_add,CHI_prison_add, method =c("jw"),p=0.1))) %>% 
  filter(ADD_SIM==0) %>% 
  select(UPI_NUMBER)

PRISONERS2<- (dbGetQuery(con, statement =glue::glue("SELECT DISTINCT PRISONS.PRISON,
                                                    PRISONS.LINE1,
                                                    PRISONS.LINE2,
                                                    PRISONS.PRISONPC,
                                                    PRISONS.HC_LINE1,
                                                    PRISONS.HC_LINE2,
                                                    PRISONS.HC_PC,
                                                    PRISONS.GPPRISON,
                                                    L_UPI_DATA.UPI_NUMBER,
                                                    L_UPI_DATA.CURRENT_LINE1,
                                                    L_UPI_DATA.CURRENT_LINE2,
                                                    L_UPI_DATA.CURRENT_LINE3,
                                                    L_UPI_DATA.CURRENT_POSTCODE,
                                                    L_UPI_DATA.CHI_STATUS,
                                                    L_UPI_DATA.DATE_OF_BIRTH,
                                                    L_UPI_DATA.SEX,
                                                    L_UPI_DATA.DATE_OF_DEATH,
                                                    L_UPI_DATA.GP_PRAC_NO,
                                                    L_UPI_DATA.INSTITUTION_CODE,
                                                    L_UPI_DATA.FIRST_FORENAME,
                                                    L_UPI_DATA.SURNAME
                                                    FROM {toupper(Sys.info()['user'])}.PRISONS PRISONS INNER JOIN UPIP.L_UPI_DATA L_UPI_DATA
                                                    ON (PRISONS.PRISONPC = L_UPI_DATA.CURRENT_POSTCODE)
                                                    WHERE ( (L_UPI_DATA.CHI_STATUS = 'C') OR (L_UPI_DATA.CHI_STATUS IS NULL))")))

PRISONERS2a<-PRISONERS2 %>% 
  mutate(LINE1=toupper(LINE1)) %>% 
  mutate(LINE2=toupper(LINE2)) %>% 
  mutate(PRISON=toupper(PRISON)) %>%
  mutate(prison_add=paste0(PRISON,LINE1,LINE2,PRISONPC)) %>% 
  mutate(prison_add=gsub("Road","RD",prison_add)) %>% 
  mutate(prison_add=gsub("HM PRISON","HMP",prison_add)) %>%
  mutate(prison_add=gsub(" ","",prison_add)) %>%
  mutate(PRISONaddno = gsub("[^0-9-]", "",(prison_add))) %>% 
  mutate(CHI_prison_add=paste0(CURRENT_LINE1,CURRENT_LINE2,CURRENT_LINE3,CURRENT_POSTCODE)) %>% 
  mutate(CHI_prison_add=gsub("Road","RD",CHI_prison_add)) %>% 
  mutate(CHI_prison_add=gsub("HM Prison","HMP",CHI_prison_add)) %>%
  mutate(CHI_prison_add=gsub("H.M.Prison","HMP",CHI_prison_add)) %>%
  mutate(CHI_prison_add=gsub("H.M. INSTITUTE","HMP",CHI_prison_add)) %>% 
  mutate(CHI_prison_add=gsub("H.M.P.","HMP",CHI_prison_add)) %>%
  mutate(CHI_prison_add=gsub("HMP","HMP",CHI_prison_add)) %>%
  mutate(CHI_prison_add=gsub("C/O","",CHI_prison_add)) %>%
  mutate(CHI_prison_add=gsub(" ","",CHI_prison_add)) %>%
  mutate(ADD_SIM = (stringdist(prison_add,CHI_prison_add, method =c("jw"),p=0.1))) %>% 
  mutate(CHIaddno = gsub("[^0-9-]", "",(CHI_prison_add))) %>% 
  filter(ADD_SIM<0.3) %>% 
  select(UPI_NUMBER)
#filter(PRISONaddno==CHIaddno)

ALL<-rbind(PRISONERS1,PRISONERS2a)
ALL<-ALL %>% 
  distinct() %>% 
  mutate(UPI_NUMBER=as.character(UPI_NUMBER))%>%
  mutate(Prisoner_Flag=1) %>% 
  dplyr::mutate(UPI_NUMBER= ifelse(nchar(UPI_NUMBER) == 9, paste0("0", UPI_NUMBER), paste0(UPI_NUMBER)))
fwrite(ALL,"/chi/(1) Project Folders/Case Control/Prisoners_June.csv")
################
#join these prison upi to case control and cases.

CHI_Prisoners<-fread("/chi/(1) Project Folders/Case Control/Prisoners_June.csv")
CHI_Prisoners<-CHI_Prisoners %>% 
  mutate(UPI_NUMBER=as.character(UPI_NUMBER))%>%
  dplyr::mutate(UPI_NUMBER= ifelse(nchar(UPI_NUMBER) == 9, paste0("0", UPI_NUMBER), paste0(UPI_NUMBER)))

final<-left_join(output,CHI_Prisoners,by=c("upi"="UPI_NUMBER"))
#check for dup cases
final_check<-final %>% 
  group_by(upi) %>% 
  mutate(n_upi=n()) %>% 
  filter(is.case=="TRUE")

fwrite(final,"/conf/linkage/output/y2k_cat_check/conf/case_and_controls_June_2.csv")
 
