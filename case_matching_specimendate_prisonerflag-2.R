#****************************************
#matching case and control records
#graeme 2020-04-17
#run on RStudio Server Cluster, v3.6.1
#****************************************

#changes:
#updated to extract FIRST_FORENAME, SURNAME for input and cases
#define input col types
#add chi_pad to input cases
#change formatting to include forename/surname in output
#added some extra checks
#added IC - not run this time, but for next time

#****************************************
#load packages
#****************************************

library(tidyverse)
library(odbc)
library(lubridate)
library(tidylog)
library(data.table)
library(tictoc)
library(data.table)
library(stringdist)
library(stringr)
library(glue)
library(readr)

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

#load cases, define coltypes
#n = 7298
# cases1 <- read_csv("/conf/linkage/output/y2k_cat_check/conf/ECOSS.csv",
#                    col_types = "ccdcc")

#updated file 14-05-20 specimen dates added
# cases <- readRDS("/conf/linkage/output/HPS/Covid19/ECOSS_deduped.rds")#not this one

cases<-fread("/chi/(1) Project Folders/Case Control/180520_combined_seeded.csv")
cases<-cases %>% 
filter(!is.na(UPI_NUMBER)&nchar(UPI_NUMBER)==10)
# cases<-cases %>% 
#   filter(result==1)#only positive cases 14117 cases on 14th may

#check chis for validity

#rename CHI column to match UPIP
cases <- cases %>% rename(SPECIMENDATE = SpecimenDate)

#n = 14117 cases
#0 blank rows in csv input file
#cases <- cases %>% filter(!is.na(UPI_NUMBER))

#pad chi to 10 digits
# cases <- cases %>% mutate(UPI_NUMBER = ifelse(nchar(CHI) == 9,
#                                       paste0("0", CHI),
#                                   CHI))

# cases1 <-cases %>% 
#   select (UPI_NUMBER,LabSpecimenNo,SPECIMENDATE) 
# #fwrite(cases1,"/chi/(1) Project Folders/Case Control/CCseed140520_20043.csv")
#seed data in indexer

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

#match to CHI to get extra info (GP_practice)
#select only some columns, filter to include current CHI and arrange
#update: don't need all these columns?
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
  mutate(AgeYear=floor((t-as_date(DATE_OF_BIRTH))/365.25)) %>% 
  mutate(AgeYear=as.numeric(AgeYear)) %>% 
  select(UPI_NUMBER,AgeYear,SPECIMENDATE,GP_PRAC_NO,CURRENT_POSTCODE,SEX,INSTITUTION_CODE) 

cases_ex<-cases_ex2

#recode sex and remove blank GP practices
#n = 7278 cases left, 19 removed
cases_ex <- cases_ex %>% 
              #mutate(Sex = recode(Sex, "M" = 1, "F" = 2)) %>%
              filter(!is.na(GP_PRAC_NO))

#age bandings - increased from last time
age_breaks <- c(0, 4, 9, 14, 19, 24, 29, 34, 39, 44, 49,
                54, 59, 64, 69, 74, 79, 84, 89, 94, 99, 104, 109, 114, 119)

age_labels <- c("0-4","5-9","10-14","15-19","20-24","25-29",
                "30-34","35-39","40-44","45-49","50-54","55-59",
                "60-64","65-69","70-74","75-79","80-84","85-89",
                "90-94", "95-99", "100-104", "105-109", "110-114","115-119")

#group ages into those groups above
#make key for matching later
cases_ex <- cases_ex  %>%
            mutate(age_group = as.character(cut(AgeYear, 
                                                breaks = age_breaks, 
                                                labels = age_labels, 
                                                include.lowest = TRUE))) %>% 
            mutate(mkey = paste(age_group, SEX, GP_PRAC_NO, sep = "_")) %>% 
  mutate(SPECIMENDATE=as.Date(SPECIMENDATE))

#how many keys are unique? 9099/13908
n_distinct(cases_ex$mkey)

#****************************************
#extract matched info from UPIP####
#****************************************

#get list of practices to use, remove blanks and duplicates
#n = 878
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
#remove those with a date of death and keep only current CHI
#keep only unique rows
#join to gp practices from cases
#extracting the postcode changes the order so need to arrange
#to help make sure reproducible in future
#COLLECT \NRS DEATHS DATA


#get data for these practices only
#7746686 records
upip_ex <- (dbGetQuery(con, statement ="SELECT L_UPI_DATA.UPI_NUMBER,
L_UPI_DATA.CHI_STATUS,
L_UPI_DATA.DATE_OF_BIRTH,
L_UPI_DATA.DATE_OF_DEATH,
L_UPI_DATA.SEX,
L_UPI_DATA.FIRST_FORENAME,
L_UPI_DATA.SURNAME,
L_UPI_DATA.CURRENT_POSTCODE,
L_UPI_DATA.GP_PRAC_NO,
L_UPI_DATA.INSTITUTION_CODE
FROM MARTIR03.GP GP INNER JOIN UPIP.L_UPI_DATA L_UPI_DATA
ON (GP.GP_PRAC_NO = L_UPI_DATA.GP_PRAC_NO)
WHERE ( (L_UPI_DATA.CHI_STATUS = 'C') OR (L_UPI_DATA.CHI_STATUS IS NULL))"))

#work out age in years for matching
#maybe want to fix today to a specific date?
upip_ex <- upip_ex %>% 
            mutate(age = interval(start = DATE_OF_BIRTH, end = "2020-04-01"),
                    age = floor(time_length(age, unit = "year")))

#add age groups
upip_ex <- upip_ex %>% 
            mutate(age_group = as.character(cut(age, 
                                                breaks = age_breaks, 
                                                labels = age_labels, 
                                                include.lowest = TRUE)))

#add key for matching - 40185 unique keys
upip_ex <- upip_ex %>%
            mutate(mkey = paste(age_group, SEX, GP_PRAC_NO, sep = "_"))
cases_temp<-cases_ex %>% 
  select(mkey,SPECIMENDATE)
upip_ex_a<- left_join(upip_ex,cases_temp,by="mkey")
upip_ex_a<-upip_ex_a %>% 
  distinct()
#NRS death info to be filtered to ensure all deaths are after specimen date
NRS_Deaths <- readRDS("/conf/linkage/output/HPS/Covid19/deaths/NRS_Deaths.rds")
NRS_Deaths <-NRS_Deaths %>% 
  select(CHI,Date.Death)
upip_ex_b<-left_join(upip_ex_a,NRS_Deaths,by=c("UPI_NUMBER"="CHI"))
upip_ex_b<-upip_ex_b %>% 
  mutate(DATE_OF_DEATH=as.Date(DATE_OF_DEATH)) %>% 
  mutate(DATE_OF_DEATH=fifelse(is.na(DATE_OF_DEATH),Date.Death,DATE_OF_DEATH)) %>% 
  mutate(SPECIMENDATE=as.Date(SPECIMENDATE))
#remove original cases from extract: removed 7,271 rows
#filter extract to include only those that match the key: removed 5,939,055 rows
#no point searching the entire thing

upip_ex_b <- upip_ex_b %>% filter(!UPI_NUMBER %in% cases_ex$UPI_NUMBER)
upip_ex_b <- upip_ex_b %>% filter(mkey %in% cases_ex$mkey)

#remove duplicated rows: no rows removed
upip_ex_b <- upip_ex_b %>% distinct()

#fix order of controls for sampling
upip_ex_b <- upip_ex_b %>% 
            arrange(UPI_NUMBER, mkey)

#n = 1800360 cases left in upip_ex

#****************************************
#select matched records
#****************************************
#install.packages("tictoc")

#make empty list to store

ml <- list()

#save new object in case it doesn't work
#don't want to run extract again!
upip_ex2 <- upip_ex_b %>% select(UPI_NUMBER,SPECIMENDATE,DATE_OF_DEATH, mkey)

#cases_test<-slice(cases_ex,1:100)

tictoc::tic()

#set seed using R version 3.6.1
RNGkind() #check type
#[1] "Mersenne-Twister" "Inversion" "Rejection"   
set.seed(2020)

#loop over each UPI
for (i in cases_ex$UPI_NUMBER) {
  
  #extract the key from this row using pull
  z <- cases_ex %>% filter(UPI_NUMBER == i) %>% pull(mkey)
 
  #filter extract to include key
  x <- upip_ex2 %>% filter(mkey == z & (SPECIMENDATE < DATE_OF_DEATH |is.na(DATE_OF_DEATH)))
  x <-x %>% 
    select(-SPECIMENDATE) %>% 
    distinct(.,
             .keep_all = T)
  
  #sample random rows from this filtered extract
  #if there are less than 10, then select the max no
  #could add a seed here to make reproducible?
  if(nrow(x) < 10) {
    y <- sample_n(x, nrow(x)) 
  } else {
    y <- sample_n(x, 10)
  }
  
  #add input CHI for reference
 #y$input_chi <- i
  
  y<-y %>% mutate(input_chi=i)
  #remove those already extracted
  upip_ex2 <- upip_ex2 %>% filter(!UPI_NUMBER %in% y$UPI_NUMBER)
  
  #assign samples to list
  ml[[i]] <- y 
  
}

tictoc::toc()

#make list a dataframe, join back to all cols, arrange
#72610 records
matched_cases <- bind_rows(ml) %>%
                  select(-mkey) %>% 
                  inner_join(upip_ex, by = "UPI_NUMBER") %>% 
                  arrange(mkey)

#****************************************
#some checks
#****************************************
#1) do all inputs have 10? if not, which don't? n = 174 inputs
matched_cases %>% count(input_chi) %>% arrange(n) %>% filter(n < 10)

#which don't have any? n = 131 inputs
filter(cases_ex, !UPI_NUMBER %in% matched_cases$input_chi)

#2) how many are unique?
n_distinct(matched_cases$UPI_NUMBER) #controls 149993
n_distinct(matched_cases$input_chi) #cases 15073

#3) - age distributions?
#package for lining up plots
#install.packages("cowplot")
library(cowplot)
library(ggplot2)

#check age distributions
input_age_plot <- cases_ex %>% 
                    ggplot(aes(AgeYear)) + 
                    geom_histogram(fill = "#1b9e77") + 
                    ggtitle("input cases - age distribution")

matched_age_plot <- matched_cases %>% 
                      ggplot(aes(age)) + 
                      geom_histogram(fill = "#d95f02") + 
                      ggtitle("matched cases - age distribution")

#plot age distributions
plot_grid(input_age_plot, matched_age_plot, align = "hv", axis = "tblr")

#4) are all chi status current or NA?
matched_cases %>% count(CHI_STATUS)

#5) gender distribution?
matched_cases %>% count(SEX)
cases_ex %>% count(SEX)

#6) do any controls have date of death?#29
matched_cases %>% count(DATE_OF_DEATH.y)

#7) are any controls a case? NO
matched_cases %>% filter(UPI_NUMBER %in% cases_ex$UPI_NUMBER)

#****************************************
#clean up
#****************************************

#input - remove columns and rename
cases_ex_output <- cases_ex1 %>% 
                    select(UPI_NUMBER, SEX, CURRENT_POSTCODE, 
                            GP_PRAC_NO, FIRST_FORENAME, SURNAME, INSTITUTION_CODE) %>% 
                    rename(input_chi = UPI_NUMBER,
                          input_chi_sex = SEX,
                          input_chi_postcode = CURRENT_POSTCODE,
                          input_fn = FIRST_FORENAME,
                          input_sn = SURNAME,
                          input_ic = INSTITUTION_CODE)

#remove columns and join to input
matched_cases_output <- matched_cases %>% 
                        select(UPI_NUMBER, SEX, age, CURRENT_POSTCODE, INSTITUTION_CODE, input_chi) %>% 
                        left_join(cases_ex_output, by = "input_chi")


matched_cases_outputONO <- matched_cases %>% 
  select(UPI_NUMBER, SEX, age, SURNAME,FIRST_FORENAME,CURRENT_POSTCODE, INSTITUTION_CODE, input_chi) %>% 
  left_join(cases_ex_output, by = "input_chi")

matched_cases_outputONO<-matched_cases_outputONO %>% 
  select(UPI_NUMBER,SURNAME,FIRST_FORENAME)

cases_ex_ono<-cases_ex_output %>% 
  select(input_chi,input_sn,input_fn) %>% 
  rename(UPI_NUMBER=input_chi,SURNAME=input_sn,FIRST_FORENAME=input_fn)

ono_file<-rbind(cases_ex_ono,matched_cases_outputONO)
#****************************************
#save output cust and onomap
#****************************************
write_csv(ono_file,
          path = "/chi/(1) Project Folders/Case Control/ONO_200520.csv")
#process in onomap
ono_result<-fread("/chi/(1) Project Folders/Case Control/Ono_Export_200520.csv")

matched_cases_output_ono<-left_join(matched_cases_output,ono_result,by="UPI_NUMBER")#JOIN TO MATCHED CASES
ono_result_input<-ono_result %>% 
  rename(input_chi_Onolytics_Type_Code=`Onolytics Type Code`,input_chi_Onolytics_Group=`Onolytics Group`
         ,input_chi_Onolytics_Subgroup=`Onolytics Subgroup`,input_chi_Onolytics_Type=`Onolytics Type`,
         input_chi_Onolytics_coding_case=`Onolytics coding case`,input_chi_Geographical_Area=`Geographical Area`) %>% 
  select(-Religion,-`Personal Score`)
matched_cases_output_ono2<-left_join(matched_cases_output_ono,ono_result_input,by=c("input_chi"="UPI_NUMBER"))#ononmap varaibles JOIN TO CASES

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
  dplyr::mutate(UPI_NUMBER= if_else(nchar(UPI_NUMBER) == 9, paste0("0", UPI_NUMBER), paste0(UPI_NUMBER)))
fwrite(ALL,"/chi/(1) Project Folders/Case Control/Prisoners_April.csv")
################
#join these prison upi to case control and cases.

CHI_Prisoners<-fread("/chi/(1) Project Folders/Case Control/Prisoners_April.csv")
CHI_Prisoners<-CHI_Prisoners %>% 
mutate(UPI_NUMBER=as.character(UPI_NUMBER))%>%
  dplyr::mutate(UPI_NUMBER= if_else(nchar(UPI_NUMBER) == 9, paste0("0", UPI_NUMBER), paste0(UPI_NUMBER)))

#matched_cases_output_ono2 <- fread("/conf/linkage/output/y2k_cat_check/conf/case_controls210520_v2.csv")
matched_cases_output_ono2 <-matched_cases_output_ono2 %>% 
mutate(UPI_NUMBER=as.character(UPI_NUMBER))%>%
  dplyr::mutate(UPI_NUMBER= if_else(nchar(UPI_NUMBER) == 9, paste0("0", UPI_NUMBER), paste0(UPI_NUMBER)))

matched_cases_output_ono3<-left_join(matched_cases_output_ono2,CHI_Prisoners,by="UPI_NUMBER") #join to controls
matched_cases_output_ono3<-matched_cases_output_ono3 %>% 
  rename(CC_Prisoner_Flag=Prisoner_Flag) %>% 
  mutate(input_chi=as.character(input_chi))%>%
  dplyr::mutate(input_chi= if_else(nchar(input_chi) == 9, paste0("0", input_chi), paste0(input_chi)))
matched_cases_output_ono4<-left_join(matched_cases_output_ono3,CHI_Prisoners,by=c("input_chi"="UPI_NUMBER"))#join to cases
matched_cases_output_ono4<-matched_cases_output_ono4 %>% 
  rename(input_chi_Prisoner_Flag=Prisoner_Flag)
write_csv(matched_cases_output_ono4,
          path = "/conf/linkage/output/y2k_cat_check/conf/case_controls210520_v4.csv")


check<- fread("/conf/linkage/output/y2k_cat_check/conf/case_controls210520_v4.csv")
