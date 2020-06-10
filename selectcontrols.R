## Martin W's code for selecting controls

cc<- SelectMatchedControls(cases_temp,upip_ex_d,first.stratum.number =1 )

cc <-cc %>% 
    select(-mkey) %>%
    inner_join(upip_ex, by = c("upi"="UPI_NUMBER")) 
cc <-cc %>%
    mutate(AgeYear=floor((t-as_date(DATE_OF_BIRTH))/365.25)) %>%
    mutate(AgeYear=as.numeric(AgeYear))


#****************************************
#clean up
#****************************************

#input - remove columns and rename

cases_ex_output <- cc %>%
    filter(is.case==TRUE) %>%
    dplyr::select(upi, stratum, SEX,AgeYear, CURRENT_POSTCODE,GP_PRAC_NO, FIRST_FORENAME, SURNAME, INSTITUTION_CODE) %>%
    dplyr::rename(input_chi = upi,input_chi_sex = SEX,input_chi_postcode = CURRENT_POSTCODE,input_fn = FIRST_FORENAME,input_sn = SURNAME,input_ic = INSTITUTION_CODE)

#remove columns and join to input

matched_cases_output <- cc %>%
    filter(is.case==FALSE) %>%
    select(upi, SEX, AgeYear,CURRENT_POSTCODE, INSTITUTION_CODE,stratum, CURRENT_POSTCODE,GP_PRAC_NO )

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

cases_ex_out_final<-left_join(cases_ex_output,CHI_Prisoners,by=c("input_chi"="UPI_NUMBER"))
matched_cases_output_final<-left_join(matched_cases_output,CHI_Prisoners,by=c("upi"="UPI_NUMBER"))

fwrite(cases_ex_out_final,"/conf/linkage/output/y2k_cat_check/conf/input_case_June.csv")

fwrite(matched_cases_output_final,"/conf/linkage/output/y2k_cat_check/conf/matched_case_controls_June.csv")
