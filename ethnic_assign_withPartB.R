######## coding ethnicity ##############################

## Source variables: 

## ethnic.smr - raw SMR categories
## OnolyticsType - raw Onomap types based on name classification
## Geographical Area - regional groupings of OnolyticsType

## Output variables: ethnic5.smr, eth, eth5

## this script

## (1) recodes ETHNIC_smr to specify Chinese as separate category, then collapses to 5 categories: White, South Asian, Chinese, Black, Other

## We cannot combine SMR coding with Onomap coding because this will introduce non-differential missclassification of ethnicity between cases and controls, especially for those categorised as Black.

## The Onolytics type "MUSLIM" is coded as South Asian - this is reasonable in Scotland where most people with Muslim names are South Asian.


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

collapseto5.ethnicsmr <- function(ethnic.smr) {
    ethnic5.smr <- character(length(ethnic.smr))
    ethnic5.smr[grep("^1[A-Z]$", ethnic.smr)] <- "White"
    ethnic5.smr[grep("^2[A-Z]$", ethnic.smr)] <- "Other"
    ethnic5.smr[grep("^3[ABCFGH]$", ethnic.smr)] <- "South Asian"
    ethnic5.smr[grep("^3[DJ]$", ethnic.smr)] <- "Chinese"
    ethnic5.smr[grep("^3[EZ]$", ethnic.smr)] <- "Other"
    ethnic5.smr[grep("^4[ABCDEY]$", ethnic.smr)] <- "Black"
    ethnic5.smr[grep("^5[ABCDY]$", ethnic.smr)] <- "Black"
    ethnic5.smr[grep("^6[AZ]$", ethnic.smr)] <- "Other"
    ethnic5.smr[grep("^9", ethnic.smr)] <- NA
    ethnic5.smr[ethnic5.smr==""] <- NA
    ethnic5.smr <- as.factor(ethnic5.smr)
    ethnic5.smr <- factor(ethnic5.smr, levels=levels(ethnic5.smr)[c(5, 4, 2, 1, 3)])
    return(ethnic5.smr)
}

group.onomap <- function(OnolyticsType, GeographicalArea) {    
    OnolyticsType <- car::recode(OnolyticsType,
                                 "'NOT FOUND'=NA; 'INTERNATIONAL'=NA; 'UNCLASSIFIED'=NA; 'VOID'=NA; 'VOID - FORENAME'=NA; 'VOID INITIAL'=NA")
    
    table(OnolyticsType[GeographicalArea=="SOUTH ASIA"])
    table(OnolyticsType[GeographicalArea=="AFRICA"])
    table(OnolyticsType[GeographicalArea=="BRITISH ISLES"])
    table(OnolyticsType[GeographicalArea=="EAST ASIA"])
    table(OnolyticsType[GeographicalArea=="MIDDLE EAST"])
    
    eth <- rep("Other", length(OnolyticsType))
    eth[is.na(OnolyticsType)] <- NA
    eth[OnolyticsType=="CHINESE" |
        OnolyticsType=="HONG KONGESE" |
        OnolyticsType=="SINGAPORESE"] <- "Chinese"
    eth[GeographicalArea=="EAST ASIA" &
        eth != "Chinese"] <- "Other Asia & Pacific"
    
    eth[GeographicalArea=="SOUTH ASIA"] <- "South Asian"
    
    eth[GeographicalArea=="BRITISH ISLES"] <- "Britain&Ireland"
    eth[GeographicalArea=="CENTRAL EUROPE" |
        GeographicalArea=="EASTERN EUROPE" |
        GeographicalArea=="NORTHERN EUROPE" |
        GeographicalArea=="SOUTHERN EUROPE" |
        OnolyticsType=="AFRIKAANS"] <- "Other Europe"
    
    eth[GeographicalArea=="AFRICA" &
        OnolyticsType != "AFRIKAANS" &
        OnolyticsType != "LIBYAN"] <- "Black African"
    eth[GeographicalArea=="MIDDLE EAST" &
        OnolyticsType != "MUSLIM"] <- "East Med"
    eth[OnolyticsType == "MUSLIM"] <- "Muslim, not localized"
    
    eth[OnolyticsType == "BLACK CARIBBEAN"] <- "Black Caribbean"
    return(eth)
}

collapseto5.onomap.group <- function(onomap.group) {
    ## reduce to 5 categories: White, South Asian, Chinese, Black, Other
    ethnic5 <- car::recode(onomap.group, "'Black African'='Black'; 'Black Caribbean'='Black';  'Britain&Ireland'='White';  'Other Europe'='White';  'East Med'='Other'; 'Other Asia & Pacific'='Other'; 'Muslim, not localized'='Other'") 
    ethnic5 <- relevel(as.factor(ethnic5), ref="White")
    return(ethnic5)
}


### DISSAGGREGATED RECODES FOR ETHNICITY REPORT

## extra codes: 1[DEFGHJ] White, 3[ABC] South Asian, 3[E] Chinese, 3D Other Asian, 4[ABCE]  African

collapseto9.ethnicsmr <- function(ethnic.smr) {
  ethnic9.smr <- character(length(ethnic.smr))
  ethnic9.smr[grep("^1[A-Z]$", ethnic.smr)] <- "White"
  ethnic9.smr[grep("^1[L]$", ethnic.smr)] <- "White Polish"
  ethnic9.smr[grep("^2[A-Z]$", ethnic.smr)] <- "Other"
  ethnic9.smr[grep("^3[BCFH]$", ethnic.smr)] <- "Pakistani/Bangladeshi"
  ethnic9.smr[grep("^3[AG]$", ethnic.smr)] <- "Indian"
  ethnic9.smr[grep("^3[DJ]$", ethnic.smr)] <- "Chinese"
  ethnic9.smr[grep("^3[EZ]$", ethnic.smr)] <- "Other"
  ethnic9.smr[grep("^4[AE]$", ethnic.smr)] <- "Caribbean"
  ethnic9.smr[grep("^4[BDY]$", ethnic.smr)] <- "African"
  ethnic9.smr[grep("^4[C]$", ethnic.smr)] <- "Black"
  #ethnic9.smr[grep("^4[F]$", ethnic.smr)] <- "Black" #4F not included in coding above
  ethnic9.smr[grep("^5[C]$", ethnic.smr)] <- "Caribbean"
  ethnic9.smr[grep("^5[ABDY]$", ethnic.smr)] <- "Black"
  ethnic9.smr[grep("^5[Z]$", ethnic.smr)] <- "Other"
  #ethnic9.smr[grep("^4[Z]$", ethnic.smr)] <- "Other" #not included in coding above
  ethnic9.smr[grep("^6[AZ]$", ethnic.smr)] <- "Other"
  ethnic9.smr[grep("^9", ethnic.smr)] <- NA
  ethnic9.smr[ethnic9.smr==""] <- NA
  ethnic9.smr <- as.factor(ethnic9.smr)
  ethnic9.smr <- factor(ethnic9.smr, levels=levels(ethnic9.smr)[c(8,9,3,1,2,4,7,5,6)])
  return(ethnic9.smr)
}

#White, White Polish, WHite Irish, Cariabbean, African, Black, Chinese, Pakistani/Bangladeshi, Indian, Other Asian, Other
collapseto13.ethnicsmr <- function(ethnic.smr) {
  ethnic13.smr <- character(length(ethnic.smr))
  ethnic13.smr[grep("^1[A-Z]$", ethnic.smr)] <- "White"
  ethnic13.smr[grep("^1[C]$", ethnic.smr)] <- "White Irish"
  ethnic13.smr[grep("^1[L]$", ethnic.smr)] <- "White Polish"
  ethnic13.smr[grep("^2[A-Z]$", ethnic.smr)] <- "Other"
  ethnic13.smr[grep("^3[BCFH]$", ethnic.smr)] <- "Pakistani/Bangladeshi"
  ethnic13.smr[grep("^3[AG]$", ethnic.smr)] <- "Indian"
  ethnic13.smr[grep("^3[DJ]$", ethnic.smr)] <- "Chinese"
  ethnic13.smr[grep("^3[EZ]$", ethnic.smr)] <- "Other Asian"
  ethnic13.smr[grep("^4[AE]$", ethnic.smr)] <- "Caribbean"
  ethnic13.smr[grep("^4[BD]$", ethnic.smr)] <- "African"
  ethnic13.smr[grep("^4[Y]$", ethnic.smr)] <- "Other African"
  ethnic13.smr[grep("^4[C]$", ethnic.smr)] <- "Black"
  #ethnic13.smr[grep("^4[F]$", ethnic.smr)] <- "Black" #4F not included in coding above
  ethnic13.smr[grep("^5[C]$", ethnic.smr)] <- "Caribbean"
  ethnic13.smr[grep("^5[D]$", ethnic.smr)] <- "Black"
  ethnic13.smr[grep("^5[Y]$", ethnic.smr)] <- "Other Caribbean or Black"
  #ethnic13.smr[grep("^4[Z]$", ethnic.smr)] <- "Other" #not included in coding above
  ethnic13.smr[grep("^6[AZ]$", ethnic.smr)] <- "Other"
  ethnic13.smr[grep("^9", ethnic.smr)] <- NA
  ethnic13.smr[ethnic13.smr==""] <- NA
  ethnic13.smr <- as.factor(ethnic13.smr)
  ethnic13.smr <- factor(ethnic13.smr, levels=levels(ethnic13.smr)[c(11,12,13,1,7,3,2,9,4,10,5,8,6)])
  return(ethnic13.smr)
}



group.onomap2 <- function(OnolyticsType, GeographicalArea) {    
  OnolyticsType <- car::recode(OnolyticsType,
                               "'NOT FOUND'=NA; 'INTERNATIONAL'=NA; 'UNCLASSIFIED'=NA; 'VOID'=NA; 'VOID - FORENAME'=NA; 'VOID INITIAL'=NA")
  
  table(OnolyticsType[GeographicalArea=="SOUTH ASIA"])
  table(OnolyticsType[GeographicalArea=="AFRICA"])
  table(OnolyticsType[GeographicalArea=="BRITISH ISLES"])
  table(OnolyticsType[GeographicalArea=="EAST ASIA"])
  table(OnolyticsType[GeographicalArea=="MIDDLE EAST"])
  
  eth8 <- rep("Other", length(OnolyticsType))
  eth8[is.na(OnolyticsType)] <- NA
  eth8[OnolyticsType=="CHINESE" |
         OnolyticsType=="HONG KONGESE" |
         OnolyticsType=="SINGAPORESE"] <- "Chinese"
  eth8[GeographicalArea=="EAST ASIA" &
         eth8 != "Chinese"] <- "Other Asia & Pacific"
  
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
  return(eth8)
}

collapseto8.onomap.group <- function(onomap.group) {
  ## reduce to 5 categories: White, South Asian, Chinese, Black, Other
  ethnic8 <- car::recode(onomap.group, "'Britain&Ireland'='White';  'Other Europe'='White';  'East Med'='Other'; 'Other Asia & Pacific'='Other'") 
  ethnic8 <- as.factor(ethnic8)
  ethnic8 <- factor(ethnic8, levels=levels(ethnic8)[c(7,8,2,1,3,4,6,5)])
  ethnic8 <- relevel(as.factor(ethnic8), ref="White")
  return(ethnic8)
}

collapseto7.onomap.group <- function(onomap.group) {
  ## reduce to 5 categories: White, South Asian, Chinese, Black, Other
  ethnic8 <- car::recode(onomap.group, "'Britain&Ireland'='White';  'Other Europe'='White';  'East Med'='Other'; 'Other Asia & Pacific'='Other'") 
  ethnic7 <- as.character(ethnic8)
  ethnic7 <- car::recode(ethnic7, "'Black Caribbean'='Black'")
  ethnic7 <- car::recode(ethnic7, "'Black African'='Black'")
  ethnic7 <- as.factor(ethnic7)
  ethnic7 <- factor(ethnic7, levels=levels(ethnic7)[c(6,7,1,2,3,5,4)])
  ethnic7 <- relevel(as.factor(ethnic7), ref="White")
  return(ethnic7)
}

collapseto6.onomap.group <- function(onomap.group) {
  ## reduce to 5 categories: White, South Asian, Chinese, Black, Other
  ethnic7 <- car::recode(onomap.group, "'Britain&Ireland'='White';  'Other Europe'='White';  'East Med'='Other'; 'Other Asia & Pacific'='Other'") 
  ethnic6 <- as.character(ethnic7)
  ethnic6 <- car::recode(ethnic6, "'Black Caribbean'='Other'")
  ethnic6 <- car::recode(ethnic6, "'Black African'='Other'")
  ethnic6 <- as.factor(ethnic6)
  #print(levels(ethnic6))
  ethnic6 <- factor(ethnic6, levels=levels(ethnic6)[c(5,6,1,2,4,3)])
  ethnic6 <- relevel(as.factor(ethnic6), ref="White")
  return(ethnic6)
}
