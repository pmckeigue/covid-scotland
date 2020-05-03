######## coding ethnicity ##############################

## dataset is named cc.all

## this script should be edited as a function that will work with any dataset, with user-specified names for source variables and output variables

## Source variables: 

## ETHNIC_GROUP_last_last - raw SMR categories
## ETHNIC_smr - collapsed SMR categories derived by Amanda
## OnolyticsType - raw Onomap types based on name classification
## Geographical Area - regional groupings of OnolyticsType

## Output variables: ETHNIC_smr, eth, eth5

## this script

## (1) recodes ETHNIC_smr to specify Chinese as separate category, then collapses to 5 categories: White, South Asian, Chinese, Black, Other

## (2) derives a new variable "eth" from Onomap types

## (3) collapses eth to eth5, which has five categories as above

## We cannot combine SMR coding with Onomap coding because this will introduce non-differential missclassification of ethnicity between cases and controls, especially for those categorised as Black.

## The Onolytics type "MUSLIM" is coded as South Asian - this is reasonable in Scotland where most people with Muslim names are South Asian.



## label the codes used in the raw SMR variable
cc.all$ethnic_smr <- recode(cc.all$ETHNIC_GROUP_last_last,
                            "''=NA; '9'=NA;
'1A'='Scottish';
'1B'='Other British';
'1C'='Irish';
'1K'='Gypsy/ Traveller';
'1L'='Polish';
'1Z'='Any other white ethnic group';
'2A'='Any mixed or multiple ethnic groups';
'3F'='Pakistani, Pakistani Scottish or Pakistani British';
'3G'='Indian, Indian Scottish or Indian British';
'3H'='Bangladeshi, Bangladeshi Scottish or Bangladeshi British';
'3J'='Chinese, Chinese Scottish or Chinese British';
'3Z'='Other Asian, Asian Scottish or Asian British';
'4D'='African, African Scottish or African British';
'4Y'='Other African';
'5C'='Caribbean, Caribbean Scottish or Caribbean British';
'5D'='Black, Black Scottish or Black British';
'5Y'='Other Caribbean or Black';
'6A'='Arab, Arab Scottish or Arab British';
'6Z'='Other ethnic group';
'98'='Refused/Not provided by patient';
'99'='Not Known'")

## ETHNIC_smr should be collapsed to 5 categories, Chinese separate
cc.all$ETHNIC_smr[cc.all$ETHNIC_GROUP_last_last=="3J"] <- "Chinese"
cc.all$ETHNIC_smr <- recode(cc.all$ETHNIC_smr,
                            "''=NA; '9'=NA; 'White Irish'='White'; 'Mixed'='Other'; 'Caribb_Black'='Black'; 'African'='Black'; 'Asian'='Other'")
cc.all <- within(cc.all, ETHNIC_smr <- relevel(as.factor(ETHNIC_smr), ref="White"))

cc.all$OnolyticsType <- recode(cc.all$OnolyticsType,
                               "'NOT FOUND'=NA; 'INTERNATIONAL'=NA; 'UNCLASSIFIED'=NA; 'VOID'=NA; 'VOID - FORENAME'=NA; 'VOID INITIAL'=NA")

table(cc.all$OnolyticsType[cc.all$GeographicalArea=="SOUTH ASIA"])
table(cc.all$OnolyticsType[cc.all$GeographicalArea=="AFRICA"])
table(cc.all$OnolyticsType[cc.all$GeographicalArea=="BRITISH ISLES"])
table(cc.all$OnolyticsType[cc.all$GeographicalArea=="EAST ASIA"])
table(cc.all$OnolyticsType[cc.all$GeographicalArea=="MIDDLE EAST"])

cc.all$eth <- rep("Other", nrow(cc.all))
cc.all$eth[is.na(cc.all$OnolyticsType)] <- NA
cc.all$eth[cc.all$OnolyticsType=="CHINESE" |
           cc.all$OnolyticsType=="HONG KONGESE" |
           cc.all$OnolyticsType=="SINGAPORESE"] <- "Chinese"
cc.all$eth[cc.all$GeographicalArea=="EAST ASIA" &
           cc.all$eth != "Chinese"] <- "Other Asia & Pacific"
           
cc.all$eth[cc.all$GeographicalArea=="SOUTH ASIA"] <- "South Asian"

cc.all$eth[cc.all$GeographicalArea=="BRITISH ISLES"] <- "Britain&Ireland"
cc.all$eth[cc.all$GeographicalArea=="CENTRAL EUROPE" |
           cc.all$GeographicalArea=="EASTERN EUROPE" |
           cc.all$GeographicalArea=="NORTHERN EUROPE" |
           cc.all$GeographicalArea=="SOUTHERN EUROPE" |
           cc.all$OnolyticsType=="AFRIKAANS"] <- "Other Europe"

cc.all$eth[cc.all$GeographicalArea=="AFRICA" &
           cc.all$OnolyticsType != "AFRIKAANS" &
           cc.all$OnolyticsType != "LIBYAN"] <- "Black African"
cc.all$eth[cc.all$GeographicalArea=="MIDDLE EAST" &
           cc.all$OnolyticsType != "MUSLIM"] <- "East Med"
cc.all$eth[cc.all$OnolyticsType == "MUSLIM"] <- "South Asian" # "Muslim, not localized"

cc.all$eth[cc.all$OnolyticsType == "BLACK CARIBBEAN"] <- "Black Caribbean"

## reduce to 5 categories: White, South Asian, Chinese, Black, Other
cc.all$eth5 <- recode(cc.all$eth, "'Black African'='Black'; 'Black Caribbean'='Black';  'Britain&Ireland'='White';  'Other Europe'='White';  'East Med'='Other'; 'Other Asia & Pacific'='Other'; 'Muslim, not localized'='Other'") 
cc.all <- within(cc.all, eth5 <- relevel(as.factor(eth5), ref="White"))

## derive variable based on Onomap only
## this will misclassify some South Asian Muslims as other, and most Black Caribbean as White
cc.all$ethnic <- as.factor(cc.all$eth5) 
