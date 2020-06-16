
## move to main script
cc.all$dispensing.days <- as.integer(cc.all$SPECDATE - as.Date("2019-06-01"))

getdaysbeforetest <- function(scrips.in, specdates) {
    ## merge SPECDATE into scrips table and add columns daysbefore, dispensingdays, intervalTWday
    ## checked that minimum value is 16 days -- the 15-day cutoff has been applied to the scrips table
    scrips <- merge(scrips.in, specdates, by="ANON_ID", all.x=TRUE) 
    scrips$daysbefore <- as.integer(scrips$SPECDATE - scrips$dispensed_date)
    scrips$dispensing.days <- as.integer(scrips$SPECDATE - 15 - 
                                         as.Date("2019-06-01"))
    scrips$intervalTWday <- ceiling((scrips$daysbefore - 15) / TW)
    scrips$intervalTWday[scrips$intervalTWday > 3] <- 3
    return(scrips)
}

getdoseTWday.wide <- function(scrips.in, varname.prefix="dose.") {
    doseTWday <-  reshape2::dcast(data=scrips.in,
                                  formula=ANON_ID + intervalTWday ~ approved_name,
                                  value.var="dose",
                                  fun.aggregate=sum)
    ## strictly, each interval should be divided by the number of days in that interval for that individual
    ## for now, use TW days
    for(j in 1:5) { # loop over approved names to divide by DDD x TW 
        doseTWday[, j + 2] <- doseTWday[, j + 2] / TW
    }
    doseTWday$doseTWday.all <- rowSums(doseTWday[, -(1:2)])
    ## drop columns for individual drugs, and cast again to get one column for each interval
    doseTWday <- doseTWday[, c("ANON_ID", "intervalTWday", "doseTWday.all")]
    doseTWday <- doseTWday[!is.na(doseTWday$intervalTWday), ]
    doseTWday.wide <-  reshape2::dcast(data=doseTWday,
                                       formula=ANON_ID ~ intervalTWday,
                                       value.var="doseTWday.all",
                                       fun.aggregate=sum)
    colnames(doseTWday.wide)[-1] <-
        paste0(varname.prefix,
               c("interval1", "interval2", "interval3"))
    return(doseTWday.wide)
}

## reading scrips

scripvars <- c("ANON_ID", "dispensed_date", "bnf_paragraph_code",
               "formulation_code",
               "item_strength", "item_strength_uom", "item_code",
               "num_items", "quantity")
## lookup table for subpara code and item code
drugvars <- c("bnf_paragraph_code", "bnf_paragraph_description",
              "item_code", "approved_name")

##############################################################

if(FALSE) { # coding antihypertensives
antihypertensive.classes <- c("vasodilator_ah6.bnf", 
                              "centrally_acting_ah6.bnf",  
                              "adrenergic_neurone_block6.bnf",  
                              "alpha_adrenoceptor6.bnf", 
                              "ace6.bnf",  
                              "angio6.bnf",  
                              "renin_angiotensin6.bnf",  
                              "thiazides6.bnf",  
                              "calcium_channel6.bnf")
## first three classes are rarely prescribed
}

#################### add fields to scrips ########################################

## for now, just keep id, paragraph code, paragraph_description

scrips$bnf_paragraph_description <- as.factor(scrips$bnf_paragraph_description)
cat("scrips object uses", object.size(scrips) * 1E-6, "MB\n")

## we have 7 digits on scrips, giving resolution to subpara level only
length(table(substr(scrips$bnf_paragraph_code, 1, 2))) # chapter
length(table(substr(scrips$bnf_paragraph_code, 1, 4))) # chapter, section
length(table(substr(scrips$bnf_paragraph_code, 1, 6)))  # chapter, section, paragraph
length(table(scrips$bnf_paragraph_code)) # 537 groups

## we need integer variables chapter, sectioncode, paracode for use with reshape2::dcast
scrips$chapternum <- as.integer(substr(scrips$bnf_paragraph_code, 1, 2))
scrips$sectioncode <- as.integer(substr(scrips$bnf_paragraph_code, 1, 4))
scrips$paracode <- as.integer(substr(scrips$bnf_paragraph_code, 1, 6))

## recode scrips$bnf.chapter values > 14 or NA to 14
scrips$chapternum[is.na(scrips$chapternum)] <- 14
scrips$chapternum[scrips$chapternum > 14] <- 14

## separate tables for prescriptions of proton pump, opioids, nonopioids

scrips.protonpump <-
    readRDS(scrips.filename) %>%
    subset(bnf_paragraph_code=="0103050")

## https://www.whocc.no/atc_ddd_index/?code=A02BC&showdescription=yes gives defined daily doses
##
##ATC code  	Name  	DDD 	 U 	Adm.R	 Note
##A02BC05 	esomeprazole 	30 	mg 	O
##A02BC03 	lansoprazole 	30 	mg 	O 
##A02BC01 	omeprazole 	20 	mg 	O 
##A02BC02 	pantoprazole 	40 	mg 	O 	
##A02BC04 	rabeprazole 	20 	mg 	O 	

DDD.protonpump <- data.frame(approved_name=toupper(c("esomeprazole",
                                             "lansoprazole",
                                             "omeprazole",
                                             "pantoprazole",
                                             "rabeprazole")), 
                             DDD=c(30, 30, 20, 40, 20))
scrips.protonpump <- merge(scrips.protonpump, DDD.protonpump, by="approved_name", all.x=TRUE)

scrips.opioid <-
    readRDS(scrips.filename) %>%
    subset(bnf_paragraph_code=="0407020")

scrips.opioid$mode <- car::recode(scrips.opioid$formulation_code,
                                  "c('CAPS', 'LOZ', 'SACH', 'SOLN', 'TABS')='Oral';
                                    'PATCH'='Patch'; 'INJ'='Injection';")

table(scrips.opioid$approved_name, scrips.opioid$mode)
## only buprenorphine and fentanyl are used as patch formulations
print(table(scrips.opioid$approved_name, scrips.opioid$mode))

## item strength missing for patch formulations
subset(scrips.opioid, subset=mode=="Patch") %>% summary()
#######################################################

opiates <- subset(bnfchemicalcodes, subset=subparacode==407020)
opiates <- as.data.frame(opiates[!duplicated(opiates$chemicalname), ])

## create column in BNF table to match name used in PIS table
opiates$PISname <- opiates$chemicalname
opiates$PISname <- gsub(" \\(Systemic\\)", "", opiates$PISname)
opiates$PISname <- gsub("Oxycodone hydrochloride/naloxone hydrochloride",
                             "Oxycodone and naloxone", opiates$PISname)
opiates$PISname <- gsub("Morphine sulfate", "Morphine", opiates$PISname)
opiates$PISname <- gsub("Meptazinol hydrochloride", "Meptazinol", opiates$PISname)
opiates$PISname <- gsub("Tapentadol hydrochloride", "Tapentadol", opiates$PISname)

## cannot match: dipipanone with cyclizine (rare), paracetamol with tramadol (70 scrips)

opiates$PISname <- toupper(opiates$PISname)

scrips.opiates.names <- data.frame(approved_name=unique(scrips.opioid$approved_name))

scrips.opiates.names <- merge(scrips.opiates.names, opiates,
                              by.x="approved_name", by.y="PISname",
                              all.x=TRUE)
scrips.opiates.names <- subset(scrips.opiates.names,
                               select=c("approved_name", "chemicalname", "chemicalcode"))

freqs <- table(scrips.opioid$approved_name)
freqs <- data.frame(approved_name=names(freqs), numscrips=as.integer(freqs))

scrips.opiates.names <- merge(scrips.opiates.names, freqs)
conversionMME <- 
"
Opioid 	Conversion Factor 
Codeine 	0.15
Fentanyl	2.4
Hydrocodone 	1
Hydromorphone 	4
Morphine 	1
Oxycodone 	1.5
Oxymorphone 	3
Methadone hydrochloride
Morphine sulphate
Pentazocine hydrochloride
Pentazocine lactate
Pethidine hydrochloride
Nalbuphine hydrochloride
"

## https://fpm.ac.uk/opioids-aware-structured-approach-opioid-prescribing/dose-equivalents-and-changing-opioids

scrips.opiates.names$factor.oral <- NA
scrips.opiates.names$factor.oral[scrips.opiates.names$approved_name=="CODEINE PHOSPHATE"] <- 0.1
scrips.opiates.names$factor.oral[scrips.opiates.names$approved_name=="DIAMORPHINE HYDROCHLORIDE"] <- 3
scrips.opiates.names$factor.oral[scrips.opiates.names$approved_name=="DIHYDROCODEINE TARTRATE"] <- 0.1
scrips.opiates.names$factor.oral[scrips.opiates.names$approved_name=="FENTANYL"] <- 2.4
scrips.opiates.names$factor.oral[scrips.opiates.names$approved_name=="HYDROMORPHONE HYDROCHLORIDE"] <- 7.5
scrips.opiates.names$factor.oral[scrips.opiates.names$approved_name=="MEPTAZINOL"] <- 6/200
scrips.opiates.names$factor.oral[scrips.opiates.names$approved_name=="METHADONE HYDROCHLORIDE"] <- 4
scrips.opiates.names$factor.oral[scrips.opiates.names$approved_name=="MORPHINE"] <- 1
scrips.opiates.names$factor.oral[scrips.opiates.names$approved_name=="OXYCODONE"] <- 2
scrips.opiates.names$factor.oral[scrips.opiates.names$approved_name=="OXYCODONE AND NALOXONE"] <- 2 # ignore the naloxone 
scrips.opiates.names$factor.oral[scrips.opiates.names$approved_name=="PARACETAMOL WITH TRAMADOL HYDROCHLORIDE"] <- 0.15
scrips.opiates.names$factor.oral[scrips.opiates.names$approved_name=="PETHIDINE HYDROCHLORIDE"] <- 0.12
scrips.opiates.names$factor.oral[scrips.opiates.names$approved_name=="TAPENTADOL"] <- 0.4
scrips.opiates.names$factor.oral[scrips.opiates.names$approved_name=="TRAMADOL HYDROCHLORIDE"] <- 0.15

## item strength of methadone hydrochloride is always 5, CDC gives conversion factor of 4 for dose up to 20 mg/day
## https://www.cdc.gov/drugoverdose/pdf/calculating_total_daily_dose-a.pdf

## oxycodone with nalogone comes in strengths 5, 10, 20, 40 mg oxyxcodone
scrips.opioid[scrips.opioid$approved_name=="OXYCODONE AND NALOXONE",
              "item_strength"] <- 20

## paracetamol with tramadol hydrochloride has strength 37.5 mg tramadol
scrips.opioid[scrips.opioid$approved_name=="PARACETAMOL WITH TRAMADOL HYDROCHLORIDE",
              "item_strength"] <- 37.5

## transdermal buprenorphine: assume patch strength 10 microgram/h changed weekly, would be equivalent to 24 mg/day morphine over 7 days.
## so conversion factor is 1 patch = 24 * 7 
## transdermal fentanyl: assume patch strength 50 microgram/h changed every 3 days, frequency not given, equivalent to 180 mg/day morphine
## so conversion factor is 1 patch = 180 * 3

scrips.opiates.names$factor.patch <- NA
scrips.opiates.names$factor.patch[scrips.opiates.names$approved_name=="BUPRENORPHINE"] <- 24 * 7
scrips.opiates.names$factor.patch[scrips.opiates.names$approved_name=="FENTANYL"] <- 180 * 3

print(subset(scrips.opiates.names,
             select=c("approved_name", "numscrips", "factor.oral", "factor.patch")))

scrips.opioid <- merge(scrips.opioid,
                       subset(scrips.opiates.names,
                              select=c("approved_name", "factor.oral", "factor.patch")),
                       by="approved_name", all.x=TRUE)

##### compound opiates  #######################################################

nonopiates <-  subset(bnfchemicalcodes, subset=subparacode==407010)
nonopiates <- nonopiates[!duplicated(nonopiates$chemicalname), ]
compound.opiates <-
    as.data.frame(subset(nonopiates,
                         subset=substr(chemicalcode, 8, 9) %in%
                             c("A0", "F0", "M0", "N0", "Q0")))

scrips.compound.opiates <-
    readRDS(scrips.filename) %>%
    subset(bnf_paragraph_code=="0407010")

scrips.compound.opiates <- subset(scrips.compound.opiates,
                                  subset=scrips.compound.opiates$approved_name %in%
                                      names(table(scrips.compound.opiates$approved_name))[2:4])
print(table(scrips.compound.opiates$approved_name))

## co-proxamol does not appear in scrips table

## Co-codamol tablets and capsules come in 3 different strengths. They contain 8mg, 15mg or 30mg of codeine. The higher strengths (15/500 and 30/500) are only available on prescription from a doctor.
## for now  assign strength as 30

## Co-codamol with buclizine hydrochloride is 8 mg strength

## Co-dydramol comes in 4 different strengths. The tablets contain either 7.46mg, 10mg, 20mg or 30mg of dihydrocodeine.
## For now assign strength as 30

scrips.compound.opiates$item_strength[scrips.compound.opiates$approved_name=="CO-CODAMOL"] <- 15
scrips.compound.opiates$item_strength[scrips.compound.opiates$approved_name=="CO-CODAMOL WITH BUCLIZINE HYDROCHLORIDE"] <- 8
scrips.compound.opiates$item_strength[scrips.compound.opiates$approved_name=="CO-DYDRAMOL"] <- 20

scrips.compound.opiates$factor.oral <- NA
scrips.compound.opiates$factor.oral[scrips.compound.opiates$approved_name=="CO-CODAMOL"] <- 0.1
scrips.compound.opiates$factor.oral[scrips.compound.opiates$approved_name=="CO-CODAMOL WITH BUCLIZINE HYDROCHLORIDE"] <- 0.1
scrips.compound.opiates$factor.oral[scrips.compound.opiates$approved_name=="CO-DYDRAMOL"] <- 0.1

table(scrips.compound.opiates$approved_name, scrips.compound.opiates$item_strength, exclude=NULL)
  
######################################################################
########## all merges with case-control dataset are after this point ############# 

TW <- 120 # time window in days

####### proton pump exposure   ##################

scrips.protonpump <- getdaysbeforetest(scrips.in=scrips.protonpump,
                                       specdates=cc.all[, c("ANON_ID", "SPECDATE")])

## standardize by converting to DDDs
scrips.protonpump$dose <- scrips.protonpump$item_strength * scrips.protonpump$quantity /
    scrips.protonpump$DDD

## calculate dose of each drug over entire period
dose.protonpump <-  reshape2::dcast(data=scrips.protonpump,
                                    formula=ANON_ID + dispensing.days ~ approved_name,
                                    value.var="dose",
                                    fun.aggregate=sum)

for(j in 3:ncol(dose.protonpump)) { # loop over approved names to divide by dispensing.days 
    dose.protonpump[, j] <- dose.protonpump[, j] / dose.protonpump$dispensing.days 
}

dose.protonpump$DDDs.all <- rowSums(dose.protonpump[, -(1:2)])

cc.all <- merge(cc.all, dose.protonpump, by="ANON_ID", all.x=TRUE)
cc.all$DDDs.all[is.na(cc.all$DDDs.all)] <- 0
cc.all$ESOMEPRAZOLE[is.na(cc.all$ESOMEPRAZOLE)] <- 0
cc.all$LANSOPRAZOLE[is.na(cc.all$LANSOPRAZOLE)] <- 0
cc.all$OMEPRAZOLE[is.na(cc.all$OMEPRAZOLE)] <- 0
cc.all$PANTOPRAZOLE[is.na(cc.all$PANTOPRAZOLE)] <- 0
cc.all$RABEPRAZOLE[is.na(cc.all$RABEPRAZOLE)] <- 0

cc.all$DDDsgr <- 0.5 * ceiling(2 * cc.all$DDDs.all)
cc.all$DDDsgr <- as.factor(car::recode(cc.all$DDDsgr, "2:hi='2 or more'"))

###############  average dose in each TW-day interval ######################## 

doseTWday.protonpump.wide <- getdoseTWday.wide(scrips.protonpump, varname.prefix="DDD.")

cc.all <- merge(cc.all, doseTWday.protonpump.wide, by="ANON_ID", all.x=TRUE)
cc.all$DDD.interval1[is.na(cc.all$DDD.interval1)] <- 0
cc.all$DDD.interval2[is.na(cc.all$DDD.interval2)] <- 0
cc.all$DDD.interval3[is.na(cc.all$DDD.interval3)] <- 0

################################################

####### opiate exposure   ##################

## merge SPECDATE to get days before test and time window
scrips.opiate <- getdaysbeforetest(scrips.in=scrips.opioid,
                                       specdates=cc.all[, c("ANON_ID", "SPECDATE")])

## calculate oral dose as morphine equivalent
scrips.opiate$dose <- scrips.opiate$item_strength * scrips.opiate$quantity *
    scrips.opiate$factor.oral
## separate calculation for patch dose
mode.patch <- scrips.opiate$mode=="Patch"
scrips.opiate$dose[mode.patch] <- scrips.opiate$quantity[mode.patch] *
    scrips.opiate$factor.patch[mode.patch]

## calculate dose of each drug over entire period
dose.opiate <-  reshape2::dcast(data=scrips.opiate,
                                    formula=ANON_ID + dispensing.days ~ approved_name,
                                    value.var="dose",
                                    fun.aggregate=sum)
for(j in 3:ncol(dose.opiate)) { # loop over approved names to divide by dispensing.days 
    dose.opiate[, j] <- dose.opiate[, j] / dose.opiate$dispensing.days
}
dose.opiate$dose.opiate.daily <- rowSums(dose.opiate[, -(1:2)])

quantile(dose.opiate$dose.opiate.daily, probs=seq(0, 1, by=0.1), na.rm=TRUE)

## https://www.cdc.gov/drugoverdose/pdf/calculating_total_daily_dose-a.pdf
## quotes a study showing that dose >= 50 MME/day doubles risk of overdose
## this is about the 90th centile in ever-exposed Scots

cc.all <- merge(cc.all, dose.opiate[, c("ANON_ID", "dose.opiate.daily")], by="ANON_ID", all.x=TRUE)
## group dose into categories
cc.all$dose.opiate.daily[is.na(cc.all$dose.opiate.daily)] <- 0
cc.all$dosegr.opiate <- 10 * ceiling(cc.all$dose.opiate.daily / 10)
cc.all$dosegr.opiate <- as.factor(car::recode(cc.all$dosegr.opiate,
                                              "10:20='1-20'; 30:50='21-50'; 60:hi='>50'"))
cc.all$dosegr.opiate <- factor(cc.all$dosegr.opiate,
                               levels=levels(cc.all$dosegr.opiate)[c(2:4, 1)])

doseTWday.opiate.wide <- getdoseTWday.wide(scrips.opiate, varname.prefix="opiateMME.")
cc.all <- merge(cc.all, doseTWday.opiate.wide, by="ANON_ID", all.x=TRUE)
cc.all$opiateMME.interval1[is.na(cc.all$opiateMME.interval1)] <- 0
cc.all$opiateMME.interval2[is.na(cc.all$opiateMME.interval2)] <- 0
cc.all$opiateMME.interval3[is.na(cc.all$opiateMME.interval3)] <- 0

####################################################################

## compound opiate exposure

## merge SPECDATE to get days before test and time window
scrips.compound.opiates <- getdaysbeforetest(scrips.in=scrips.compound.opiates,
                                       specdates=cc.all[, c("ANON_ID", "SPECDATE")])

## calculate oral dose as morphine equivalent
scrips.compound.opiates$dose <- scrips.compound.opiates$item_strength *
    scrips.compound.opiates$quantity *
    scrips.compound.opiates$factor.oral

## calculate dose of each drug over entire period
dose.compound.opiates <-  reshape2::dcast(data=scrips.compound.opiates,
                                    formula=ANON_ID + dispensing.days ~ approved_name,
                                    value.var="dose",
                                    fun.aggregate=sum)
for(j in 3:ncol(dose.compound.opiates)) { # loop over approved names to divide by dispensing.days 
    dose.compound.opiates[, j] <- dose.compound.opiates[, j] / dose.compound.opiates$dispensing.days
}
dose.compound.opiates$dose.compound.opiates.daily <- rowSums(dose.compound.opiates[, -(1:2)])

quantile(dose.compound.opiates$dose.compound.opiates.daily, probs=seq(0, 1, by=0.1), na.rm=TRUE)

cc.all <- merge(cc.all, dose.compound.opiates[, c("ANON_ID", "dose.compound.opiates.daily")], by="ANON_ID", all.x=TRUE)
cc.all$dose.compound.opiates.daily[is.na(cc.all$dose.compound.opiates.daily)] <- 0
cc.all$dosegr.compound.opiates <- 5 * ceiling(cc.all$dose.compound.opiates.daily / 5)
cc.all$dosegr.compound.opiates <- as.factor(car::recode(cc.all$dosegr.compound.opiates,
                                              "5:10='1-10'; 15:hi='>10'"))
cc.all$dosegr.compound.opiates <- factor(cc.all$dosegr.compound.opiates,
                               levels=levels(cc.all$dosegr.compound.opiates)[c(2:4, 1)])
table(cc.all$dosegr.compound.opiates)
                                        #
########### other drugs of interest coded as binary ####################

ids.antiplatelet <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 209])
cc.all$antiplatelet <- as.factor(as.integer(cc.all$ANON_ID %in% ids.antiplatelet))

ids.nsaid <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 100101])
cc.all$nsaid <- as.factor(as.integer(cc.all$ANON_ID %in% ids.nsaid))

ids.opioid.analgesic <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 40702])
cc.all$opioid.analgesic <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                      ids.opioid.analgesic))

ids.nonopioid.analgesic <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 40701])
cc.all$nonopioid.analgesic <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                   ids.nonopioid.analgesic))

ids.antipsychotic <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 40201])
cc.all$antipsychotic <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                      ids.antipsychotic))

ids.osmotic.laxative <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 10604])
cc.all$osmotic.laxative <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                ids.osmotic.laxative))

ids.anticoagulant.any <-
    unique(scrips$ANON_ID[scrips$bnf_paragraph_code == "0208010" |
                          scrips$bnf_paragraph_code == "0208020"])
cc.all$anticoagulant.any <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                ids.anticoagulant.any))

ids.hydroxychloroquine <- readRDS(scrips.filename) %>%
    subset(subset=approved_name=="HYDROXYCHLOROQUINE SULFATE", select=ANON_ID)
ids.hydroxychloroquine <- unique(ids.hydroxychloroquine$ANON_ID)
cc.all$hydroxychloroquine <- as.factor(as.integer(cc.all$ANON_ID %in% ids.hydroxychloroquine))

ids.protonpump <- unique(scrips$ANON_ID[scrips$bnf_paragraph_code == "0103050"])
cc.all$protonpump <- as.factor(as.integer(cc.all$ANON_ID %in% ids.protonpump))
cc.all$y.protonpump <- as.integer(cc.all$protonpump =="1")

###############################################################
