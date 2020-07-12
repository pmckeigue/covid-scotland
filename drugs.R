
get.timewindow <- function(scrips.in) {
    ## add columns daysbefore, intervalTWday, dispensing.days
    ## the 15-day cutoff has been applied to the scrips table
    scrips.in$intervalTWday <- ceiling((scrips.in$daysbefore - 15) / TW)
    scrips.in$intervalTWday[scrips.in$intervalTWday > 2] <- 2
    scrips.in$dispensing.days <- 240  # same for all individuals
    return(scrips.in)
}

getdoseTWday.wide <- function(scrips.in, varname.prefix="dose.") {
    doseTWday <-  data.table::dcast(data=scrips.in,
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
    doseTWday.wide <-  data.table::dcast(data=doseTWday,
                                       formula=ANON_ID ~ intervalTWday,
                                       value.var="doseTWday.all",
                                       fun.aggregate=sum)
    colnames(doseTWday.wide)[-1] <-
        paste0(varname.prefix,
               c("interval1", "interval2"))
    return(doseTWday.wide)
}

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

## we need integer variables chapter, sectioncode, paracode for use with data.table::dcast
scrips[, chapternum := as.integer(substr(bnf_paragraph_code, 1, 2))]
scrips[, sectioncode := as.integer(substr(bnf_paragraph_code, 1, 4))]
scrips[, paracode := as.integer(substr(bnf_paragraph_code, 1, 6))]

## recode scrips$bnf.chapter values > 14 or NA to 14
scrips[is.na(chapternum), chapternum := 14]
scrips[chapternum > 14, chapternum := 14] 

if(!old) {
    
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
#print(table(scrips.opioid$approved_name, scrips.opioid$mode))

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

    scrips.opioid <- merge(scrips.opioid,
                       subset(scrips.opiates.names,
                              select=c("approved_name", "factor.oral", "factor.patch")),
                       by="approved_name", all.x=TRUE)

##### compound opiates  #######################################################

    nonopiates.orcompound <-  subset(bnfchemicalcodes, subset=subparacode==407010)
    nonopiates.orcompound <-
        nonopiates.orcompound[!duplicated(nonopiates.orcompound$chemicalname), ]
    compound.opiates <-
        as.data.frame(subset(nonopiates.orcompound,
                         subset=substr(chemicalcode, 8, 9) %in%
                             c("A0", "F0", "M0", "N0", "Q0")))
    nonopiates <-
        as.data.frame(subset(nonopiates.orcompound,
                             subset=!(substr(chemicalcode, 8, 9) %in%
                                      c("A0", "F0", "M0", "N0", "Q0"))))
    scrips.compound.opiates <-
        readRDS(scrips.filename) %>%
        subset(bnf_paragraph_code=="0407010") %>% 
        subset(substring(bnf_item_code, 1, 9) %in% compound.opiates$chemicalcode)

    scrips.nonopiates <-
        readRDS(scrips.filename) %>%
        subset(bnf_paragraph_code=="0407010") %>% 
        subset(!(substring(bnf_item_code, 1, 9) %in% compound.opiates$chemicalcode))
    
    scrips.compound.opiates <- as.data.table(scrips.compound.opiates)
    scrips.compound.opiates <-
        subset(scrips.compound.opiates,
               subset=scrips.compound.opiates$approved_name %in%
                   names(table(scrips.compound.opiates$approved_name))[2:4])
    
item.desc <- subset(scrips.compound.opiates, approved_name=="CO-DYDRAMOL",
                    select="bnf_item_description") %>% table() %>% names()
     
#print(table(scrips.compound.opiates$approved_name))

## co-proxamol does not appear in scrips table
## Co-codamol tablets and capsules come in 3 different strengths. They contain 8mg, 15mg or 30mg of codeine. The higher strengths (15/500 and 30/500) are only available on prescription from a doctor.
## Co-codamol with buclizine hydrochloride is 8 mg strength
## Co-dydramol comes in 4 different strengths. The tablets contain either 7.46mg, 10mg, 20mg or 30mg of dihydrocodeine.
## REMEDEINE FTE_TAB is 30 mg, REMEDEINE_TAB is 20 mg            

scrips.compound.opiates[approved_name=="CO-CODAMOL",
                        item_strength := as.numeric(gsub("(.+ )([0-9.]+)(MG\\/.+)",
                                                         "\\2", bnf_item_description))]

scrips.compound.opiates[approved_name=="CO-CODAMOL WITH BUCLIZINE HYDROCHLORIDE",
                        item_strength := 8]
scrips.compound.opiates[approved_name=="CO-DYDRAMOL",
                        item_strength := as.numeric(gsub("(.+ )([0-9.]+)(MG\\/.+)",
                                                         "\\2", bnf_item_description))]
scrips.compound.opiates[bnf_item_description=="REMEDEINE_FTE_TAB",
                        item_strength :=30]
scrips.compound.opiates[bnf_item_description=="REMEDEINE_TAB",
                        item_strength :=20]
 
mutate(scrips.compound.opiates, factor.oral = NA)
    scrips.compound.opiates[approved_name=="CO-CODAMOL", factor.oral := 0.1]
scrips.compound.opiates[approved_name=="CO-CODAMOL WITH BUCLIZINE HYDROCHLORIDE", factor.oral := 0.1]
scrips.compound.opiates[approved_name=="CO-DYDRAMOL", factor.oral := 0.1]

with(scrips.compound.opiates, table(approved_name, item_strength, exclude=NULL))

######################################################################
########## all merges with case-control dataset are after this point ############# 
cat("Merging with case-control table ...\n")

TW <- 120 # time window in days

subparas.laporte <-
    c("subpara.103050.Proton pump inhibitors",
      "subpara.102000.Antispasmodic and other drugs altering gut motility",       
      "subpara.304010.Antihistamines",                                    
      "subpara.401010.Hypnotics",                                  
      "subpara.401020.Anxiolytics",                                
      "subpara.402010.Antipsychotic drugs",                               
      "subpara.402030.Drugs used for mania and hypomania",                         
      "subpara.403010.Tricyclic and related antidepressant drugs",                 
      "subpara.403030.Selective serotonin re-uptake inhibitors",              
      "subpara.403040.Other antidepressant drugs"   ,              
      "subpara.406000.Drugs used in nausea and vertigo",             
      "subpara.407020.Opioid analgesics",            
      "gabapentinoids", #"subpara.408010.Control of epilepsy",           
      "subpara.409020.Antimuscarinic drugs used in parkinsonism",

      ## mirabegron is not anticholinergic
      ## trospium, darifenacin are claimed to be selective for the bladder
                                        # "urinary.antispasmodics",
      "subpara.704020.Drugs for urinary frequency enuresis and incontinence",      
      "subpara.1001010.Non-steroidal anti-inflammatory drugs"                     
      )

    
####### proton pump exposure   ##################

scrips.protonpump <- get.timewindow(scrips.in=scrips.protonpump)

## standardize by converting to DDDs
scrips.protonpump$dose <- scrips.protonpump$item_strength * scrips.protonpump$quantity /
    scrips.protonpump$DDD

## calculate dose of each drug over entire period
dose.protonpump <-  data.table::dcast(data=scrips.protonpump,
                                    formula=ANON_ID + dispensing.days ~ approved_name,
                                    value.var="dose",
                                    fun.aggregate=sum)

for(j in 3:ncol(dose.protonpump)) { # loop over approved names to divide by dispensing.days 
    dose.protonpump[[j]] <- dose.protonpump[[j]] / dose.protonpump$dispensing.days 
}
dose.protonpump$DDDs.all <- rowSums(dose.protonpump[, -(1:2)])
dose.protonpump <- as.data.table(dose.protonpump, key="ANON_ID")
## shouldn't need this line if left join works correctly
dose.protonpump <- subset(dose.protonpump, ANON_ID %in% cc.all$ANON_ID)
dose.protonpump[, dispensing.days:= NULL]

setkey(cc.all, ANON_ID)
cc.all <- dose.protonpump[cc.all]

cc.all[, `:=`(DDDs.all = coalesce(DDDs.all, 0), 
              ESOMEPRAZOLE = coalesce(ESOMEPRAZOLE, 0),
              LANSOPRAZOLE = coalesce(LANSOPRAZOLE, 0),
              OMEPRAZOLE = coalesce(OMEPRAZOLE, 0),
              PANTOPRAZOLE = coalesce(PANTOPRAZOLE, 0),
              RABEPRAZOLE = coalesce(RABEPRAZOLE, 0))]

cc.all[, DDDsgr := 0.5 * ceiling(2 * DDDs.all)]
cc.all[, DDDsgr := as.factor(car::recode(DDDsgr, "2:hi='2 or more'"))]

###############  average dose in each TW-day interval ######################## 

doseTWday.protonpump.wide <- getdoseTWday.wide(scrips.protonpump, varname.prefix="DDD.")
doseTWday.protonpump.wide <- as.data.table(doseTWday.protonpump.wide, key="ANON_ID")
doseTWday.protonpump.wide <- subset(doseTWday.protonpump.wide, ANON_ID %in% cc.all$ANON_ID)

setkey(cc.all, ANON_ID)
cc.all <- doseTWday.protonpump.wide[cc.all] 

cc.all[, `:=`(DDD.interval1 = coalesce(DDD.interval1, 0),
              DDD.interval2 = coalesce(DDD.interval2, 0))]
           
protonpump.exposure.nonrecent <- as.integer(cc.all$DDD.interval2 > 0)
protonpump.exposure.recent <- as.integer(cc.all$DDD.interval1 > 0)

cc.all[, protonpump.exposuregr := as.integer(DDDs.all > 0)]
cc.all[protonpump.exposure.nonrecent==1 & protonpump.exposure.recent==0, protonpump.exposuregr := 1]
cc.all[protonpump.exposure.nonrecent==0 & protonpump.exposure.recent==1, protonpump.exposuregr := 2]
cc.all[protonpump.exposure.nonrecent==1 & protonpump.exposure.recent==1, protonpump.exposuregr := 3]
protonpump.exposure.cat <- as.factor(car::recode(cc.all[["protonpump.exposuregr"]],
                               "0='No prescriptions'; 1='Non-recent only';
                                2='Recent only';
                                3='Prescriptions in both time windows'"))
protonpump.exposure.cat <- factor(protonpump.exposure.cat,
                                  levels=levels(protonpump.exposure.cat)[c(2, 4, 3, 1)])
## problem with changing the type of an existing column
cc.all[, protonpump.exposurecat := protonpump.exposure.cat]

####################################################################

## compound opiate exposure

## merge SPECDATE to get days before test and time window
scrips.compound.opiates <- get.timewindow(scrips.in=scrips.compound.opiates)

## calculate oral dose as morphine equivalent
scrips.compound.opiates$dose <- scrips.compound.opiates$item_strength *
    scrips.compound.opiates$quantity *
    scrips.compound.opiates$factor.oral

## calculate dose of each drug over entire period
dose.compound.opiates <-  data.table::dcast(data=scrips.compound.opiates,
                                    formula=ANON_ID + dispensing.days ~ approved_name,
                                    value.var="dose",
                                    fun.aggregate=sum)
for(j in 3:ncol(dose.compound.opiates)) { # loop over approved names to divide by dispensing.days 
    dose.compound.opiates[[j]] <- dose.compound.opiates[[j]] / dose.compound.opiates$dispensing.days
}
dose.compound.opiates$dose.compound.opiates.daily <- rowSums(dose.compound.opiates[, -(1:2)])

quantile(dose.compound.opiates$dose.compound.opiates.daily, probs=seq(0, 1, by=0.1), na.rm=TRUE)

dose.compound.opiates <- subset(dose.compound.opiates,
                                select=c(ANON_ID, dose.compound.opiates.daily))
dose.compound.opiates <- as.data.table(dose.compound.opiates, key="ANON_ID") 
cc.all <- dose.compound.opiates[cc.all]

cc.all[, dose.compound.opiates.daily := coalesce(dose.compound.opiates.daily, 0)]

cc.all$dosegr.compound.opiates <- 5 * ceiling(cc.all$dose.compound.opiates.daily / 5)
cc.all$dosegr.compound.opiates <- as.factor(car::recode(cc.all$dosegr.compound.opiates,
                                              "5:10='1-10'; 15:hi='>10'"))
cc.all[, dosegr.compound.opiates := factor(dosegr.compound.opiates,
                                           levels=levels(dosegr.compound.opiates)[c(2:4, 1)])]
table(cc.all$dosegr.compound.opiates)

####### opiate exposure   ##################

scrips.opiate <- get.timewindow(scrips.in=scrips.opioid)

## calculate oral dose as morphine equivalent
scrips.opiate$dose <- scrips.opiate$item_strength * scrips.opiate$quantity *
    scrips.opiate$factor.oral
## separate calculation for patch dose
mode.patch <- scrips.opiate$mode=="Patch"
scrips.opiate$dose[mode.patch] <- scrips.opiate$quantity[mode.patch] *
    scrips.opiate$factor.patch[mode.patch]

## calculate dose of each drug over entire period
dose.opiate <-  data.table::dcast(data=scrips.opiate,
                                    formula=ANON_ID + dispensing.days ~ approved_name,
                                    value.var="dose",
                                    fun.aggregate=sum)
for(j in 3:ncol(dose.opiate)) { # loop over approved names to divide by dispensing.days 
    dose.opiate[[j]] <- dose.opiate[[j]] / dose.opiate$dispensing.days
}
dose.opiate$dose.opiate.daily <- rowSums(dose.opiate[, -(1:2)])

quantile(dose.opiate$dose.opiate.daily, probs=seq(0, 1, by=0.1), na.rm=TRUE)

## https://www.cdc.gov/drugoverdose/pdf/calculating_total_daily_dose-a.pdf
## quotes a study showing that dose >= 50 MME/day doubles risk of overdose
## this is about the 90th centile in ever-exposed Scots

dose.opiate <- as.data.table(dose.opiate[, c("ANON_ID", "dose.opiate.daily")], key="ANON_ID")
cc.all <- dose.opiate[cc.all]
cc.all[, dose.opiate.daily := coalesce(dose.opiate.daily, 0)]

##### add opiate dose from compound analgesics
cc.all[, dose.opiate.daily := dose.opiate.daily + dose.compound.opiates.daily]

#####################################################

## group dose into categories
cc.all$dosegr.opiate <- 10 * ceiling(cc.all$dose.opiate.daily / 10)
cc.all$dosegr.opiate <- as.factor(car::recode(cc.all$dosegr.opiate,
                                              "10:20='1-20'; 30:50='21-50'; 60:hi='>50'"))
cc.all[, dosegr.opiate := factor(dosegr.opiate,
                                 levels=levels(dosegr.opiate)[c(2:4, 1)])]

doseTWday.opiate.wide <- getdoseTWday.wide(scrips.opiate, varname.prefix="opiateMME.")
cc.all <- merge(cc.all, doseTWday.opiate.wide, by="ANON_ID", all.x=TRUE)

cc.all[, `:=`(opiateMME.interval1 = coalesce(opiateMME.interval1, 0),
              opiateMME.interval2 = coalesce(opiateMME.interval2, 0))]
           
opiate.exposure.nonrecent <- as.integer(cc.all$opiateMME.interval2 > 0) 
opiate.exposure.recent <- as.integer(cc.all$opiateMME.interval1 > 0)

cc.all <- mutate(cc.all, opiate.exposuregr = as.integer(dose.opiate.daily > 0))
cc.all[opiate.exposure.nonrecent==1 & opiate.exposure.recent==0, opiate.exposuregr := 1]
cc.all[opiate.exposure.nonrecent==0 & opiate.exposure.recent==1, opiate.exposuregr := 2]
cc.all[opiate.exposure.nonrecent==1 & opiate.exposure.recent==1, opiate.exposuregr := 3]
opiate.exposure.cat <- as.factor(car::recode(cc.all[["opiate.exposuregr"]],
                               "0='No prescriptions'; 1='Non-recent only';
                                2='Recent only';
                                3='Prescriptions in both time windows'"))
opiate.exposure.cat <- factor(opiate.exposure.cat,
                             levels=levels(opiate.exposure.cat)[c(2, 4, 3, 1)])

## data.table behaves strangely when type of an existing column is changed
    cc.all <- mutate(cc.all, opiate.exposurecat = opiate.exposure.cat)
    
    ids.opioid.analgesic <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 40702])
    cc.all$opioid.analgesic <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                    ids.opioid.analgesic))
    
    ids.nonopioid.analgesic <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 40701])
    cc.all$nonopioid.analgesic <- as.factor(as.integer(cc.all$ANON_ID %in%
                                                       ids.nonopioid.analgesic))
    
    ## this variable must be class integer to be included in numdrugs
    cc.all <- mutate(cc.all,
                     anyopiate = as.integer(ANON_ID %in%
                                            unique(scrips.compound.opiates$ANON_ID) |
                                            ANON_ID %in% ids.opioid.analgesic))                 
}

########### other drugs of interest coded as binary ####################

ids.antiplatelet <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 209])
cc.all[, antiplatelet := as.factor(as.integer(ANON_ID %in% ids.antiplatelet))]

ids.nsaid <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 100101])
cc.all[, nsaid := as.factor(as.integer(ANON_ID %in% ids.nsaid))]


ids.antipsychotic <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 40201])
cc.all[, antipsychotic := as.factor(as.integer(ANON_ID %in%
                                               ids.antipsychotic))]

ids.osmotic.laxative <- unique(scrips$ANON_ID[as.integer(substr(scrips$bnf_paragraph_code, 1, 6)) == 10604])
cc.all[, osmotic.laxative := as.factor(as.integer(ANON_ID %in%
                                                  ids.osmotic.laxative))]

ids.anticoagulant.any <-
    unique(scrips$ANON_ID[scrips$bnf_paragraph_code == "0208010" |
                          scrips$bnf_paragraph_code == "0208020"])
cc.all[, anticoagulant.any := as.factor(as.integer(ANON_ID %in%
                                                   ids.anticoagulant.any))]

ids.hydroxychloroquine <- unique(scrips[approved_name=="HYDROXYCHLOROQUINE SULFATE", ANON_ID])
cc.all[, hydroxychloroquine := as.factor(as.integer(ANON_ID %in% ids.hydroxychloroquine))]

ids.protonpump <- unique(scrips[bnf_paragraph_code == "0103050", ANON_ID])
cc.all[, protonpump := as.factor(as.integer(ANON_ID %in% ids.protonpump))]
cc.all[, y.protonpump := as.integer(protonpump =="1")]

cc.all[, compound.opiate := as.factor(as.integer(ANON_ID %in%
                                                 scrips.compound.opiates$ANON_ID))]

ids.gabapentinoids <- unique(scrips[approved_name=="GABAPENTIN" |
                                    approved_name=="PREGABALIN", ANON_ID])
cc.all[, gabapentinoids := as.factor(as.integer(ANON_ID %in% ids.gabapentinoids))]

ids.urinary.antispasmodics <- unique(scrips[bnf_paragraph_code=="0704020" & 
                                    approved_name !="MIRABEGRON", ANON_ID])
cc.all[, urinary.antispasmodics := as.factor(as.integer(ANON_ID %in% ids.urinary.antispasmodics))]

###############################################################

## merge drugs, one variable per chapter
scrips.wide <- data.table::dcast(scrips, ANON_ID ~ chapternum, fun.aggregate=length, 
                               value.var="chapternum")
shortnames.cols <-  bnfchapters$shortname[match(as.integer(colnames(scrips.wide)[-1]),
                                                as.integer(bnfchapters$chapternum))]
colnames(scrips.wide)[-1] <- paste("BNF", colnames(scrips.wide)[-1], shortnames.cols,
                                   sep="_")
scrips.wide <- as.data.table(scrips.wide, key="ANON_ID")
cc.all <- scrips.wide[cc.all]

bnfcols <- as.integer(grep("^BNF", colnames(cc.all)))
## recode indicator variables
cc.all[, (bnfcols) := lapply(.SD, recode.indicator), .SDcols = bnfcols]

## merge BNF chapters, one variable per subpara
cat("Merging BNF subparagraph codes ...")
chnums = 1:13
cc.all <- merge.bnfsubparas(chnums=chnums, data=cc.all)
cat("done\n")

subparacols <- grep("^subpara\\.", names(cc.all))
x <- cc.all[,  ..subparacols]
for(j in 1:ncol(x)) set(x, j=j, value=as.integer(x[[j]]) - 1)
head(sapply(x[, 1:5], class))
cc.all[, numdrugs.subpara := rowSums(x)]
cc.all[, numdrugsgr := 3 * ceiling(cc.all$numdrugs.subpara / 3)]
cc.all[, numdrugsgr := as.factor(car::recode(cc.all$numdrugsgr,
                                              "3='1 to 3'; 6='4 to 6'; 9='7 to 9';
                                              12='10 to 12'; 15:hi='>12'"))]
cc.all[, numdrugsgr := factor(numdrugsgr,
                                 levels=levels(numdrugsgr)[c(2, 3, 5, 6, 4, 1)])]

cardiovasc.subparacols <- grep("^subpara\\.2", names(cc.all))
x <- cc.all[,  ..cardiovasc.subparacols]
for(j in names(x)) set(x, j=j, value=as.integer(x[[j]]) - 1) 

cc.all[, numdrugs.cardiovasc := rowSums(x)]
cc.all[, numdrugs.notcardiovasc := numdrugs.subpara - numdrugs.cardiovasc]


x <- cc.all[,  ..subparas.laporte]
for(j in names(x)) set(x, j=j, value=as.integer(x[[j]]) - 1) 
cc.all[, numdrugs.laporte := rowSums(x)]

