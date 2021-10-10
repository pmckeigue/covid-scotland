library(data.table)
source("helperfunctions.R")

batches <- as.data.table(read.table("data/BatchReference.csv", sep="\t", header=TRUE))
setnames(batches, "Batch", "shield.batch")
batches$shield.batch <- as.integer(gsub("Batch", "", batches$shield.batch))
setkey(batches, shield.batch)
batches$Date.Sent <- gsub("Friday ", "", batches$Date.Sent)
batches$Date.Sent <- gsub("^([[:digit:]]+)([stndrdth]{2} )", "\\1 ", batches$Date.Sent)
batches$Date.Sent <- as.Date(batches$Date.Sent, "%d %B %Y")

linkdate <- "sep22"

if(linkdate=="jul28") {
    datadir <- "./data/2021-07-28/"
    shielding.full.filename <- paste0(datadir,
                                  "CC_shielding_patients_anon_2021-07-28.csv")
    shielded.full <- fread(shielding.full.filename)
} else if(linkdate=="sep02") {
    datadir <- "./data/2021-09-02/"
    shielding.full.filename <- paste0(datadir,
                                  "CC_shielding_patients_anon_2021-09-02.rds")
    shielding.deaths.filename <- paste0(datadir,
                                  "shielding_deaths_anon_2021-09-02.rds")
    shielded.full <- RDStodt(shielding.full.filename)
    shielded.deaths <- RDStodt(shielding.deaths.filename)
    vaxshield.filename <- paste0(datadir, "shielding_vacc_anon_2021-09-02.rds")
    vax.shield <- RDStodt(vaxshield.filename)
} else if(linkdate=="sep22") {
    datadir <- "./data/2021-09-22/"
    shielding.full.filename <- paste0(datadir,
                                  "CC_shielding_patients_anon_2021-09-22.rds")
    shielding.deaths.filename <- paste0(datadir,
                                  "shielding_deaths_anon_2021-09-22.rds")
    shielded.full <- RDStodt(shielding.full.filename)
    shielded.deaths <- RDStodt(shielding.deaths.filename)
    ## FIXME: no current vaxshield table
    vaxshield.filename <- "data/2021-09-02/shielding_vacc_anon_2021-09-02.rds"
    vax.shield <- RDStodt(vaxshield.filename)
}    

setnames(shielded.full, "group", "shield.group", skip_absent=TRUE)
setnames(shielded.full, "batch", "shield.batch", skip_absent=TRUE)
setnames(shielded.full, "age", "age_years", skip_absent=TRUE)
shielded.full[, sex := as.factor(car::recode(sex, "1='Male'; 2='Female'"))]
setkey(shielded.full, shield.batch)
shielded.full <- batches[shielded.full]
shielded.full[, removal := nchar(as.character(removal))]
shielded.full <- unique(shielded.full, by="anon_id") # only 2 duplicated ids

shielded.full[, shield.group := car::recode(shield.group,
                                                "1='Solid organ transplant';
                                      2='Specific cancers';
                                      3='Severe respiratory';
                                      4='Rare diseases';
                                      5='On immunosuppressants';
                                      6='Pregnant with heart disease';
                                      7='Additional conditions'",
                                      as.factor=TRUE,
                                      levels=c(
                                          "Solid organ transplant",
                                         "Specific cancers",
                                         "Severe respiratory",
                                         "Rare diseases",
                                         "On immunosuppressants",
                                         "Pregnant with heart disease",
                                         "Additional conditions"
                                     ))]
## remove obvious wrong assignments
shielded.full <- shielded.full[!(shield.group == "Pregnant with heart disease" &
                                 (sex=="Male" | age_years > 55))]

print(nrow(shielded.full)) # 202510
      
save(shielded.full, file=paste0(datadir, "shielded.full.RData")) # this will be used for a left join of cc.all with shielding cohort

######################################################################

## import deaths into shielded.full
shielded.deaths[, date_of_death := as.Date(date_of_death)]
shielded.deaths <- shielded.deaths[, .(anon_id, date_of_death)]
setkey(shielded.deaths, anon_id)

setkey(shielded.full, anon_id)
shielded.full <- shielded.deaths[shielded.full]

if(!exists("cc.all")) {
    load(paste0(datadir, "cc.all.RData"))
}

## import case status into shielded.full
cases <- cc.all[CASE==1, .(anon_id, specimen_date, casegr, casegr2, casegr3)]
setorder(cases, specimen_date)
cases <- unique(cases, by="anon_id") ## there shouldn't be any duplicates
setkey(cases, anon_id)
setkey(shielded.full, anon_id)
shielded.full <- cases[shielded.full] ## this merge adds 2 more rows
rm(cases)

shielded.full[, casegr := as.character(casegr)]
shielded.full[is.na(casegr), casegr := "Not diagnosed as case"]
shielded.full[, casegr := factor(casegr, levels=c("Not diagnosed as case",
                                                  levels(cc.all$casegr)))]
shielded.full[, casegr2 := as.character(casegr2)]
shielded.full[is.na(casegr2), casegr2 := "Not diagnosed as case"]
shielded.full[, casegr2 := factor(casegr2, levels=c("Not diagnosed as case",
                                                    levels(cc.all$casegr2)))]
shielded.full[, casegr3 := as.character(casegr3)]
shielded.full[is.na(casegr3), casegr3 := "Not diagnosed as case"]
shielded.full[, casegr3 := factor(casegr3, levels=c("Not diagnosed as case",
                                                    levels(cc.all$casegr3)))]
shielded.full[, agegr20 := as.factor(car::recode(as.integer(age_years),
                                                     "0:39='0-39'; 40:59='40-59';
                               60:74='60-74'; 75:hi='75 or more'"))]

shielded.full[, severecase := as.integer(casegr2=="Severe")]
shielded.full[, shield.group := car::recode(shield.group, recode="'Pregnant with heart disease'='Additional conditions'")]
shielded.full[, shield.group := relevel(shield.group, "Severe respiratory")]

## assign entry and exit dates
shielded.full[, exitdate := specimen_date]
shielded.full[is.na(exitdate) & !is.na(date_of_death),
              exitdate := date_of_death]
shielded.full[is.na(exitdate), exitdate := max(shielded.full$specimen_date, na.rm=TRUE)]
shielded.full[, entrydate := as.Date("2020-03-01")]
shielded.full[entrydate >= exitdate, entrydate := NA]

table.casegr.shielded <- with(shielded.full, table(casegr, removal, exclude=NULL))
table.casegr.shielded <- rbind(table.casegr.shielded, colSums(table.casegr.shielded))
rowsums.table.casegr.shielded <- rowSums(table.casegr.shielded)
table.casegr.shielded <- paste.rowpercent(table.casegr.shielded)
table.casegr.shielded <- cbind(table.casegr.shielded, rowsums.table.casegr.shielded)
print(table.casegr.shielded)

## import vaccine data into shielded.full
vax.shield <- vax.shield[vacc_status == "completed"]
setnames(vax.shield, "vacc_occurence_time", "vaxdate")
vax.shield[, vaxdate := as.Date(vaxdate)]
vax.shield[vaxdate < as.Date("2020-12-01"), vaxdate := as.Date(vacc_event_created_at)]

vax.shield <- vax.shield[, .(anon_id, vacc_dose_number, vaxdate, vacc_product_name,
                     vacc_batch_number)]
setkey(vax.shield, anon_id, vacc_dose_number) 
vax.shield <- unique(vax.shield, by=key(vax.shield))
vacc.wide <- dcast(vax.shield, anon_id ~ vacc_dose_number,
                   value.var=c("vaxdate", "vacc_product_name",
                               "vacc_batch_number"))
rm(vax.shield)

vacc.wide <- unique(vacc.wide)
vacc.wide[, weekdose1 := floor(as.integer(vaxdate_1 - as.Date("2020-12-01")) / 7)]
vacc.wide[, weekdose1 := as.Date("2010-12-01") + 7 * weekdose1]
setkey(vacc.wide, anon_id)

shielded.full <- vacc.wide[shielded.full]
rm(vacc.wide)

save(shielded.full, file=paste0(datadir, "shielded.linked.RData"))

rm(shielded.full)

####### rheumatology ######################################

datadir <- "data/2021-09-22/"

glasgow <- RDStodt(paste0(datadir, "Glasgow_data_rheum_2021-09-22.rds"))
setnames(glasgow, "main_diagnosis_text", "diagnosis")
setnames(glasgow, "umbrella", "diagnosis_recoded")
setnames(glasgow, "drug_class", "drug_recoded")
grampian1 <- RDStodt(paste0(datadir, "Grampian_data1_rheum_2021-09-22.rds"))
grampian2 <- RDStodt(paste0(datadir, "Grampian_data2_rheum_2021-09-22.rds")) # 39 of 198 in grampian1 are in grampian2
lothian <- RDStodt(paste0(datadir, "Lothian_data_rheum_2021-09-22.rds"))
setnames(lothian, "diseases", "diagnosis_recoded", skip_absent=TRUE)
setnames(lothian, "drugs", "drug_recoded")

rheumatol <- rbind(glasgow, grampian1, grampian2, lothian, fill=TRUE)
rheumatol <- rheumatol[, .(anon_id, diagnosis, diagnosis_recoded, drug, drug_recoded)]
names(rheumatol) <- c("anon_id", "diagnosis", "diagnosis.group", "drug", "drug.class")
rheumatol[is.na(drug.class), drug.class := drug]
with(rheumatol, table(drug.class, exclude=NULL))
rheumatol[drug=="Abatacept", drug.class := "T cell co-stim"]
## Ustekinumab is an IL12/23 inhibitor
rheumatol[, drug.class := car::recode(drug.class,
                                      "c('ABA', 'RTX')='B cell depletion';
                                       c('BEL', 'Belimumab')='BLySi';  
                                       'TNF'='TNFi';  
                                       c('IL12', 'IL12/23')='IL12/23i';  
                                       'IL17'='IL17i';  
                                       'IL6'='IL6i';  
                                       'IL1'='IL1i';  
                                       'JAK'='JAKi'")]  
rheumatol[, drugclassgr7 := drug.class]
rheumatol[, drugclassgr7:= car::recode(drugclassgr7,
                                        "c('BLySi', 'CYC', 'IL12/23i', 'IL1i', 'T cell co-stim')='Other biologic'",
                                        as.factor=TRUE,
                                       levels=c("TNFi", "IL6i", "IL17i", "B cell depletion",
                                                "Other biologic", "JAKi", "Methotrexate"))]
table(rheumatol$drugclassgr7)

with(rheumatol, table(diagnosis.group, exclude=NULL))
     
rheumatol[, diagnosis.group := car::recode(diagnosis.group,
                                      "c('AXPSA', 'AxSpA/Psoriatic', 'PsA', 'Sacroiliitis')='AXSPA';
                                       c('Connective Tissue disease', 'CTD', 'CYD')='CTD';  
                                       c('IEOSJA', 'SAPHO', 'OTH')='Other'")]
rheumatol[grep("Behcet|Bechet", diagnosis.group, ignore.case=TRUE), diagnosis.group := "Other"]
with(rheumatol, table(diagnosis.group, exclude=NULL))
        
rheumatol[grep("rheumatoid|sero.?positive|^RA|^Palindromic|^juvenile|Still'?s|^JIA", diagnosis, ignore.case=TRUE),
          diagnosis.group := "RA"]
rheumatol[grep("^mixed|lupus|SLE|scleroderma|polym[iy]ositis|dermatomyositis|CTD|polyarteritis|polymyalgia|sjogren", diagnosis, ignore.case=TRUE),
          diagnosis.group := "CTD"]
rheumatol[grep("spondyl|psor|^PSA|Enteropathic|Crohn|bowel", diagnosis, ignore.case=TRUE),
          diagnosis.group := "AXSPA"]
rheumatol[grep("Behcet|Bechet", diagnosis, ignore.case=TRUE), diagnosis.group := "Other"]
rheumatol[is.na(diagnosis.group) & grepl("inflammatory arthritis", diagnosis, ignore.case=TRUE),
          diagnosis.group := "RA"]
rheumatol[is.na(diagnosis.group) & grepl("sero[ -]?negative arthritis", diagnosis, ignore.case=TRUE), diagnosis.group := "RA"]

with(rheumatol[is.na(diagnosis.group)], table(diagnosis, exclude=NULL))
with(rheumatol, table(diagnosis.group, exclude=NULL))

rheumatol[is.na(diagnosis.group) & !is.na(diagnosis), diagnosis.group := "Other"]
rheumatol[, diagnosis.group := car::recode(diagnosis.group,
                                           "c('Other', NA)='Other/not recorded'",
                                           as.factor=TRUE,
                                           levels=c("RA", "CTD", "AXSPA", "Other/not recorded"))]

with(rheumatol, table(diagnosis.group, exclude=NULL))
rheumatol <- unique(rheumatol[, .(anon_id, diagnosis.group, drugclassgr7)])
setkey(rheumatol, anon_id)

## dcast diagnoses in wide format
rheumatol.diagnoses <- unique(rheumatol[, .(anon_id, diagnosis.group)])
rheumatol.diagnoses <- dcast(data=rheumatol.diagnoses, formula=anon_id ~ diagnosis.group,
                             fun.aggregate=length)
rheumatol.diagnoses[`Other/not recorded`==1 & (AXSPA==1 | RA==1 | CTD==1), `Other/not recorded` := 0]
rheumatol.diagnoses[AXSPA==1 & (RA==1 | CTD==1), AXSPA := 0] 
rheumatol.diagnoses[CTD==1 & RA==1, CTD := 0] 
rheumatol.diagnoses <- melt(rheumatol.diagnoses, id="anon_id", variable.name="diagnosis.group")
rheumatol.diagnoses <- rheumatol.diagnoses[value==1]
rheumatol.diagnoses[, value := NULL]

## dcast drug classes in wide fromat
rheumatol.drugclasses <- unique(rheumatol[, .(anon_id, drugclassgr7)])
rheumatol.drugclasses <- dcast(data=rheumatol.drugclasses, formula=anon_id ~ drugclassgr7,
                             fun.aggregate=length)

setkey(rheumatol.diagnoses, anon_id)
setkey(rheumatol.drugclasses, anon_id)
rheumatol.wide <- rheumatol.drugclasses[rheumatol.diagnoses]
setkey(rheumatol.wide, anon_id)

save(rheumatol.wide, file=paste0(datadir, "rheumatol.wide.RData"))

