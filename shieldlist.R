
datadir <- "./data/2021-07-28/"
shielding.full.filename <- paste0(datadir,
                                         "CC_shielding_patients_anon_2021-07-28.csv")

shielded.full <- fread(shielding.full.filename)
setkey(shielded.full, anon_id)
setnames(shielded.full, "group", "shield.group", skip_absent=TRUE)
setnames(shielded.full, "batch", "shield.batch", skip_absent=TRUE)
setnames(shielded.full, "age", "age_years", skip_absent=TRUE)
shielded.full[, sex := as.factor(car::recode(sex, "1='Male'; 2='Female'"))]
setkey(shielded.full, shield.batch)
shielded.full <- batches[shielded.full]
shielded.full[, removal := nchar(removal)]

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


## import case status into shielded.full
cases <- cc.all[CASE==1, .(anon_id, specimen_date, casegr, casegr2, casegr3)]
setkey(cases, anon_id)
setkey(shielded.full, anon_id)
shielded.full <- cases[shielded.full]
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

shielded.full[, exitdate := specimen_date]
shielded.full[is.na(exitdate), exitdate := max(shielded.full$specimen_date, na.rm=TRUE)]
shielded.full[, entrydate := as.Date("2020-03-01")]
shielded.full[entrydate >= exitdate, entrydate := NA]

save(shielded.full, file=paste0(datadir, "shielded.full.RData"))
rm(shielded.full)

####### rheumatology ######################################

glasgow <- fread(paste0(datadir, "Glasgow_data_rheum_2021-07-28.csv"))
setnames(glasgow, "main_diagnosis_text", "diagnosis")
setnames(glasgow, "umbrella", "diagnosis_recoded")
setnames(glasgow, "drug_class", "drug_recoded")
grampian1 <- fread(paste0(datadir, "Grampian_data1_rheum_2021-07-28.csv"))
grampian2 <- fread(paste0(datadir, "Grampian_data2_rheum_2021-07-28.csv")) # 39 of 198 in grampian1 are in grampian2
lothian <- fread(paste0(datadir, "Lothian_data_rheum_2021-07-28.csv"))
setnames(lothian, "diseases", "diagnosis_recoded")
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
#with(rheumatol, table(drug.class, exclude=NULL))
#with(rheumatol, table(drug, drug.class, exclude=NULL))

rheumatol[grep("rheumatoid|sero.?positive", diagnosis, ignore.case=TRUE), diagnosis.group := "RA"]
rheumatol[grep("^mixed|lupus|SLE|scleroderma|polymyositis", diagnosis, ignore.case=TRUE), diagnosis.group := "Lupus_scleroderma_mixed"]
rheumatol[grep("spondyl|psor", diagnosis, ignore.case=TRUE), diagnosis.group := "AXSPA"]
rheumatol[, diagnosis.group := car::recode(diagnosis.group,
                                      "c('AXPSA', 'AxSpA/Psoriatic', 'PsA', 'Sacroiliitis')='AXSPA';
                                       c('Connective Tissue disease', 'CYD')='CTD';  
                                       c('IEOSJA', 'SAPHO', 'OTH')='Other'")]
rheumatol[diagnosis.group=="Behcet's Disease", diagnosis.group := "Other"]
with(rheumatol, table(diagnosis.group, exclude=NULL))
#with(rheumatol, table(diagnosis, diagnosis.group, exclude=NULL))
with(rheumatol, table(diagnosis.group, exclude=NULL))
rheumatol[diagnosis.group=="CTD" | is.na(diagnosis.group), diagnosis.group := "Other"]
rheumatol[diagnosis.group=="RA" | diagnosis.group=="Lupus_scleroderma_mixed", diagnosis.group := "RArelated"]
with(rheumatol, table(diagnosis.group, exclude=NULL))
rheumatol <- unique(rheumatol[, .(anon_id, diagnosis.group, drug.class)])
rheumatol[, drug.class := car::recode(drug.class,
                       "c('BLySi', 'CYC', 'IL12/23i', 'IL1i', 'T cell co-stim')='Other'")]
table(rheumatol$drug.class)

rheumatol[, diagnosis.group := as.factor(diagnosis.group)]
rheumatol[, drug.class := as.factor(drug.class)] # 4743 records
setkey(rheumatol, anon_id)

rheumatol.diagnoses <- unique(rheumatol[, .(anon_id, diagnosis.group)])
rheumatol.diagnoses <- dcast(data=rheumatol.diagnoses, formula=anon_id ~ diagnosis.group,
                             fun.aggregate=length)
rheumatol.diagnoses[Other==1 & (AXSPA==1 | RArelated==1), Other := 0]
rheumatol.diagnoses[AXSPA==1 & RArelated==1, AXSPA := 0] # only 4 of these
rheumatol.diagnoses <- melt(rheumatol.diagnoses, id="anon_id", variable.name="diagnosis.group")
rheumatol.diagnoses <- rheumatol.diagnoses[value==1]
rheumatol.diagnoses[, value := NULL]

rheumatol.drugclasses <- unique(rheumatol[, .(anon_id, drug.class)])
rheumatol.drugclasses <- dcast(data=rheumatol.drugclasses, formula=anon_id ~ drug.class,
                             fun.aggregate=length)

setkey(rheumatol.diagnoses, anon_id)
setkey(rheumatol.drugclasses, anon_id)
rheumatol.wide <- rheumatol.drugclasses[rheumatol.diagnoses]
setkey(rheumatol.wide, anon_id)

