library(data.table)
library(survival)

library(data.table)
source("helperfunctions.R")

datadir <- "data/2021-07-28/"

if(!exists("cc.all")) {
    load(paste0(datadir, "cc.all.RData"))
}
load(paste0(datadir, "rheumatol.wide.RData"))
load(paste0(datadir, "shielded.full.RData"))
load(paste0(datadir, "drugdoses.RData"))

smr00.filename <- paste0(datadir, "CC_SMR00_2021-07-28.rds")
smr00 <- RDStodt(smr00.filename, keyname="anon_id") # outpatients
smr00[, clinic_date := as.Date(clinic_date)] # convert PosixCT
smr00.rheumatol <- smr00[specialty=="AR", .(anon_id, clinic_date)]


setkey(cc.all, anon_id, specimen_date)
cc.all <- drugdoses[cc.all]

setkey(cc.all, anon_id)
cc.all <- rheumatol.wide[cc.all]
cc.all[, diagnosis.group := as.character(diagnosis.group)]

# biologics list covers Lothian, Grampian, part of GGC, Highland and islands 
print(with(cc.all, table(hb2019name, is.na(diagnosis.group))))

## almost all those on biologics list had a rheumatology outpatient consultation
table.diag <- with(cc.all[!is.na(diagnosis.group)],
                   table(rheumatol.diag,  anon_id %in% smr00.rheumatol$anon_id, exclude=NULL))
colsums.table.diag <- colSums(table.diag)
print(rbind(table.diag, colsums.table.diag))


cc.all[is.na(diagnosis.group) & rheumatol.diag=="No rheumatologic diagnosis",
       diagnosis.group := "No rheumatologic diagnosis"] 
cc.all[is.na(diagnosis.group) & rheumatol.diag !="No rheumatologic diagnosis",
       diagnosis.group := rheumatol.diag]
cc.all[diagnosis.group=="Other/not recorded" &
        (rheumatol.diag =="AXSPA" | rheumatol.diag=="CTD"),
    diagnosis.group := rheumatol.diag]
cc.all[, diagnosis.group := factor(diagnosis.group,
                                   levels=c("No rheumatologic diagnosis",
                                            "RA",
                                            "CTD",
                                            "AXSPA",
                                            "Other/not recorded"))]
setnafill(cc.all, fill=0,
          cols=c("B cell depletion", "IL17i", "IL6i", "JAKi", "Methotrexate", "Other biologic", "TNFi"))

cc.all[TNFi==1, biologic3 := "TNFi only"]
cc.all[`B cell depletion`==1 | IL17i==1 | IL6i==1 | JAKi==1 | `Other biologic`==1,
       biologic3 := "Other biologic or JAKi"]
cc.all[is.na(biologic3), biologic3 := "No biologic"]
cc.all[, biologic3 := factor(biologic3, levels=c("No biologic", "TNFi only", "Other biologic or JAKi"))]
cc.all[, ar.outpatient := as.integer(anon_id %in% smr00.rheumatol$anon_id)]

rheumatol.wide[, ar.outpatient := as.integer(anon_id %in% smr00.rheumatol$anon_id)]
with(rheumatol.wide, table(ar.outpatient))


cases <- cc.all[CASE==1, .(anon_id, specimen_date, casegroup, fatal.casegroup, fatalcase, casegr, casegr2, casegr3)]
setkey(cases, anon_id)
setkey(rheumatol.wide, anon_id)
rheumatol.wide <- cases[rheumatol.wide]
rheumatol.wide[, casegroup := car::recode(casegroup,
                                          "NA='Not diagnosed as case';
                                          'A'='Severe';
                                          'B'='Hospitalised not severe';
                                   c('C', 'Unclassified')='Not hospitalised not severe'",
                                   as.factor=TRUE, 
                                   levels=c("Not diagnosed as case",
                                            "Not hospitalised not severe",
                                            "Hospitalised not severe", "Severe"))]

rheumatol.wide[, shielded := anon_id %in% shielded.full$anon_id]

table.rheumatol.shielded <- with(rheumatol.wide, table(casegroup, shielded))
table.rheumatol.shielded <- cbind(table.rheumatol.shielded, rowSums(table.rheumatol.shielded))
colsums.table.rheumatol.shielded <- colSums(table.rheumatol.shielded)
table.rheumatol.shielded <- paste.colpercent(table.rheumatol.shielded, digits=1)
table.rheumatol.shielded <- rbind(table.rheumatol.shielded, colsums.table.rheumatol.shielded)
colnames(table.rheumatol.shielded) <- c("Not added to shielding list", "On shielding list", "All")
rownames(table.rheumatol.shielded)[5] <- "All"
table.rheumatol.shielded

## 2. case-control analysis of those with rheumatology diagnoses, on or off shielding list
table.rheumatol.diagnosis <-
    tabulate.freqs.regressions(varnames=c("vax14.factor",
                                          "diagnosis.group"),
                               data=cc.all[(casegroup=="A" | casegroup=="B")])

## 3. case-control analysis of associations with DMARDS, biologic and otherwise.  
# exclude ibd and specific cancers from this analysis
table.rheumatol.drugs <-
    tabulate.freqs.regressions(varnames=c("vax14.factor", "methotrexate",
                                          "hydroxychloroquine", "sulfasalazine",
                                          "leflunomide", "prednisolone", 
                                          "TNFi", "B cell depletion", "IL6i", "IL17i", "JAKi"),
                               data=cc.all[(casegroup=="A" | casegroup=="B") # &
                                       #    ibd==0 & 
                                       #    shield.group != "Specific cancers" &
                                       #    shield.group != "Severe respiratory" &
                                       #    chronresp==0
                                           ])

table.rheumatoldiag.drugs <-
    tabulate.freqs.regressions(varnames=c("vax14.factor", "methotrexate",
                                          "hydroxychloroquine", "sulfasalazine",
                                          "leflunomide", "prednisolone", 
                                          "TNFi", "B cell depletion", "IL6i", "IL17i", "JAKi"),
                               data=cc.all[(casegroup=="A" | casegroup=="B") &
                                           (diagnosis.group != "No rheumatologic diagnosis" |
                                            ar.outpatient==1)])

 
                               
###############################################

rmarkdown::render("rheumatol.Rmd")
