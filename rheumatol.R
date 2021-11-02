library(data.table)
library(survival)

library(data.table)
source("helperfunctions.R")

datadir <- "data/2021-09-22/"

if(!exists("cc.all")) {
    load(paste0(datadir, "cc.all.RData"))
}
load(paste0(datadir, "rheumatol.wide.RData"))
load(paste0(datadir, "shielded.full.RData"))
load(paste0(datadir, "drugdoses.RData"))

smr00.filename <- paste0(datadir, "CC_SMR00_2021-07-28.rds")
smr00.filename <- paste0(datadir, "CC_SMR00_2021-09-22.rds")
smr00 <- RDStodt(smr00.filename, keyname="anon_id") # outpatients
smr00[, clinic_date := as.Date(clinic_date)] # convert PosixCT
smr00.rheumatol <- smr00[specialty=="AR", .(anon_id, clinic_date)]


setkey(cc.all, anon_id, specimen_date)
cc.all <- drugdoses[cc.all]
cc.all[, prednisolone.gt5mgday := as.integer(prednisolone.equiv / 7 > 5)]
cc.all[, csdmard := as.integer(prednisolone.equiv +
                               sulfasalazine +
                               leflunomide +
                               hydroxychloroquine +
                               methotrexate > 0)]
cc.all[, rheumatol.opd := as.integer(anon_id %in% smr00.rheumatol$anon_id)]
cc.all[, rheumatol.csdmard := as.integer(csdmard==1 & rheumatol.opd==1)]
setkey(cc.all, anon_id)

cc.all <- rheumatol.wide[cc.all]
cc.all[, diagnosis.group := as.character(diagnosis.group)]

# biologics list covers Lothian, Grampian, part of GGC, Highland and islands 
print(with(cc.all, table(hb2019name, is.na(diagnosis.group))))

table.hospdiag <- with(cc.all,
                   table(rheumatol.diag,  rheumatol.opd==1 & csdmard==1, exclude=NULL))
colsums.table.hospdiag <- colSums(table.hospdiag)
print(rbind(table.hospdiag, colsums.table.hospdiag))

## almost all those on biologics list had a rheumatology outpatient consultation in last 5 years
table.diag <- with(cc.all[!is.na(diagnosis.group)],
                   table(diagnosis.group,  rheumatol.opd==1, exclude=NULL))
colsums.table.diag <- colSums(table.diag)
print(rbind(table.diag, colsums.table.diag))

## but only half of these had a csdmard prescription in last 240 days
table.csdmard <- with(cc.all[!is.na(diagnosis.group)],
                   table(diagnosis.group,  rheumatol.opd==1 & csdmard==1, exclude=NULL))
colsums.table.csdmard <- colSums(table.csdmard)
print(rbind(table.csdmard, colsums.table.csdmard))

cc.all[rheumatol.opd==1 & csdmard==1 & 
       rheumatol.diag=="No rheumatologic diagnosis",
       rheumatol.diag := "Rheumatology OP on csDMARDs, no discharge.diagnosis"]

with(cc.all, table(rheumatol.diag))

cc.all[, rheumatol.diagf := factor(rheumatol.diag)]
with(cc.all, table(rheumatol.diagf))
cc.all[, rheumatol.diagf := factor(rheumatol.diagf,
                                   levels=levels(rheumatol.diagf)[c(3, 5, 4, 1, 2)])]
with(cc.all, table(rheumatol.diagf))

 
## fill in missing diagnoses on biologics list
cc.all[diagnosis.group=="Other/not recorded" &
        (rheumatol.diag=="RA" | rheumatol.diag=="AXSPA" | rheumatol.diag=="CTD"),
       diagnosis.group := rheumatol.diag]
    
setnafill(cc.all, fill=0,
          cols=c("B cell depletion", "IL17i", "IL6i", "JAKi", "Methotrexate", "Other biologic", "TNFi"))

cc.all[TNFi==1, biologic3 := "TNFi only"]
cc.all[`B cell depletion`==1 | IL17i==1 | IL6i==1 | JAKi==1 | `Other biologic`==1,
       biologic3 := "Other biologic or JAKi"]
cc.all[is.na(biologic3), biologic3 := "No biologic"]
cc.all[, biologic3 := factor(biologic3, levels=c("No biologic", "TNFi only", "Other biologic or JAKi"))]
# cc.all[, ar.outpatient := as.integer(anon_id %in% smr00.rheumatol$anon_id)]

cc.all[, methotrexate := methotrexate * rheumatol.opd]
cc.all[, hydroxychloroquine := hydroxychloroquine * rheumatol.opd]
cc.all[, sulfasalazine := sulfasalazine * rheumatol.opd]
cc.all[, leflunomide := leflunomide * rheumatol.opd]
cc.all[, prednisolone.gt5mgday := prednisolone.gt5mgday * rheumatol.opd]
cc.all[, prednisolone.dosegr := cut(prednisolone.equiv / 7 * rheumatol.opd, 
                                    breaks=c(-1, 0, 5, 10,
                                             ceiling(max(prednisolone.equiv/7))),
                                    labels=c("0", ">0 to 5", ">5 to 10", ">10"))]
with(cc.all[CASE==1 & (casegroup=="A" | casegroup=="B") &
            rheumatol.opd==1 & csdmard==1], table(prednisolone.dosegr, exclude=NULL))
     
rheumatol.wide[, ar.outpatient := as.integer(anon_id %in% smr00.rheumatol$anon_id)]
with(rheumatol.wide, table(ar.outpatient))

rheumatol.wide[, sampled := as.integer(anon_id %in% cc.all$anon_id)]
rheumatol.wide[, csdmard := as.integer(anon_id %in% cc.all[csdmard==1, anon_id])]

rheumatol.wide[, hosp.diag := as.integer(anon_id %in% cc.all[rheumatol.diag != "No rheumatologic diagnosis", anon_id])]

cases <- cc.all[CASE==1, .(anon_id, specimen_date, casegroup, fatal.casegroup, fatalcase, casegr, casegr2, casegr3, adm.within14, critical.within21)]
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
    tabulate.freqs.regressions(varnames=c("vax14.factor", "care.home", "inpat.recent",
                                          "rheumatol.diagf"),
                               data=cc.all[(casegroup=="A" | casegroup=="B")])
#rownames(table.rheumatol.diagnosis)[7:9] <- c("Rheumatoid", "Connective tissue disease", #"Psoriasis / other seronegative")

## 3. case-control analysis of associations with DMARDS, biologic and otherwise.  
# exclude ibd and specific cancers from this analysis
table.rheumatol.drugs <-
    tabulate.freqs.regressions(varnames=c("vax14.factor", "care.home", "inpat.recent",
                                          "methotrexate",
                                          "hydroxychloroquine", "sulfasalazine",
                                          "leflunomide", "prednisolone.dosegr", 
                                          "TNFi", "B cell depletion", "IL6i", "IL17i", "JAKi"),
                               data=cc.all[(casegroup=="A" | casegroup=="B")])
rownames(table.rheumatol.drugs) <- firstup(rownames(table.rheumatol.drugs))
rownames(table.rheumatol.drugs)[1:3] <- c("Unvaccinated", "1 dose", "2 doses")

table.rheumatol.drugs3 <-
    tabulate.freqs.regressions(varnames=c("vax14.factor", "care.home", "inpat.recent",
                                          "rheumatol.csdmard","biologic3"),
                               data=cc.all[(casegroup=="A" | casegroup=="B")])

table.shielded.hosp <- paste.rowpercent(with(rheumatol.wide,
                                             table(shielded,
                                                   adm.within14, exclude=NULL)), digits=1)
 
rmarkdown::render("rheumatol.Rmd")
rmarkdown::render("ssr.Rmd")
