library(data.table)
library(survival)
library(ggplot2)

source("helperfunctions.R")

#datadir <- "data/2021-07-28/"
#datadir <- "data/2021-09-02/"
datadir <- "data/2021-09-22/"

if(!exists("cc.all")) {
    load(paste0(datadir, "cc.all.RData"))
}

source("rrtimewin.R", verbose=TRUE)

gc()
objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

interval.length <- 28  # for poisson regression of shielding cohort
lastdate.cases <- as.Date("2021-09-19")

## tabulate case groups by vax status in shielding cohort
paste.colpercent(with(cc.all[CASE==1 & specimen_date >= as.Date("2020-12-01") &
                             listedgr3=="Eligible for shielding"],
                      table(casegr, vax14.factor)))

paste.rowpercent(with(cc.all[CASE==1 & specimen_date >= as.Date("2020-12-01") &
                             vax14.factor=="2"], 
                      table(shield.group, casegroup)))

if(FALSE) {
## tabulate severe cases after 1 April
 
    severe.byagevax <- with(cc.all[CASE==1 & casegroup=="A" &
                             specimen_date >= as.Date("2021-04-01") & specimen_date <= lastdate.cases],
                      table(agegr20plus5, vax14.factor))
severe.byagevax.rowsums <- rowSums(severe.byagevax)
severe.byagevax <- cbind(severe.byagevax, severe.byagevax.rowsums)
severe.byagevax <- paste.colpercent(severe.byagevax)
colnames(severe.byagevax)[4] <- "All"

severe.byagerisk <- with(cc.all[CASE==1 & casegroup=="A" &
                             specimen_date >= as.Date("2021-04-01") & specimen_date <= lastdate.cases],
                      table(agegr20plus5, listedgr3))
severe.byagerisk.colsums <- colSums(severe.byagerisk)
severe.byagerisk <- rbind(severe.byagerisk, severe.byagerisk.colsums)
severe.byagerisk <- paste.rowpercent(severe.byagerisk)
rownames(severe.byagerisk)[6] <- "All"

casesnorisk <- cc.all[casegroup=="A" & listedgr3=="No risk condition" & 
                      specimen_date >= as.Date("2021-04-01") &
                      specimen_date <= lastdate.cases &
                      age_years < 50] # 91 cases

setkey(casesnorisk, anon_id, specimen_date)
setkey(cc.diagnoses, anon_id, specimen_date)
casesnorisk.diagnoses <- cc.diagnoses[casesnorisk]
casesnorisk.chnums <- unique(casesnorisk.diagnoses,
                                by=c("anon_id", "specimen_date", "chnum"))

chnum.table <- with(casesnorisk.chnums, table(chnum, CASE))
colnames(chnum.table) <- c("Controls", "Cases")
chnum.table <- data.table(chnum=as.integer(rownames(chnum.table)),
                             as.data.frame.matrix(chnum.table))
chnum.table[, descr := icdchapters$name[as.integer(chnum)]]
chnum.table[, total := Controls + Cases]
chnum.table[, Controls := paste0(Controls, " (",
                                 round(100 * Controls / casesnorisk[CASE==0, .N]),
                                 "\\%)")]
chnum.table[, Cases := paste0(Cases, " (",
                                 round(100 * Cases / casesnorisk[CASE==1, .N]),
                                 "\\%)")]
chnum.table[, chnum := as.roman(chnum)]
chnum.table <- chnum.table[total >= 20, .(chnum, descr, Controls, Cases)]

casesnorisk.subchnums <- unique(casesnorisk.diagnoses,
                                by=c("anon_id", "specimen_date", "chnum", "subchnum"))
subchnum.table <- with(casesnorisk.subchnums, table(subchnum, CASE))
colnames(subchnum.table) <- c("Controls", "Cases")
subchnum.table <- data.table(subchnum=as.integer(rownames(subchnum.table)),
                             as.data.frame.matrix(subchnum.table))
subchnum.table[, descr := icdsubchapters$name[as.integer(subchnum)]]
subchnum.table[, total := Controls + Cases]
subchnum.table[, Controls := paste0(Controls, " (",
                                 round(100 * Controls / casesnorisk[CASE==0, .N]),
                                 "\\%)")]
subchnum.table[, Cases := paste0(Cases, " (",
                                 round(100 * Cases / casesnorisk[CASE==1, .N]),
                                 "\\%)")]
subchnum.table <- subchnum.table[total >= 20, .(subchnum, descr, Controls, Cases)]

## subch 66 is alcohol-related diagnoses F00

casesnorisk.scrips <- cc.scripsbnf[casesnorisk]
casesnorisk.scrips[, bnfchnum := substr(bnf_item_code, 1, 2)]
casesnorisk.bnfchnums <- unique(casesnorisk.scrips, 
                                by=c("anon_id", "specimen_date", "bnfchnum"))
bnfchnum.table <- with(casesnorisk.bnfchnums, table(bnfchnum, CASE))
colnames(bnfchnum.table) <- c("Controls", "Cases")
bnfchnum.table <- data.table(bnfchnum=as.integer(rownames(bnfchnum.table)),
                             as.data.frame.matrix(bnfchnum.table))
bnfchnum.table[, total := Controls + Cases]
bnfchnum.table[, Controls := paste0(Controls, " (",
                                 round(100 * Controls / casesnorisk[CASE==0, .N]),
                                 "\\%)")]
bnfchnum.table[, Cases := paste0(Cases, " (",
                                 round(100 * Cases / casesnorisk[CASE==1, .N]),
                                 "\\%)")]
bnfchnum.table[, descr := bnfchapters$chaptername[bnfchnum]]
bnfchnum.table <- bnfchnum.table[total >= 20, .(bnfchnum, descr, Controls, Cases)]

table.numdrugs <- cbind(summary(casesnorisk[CASE==0]$numdrugs),
                        summary(casesnorisk[CASE==1]$numdrugs))
colnames(table.numdrugs) <- c("Controls", "Cases")
}
########################################################
## tabulate freqs of risk factors including shielding conditions in cases and controls

freqs1.riskgroups <- univariate.tabulate(varnames=c("care.home", "qSIMD.integer", "hh.over18", 
                                                     "numdrugs", "inpat.recent",
                                                     "listedgr3", 
                                                     "vax14.factor"),
                                        data=cc.all[casegroup=="A" & age_years < 50 & 
                                                #    care.home=="Independent" & 
                                                    specimen_date >= as.Date("2021-04-01") & specimen_date <= lastdate.cases],
                                        drop.reflevel=FALSE)
freqs2.riskgroups <- univariate.tabulate(varnames=c("care.home", "qSIMD.integer", "hh.over18", 
                                                     "numdrugs", "inpat.recent",
                                                     "listedgr3", 
                                                     "vax14.factor"),
                                        data=cc.all[casegroup=="A" & age_years >= 50 & 
                                                #    care.home=="Independent" & 
                                                    specimen_date >= as.Date("2021-04-01") & specimen_date <= lastdate.cases],
                                        drop.reflevel=FALSE)
freqs.riskgroups <- cbind(freqs1.riskgroups, freqs2.riskgroups)
rownames(freqs.riskgroups) <- replace.names(rownames(freqs.riskgroups))

freqs.vaxgr <- univariate.tabulate(varnames=c("care.home", "qSIMD.integer", "hh.over18", 
                                                   "numdrugs", "inpat.recent", "shield.group", 
                                                   "vaxgr"),
                                        data=cc.all[casegroup=="A" &
                                               #     care.home=="Independent" & 
                                                    specimen_date >= as.Date("2020-12-01")],
                                        drop.reflevel=FALSE)
rownames(freqs.vaxgr) <- replace.names(rownames(freqs.vaxgr))

#####################################
## rate ratios by risk group x dose

table.riskgr <- summary(clogit(data=cc.all[casegroup=="A" & #care.home=="Independent" & 
                           specimen_date >= as.Date("2020-12-01")],
                           formula=CASE ~ care.home + hh.over18 + numdrugs + inpat.recent +
                               listedgr3/vax14.factor + strata(stratum)))
table.riskgr <- table.riskgr$coefficients
effect <- rownames(table.riskgr)
table.riskgr <- as.data.table(table.riskgr)
table.riskgr[, rateratio := or.ci(coef, `se(coef)`)]
table.riskgr[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
table.riskgr <- table.riskgr[, .(rateratio, pvalue)]
table.riskgr <- data.table(effect=effect, # freqs.riskgr,
                               table.riskgr)
table.riskgr[, effect := gsub("^care.home", "", effect)]
table.riskgr[, effect := gsub("^listedgr3", "", effect)]
table.riskgr[, effect := gsub("TRUE$", "", effect)]
table.riskgr[, effect := gsub(":vax14\\.factor", ": vax dose ", effect)]
table.riskgr[, effect := replace.names(effect)]
table.riskgr <- table.riskgr[c(1:6, 7, 10, 8, 11, 9, 12), ]
table.riskgr <- table.riskgr[7:12, ]
table.riskgr[, dose := rep(c(1, 2), 3)]
table.riskgr[, effect := gsub(":.*", "", effect)]
table.riskgr <- dcast(table.riskgr, effect ~ dose, value.var=c("rateratio", "pvalue"))
table.riskgr <- table.riskgr[c(3, 2, 1), ]
table.riskgr <- table.riskgr[, c(1, 2, 4, 3, 5)]

# repeat rate ratios by risk group x dose excluding fatal cases without covid_ucod==1
table.nonfatalorucod <- summary(clogit(data=cc.all[casegroup=="A" & #care.home=="Independent" & 
                                           specimen_date >= as.Date("2020-12-01") &
                                           (!fatalcase | covid_ucod==1)],
                           formula=CASE ~ care.home + hh.over18 + numdrugs + inpat.recent +
                               listedgr3/vax14.factor + strata(stratum)))
table.nonfatalorucod <- table.nonfatalorucod$coefficients
effect <- rownames(table.nonfatalorucod)
table.nonfatalorucod <- as.data.table(table.nonfatalorucod)
table.nonfatalorucod[, rateratio := or.ci(coef, `se(coef)`)]
table.nonfatalorucod[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
table.nonfatalorucod <- table.nonfatalorucod[, .(rateratio, pvalue)]
table.nonfatalorucod <- data.table(effect=effect, # freqs.riskgr,
                               table.nonfatalorucod)
table.nonfatalorucod[, effect := gsub("^care.home", "", effect)]
table.nonfatalorucod[, effect := gsub("^listedgr3", "", effect)]
table.nonfatalorucod[, effect := gsub("TRUE$", "", effect)]
table.nonfatalorucod[, effect := gsub(":vax14\\.factor", ": vax dose ", effect)]
table.nonfatalorucod[, effect := replace.names(effect)]
table.nonfatalorucod <- table.nonfatalorucod[c(1:6, 7, 10, 8, 11, 9, 12), ]

# rate ratios by shield group x dose
table.shieldgr <- summary(clogit(data=cc.all[casegroup=="A" & #care.home=="Independent" & 
                           #specimen_date <= as.Date("2021-03-16") &
                           specimen_date >= as.Date("2020-12-01")],
                           formula=CASE ~ care.home + hh.over18 + numdrugs + inpat.recent +
                               shield.group/vax14.factor + strata(stratum)))
table.shieldgr <- table.shieldgr$coefficients
effect <- rownames(table.shieldgr)
table.shieldgr <- as.data.table(table.shieldgr)
table.shieldgr[, rateratio := or.ci(coef, `se(coef)`)]
table.shieldgr[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
table.shieldgr <- table.shieldgr[, .(rateratio, pvalue)]
table.shieldgr <- data.table(effect=effect, # freqs.riskgr,
                               table.shieldgr)
table.shieldgr[, effect := gsub("^care.home", "", effect)]
table.shieldgr[, effect := gsub("^shield\\.group", "", effect)]
table.shieldgr[, effect := gsub("TRUE$", "", effect)]
table.shieldgr[, effect := gsub(":vax14\\.dose", ": vax dose", effect)]
table.shieldgr[, effect := replace.names(effect)]
table.shieldgr <- table.shieldgr[c(14:19, 22:27), ]
table.shieldgr[, dose := rep(c(1, 2), each=6)]
table.shieldgr[, effect := gsub(":.*", "", effect)]
table.shieldgr <- dcast(table.shieldgr, effect ~ dose, value.var=c("rateratio", "pvalue"))
table.shieldgr <- table.shieldgr[c(5, 6, 4, 3, 2, 1), ]
table.shieldgr <- table.shieldgr[, c(1, 2, 4, 3, 5)]

cc.all[, shieldgr9 := as.character(shield.group)]
cc.all[shield.group=="Specific cancers" & bloodcancer==1 , shieldgr9 := "Blood cancer"]
cc.all[shieldgr9=="Specific cancers" & bloodcancer==0 , shieldgr9 := "Other specific cancers"]
cc.all[, shieldgr9 := factor(shieldgr9, levels=c("No risk condition",
                                                 "Moderate risk condition",
                                                 "Solid organ transplant",
                                                 "Blood cancer",
                                                 "Other specific cancers",
                                                 "Severe respiratory",
                                                 "Rare diseases",
                                                 "On immunosuppressants",
                                                 "Additional conditions"))]
table.shieldgr9 <- summary(clogit(data=cc.all[casegroup=="A" & #care.home=="Independent" & 
                           #specimen_date <= as.Date("2021-03-16") &
                           specimen_date >= as.Date("2020-12-01")],
                           formula=CASE ~ care.home + hh.over18 + numdrugs + inpat.recent +
                               shieldgr9/vax14.factor + strata(stratum)))
table.shieldgr9 <- table.shieldgr9$coefficients
effect <- rownames(table.shieldgr9)
table.shieldgr9 <- as.data.table(table.shieldgr9)
table.shieldgr9[, rateratio := or.ci(coef, `se(coef)`)]
table.shieldgr9[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
table.shieldgr9 <- table.shieldgr9[, .(rateratio, pvalue)]
table.shieldgr9 <- data.table(effect=effect, # freqs.riskgr,
                               table.shieldgr9)


##  rate ratios by risk group x dose-product
table.vaxgr <- summary(clogit(data=cc.all[casegroup=="A" & #care.home=="Independent" & 
                           specimen_date >= as.Date("2020-12-01")],
                           formula=CASE ~ care.home + hh.over18 + numdrugs + inpat.recent +
                               listedgr3/vaxgr + strata(stratum)))
table.vaxgr <- table.vaxgr$coefficients
effect <- rownames(table.vaxgr)
table.vaxgr <- as.data.table(table.vaxgr)
table.vaxgr[, rateratio := or.ci(coef, `se(coef)`)]
table.vaxgr[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
table.vaxgr <- table.vaxgr[, .(rateratio, pvalue)]
table.vaxgr <- data.table(effect=effect, # freqs.vaxgr,
                               table.vaxgr)
table.vaxgr[, effect := gsub("^care.home", "", effect)]
table.vaxgr[, effect := gsub("^listedgr3", "", effect)]
table.vaxgr[, effect := gsub("TRUE$", "", effect)]
table.vaxgr[, effect := gsub(":vaxgr", ": ", effect)]
table.vaxgr[, effect := replace.names(effect)]
table.vaxgr <- table.vaxgr[c(1:6,
                             7, 13, 10, 16,
                             8, 14, 11, 17,
                             9, 15, 12, 18), ]
table.vaxgr <- table.vaxgr[7:18, ]
table.vaxgr[, dose := rep(c(1, 2), 6)]
table.vaxgr[, effect := gsub(". doses? ", "", effect)]
table.vaxgr <- dcast(table.vaxgr, effect ~ dose, value.var=c("rateratio", "pvalue"))
table.vaxgr <- table.vaxgr[c(6, 5, 4, 3, 2, 1), ]
table.vaxgr <- table.vaxgr[, c(1, 2, 4, 3, 5)]

vaxshield.effect.model <- clogit(data=cc.all[casegroup=="A" & #care.home=="Independent" & 
                           specimen_date >= as.Date("2020-12-01")],
                           formula=CASE ~ care.home + hh.over18 + numdrugs + inpat.recent +
                               listedgr3 * vax14.factor + 
                               strata(stratum))
table.vaxshield <- summary(vaxshield.effect.model)$coefficients

vaxriskgr.interaction <- or.ci(table.vaxshield[12, 1], table.vaxshield[12, 3], ndigits=1)

cc.all[, vaxproduct := car::recode(vacc_product_name_1,
                                   "c('Pfizer', 'Moderna')='mRNA'",
                                   as.factor=TRUE,
                                   levels=c("mRNA", "AstraZeneca"))]
    
cc.shielded <- cc.all[(casegroup=="A" | casegroup=="B") &
                      listedgr3=="Eligible for shielding" & 
                      specimen_date >= as.Date("2020-12-01"),
                      .(CASE, care.home, hh.over18, numdrugs, inpat.recent, 
                        vax14.dose, vaxproduct, stratum)]
cc.shielded <- na.omit(cc.shielded)
cc.shielded[, vax14.dose := scale(vax14.dose, scale=FALSE)]
cc.shielded[, vaxproduct := scale(as.integer(vaxproduct), scale=FALSE)]
            
vaxshield.product.model <- clogit(data=cc.shielded,
                        formula=CASE ~ care.home + hh.over18 + numdrugs + inpat.recent +
                                      vax14.dose + vaxproduct + strata(stratum))
summary.vaxshield.product <- summary(vaxshield.product.model)
summary.vaxshield.product

########################################################
freqs1.hosp.riskgroups <- univariate.tabulate(varnames=c("care.home", "qSIMD.integer", "hh.over18", 
                                                     "numdrugs", "inpat.recent",
                                                     "listedgr3", 
                                                     "vax14.factor"),
                                        data=cc.all[(casegroup=="A" | casegroup=="B") & age_years < 50 & 
                                                #    care.home=="Independent" & 
                                                    specimen_date >= as.Date("2021-04-01") & specimen_date <= lastdate.cases],
                                        drop.reflevel=FALSE)
freqs2.hosp.riskgroups <- univariate.tabulate(varnames=c("care.home", "qSIMD.integer", "hh.over18", 
                                                     "numdrugs", "inpat.recent",
                                                     "listedgr3", 
                                                     "vax14.factor"),
                                        data=cc.all[(casegroup=="A" | casegroup=="B") & age_years >= 50 & 
                                                #    care.home=="Independent" & 
                                                    specimen_date >= as.Date("2021-04-01") & specimen_date <= lastdate.cases],
                                        drop.reflevel=FALSE)

freqs.hosp.riskgroups <- cbind(freqs1.hosp.riskgroups, freqs2.hosp.riskgroups)
rownames(freqs.hosp.riskgroups) <- replace.names(rownames(freqs.hosp.riskgroups))

freqs.hosp.vaxgr <- univariate.tabulate(varnames=c("care.home", "qSIMD.integer", "hh.over18", 
                                                   "numdrugs", "inpat.recent", "shield.group", 
                                                   "vaxgr"),
                                        data=cc.all[(casegroup=="A" | casegroup=="B") &
                                               #     care.home=="Independent" & 
                                                    specimen_date >= as.Date("2020-12-01")],
                                        drop.reflevel=FALSE)
rownames(freqs.hosp.vaxgr) <- replace.names(rownames(freqs.hosp.vaxgr))


###############################
# rate ratios by risk group x dose
table.hosp.riskgr <- summary(clogit(data=cc.all[(casegroup=="A" | casegroup=="B") & #care.home=="Independent" & 
                           specimen_date >= as.Date("2020-12-01")],
                           formula=CASE ~ care.home + hh.over18 + numdrugs + inpat.recent +
                               listedgr3/vax14.factor + strata(stratum)))
table.hosp.riskgr <- table.hosp.riskgr$coefficients
effect <- rownames(table.hosp.riskgr)
table.hosp.riskgr <- as.data.table(table.hosp.riskgr)
table.hosp.riskgr[, rateratio := or.ci(coef, `se(coef)`)]
table.hosp.riskgr[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
table.hosp.riskgr <- table.hosp.riskgr[, .(rateratio, pvalue)]
table.hosp.riskgr <- data.table(effect=effect, # freqs.riskgr,
                               table.hosp.riskgr)
table.hosp.riskgr[, effect := gsub("^care.home", "", effect)]
table.hosp.riskgr[, effect := gsub("^listedgr3", "", effect)]
table.hosp.riskgr[, effect := gsub("TRUE$", "", effect)]
table.hosp.riskgr[, effect := gsub(":vax14\\.factor", ": vax dose ", effect)]
table.hosp.riskgr[, effect := replace.names(effect)]
table.hosp.riskgr <- table.hosp.riskgr[c(1:6, 7, 10, 8, 11, 9, 12), ]
table.hosp.riskgr <- table.hosp.riskgr[7:12, ]
table.hosp.riskgr[, dose := rep(c(1, 2), 3)]
table.hosp.riskgr[, effect := gsub(":.*", "", effect)]
table.hosp.riskgr <- dcast(table.hosp.riskgr, effect ~ dose, value.var=c("rateratio", "pvalue"))
table.hosp.riskgr <- table.hosp.riskgr[c(3, 2, 1), ]
table.hosp.riskgr <- table.hosp.riskgr[, c(1, 2, 4, 3, 5)]

# rate ratios by shield group x dose
table.hosp.shieldgr <- summary(clogit(data=cc.all[(casegroup=="A" | casegroup=="B") &  
                           #specimen_date <= as.Date("2021-03-16") &
                           specimen_date >= as.Date("2020-12-01")],
                           formula=CASE ~ care.home + hh.over18 + numdrugs + inpat.recent +
                               shield.group/vax14.factor + strata(stratum)))
table.hosp.shieldgr <- table.hosp.shieldgr$coefficients
effect <- rownames(table.hosp.shieldgr)
table.hosp.shieldgr <- as.data.table(table.hosp.shieldgr)
table.hosp.shieldgr[, rateratio := or.ci(coef, `se(coef)`)]
table.hosp.shieldgr[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
table.hosp.shieldgr <- table.hosp.shieldgr[, .(rateratio, pvalue)]
table.hosp.shieldgr <- data.table(effect=effect, # freqs.riskgr,
                               table.hosp.shieldgr)
table.hosp.shieldgr[, effect := gsub("^care.home", "", effect)]
table.hosp.shieldgr[, effect := gsub("^shield\\.group", "", effect)]
table.hosp.shieldgr[, effect := gsub("TRUE$", "", effect)]
table.hosp.shieldgr[, effect := gsub(":vax14\\.dose", ": vax dose", effect)]
table.hosp.shieldgr[, effect := replace.names(effect)]
table.hosp.shieldgr <- table.hosp.shieldgr[c(14:19, 22:27), ]
table.hosp.shieldgr[, dose := rep(c(1, 2), each=6)]
table.hosp.shieldgr[, effect := gsub(":.*", "", effect)]
table.hosp.shieldgr <- dcast(table.hosp.shieldgr, effect ~ dose, value.var=c("rateratio", "pvalue"))
table.hosp.shieldgr <- table.hosp.shieldgr[c(5, 6, 4, 3, 2, 1), ]
table.hosp.shieldgr <- table.hosp.shieldgr[, c(1, 2, 4, 3, 5)]

##  rate ratios by risk group x dose-product
table.hosp.vaxgr <- summary(clogit(data=cc.all[(casegroup=="A" | casegroup=="B") & #care.home=="Independent" & 
                           specimen_date >= as.Date("2020-12-01")],
                           formula=CASE ~ care.home + hh.over18 + numdrugs + inpat.recent +
                               listedgr3/vaxgr + strata(stratum)))
table.hosp.vaxgr <- table.hosp.vaxgr$coefficients
effect <- rownames(table.hosp.vaxgr)
table.hosp.vaxgr <- as.data.table(table.hosp.vaxgr)
table.hosp.vaxgr[, rateratio := or.ci(coef, `se(coef)`)]
table.hosp.vaxgr[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
table.hosp.vaxgr <- table.hosp.vaxgr[, .(rateratio, pvalue)]
table.hosp.vaxgr <- data.table(effect=effect, # freqs.vaxgr,
                               table.hosp.vaxgr)
table.hosp.vaxgr[, effect := gsub("^care.home", "", effect)]
table.hosp.vaxgr[, effect := gsub("^listedgr3", "", effect)]
table.hosp.vaxgr[, effect := gsub("TRUE$", "", effect)]
table.hosp.vaxgr[, effect := gsub(":vaxgr", ": ", effect)]
table.hosp.vaxgr[, effect := replace.names(effect)]
table.hosp.vaxgr <- table.hosp.vaxgr[c(1:6,
                             7, 13, 10, 16,
                             8, 14, 11, 17,
                             9, 15, 12, 18), ]
table.hosp.vaxgr <- table.hosp.vaxgr[7:18, ]
table.hosp.vaxgr[, dose := rep(c(1, 2), 6)]
table.hosp.vaxgr[, effect := gsub(". doses? ", "", effect)]
table.hosp.vaxgr <- dcast(table.hosp.vaxgr, effect ~ dose, value.var=c("rateratio", "pvalue"))
table.hosp.vaxgr <- table.hosp.vaxgr[c(6, 5, 4, 3, 2, 1), ]
table.hosp.vaxgr <- table.hosp.vaxgr[, c(1, 2, 4, 3, 5)]

###########################################################
riskgr.severe.vax2 <- with(cc.all[CASE==1 & casegroup=="A" & vax14.dose==1 & 
                                  specimen_date >= as.Date("2020-12-01")],
                           table(listedgr3))
pct.riskgr.severe.vax2 <- round(100 * (1 - (riskgr.severe.vax2[1]) / sum(riskgr.severe.vax2)))
pct.noriskgr.severe.vax2 <- round(100 * (riskgr.severe.vax2[1]) / sum(riskgr.severe.vax2))
pct.modriskgr.severe.vax2 <- round(100 * (riskgr.severe.vax2[2]) / sum(riskgr.severe.vax2))
pct.cev.severe.vax2 <- round(100 * (riskgr.severe.vax2[3]) / sum(riskgr.severe.vax2))

riskgr.hosp.vax2 <- with(cc.all[CASE==1 & (casegroup=="A" | casegroup=="B")
                                  & vax14.dose==1 & 
                                  specimen_date >= as.Date("2020-12-01")],
                           table(listedgr3))
pct.riskgr.hosp.vax2 <- round(100 * (1 - (riskgr.hosp.vax2[1]) / sum(riskgr.hosp.vax2)))

cc.all[, months.sincetransplant := as.integer(specimen_date - date.lasttransplant) * 12 / 365.25]
    
table.severe.vax2 <-
    tabulate.freqs.regressions(data=cc.all[casegroup=="A" & vax14.dose==1 & 
                                           specimen_date >= as.Date("2020-12-01")],
                               varnames=c("care.home", "qSIMD.integer", "hh.over18",
                                          "numicdchapters", "numdrugs.notcv", "inpat.recent",
                                          "shield.group"))

table.hosp.vax2 <-
    tabulate.freqs.regressions(data=cc.all[(casegroup=="A" | casegroup=="B") &
                                           vax14.dose==1 & 
                                           specimen_date >= as.Date("2020-12-01")],
                               varnames=c("care.home", "qSIMD.integer", "hh.over18",
                                          "numicdchapters", "numdrugs.notcv", "inpat.recent",
                                          "shield.group"))

model.hosp.vax2 <-
    summary(clogit(data=cc.all[(casegroup=="A" | casegroup=="B") &
                                           vax14.dose==1 & 
                               specimen_date >= as.Date("2020-12-01")],
                   formula=CASE ~ hh.over18 + 
                       numicdchapters + numdrugs.notcv + inpat.recent + 
                       shield.group + strata(stratum)))

model.hosp.shieldonly.vax2 <-
    summary(clogit(data=cc.all[(casegroup=="A" | casegroup=="B") &
                                           vax14.dose==1 & 
                               specimen_date >= as.Date("2020-12-01")],
                   formula=CASE ~ shield.any + strata(stratum)))

########### plots ##################################

load(paste0(datadir, "shieldcohort.models.RData"))
## this should load data shield.all

p.anycase <- ggplot(data=shield.all, aes(x=date, y=probmonth, color=casegr)) +
    geom_line() +
    labs(x=paste0("Presentation date: start of ", interval.length, "-day interval"),
         y="Incidence rate per month") +
    scale_y_continuous(breaks=seq(0, 0.012, by=0.002), limits=c(0, 0.011), expand=c(0, 0)) +  
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = as.Date("2020-03-01"),
                                                  to = lastdate, by = "month"),
                                         "%d %b")
                             ),
                 limits=c(as.Date("2020-12-01"), lastdate))  +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.6, 0.7)) 

p.anycase

cc.severe.num <- cc.all[casegroup=="A"  & 
                        specimen_date >= as.Date("2020-12-01"), .N, by=CASE]
cc.severe.vax2.num <- cc.all[casegroup=="A" & vax14.dose==1 & 
                             specimen_date >= as.Date("2020-12-01"), .N, by=CASE]
cc.hosp.vax2.num <- cc.all[(casegroup=="A" | casegroup=="B") & vax14.dose==1 & 
                             specimen_date >= as.Date("2020-12-01"), .N, by=CASE]

severe.pct <- paste0(round(100 * prop.table(with(cc.all[casegroup=="A"  &
                                                        CASE==1 &  
                                   specimen_date >= as.Date("2020-12-01")],
                            table(listedgr3)))), "%")
severe.byvax.pct <- paste0(round(100 * prop.table(with(cc.all[casegroup=="A"  &
                                                        CASE==1 &  
                                   specimen_date >= as.Date("2020-12-01")],
                            table(vax14.factor)))), "%")

severe.vax2.pct <- paste0(round(100 * prop.table(with(cc.all[casegroup=="A"  &
                                                             vax14.dose==1 & CASE==1 &  
                                   specimen_date >= as.Date("2020-12-01")],
                            table(listedgr3)))), "%")
severe.pct
severe.vax2.pct 
severe.byvax.pct

load(paste0(datadir, "severe2vax.RData")) 

rmarkdown::render("vaxshieldupdate.Rmd")
rmarkdown::render("vaxtrend.Rmd")
rmarkdown::render("riskfactorsupdate.Rmd")


rowSums(with(cc.all[CASE==1 & vax14.factor=="2"], table(floor((specimen_date - date.lasttransplant) / 180), casegroup)))

########################################################

cc.all[, vax14.1dose := as.integer(vax14.dose > 0)]

table.riskgr16March <- summary(clogit(data=cc.all[casegroup=="A" & #care.home=="Independent" & 
                           specimen_date <= as.Date("2021-03-16") & 
                           specimen_date >= as.Date("2020-12-01") &
                           vax14.factor != "2"],
                           formula=CASE ~ care.home + qSIMD.integer + 
                                  #       hh.over18gr + hh.schoolagegr + 
                               #numdrugs.notcv +
                               #inpat.recent +
                               shield.any/vax14.1dose + strata(stratum)))
table.riskgr16March <- table.riskgr16March$coefficients
effect <- rownames(table.riskgr16March)
table.riskgr16March <- as.data.table(table.riskgr16March)
table.riskgr16March[, rateratio := or.ci(coef, `se(coef)`)]
table.riskgr16March[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
table.riskgr16March <- table.riskgr16March[, .(rateratio, pvalue)]
table.riskgr16March <- data.table(effect=effect, # freqs.riskgr,
                               table.riskgr16March)
table.riskgr16March[, effect := gsub("^care.home", "", effect)]
table.riskgr16March[, effect := gsub("^listedgr3", "", effect)]
table.riskgr16March[, effect := gsub("TRUE$", "", effect)]
table.riskgr16March[, effect := gsub(":vax14\\.factor", ": vax dose ", effect)]
table.riskgr16March[, effect := replace.names(effect)]
table.riskgr16March <- table.riskgr16March[grep(":", effect)]
table.riskgr16March

png("p,rateratio.png")
p.rateratio
dev.off()

table.shield.schoolage <- with(cc.all, table(shieldedonly.group, hh.schoolage.any))
table.shield.schoolage <- rbind(table.shield.schoolage, colSums(table.shield.schoolage))
table.shield.schoolage <- paste.rowpercent(table.shield.schoolage)

