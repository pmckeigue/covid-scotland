library(data.table)
library(survival)
library(ggplot2)

source("helperfunctions.R")

rollmean.dtN <- function(dt, k) {
    ## FIXME: add a by= argument
    N.dates <- dt[, .N, by=SPECIMENDATE]
    setkey(N.dates, SPECIMENDATE)
    all.dates <- data.table(date=seq(from=min(N.dates$SPECIMENDATE),
                                     to=max(N.dates$SPECIMENDATE),
                                     by=1))
    setkey(all.dates, date)

    ## add rows for missing dates by a left join
    all.dates <- N.dates[all.dates]
    all.dates[is.na(N), N := 0]
    return(data.table(date=dt$SPECIMENDATE,
                      avdaily=zoo::rollmean(dt$N, k=k, fill=NA)))
}

rollsum.datewin <- function(dt, k, datevar) {
    ## returns a table of rolling sums of width k centred on date
    ## dt is a dataset of totals (N) by date, which may have no rows for some dates in the
    ## range of interest
    ## add rows for missing dates by a left join of all.dates with dt
    all.dates <- data.table(date=seq(from=min(dt[[datevar]]),
                                     to=max(dt[[datevar]]),
                                     by=1))
    setkey(all.dates, date)
    # setkeyv(dt, datevar) can't set physical key here
    all.dates <- dt[all.dates]

    return(data.table(date=all.dates[[datevar]],
                      rsum=zoo::rollsum(all.dates[, N], k=k, fill=0, align="center")))
}

format.estcipv <- function(x.ci, x.pval) {
    x <- gsub(" \\(", " \\(95\\% CI ", as.character(x.ci))
    x <- gsub(", ", " to ", x)
    x <- gsub("\\)", paste0(", _p_=\\", x.pval, "\\)"), x)
    x <- gsub("times", "\\\\times", x)
    return(x)
}

format.estci <- function(x.ci) {
    x <- gsub(" \\(", " \\(95\\% CI ", as.character(x.ci))
    x <- gsub(", ", " to ", x)
    #x <- gsub("\\)", paste0(", _p_=\\", x.pval, "\\)"), x)
    #x <- gsub("times", "\\\\times", x)
    return(x)
}

datadir <- "data/2021-07-28/"
if(!exists("cc.all")) {
    load(paste0(datadir, "cc.all.RData"))
}

source("rrtimewin.R")

interval.length <- 28  # for poisson regression of shielding cohort
lastdate.cases <- as.Date("2021-07-12")

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
# rate ratios by risk group x dose
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
                               shield.group/vax14.dose + strata(stratum)))
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

# rate ratios by shield group x dose
table.hosp.shieldgr <- summary(clogit(data=cc.all[(casegroup=="A" | casegroup=="B") & #care.home=="Independent" & 
                           #specimen_date <= as.Date("2021-03-16") &
                           specimen_date >= as.Date("2020-12-01")],
                           formula=CASE ~ care.home + hh.over18 + numdrugs + inpat.recent +
                               shield.group/vax14.dose + strata(stratum)))
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
###########################################################

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

p.anycase <- ggplot(data=shield.all, aes(x=date, y=probmonth, color=casegr)) +
    geom_line() +
    labs(x=paste0("Presentation date: start of ", interval.length, "-day interval"),
         y="Incidence rate per month") +
    scale_y_continuous(breaks=seq(0, 0.01, by=0.002), limits=c(0, 0.01), expand=c(0, 0)) +  
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = as.Date("2021-07-25"), by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = as.Date("2020-03-01"),
                                                  to = as.Date("2021-07-25"), by = "month"),
                                         "%d %b")
                             ),
                 limits=c(as.Date("2020-12-01"), as.Date("2021-07-25")))  +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.6, 0.7)) 

p.anycase

load(paste0(datadir, "severe2vax.RData")) 

rmarkdown::render("vaxshield230821.Rmd")
rmarkdown::render("vaxtrend.Rmd")

########################################################

table.riskgr16March <- summary(clogit(data=cc.all[casegroup=="A" & #care.home=="Independent" & 
                           specimen_date <= as.Date("2021-03-16")& 
                           specimen_date >= as.Date("2020-12-01")],
                           formula=CASE ~ care.home + hh.over18 + numdrugs + inpat.recent +
                               listedgr3/vax14.dose + strata(stratum)))
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
table.riskgr16March <- table.riskgr16March[c(1:6, 7, 10, 8, 11, 9, 12), ]

png("p,rateratio.png")
p.rateratio
dev.off()
