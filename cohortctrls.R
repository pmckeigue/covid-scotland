library(survival)
library(data.table)
library(ggplot2)

source("helperfunctions.R")

lastdate <- as.Date("2021-09-22")

datadir <- "./data/2021-09-22/"
if(!exists("cc.all")) {
    load(paste0(datadir, "cc.all.RData"))
}

## select cohort of controls
cohort2 <- cc.all[CASE==0 & care.home=="Independent",
                  .(anon_id, specimen_date, 
                    age_years, sex, listedgr3, 
                    vaxdate_1, date_of_death, underlying_cause_of_death)]
rm(cc.all)

## don't exclude those in cohort 1: on first testing positive they transfer after 90 days from cohort 2 to cohort 1 exit 
setnames(cohort2, "specimen_date", "controlsampling_date")
## select unique anon_ids
setorder(cohort2, "controlsampling_date")
cohort2 <- unique(cohort2, by="anon_id")
setkey(cohort2, anon_id)
## restrict to those who had not received first dose on or before controlsampling_date
cohort2 <- cohort2[is.na(vaxdate_1) | vaxdate_1 > controlsampling_date]

###################################################################################

cohort2[, event := as.integer(!is.na(date_of_death))]
cohort2[, cancer := !is.na(date_of_death) &
              grepl("^[CD]", underlying_cause_of_death)]
cohort2[, circulatory := !is.na(date_of_death) &
              grepl("^I", underlying_cause_of_death)]
cohort2[, external := !is.na(date_of_death) &
              grepl("^[VWX]", underlying_cause_of_death)]

cohort2[, tstart := as.integer(controlsampling_date)]
cohort2[, tstop := as.integer(pmin(date_of_death,
                                   max(date_of_death, na.rm=TRUE), na.rm=TRUE))]
cohort2 <- cohort2[tstart < tstop]
cohort2 <- cohort2[age_years >= 60]
        
interval.days <- 7

## split cohort2 into person-time intervals by calendar time to allow time-updated recent vax datus as covariate
min.tstart <- min(cohort2$tstart)
cohort2.split <- as.data.table(survSplit(Surv(tstart, tstop, event) ~ .,
                                         data=cohort2[, .(tstart, tstop, event, cancer,
                                                          circulatory, external, 
                                                          age_years, sex,
                                                          listedgr3, vaxdate_1)],
                                          cut=min.tstart + interval.days * (1:100),
                                          episode=paste0("interval", interval.days, "day")))

cohort2.split[, datestart := as.Date(tstart, origin=as.Date("1970-01-01"))]
cohort2.split[, vax1.last28 := as.integer(!is.na(vaxdate_1) &
                                          datestart - vaxdate_1 <= 28 &
                                          datestart - vaxdate_1 > 0)]
cohort2.split[, vaxdate_1 := NULL]

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))
gc()


cohort2.split[, cancerdeath := as.integer(event==1 & cancer)]
cohort2.split[, circulatorydeath := as.integer(event==1 & circulatory)]
cohort2.split[, externaldeath := as.integer(event==1 & external)]
cohort2.split[, otherdeath := as.integer(event==1 & !(cancer | circulatory | external))]

cox.varnames <- c("age_years", "sex",  
                  "listedgr3", "vax1.last28")

######################################################################
formulastring <- paste0("Surv(tstart, tstop, cancerdeath) ~ ", paste(cox.varnames, collapse=" + "))
coeffs.cancer <- summary(coxph(as.formula(formulastring),
                              data=cohort2.split))$coefficients
######################################################################
formulastring <- paste0("Surv(tstart, tstop, otherdeath) ~ ", paste(cox.varnames, collapse=" + "))
coeffs.other <- summary(coxph(as.formula(formulastring),
                              data=cohort2.split))$coefficients
########################################################################
formulastring <- paste0("Surv(tstart, tstop, externaldeath) ~ ", paste(cox.varnames, collapse=" + "))
coeffs.external <- summary(coxph(as.formula(formulastring),
                              data=cohort2.split))$coefficients
coeffs.external
#######################################################################
formulastring <- paste0("Surv(tstart, tstop, circulatorydeath) ~ ", paste(cox.varnames, collapse=" + "))
coeffs.circ <- summary(coxph(as.formula(formulastring),
                              data=cohort2.split))$coefficients


coeffs <- data.table(cause=c("Cancer", "Circulatory", "External", "Other"),
                     rbind(coeffs.cancer[5, ],
                           coeffs.circ[5, ],
                           coeffs.external[5, ],
                           coeffs.other[5, ]))
coeffs[, oddsratio := or.ci(coef, `se(coef)`)]
