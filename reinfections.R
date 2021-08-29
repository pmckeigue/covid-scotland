library(survival)
library(data.table)
library(ggplot2)

source("helperfunctions.R")

datadir <- "./data/2021-07-28/"

if(!exists("cc.all")) {
    load(paste0(datadir, "cc.all.RData"))
}

load(paste0(datadir, "reinfections.RData"))
load(paste0(datadir, "ecoss.firstpos.RData"))

interval.days <- 7

## merge reinfections with first infections to set up a time to event model for reinfection, conditioning on baseline calendar time and time since first infection
## reinfections should have specimen dates of 1st and 2nd infection
setkey(cc.all, anon_id)
setkey(reinfections, anon_id)
setnames(reinfections, "specimen_date", "reinfection_date", skip_absent=TRUE)
reinfections <- reinfections[, .(anon_id, reinfection_date)]


## use rapid to derive time from specimen date to hospitalisation 
rapid.filename <- paste0(datadir, "CC_RAPID_ANON_2021-07-28.rds")
rapid <- RDStodt(rapid.filename, keyname="anon_id")
#rapid[, admit_time := NULL]
#rapid[, disc_time := NULL]
setnames(rapid, "admission_date", "AdmissionDate.rapid", skip_absent=TRUE)
setnames(rapid, "discharge_date", "DischargeDate.rapid", skip_absent=TRUE)
lastdate.rapid <- max(c(rapid$AdmissionDate.rapid, rapid$DischargeDate.rapid), na.rm=TRUE) 
rapid[is.na(DischargeDate.rapid), DischargeDate.rapid := lastdate.rapid]
rapid[, joindate := AdmissionDate.rapid]
setkey(rapid, anon_id, joindate)
reinfections[, joindate := reinfection_date]
setkey(reinfections, anon_id, joindate)

## left join of reinfections with next admission date in rapid
## so roll backwards up to 28 days from admission date
## variable AdmissionDate.rapid is the first admission date within 28 days
## variable daystoadmission is the number of days from specimen date to admission 
## excludes admission dates before specimen date 
cc.nexthosp <- rapid[reinfections, roll=-28] # roll backwards from admission date
## retain only records for which an admission date after specimendate was found in rapid
cc.nexthosp <- cc.nexthosp[!is.na(AdmissionDate.rapid)]
cc.nexthosp[, daystoadmission := as.integer(AdmissionDate.rapid - reinfection_date)]
setorder(cc.nexthosp, anon_id, reinfection_date, daystoadmission)
setkey(cc.nexthosp, anon_id, reinfection_date)
cc.nexthosp <- unique(cc.nexthosp, by=key(cc.nexthosp))
setkey(reinfections, anon_id, reinfection_date)
reinfections <- cc.nexthosp[, .(anon_id, AdmissionDate.rapid, reinfection_date,
                                daystoadmission)][reinfections]

atleast90day.firstpos <- 
    reinfections[cc.all[testpositive.case==TRUE,
                        .(anon_id, specimen_date, date_of_death, covid_ucod,
                          age_years, agegr20, sex, care.home, hh.over18,
                          hh.schoolage, occup, listedgr3, shield.any,
                          adm.within28,
                          vaxdate_1, vaxdate_2)]]
atleast90day.firstpos[, adm.within28.reinf := as.integer(!is.na(daystoadmission) &
                                                         daystoadmission <= 28)]
setnames(atleast90day.firstpos, "specimen_date", "firstinfection_date", skip_absent=TRUE)


## 24 reinfections fatal within 28 days
atleast90day.firstpos[, fatal.reinf := as.integer(!is.na(date_of_death) &
                                                  as.integer(date_of_death -
                                                             reinfection_date) <= 28)]
atleast90day.firstpos[is.na(reinfection_date), fatal.reinf := NA]

atleast90day.firstpos[covid_ucod==1, fatal.reinf := 2]
atleast90day.firstpos[, fatal.reinf := car::recode(fatal.reinf,
          recode="0='Non-fatal';
                  1='Fatal, other underlying cause';
                  2='Fatal, COVID-19 as underlying cause'", 
          as.factor=TRUE,
          levels=c("Non-fatal",
                   "Fatal, other underlying cause",
                   "Fatal, COVID-19 as underlying cause"))]

table.fatal <- with(atleast90day.firstpos[!is.na(reinfection_date)],
                    table(agegr20, fatal.reinf, exclude=NULL))
table.fatal <- rbind(table.fatal, colSums(table.fatal))
table.fatal <- paste.rowpercent(table.fatal, digits=1)
rownames(table.fatal)[5] <- "All"

## merge reinfections with deaths, rapid, sicsag to code as hospitalized, critical care, fatal

## set up dataset as cohort: enter at 90 days after first infection,
atleast90day.firstpos[, entrydate := firstinfection_date + 90]
atleast90day.firstpos <- atleast90day.firstpos[is.na(date_of_death) |
                                               date_of_death > entrydate]
atleast90day.firstpos[, event := as.integer(!is.na(reinfection_date))]
atleast90day.firstpos[, exitdate := as.Date("2021-07-28")]
atleast90day.firstpos[!is.na(date_of_death), exitdate := pmin(exitdate, date_of_death)]
atleast90day.firstpos[!is.na(reinfection_date), exitdate := pmin(exitdate, reinfection_date)]
atleast90day.firstpos <- atleast90day.firstpos[entrydate < exitdate]
atleast90day.firstpos[, tstart := as.integer(entrydate)]
atleast90day.firstpos[, tstop := as.integer(exitdate)]
atleast90day.firstpos[, tobs := tstop - tstart]

## split into person-time intervals
min.tstart <- min(atleast90day.firstpos$tstart)
firstpos.split <- as.data.table(survSplit(Surv(tstart, tstop, event) ~ .,
                                          data=atleast90day.firstpos,
                                          cut=min.tstart + interval.days * (1:100),
                                          episode=paste0("interval", interval.days, "day")))
firstpos.split[, datestart := as.Date(tstart, origin=as.Date("1970-01-01"))]
firstpos.split[, datestop := as.Date(tstop, origin=as.Date("1970-01-01"))]
firstpos.split[, tobs := tstop - tstart]
firstpos.split[, months.sincefirstinfection := (tstart - as.integer(firstinfection_date)) * 12 / 365.25]
summary(firstpos.split[, .(anon_id, entrydate, exitdate, days.sincefirstinfection, tobs, event)])

# time-updated vaccination status
firstpos.split[, vax14.1dose := as.integer(!is.na(vaxdate_1) &
                                           as.integer(datestart - vaxdate_1) >= 14)]
firstpos.split[, vax14.2dose := as.integer(!is.na(vaxdate_2) &
                                           as.integer(datestart - vaxdate_2) >= 14)]
firstpos.split[, vax14.factor := as.factor(vax14.1dose + vax14.2dose)]


## Cox model for time to reinfection
cox.varnames <- c("age_years", "sex", "care.home", "occup",
                  "hh.over18", "hh.schoolage",
                  "listedgr3", "months.sincefirstinfection", "vax14.factor")
formulastring <- paste0("Surv(tstart, tstop, event) ~ ", paste(cox.varnames, collapse=" + "))
cox.formula <- with(firstpos.split, as.formula(formulastring))
reinf.coeffs <- summary(coxph(cox.formula,
                              data=firstpos.split))$coefficients

reinf.coeffs <- data.table(effect=rownames(reinf.coeffs), reinf.coeffs)
reinf.coeffs[, rateratio := or.ci(coef, `se(coef)`)]
reinf.coeffs[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
reinf.coeffs <- reinf.coeffs[, .(effect, rateratio, pvalue)]

levels.coeffs <- levels.toshow(varnames=cox.varnames,
                               dt=firstpos.split)

reinf.coeffs <- pad.coeffs(reinf.coeffs, levels.coeffs)

reinf.freqs <- univariate.tabulate(varnames=cox.varnames, 
                    outcome="event",
                    data=firstpos.split, drop.reflevel=FALSE)
colnames(reinf.freqs) <- gsub("^0", "No event in interval", colnames(reinf.freqs))
colnames(reinf.freqs) <- gsub("^1", "Reinfections", colnames(reinf.freqs))

table.reinf.coeffs <- cbind(reinf.coeffs[, 1], reinf.freqs, reinf.coeffs[, -1])

table.reinf.coeffs[, effect := gsub("age_years", "Age", effect)] 
gc()
############################################################

reinf.poissonmodel <- glm(event ~ offset(log(tobs)) + age_years + sex + care.home + occup +
                              shield.any + vax14.factor +
                              splines::ns(tstart, 6),
                          family="poisson", data=firstpos.split)
 
firstpos.split[, reinf.loglambda := as.numeric(predict(reinf.poissonmodel,
                                                       newdata=firstpos.split,
                                                      type="link")) - log(tobs)]
rm(reinf.poissonmodel)
gc()

firstpos.split[, reinf.probmonth := 1 - exp(-exp(reinf.loglambda) * 365.25 / 12)]

reinf.anycase  <- firstpos.split[, .(probmonth=mean(reinf.probmonth, na.rm=TRUE)),
                                 by=c("datestart", "vax14.factor")]

k <- 5
reinf.rmeans <- rollmean.bydate.grouped(reinf.anycase, datevar="datestart",
                                        groupvar="vax14.factor", xvar="probmonth", k=k)


p.reinf <- ggplot(data=reinf.rmeans, aes(x=date, y=avdaily, color=vax14.factor)) +
    geom_line() +
    labs(x=paste0("Presentation date: start of ", interval.days, "-day interval, rolling mean of ", k, " intervals"),
         y="Incidence rate per month") +
    scale_y_continuous(breaks=seq(0, 0.008, by=0.001), limits=c(0, 0.008), expand=c(0, 0)) +  
    scale_x_date(breaks = seq.Date(from = as.Date("2020-05-01"),
                                   to = as.Date("2021-07-28"), by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = as.Date("2020-05-01"),
                                                  to = as.Date("2021-07-28"), by = "month"),
                                         "%d %b")
                             ),
                 limits=c(as.Date("2020-05-01"), as.Date("2021-07-28"))) +
    labs(color="Vaccine doses") +
    theme(legend.position = c(0.2, 0.8)) 

p.reinf

atleast90day.firstpos[, mean(tobs)] * 12 / 365.25

rmarkdown::render("reinfections.Rmd")
