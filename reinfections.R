library(survival)
library(data.table)
library(ggplot2)

source("helperfunctions.R")

lastdate <- as.Date("2021-09-22")

datadir <- "./data/2021-09-22/"
if(!exists("cc.all")) {
    load(paste0(datadir, "cc.all.RData"))
}

covid.specialtycodes <- "^A[16BQ]"
noncovid.specialtycodes <- "^A[289DGHMR]|^[CEFHJ]"

load(paste0(datadir,"ecoss.alltests.RData"))
ecoss <- unique(ecoss)
setkey(ecoss, SpecimenDate)

## select definite positive records in ecoss based on > 1 test or Ct < 30
ecoss.pos <- ecoss[ecoss.result=="Positive"]
ecoss.firstpos <- ecoss.pos[!duplicated(anon_id)]
setkey(ecoss.firstpos, anon_id, SpecimenDate)

diff0 <- function(x) c(0, diff(x))
ecoss.pos[, interval := diff0(SpecimenDate), by=anon_id]
ecoss.pos[, testnum := 1:.N, by=anon_id]
ecoss.pos[, numtests := .N, by=anon_id]
## select tests 1 and 2 on those with > 1 test
ecoss.first2tests <- ecoss.pos[numtests > 1 & testnum <= 2]
## reshape to wide
ecoss.first2tests.wide <- dcast(ecoss.first2tests, anon_id ~ testnum,
                                value.var="interval")
colnames(ecoss.first2tests.wide)[3] <- "intervalto2nd"
## keep records with 2nd test within 14 days of 1st test
ids.twicepos <- ecoss.first2tests.wide[intervalto2nd <= 14, anon_id]
ids.ctlt30 <- ecoss.firstpos[!is.na(ct.lt30.2genes) & ct.lt30.2genes == 1, anon_id]

## restrict to ids.twicepos or Ct < 30
ids.definite <- unique(c(ids.twicepos, ids.ctlt30)) 


## now restrict ecoss.firstpos to those in ids.definite

## create table of individuals with definite first positive result
ecoss.firstpos.definite <- ecoss.firstpos[anon_id %in% ids.definite]
setkey(ecoss.firstpos.definite, anon_id, SpecimenDate)
setkey(cc.all, anon_id, specimen_date)
ecoss.firstpos.definite <-
    cc.all[, .(anon_id, specimen_date, date_of_death)][ecoss.firstpos.definite]

# restrict to those who survived 90 days
num.firstpos <- length(unique(ecoss.firstpos.definite$anon_id))
ecoss.firstpos.definite <-
    ecoss.firstpos.definite[is.na(date_of_death) |
                            as.integer(date_of_death - specimen_date) >= 90]
num.firstpos.gt90 <- length(unique(ecoss.firstpos.definite$anon_id))

## generate indicators for definite firstpos in ecoss table and ecoss.pos
ecoss[, firstpos.definite := as.integer(anon_id %in% ecoss.firstpos.definite$anon_id)]
ecoss.pos[, firstpos.definite := as.integer(anon_id %in% ecoss.firstpos.definite$anon_id)]

## CDC criterion - reinfections
## select reinfections based on 90-day limit
reinfections <- ecoss.pos[firstpos.definite==1 & interval >=90]
setnames(reinfections, "SpecimenDate", "specimen_date", skip_absent=TRUE)
save(reinfections, file=paste0(datadir, "reinfections.RData"))

cat(nrow(ecoss.firstpos.definite), "definite first positives, ",
    nrow(reinfections), "reinfections\n")

setnames(reinfections, "specimen_date", "reinfection_date", skip_absent=TRUE)
reinfections <- reinfections[, .(anon_id, reinfection_date)]

#########################################################################

## use rapid to derive time from specimen date to hospitalisation 
load(paste0(datadir, "rapid.RData"))

rapid[, joindate := AdmissionDate.rapid]
setkey(rapid, anon_id, joindate)
reinfections[, joindate := reinfection_date]
setkey(reinfections, anon_id, joindate)

## left join of reinfections with next admission date in rapid
## so roll backwards up to 28 days from admission date
## variable AdmissionDate.rapid is the first admission date within 28 days
## variable daystoadmission is the number of days from specimen date to admission 
## excludes admission dates before specimen date 
cc.nexthosp <- rapid[reinfections, roll=-14] # roll backwards from admission date
## retain only records for which an admission date after specimendate was found in rapid
cc.nexthosp <- cc.nexthosp[!is.na(AdmissionDate.rapid)]
cc.nexthosp[, daystoadmission := as.integer(AdmissionDate.rapid - reinfection_date)]
setorder(cc.nexthosp, anon_id, reinfection_date, daystoadmission)
setkey(cc.nexthosp, anon_id, reinfection_date)
cc.nexthosp <- unique(cc.nexthosp, by=key(cc.nexthosp))
setkey(reinfections, anon_id, reinfection_date)
reinfections <- cc.nexthosp[, .(anon_id, AdmissionDate.rapid, reinfection_date,
                                daystoadmission, specialty)][reinfections]

## 2633 reinfections, 61 admitted within 14 days of positive test

## left join all test-positive cases in cc.all with reinfections to set up a cohort.defpos of all those at risk for reinfection
setkey(cc.all, anon_id)
setkey(reinfections, anon_id)
cohort.defpos <- 
    reinfections[cc.all[testpositive.case==TRUE,
                        .(anon_id, specimen_date, stratum, date_of_death, underlying_cause_of_death,
                          covid_ucod,
                          age_years, agegr3, agegr20, sex, care.home, hh.over18gr3,
                          hh.schoolage, occup, sector, listedgr3, shield.any,
                          adm.within14,
                          vaxdate_1, vaxdate_2)]]

## 508621 test-positive cases in cc.all
## restrict to those with definite positive test first time who survived at least 90 days
cohort.defpos <- cohort.defpos[anon_id %in% ecoss.firstpos.definite$anon_id]
setnames(cohort.defpos, "specimen_date", "firstinfection_date", skip_absent=TRUE)
cohort.defpos[, case.reinf := as.integer(!is.na(reinfection_date))]
cohort.defpos[, sexf := relevel(sex, ref="Male")]
table(cohort.defpos$case.reinf)
numtestpos.definite <- nrow(cohort.defpos)
## 362130 definite test-positive cases, 2598 reinfections

## create variable for vax status at first infection 
cohort.defpos[, vaxedfirst := "1 dose"]
cohort.defpos[is.na(vaxdate_1) | vaxdate_1 > firstinfection_date, vaxedfirst := "Not vaccinated"]
cohort.defpos[!is.na(vaxdate_2) & vaxdate_2 < firstinfection_date, vaxedfirst := "2 doses"]
cohort.defpos[, vaxedfirst := factor(vaxedfirst,
                                     levels=c("Not vaccinated", "1 dose", "2 doses"))]
table(cohort.defpos$vaxedfirst)

## assign variable for admission within 14 days of reinfection
cohort.defpos[, adm.within14.reinf := as.integer(!is.na(daystoadmission) &
                                                         daystoadmission <= 14)]
table(cohort.defpos$adm.within14.reinf, exclude=NULL)

## 6 reinfections fatal within 28 days
cohort.defpos[!is.na(date_of_death),
                             fatal28.reinf := as.integer(as.integer(date_of_death -
                                                                    reinfection_date) <= 28)]
cohort.defpos[is.na(date_of_death) | is.na(reinfection_date), fatal28.reinf := 0]
                                                                    

## exclude 16 deaths certified with COVID beyond 28 days
with(cohort.defpos, table(covid_ucod, fatal28.reinf, exclude=NULL))

## 71 hosp or fatal
cohort.defpos[, hosporfatal := as.integer(adm.within14.reinf==1 | fatal28.reinf==1)]

     
table(cohort.defpos$case.reinf, exclude=NULL)
table(cohort.defpos$adm.within14.reinf, exclude=NULL)
cohort.defpos[, outcome := case.reinf + adm.within14.reinf]
cohort.defpos[fatal28.reinf==1, outcome := 3]
cohort.defpos[fatal28.reinf==1 & covid_ucod==0, outcome := 4]

cohort.defpos[, outcome := car::recode(outcome,
          recode="0='No detected reinfection';
                  1='Test-positive, not hospitalised';
                  2='Hospitalised, non-fatal';
                  3='Fatal, COVID-19 as underlying cause';
                  4='Fatal, other underlying cause'",
          as.factor=TRUE,
          levels=c("No detected reinfection",
                   "Test-positive, not hospitalised",
                   "Hospitalised, non-fatal",
                   "Fatal, COVID-19 as underlying cause",
                   "Fatal, other underlying cause"))]

## set up entry date, event, exitdate: enter at 90 days after first infection or control sampling
cohort.defpos[, entrydate := firstinfection_date + 90]
cohort.defpos <- cohort.defpos[is.na(date_of_death) | date_of_death > entrydate] 
print(nrow(cohort.defpos))
with(cohort.defpos, table(vaxedfirst, exclude=NULL))
## 362123 still under observation


## event defined as detected reinfection
cohort.defpos[, event := as.integer(!is.na(reinfection_date))]
## 120 hospitalised / fatal reinfections
cohort.defpos[, eventhosp := as.integer(!is.na(reinfection_date) & (adm.within14.reinf==1 | fatal28.reinf==1))]

cohort.defpos[, exitdate := lastdate]
cohort.defpos[!is.na(date_of_death), exitdate := pmin(exitdate, date_of_death)]
cohort.defpos[event==1, exitdate := pmin(exitdate, reinfection_date)]
with(cohort.defpos, table(vaxedfirst, exclude=NULL))

## this line loses about 95% of those vaxed by first infection
## because their entrydate (firstinfection_date + 90) was later than lastdate ??
## i.e. vaccinated twice then infected after about 25 June 2021. 

cohort.defpos <- cohort.defpos[entrydate < exitdate]
with(cohort.defpos, table(vaxedfirst, exclude=NULL))

cohort.defpos[, tstart := as.integer(entrydate)] ## this has to be an integer for survsplit 
cohort.defpos[, tstop := as.integer(exitdate)]
cohort.defpos[, tobs := tstop - tstart] 
cohort.defpos[, tobs.months := tobs * 12 / 365.25] 

cohort.defpos.rates <- cohort.defpos[, .(events=sum(event),
                                         eventshosp=sum(eventhosp), tobs=sum(tobs.months)),
                                     by=vaxedfirst]
cohort.defpos.rates[, rateper1000 := 1000 * events / tobs]
cohort.defpos.rates[, hosp.per1000 := 1000 * eventshosp / tobs]
cohort.defpos.rates

cohort.defpos[, vaxstatus.atexit := 0]
cohort.defpos[!is.na(vaxdate_1) & vaxdate_1 < exitdate &
              (is.na(vaxdate_2) | vaxdate_2 > exitdate), vaxstatus.atexit := 1]
cohort.defpos[!is.na(vaxdate_2) & vaxdate_2 < exitdate, vaxstatus.atexit := 2]

## cohort.defpos: calculate observed time and testing rate for each individuals by vax status
## tobs.unvax is assigned a value < tobs if vaxdate_1 is nonmissing and before exit date
##   if vaxdate_1 > entry date, assign as vaxdate_1 - entry date, otherwise as 0
## if vaxdate_1 is nonmissing and before exit date, then assign tobs.vax1dose
##   starts at pmax(vaxdate_1, entry date), ends at pmin(vaxdate_2, exitdate)
## if vaxdate_2 is nonmissing and before exit date, then assign tobs.vax2dose
##   starts at pmax(vaxdate_2, entrydate), ends at exitdate

cohort.defpos[!is.na(vaxdate_1) & vaxdate_1 <= entrydate, tobs.unvax := 0]
cohort.defpos[!is.na(vaxdate_1) & vaxdate_1 > entrydate, tobs.unvax := vaxdate_1 - entrydate]
cohort.defpos[is.na(vaxdate_1) | vaxdate_1 > exitdate, tobs.unvax := tobs]

cohort.defpos[, tobs.vax1dose := 0]
cohort.defpos[!is.na(vaxdate_1) & vaxdate_1 <= exitdate &
              (is.na(vaxdate_2) | vaxdate_2 > entrydate),
              tobs.vax1dose := pmin(vaxdate_2, exitdate, na.rm=TRUE) -
                  pmax(vaxdate_1, entrydate)]

cohort.defpos[, tobs.vax2dose := 0]
cohort.defpos[!is.na(vaxdate_2) & vaxdate_2 <= exitdate,
              tobs.vax2dose := exitdate - pmax(vaxdate_2, entrydate)]

## calculate observation time in months by vax status
cohort.defpos[, tobs.months := tobs * 12 / 365.25]
cohort.defpos[, tobs.months.unvax := tobs.unvax * 12 / 365.25]
cohort.defpos[, tobs.months.vax1dose := tobs.vax1dose * 12 / 365.25]
cohort.defpos[, tobs.months.vax2dose := tobs.vax2dose * 12 / 365.25]

tobs.cohort.defpos <- cohort.defpos[, .(tobs.months.unvax=sum(tobs.months.unvax),
                                        tobs.months.vax1dose=sum(tobs.months.vax1dose), 
                                        tobs.months.vax2dose=sum(tobs.months.vax2dose),
                                        tobs.months=sum(tobs.months),
                                        events=sum(event), eventshosp=sum(eventhosp)),
                                    by=list(vaxedfirst, vaxstatus.atexit)]
#tobs.cohort.defpos[, sumtobs := tobs.unvax + tobs.vax1dose + tobs.vax2dose]
tobs.cohort.defpos[, rateper1000.unvax := ifelse(vaxstatus.atexit==0,
                                                 1000 * events / tobs.months.unvax, 0)]
tobs.cohort.defpos[, rateper1000.vax1dose :=  ifelse(vaxstatus.atexit==1,
                                                     1000 * events / tobs.months.vax1dose, 0)]
tobs.cohort.defpos[, rateper1000.vax2dose := ifelse(vaxstatus.atexit==2,
                                                    1000 * events / tobs.months.vax2dose, 0)]
tobs.cohort.defpos[, rateper1000 := rateper1000.unvax + rateper1000.vax1dose + rateper1000.vax2dose]

tobs.cohort.defpos[, hosprateper1000.unvax := ifelse(vaxstatus.atexit==0,
                                                 1000 * eventshosp / tobs.months.unvax, 0)]
tobs.cohort.defpos[, hosprateper1000.vax1dose :=  ifelse(vaxstatus.atexit==1,
                                                     1000 * eventshosp / tobs.months.vax1dose, 0)]
tobs.cohort.defpos[, hosprateper1000.vax2dose := ifelse(vaxstatus.atexit==2,
                                                    1000 * eventshosp / tobs.months.vax2dose, 0)]
tobs.cohort.defpos[, hosprateper1000 := hosprateper1000.unvax + hosprateper1000.vax1dose + hosprateper1000.vax2dose]

setorder(tobs.cohort.defpos, vaxstatus.atexit, vaxedfirst)
tobs.cohort.defpos[, .(vaxedfirst, vaxstatus.atexit, events, rateper1000)]

rate.ci <-function(events, rate, digits=1) {
    ci.lower <- qgamma(0.025, events, events / rate)
    ci.upper <- qgamma(0.975, events, events / rate)
    rate.ci <- ifelse(events > 0,
                      paste0(round(rate, digits), " (95% CI ",
                             round(ci.lower, digits), " to ",
                             round(ci.upper, digits), ")"), 0)
    return(rate.ci)
}

tobs.cohort.defpos[, rate.ci := rate.ci(events, rateper1000)]
tobs.cohort.defpos[, hosprate.ci := rate.ci(eventshosp, hosprateper1000, digits=2)]
tobs.cohort.defpos[, .(vaxedfirst, vaxstatus.atexit, events, rate.ci, eventshosp, hosprate.ci)]

with(cohort.defpos, table(vaxedfirst, exclude=NULL))
with(cohort.defpos, table(outcome, vaxedfirst, exclude=NULL))

#############################

min.tstart <- min(cohort.defpos$tstart)
defpos.split <- as.data.table(survSplit(Surv(tstart, tstop, event) ~ .,
                                          data=cohort.defpos,
                                          cut=min.tstart + interval.days * (1:100),
                                          episode=paste0("interval", interval.days, "day")))
defpos.split[event==0, eventhosp := 0]
defpos.split[, datestart := as.Date(tstart, origin=as.Date("1970-01-01"))]
defpos.split[, datestop := as.Date(tstop, origin=as.Date("1970-01-01"))]
defpos.split[, tobs := tstop - tstart]
defpos.split[, months.sincefirstinfection := (tstart - as.integer(firstinfection_date)) * 12 / 365.25]
summary(defpos.split[, .(anon_id, entrydate, exitdate, months.sincefirstinfection, tobs, event)])

# time-updated vaccination status
defpos.split[, vax14.1dose := as.integer(!is.na(vaxdate_1) &
                                           as.integer(datestart - vaxdate_1) >= 14)]
defpos.split[, vax14.2dose := as.integer(!is.na(vaxdate_2) &
                                           as.integer(datestart - vaxdate_2) >= 14)]
defpos.split[, vax14.dose := vax14.1dose + vax14.2dose]
defpos.split[, vax14.factor := car::recode(vax14.dose,
                                             "0='Unvaccinated'; 1='1 dose'; 2='2 doses'",
                                             as.factor=TRUE,
                                             levels=c("Unvaccinated", "1 dose", "2 doses"))]
defpos.split[, weeks.vaxdose2 := as.integer(datestart - vaxdate_2 - 14) * vax14.2dose / 7]
defpos.split[is.na(weeks.vaxdose2), weeks.vaxdose2 := 0]
summary(defpos.split[vax14.factor=="2 doses", weeks.vaxdose2])

coeffs.vaxedbefore <-
    summary(coxph(formula=Surv(tstart, tstop, event) ~ age_years + sexf + occup + listedgr3 +
          vaxedfirst,
      data=defpos.split[vax14.factor=="2 doses"]))$coefficients
coeffs.vaxedbefore <- data.table(effect=rownames(coeffs.vaxedbefore), coeffs.vaxedbefore)
coeffs.vaxedbefore[, rateratio := or.ci(coef, `se(coef)`)]
coeffs.vaxedbefore[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
coeffs.vaxedbefore

################################################################################
## split cohort 1 by vax status when first infected
cohort3 <- cohort.defpos[vaxedfirst=="2 doses"]
## restrict cohort1 to those who had not received first dose on or before specimen_date
cohort1 <- cohort.defpos[vaxedfirst=="Not vaccinated"]

table(cohort1$outcome, exclude=NULL)
table(cohort3$outcome, exclude=NULL)

########################################################################

## select matched controls of the unvaxed cases who ended up in cohort.defpos
cohort2 <- cc.all[CASE==0 & stratum %in%
                  cohort.defpos[vaxedfirst=="Not vaccinated", stratum],
                  .(anon_id, specimen_date, 
                    age_years, agegr3, agegr20, sex, care.home, hh.over18gr3,
                    hh.schoolage, occup, sector, listedgr3, shield.any,
                    vaxdate_1, vaxdate_2)]
## don't exclude those in cohort 1: on first testing positive they transfer after 90 days from cohort 2 to cohort 1  
setnames(cohort2, "specimen_date", "controlsampling_date")
cohort2[, sexf := relevel(sex, ref="Male")]
## select unique anon_ids
setorder(cohort2, "controlsampling_date")
cohort2 <- unique(cohort2, by="anon_id")
setkey(cohort2, anon_id)
## now merge again on anon_id with cc.all[CASE==1] to get case, specimen_date, casegroup etc. 
setkey(cc.all, anon_id)
cohort2 <- cc.all[CASE==1, .(anon_id, specimen_date, specialty, fatalcase, casegroup, adm.within14, date_of_death, underlying_cause_of_death)][cohort2]
table(cohort2$casegroup, exclude=NULL)

## restrict cohort2 to those not vaxed by control sampling date
cohort2 <- cohort2[is.na(vaxdate_1) | vaxdate_1 > controlsampling_date]

## restrict cohort2 to those who had not tested positive by 90 days after controlsamplingdate
cohort2 <- cohort2[is.na(specimen_date) |
                   as.integer(specimen_date - controlsampling_date) >= 90]



cohort2[, hosporfatal := as.integer(adm.within14==1 | fatalcase==1)]
cohort2[is.na(hosporfatal), hosporfatal := 0]
with(cohort2, table(is.na(underlying_cause_of_death), fatalcase, exclude=NULL))
cohort2[, covid_ucod := as.integer(fatalcase==1 &
                                   grepl("^U07[12]$", underlying_cause_of_death))]


cohort2[, entrydate := controlsampling_date + 90]
cohort2 <- cohort2[is.na(date_of_death) | date_of_death > entrydate] 
cohort2[, event := as.integer(!is.na(casegroup))]

cohort2[, eventhosp := as.integer(casegroup=="A"| casegroup=="B")]
cohort2[is.na(casegroup), eventhosp := 0]

cohort2[, exitdate := lastdate]
cohort2[!is.na(date_of_death), exitdate := pmin(exitdate, date_of_death)]
cohort2[event==1, exitdate := pmin(exitdate, specimen_date)]
cohort2 <- cohort2[entrydate < exitdate]
cohort2[, tstart := as.integer(entrydate)] ## this has to be an integer for survsplit 
cohort2[, tstop := as.integer(exitdate)]
cohort2[, tobs := tstop - tstart] 
setkey(cohort2, anon_id)

cohort2[, vaxstatus.atexit := 0]
cohort2[!is.na(vaxdate_1) & vaxdate_1 < exitdate & 
              (is.na(vaxdate_2) | vaxdate_2 > exitdate), vaxstatus.atexit := 1]
cohort2[!is.na(vaxdate_2) & vaxdate_2 < exitdate, vaxstatus.atexit := 2]

cohort1.ucod <- with(cohort1[vaxstatus.atexit==0 & fatal28.reinf==1], table(covid_ucod))
cohort2.ucod <- with(cohort2[vaxstatus.atexit==0 & fatalcase==1], table(covid_ucod))

## cohort2: calculate observed time and testing rate for each individuals by vax status
## tobs.unvax is assigned a value < tobs if vaxdate_1 is nonmissing and before exit date
##   if vaxdate_1 > entry date, assign as vaxdate_1 - entry date, otherwise as 0
## if vaxdate_1 is nonmissing and before exit date, then assign tobs.vax1dose
##   starts at pmax(vaxdate_1, entry date), ends at pmin(vaxdate_2, exitdate)
## if vaxdate_2 is nonmissing and before exit date, then assign tobs.vax2dose
##   starts at pmax(vaxdate_2, entrydate), ends at exitdate

cohort2[!is.na(vaxdate_1) & vaxdate_1 <= entrydate, tobs.unvax := 0]
cohort2[!is.na(vaxdate_1) & vaxdate_1 > entrydate, tobs.unvax := vaxdate_1 - entrydate]
cohort2[is.na(vaxdate_1) | vaxdate_1 > exitdate, tobs.unvax := tobs]

cohort2[, tobs.vax1dose := 0]
cohort2[!is.na(vaxdate_1) & vaxdate_1 <= exitdate &
              (is.na(vaxdate_2) | vaxdate_2 > entrydate),
              tobs.vax1dose := pmin(vaxdate_2, exitdate, na.rm=TRUE) -
                  pmax(vaxdate_1, entrydate)]

cohort2[, tobs.vax2dose := 0]
cohort2[!is.na(vaxdate_2) & vaxdate_2 <= exitdate,
              tobs.vax2dose := exitdate - pmax(vaxdate_2, entrydate)]

## calculate observation time in months by vax status
cohort2[, tobs.months := tobs * 12 / 365.25]
cohort2[, tobs.months.unvax := tobs.unvax * 12 / 365.25]
cohort2[, tobs.months.vax1dose := tobs.vax1dose * 12 / 365.25]
cohort2[, tobs.months.vax2dose := tobs.vax2dose * 12 / 365.25]

tobs.cohort2 <- cohort2[, .(tobs.months.unvax=sum(tobs.months.unvax),
                                        tobs.months.vax1dose=sum(tobs.months.vax1dose), 
                                        tobs.months.vax2dose=sum(tobs.months.vax2dose),
                                        tobs.months=sum(tobs.months),
                                        events=sum(event), eventshosp=sum(eventhosp)),
                                    by=vaxstatus.atexit]
#tobs.cohort2[, sumtobs := tobs.unvax + tobs.vax1dose + tobs.vax2dose]
tobs.cohort2[, rateper1000.unvax := ifelse(vaxstatus.atexit==0,
                                                 1000 * events / tobs.months.unvax, 0)]
tobs.cohort2[, rateper1000.vax1dose :=  ifelse(vaxstatus.atexit==1,
                                                     1000 * events / tobs.months.vax1dose, 0)]
tobs.cohort2[, rateper1000.vax2dose := ifelse(vaxstatus.atexit==2,
                                                    1000 * events / tobs.months.vax2dose, 0)]
tobs.cohort2[, rateper1000 := rateper1000.unvax + rateper1000.vax1dose + rateper1000.vax2dose]

tobs.cohort2[, hosprateper1000.unvax := ifelse(vaxstatus.atexit==0,
                                                 1000 * eventshosp / tobs.months.unvax, 0)]
tobs.cohort2[, hosprateper1000.vax1dose :=  ifelse(vaxstatus.atexit==1,
                                                     1000 * eventshosp / tobs.months.vax1dose, 0)]
tobs.cohort2[, hosprateper1000.vax2dose := ifelse(vaxstatus.atexit==2,
                                                    1000 * eventshosp / tobs.months.vax2dose, 0)]
tobs.cohort2[, hosprateper1000 := hosprateper1000.unvax + hosprateper1000.vax1dose + hosprateper1000.vax2dose]

setorder(tobs.cohort2, vaxstatus.atexit)
tobs.cohort2[, .(vaxstatus.atexit, events, rateper1000)]
tobs.cohort2[, rate.ci := rate.ci(events, rateper1000)]
tobs.cohort2[, hosprate.ci := rate.ci(eventshosp, hosprateper1000, digits=2)]
tobs.cohort2[, .(vaxstatus.atexit, events, rate.ci, eventshosp, hosprate.ci)]

tobs.cohort2[, .(vaxstatus.atexit, events, eventshosp, rate.ci, hosprate.ci)]

paste.rowpercent(with(cohort1, table(agegr3, event)), digits=1)
paste.rowpercent(with(cohort2, table(agegr3, event)), digits=1)
paste.rowpercent(with(cohort1, table(agegr3, eventhosp)), digits=2)
paste.rowpercent(with(cohort2, table(agegr3, eventhosp)), digits=2)

with(cohort1[vaxstatus.atexit==0 & !is.na(specialty)], table(grepl(noncovid.specialtycodes, specialty)))
with(cohort2[vaxstatus.atexit==0 & !is.na(specialty)], table(grepl(noncovid.specialtycodes, specialty)))
with(cohort2[vaxstatus.atexit==0 & !is.na(specialty)], table(grepl(covid.specialtycodes, specialty)))

###################################################################

## merge cohort1 entry, exit, vax dates into ecoss 
setkey(ecoss, anon_id)
setkey(cohort1, anon_id)
## left join ecoss with cohort1 at risk of reinfection, then drop all records not matched in cohort1
cohort1.astests <- cohort1[, .(anon_id, firstinfection_date, entrydate, exitdate, 
                               vaxdate_1, vaxdate_2,
                     age_years, sex, care.home,
                     sector, occup,  listedgr3)][ecoss]
cohort1.astests <- cohort1.astests[!is.na(entrydate)]
## drop tests before entrydate
cohort1.astests <- cohort1.astests[SpecimenDate >= entrydate]
## censor cohort1.astests at first positive test
setorder(cohort1.astests, SpecimenDate)
cohort1.astests[, cumpos := cumsum(as.integer(ecoss.result=="Positive")), by=anon_id]
cohort1.astests <- cohort1.astests[cumpos <= 1]
## this keeps negative results after first positive
## so we need to drop them too
cohort1.astests <- cohort1.astests[cumpos==0 | ecoss.result=="Positive"]

setnames(cohort1.astests, "exitdate", "exitdate_followup")
## recode entrydate and exitdate for each time interval up to a test 
cohort1.astests[, exitdate := SpecimenDate]
cohort1.astests[, diff.exitdate := c(-1, diff(exitdate)), by=anon_id]
cohort1.astests[diff.exitdate > -1, entrydate := exitdate - diff.exitdate]
cohort1.astests <- cohort1.astests[exitdate > entrydate]
cohort1.astests[, testnum := seq_len(.N), by=anon_id]

## code months since first infection at time of each test
cohort1.astests[, months.sincefirstinfection := as.integer(SpecimenDate -
                                                           firstinfection_date) * 12 / 365.25]


## code vax status at time of each test
cohort1.astests[, vax14.1dose := as.integer(!is.na(vaxdate_1) &
                                  (as.integer(SpecimenDate - vaxdate_1) >= 14))]
cohort1.astests[is.na(vax14.1dose), vax14.1dose := 0]
cohort1.astests[, vax14.2dose := as.integer(!is.na(vaxdate_2) &
                                  (as.integer(SpecimenDate - vaxdate_2) >= 14))]
cohort1.astests[is.na(vax14.2dose), vax14.2dose := 0]
cohort1.astests[, vax14.dose := vax14.1dose + vax14.2dose]
cohort1.astests[, vax14.factor := car::recode(vax14.dose,
                                             "0='Unvaccinated'; 1='1 dose'; 2='2 doses'",
                                             as.factor=TRUE,
                                             levels=c("Unvaccinated", "1 dose", "2 doses"))]
setorder(cohort1.astests, anon_id, SpecimenDate)
cohort1.astests[, sexf := relevel(sex, ref="Male")]

#############################################################

## left join ecoss with cohort2 at risk of reinfection, then drop all records not matched in cohort2
cohort2.astests <- cohort2[, .(anon_id, controlsampling_date, entrydate, exitdate, vaxdate_1, vaxdate_2,
                     age_years, sex, care.home,
                     sector, occup,  listedgr3)][ecoss]
cohort2.astests <- cohort2.astests[!is.na(entrydate)]
## drop tests before entrydate
cohort2.astests <- cohort2.astests[SpecimenDate >= entrydate]
## censor cohort2.astests at first positive test
setorder(cohort2.astests, SpecimenDate)
cohort2.astests[, cumpos := cumsum(as.integer(ecoss.result=="Positive")), by=anon_id]
cohort2.astests <- cohort2.astests[cumpos <= 1]
## this keeps negative results after first positive
## so we need to drop them too
cohort2.astests <- cohort2.astests[cumpos==0 | ecoss.result=="Positive"]

setnames(cohort2.astests, "exitdate", "exitdate_followup")
## recode entrydate and exitdate for each time interval up to a test 
cohort2.astests[, exitdate := SpecimenDate]
cohort2.astests[, diff.exitdate := c(-1, diff(exitdate)), by=anon_id]
cohort2.astests[diff.exitdate > -1, entrydate := exitdate - diff.exitdate]
cohort2.astests <- cohort2.astests[exitdate > entrydate]
cohort2.astests[, testnum := seq_len(.N), by=anon_id]

## code months since first infection at time of each test
cohort2.astests[, months.sincefirstsamplingdate := as.integer(SpecimenDate - controlsampling_date) *
            12 / 365.25]

## code vax status at time of each test
cohort2.astests[, vax14.1dose := as.integer(!is.na(vaxdate_1) &
                                  (as.integer(SpecimenDate - vaxdate_1) >= 14))]
cohort2.astests[is.na(vax14.1dose), vax14.1dose := 0]
cohort2.astests[, vax14.2dose := as.integer(!is.na(vaxdate_2) &
                                  (as.integer(SpecimenDate - vaxdate_2) >= 14))]
cohort2.astests[is.na(vax14.2dose), vax14.2dose := 0]
cohort2.astests[, vax14.dose := vax14.1dose + vax14.2dose]
cohort2.astests[, vax14.factor := car::recode(vax14.dose,
                                             "0='Unvaccinated'; 1='1 dose'; 2='2 doses'",
                                             as.factor=TRUE,
                                             levels=c("Unvaccinated", "1 dose", "2 doses"))]
setorder(cohort2.astests, anon_id, SpecimenDate)
cohort2.astests[, sexf := relevel(sex, ref="Male")]

#################################################################

## generate table of number of tests for each individual by vax status
## we've already dropped all cohort1.astests records before entry date
numtests1 <- cohort1.astests[, .N, by=c("anon_id", "vax14.2dose")]
numtests1[, vax14.2dose := car::recode(vax14.2dose, "0='Ntests.unvax'; 1='Ntests.vax2dose'",
                                      as.factor=TRUE)]
numtests1.wide <- dcast(numtests1, anon_id ~ vax14.2dose, value.var="N")
setnafill(numtests1.wide, cols=c("Ntests.unvax", "Ntests.vax2dose"), fill=0)
setkey(numtests1.wide, anon_id)

numtests2 <- cohort2.astests[, .N, by=c("anon_id", "vax14.2dose")]
numtests2[, vax14.2dose := car::recode(vax14.2dose, "0='Ntests.unvax'; 1='Ntests.vax2dose'",
                                      as.factor=TRUE)]
numtests2.wide <- dcast(numtests2, anon_id ~ vax14.2dose, value.var="N")
setnafill(numtests2.wide, cols=c("Ntests.unvax", "Ntests.vax2dose"), fill=0)
setkey(numtests2.wide, anon_id)

######################################################################################
## left join cohort1 with numtests.wide
cohort1 <- numtests1.wide[cohort1] 
setnafill(cohort1, cols=c("Ntests.unvax", "Ntests.vax2dose"), fill=0)

cohort2 <- numtests2.wide[cohort2] 
setnafill(cohort2, cols=c("Ntests.unvax", "Ntests.vax2dose"), fill=0)

cohort1[, Ntests := Ntests.unvax + Ntests.vax2dose]
#cohort1[, testrate := Ntests / (tobs * 12 / 365.25)]
#cohort1[, testrate.unvax := Ntests.unvax / (tobs.unvax * 12 / 365.25)]

testrate1.unvax.permonth <- cohort1[, sum(Ntests.unvax) / sum(tobs.unvax)] * 365.25 / 12

############################################
## cohort2: calculate observed time and testing rate for each individuals by vax status
cohort2[, tobs.unvax := as.integer(pmin(vaxdate_2, exitdate) - entrydate)]
## set tobs.unvax to tobs for those with no second dose
cohort2[is.na(tobs.unvax), tobs.unvax := tobs]
cohort2[tobs.unvax < 0, tobs.unvax := 0]
cohort2[, tobs.vax2dose := tobs - tobs.unvax]
cohort2[is.na(vaxdate_2), tobs.vax := 0]

## calculate observation time in months
cohort2[, tobs.months.unvax := tobs.unvax * 12 / 365.25]
cohort2[, tobs.months.vax2dose := tobs.vax2dose * 12 / 365.25]
cohort2[, Ntests := Ntests.unvax + Ntests.vax2dose]
cohort2[, testrate := Ntests / (tobs * 12 / 365.25)]
testrate2.unvax.permonth <- cohort2[, sum(Ntests.unvax) / sum(tobs.unvax)] * 365.25 / 12

## tabulate testing rates by risk factors 
tests.care <- cohort1[, list(Ntests.unvax=sum(Ntests.unvax),
                                           tobs.months.unvax=sum(tobs.months.unvax),
                                           Ntests.vax2dose=sum(Ntests.vax2dose),
                                           tobs.months.vax2dose=sum(tobs.months.vax2dose)),
                                    by=care.home]
tests.occup <- cohort1[, list(Ntests.unvax=sum(Ntests.unvax),
                                            tobs.months.unvax=sum(tobs.months.unvax),
                                            Ntests.vax2dose=sum(Ntests.vax2dose),
                                            tobs.months.vax2dose=sum(tobs.months.vax2dose)),
                                     by=occup]
tests.listedgr3 <- cohort1[, list(Ntests.unvax=sum(Ntests.unvax),
                                            tobs.months.unvax=sum(tobs.months.unvax),
                                            Ntests.vax2dose=sum(Ntests.vax2dose),
                                            tobs.months.vax2dose=sum(tobs.months.vax2dose)),
                                     by=listedgr3]
#tests.listedgr3 <- tests.listedgr3[c(2, 1, 3), ]

tests.agegr <- cohort1[, list(Ntests.unvax=sum(Ntests.unvax),
                                            tobs.months.unvax=sum(tobs.months.unvax),
                                            Ntests.vax2dose=sum(Ntests.vax2dose),
                                            tobs.months.vax2dose=sum(tobs.months.vax2dose)),
                                     by=agegr3]

#tests.agegr <- tests.agegr[c(2, 3, 1), ]

tests.sex <- cohort1[, list(Ntests.unvax=sum(Ntests.unvax),
                                          tobs.months.unvax=sum(tobs.months.unvax),
                                          Ntests.vax2dose=sum(Ntests.vax2dose),
                                          tobs.months.vax2dose=sum(tobs.months.vax2dose)),
                     by=sex]

tests.all <- cohort1[, list(Ntests.unvax=sum(Ntests.unvax),
                                          tobs.months.unvax=sum(tobs.months.unvax),
                                          Ntests.vax2dose=sum(Ntests.vax2dose),
                                          tobs.months.vax2dose=sum(tobs.months.vax2dose))]
tests.all <- cbind(Variable="All", tests.all)



tests <- rbind(tests.all, tests.sex, tests.agegr, tests.care, tests.occup, tests.listedgr3, use.names=FALSE)

tests[, unvax.freq := round(Ntests.unvax / tobs.months.unvax, 2)]
tests[, vax.freq := round(Ntests.vax2dose / tobs.months.vax2dose, 2)]
tests[, ratio.freq := round(unvax.freq / vax.freq, 2)]

############################################################

interval.days <- 14

tests.bymonth <- cohort1.astests[, .N, by=floor(months.sincefirstinfection)]
colnames(tests.bymonth) <- c("Month", "numtests")
setorder(tests.bymonth, Month)

cohort1[, tobs.sincefirstinfection := as.integer(exitdate - entrydate)]
cohort1.personmonths <- as.data.table(survSplit(Surv(tobs.sincefirstinfection, event) ~ .,
                      data=cohort1[, .(tobs.sincefirstinfection, event)], 
                      cut=seq(1, max(cohort1$tobs.sincefirstinfection), by=365.25/12),
                      episode="Month"))
cohort1.personmonths <- cohort1.personmonths[, .N, by=Month]
cohort1.personmonths[, Month := Month + 2]

setkey(cohort1.personmonths, Month)
setkey(tests.bymonth, Month)
cohort1.personmonths <- tests.bymonth[cohort1.personmonths]
cohort1.personmonths[, rate := numtests/N]
cohort1.personmonths


## split cohort1 into person-time intervals by calendar time, censored at firstpos
min.tstart <- min(cohort1$tstart)
firstpos.split <- as.data.table(survSplit(Surv(tstart, tstop, event) ~ .,
                                          data=cohort1,
                                          cut=min.tstart + interval.days * (1:100),
                                          episode=paste0("interval", interval.days, "day")))
firstpos.split[event==0, eventhosp := 0]
firstpos.split[, datestart := as.Date(tstart, origin=as.Date("1970-01-01"))]
firstpos.split[, datestop := as.Date(tstop, origin=as.Date("1970-01-01"))]
firstpos.split[, tobs := tstop - tstart]
firstpos.split[, months.sincefirstinfection := (tstart - as.integer(firstinfection_date)) * 12 / 365.25]
summary(firstpos.split[, .(anon_id, entrydate, exitdate, months.sincefirstinfection, tobs, event)])

# time-updated vaccination status
firstpos.split[, vax14.1dose := as.integer(!is.na(vaxdate_1) &
                                           as.integer(datestart - vaxdate_1) >= 14)]
firstpos.split[, vax14.2dose := as.integer(!is.na(vaxdate_2) &
                                           as.integer(datestart - vaxdate_2) >= 14)]
firstpos.split[, vax14.dose := vax14.1dose + vax14.2dose]
firstpos.split[, vax14.factor := car::recode(vax14.dose,
                                             "0='Unvaccinated'; 1='1 dose'; 2='2 doses'",
                                             as.factor=TRUE,
                                             levels=c("Unvaccinated", "1 dose", "2 doses"))]
firstpos.split[, weeks.vaxdose2 := as.integer(datestart - vaxdate_2 - 14) * vax14.2dose / 7]
firstpos.split[is.na(weeks.vaxdose2), weeks.vaxdose2 := 0]
summary(firstpos.split[vax14.factor=="2 doses", weeks.vaxdose2])

#######################################################################

## cohort1.astests is already split into person-test observations with time-updated vaccination status
## 
cohort1.astests[, calendar14day := floor(as.integer(SpecimenDate) / 14)]

##################################################

cox.varnames <- c("age_years", "sexf", "care.home", "occup",
                  "listedgr3", "months.sincefirstinfection", "vax14.factor")

## Cox model for test-positive events with tests as timescale
formulastring.tests <- paste0("Surv(time=testnum, event=cumpos) ~ ",
                        paste(c(cox.varnames, "strata(calendar14day)"), collapse=" + "))
cox.formula.tests <- with(cohort1.astests, as.formula(formulastring.tests))
reinf.tests.coxmodel <- coxph(formula=cox.formula.tests,
                                    data=cohort1.astests)
reinf.tests.coeffs <- summary(reinf.tests.coxmodel)$coefficients
reinf.tests.coeffs <- data.table(effect=rownames(reinf.tests.coeffs), reinf.tests.coeffs)
reinf.tests.plot.coeffs <- reinf.tests.coeffs
reinf.tests.coeffs[, rateratio := or.ci(coef, `se(coef)`)]
reinf.tests.coeffs[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
reinf.tests.coeffs <- reinf.tests.coeffs[, .(effect, rateratio, pvalue)]
reinf.tests.coeffs <- pad.coeffs(reinf.tests.coeffs, varnames=cox.varnames, dt=firstpos.split)

######################################################

formulastring <- paste0("Surv(tstart, tstop, event) ~ ", paste(cox.varnames, collapse=" + "))
cox.formula <- with(firstpos.split, as.formula(formulastring))
reinf.time.coeffs <- summary(coxph(cox.formula,
                              data=firstpos.split))$coefficients
reinf.time.coeffs <- data.table(effect=rownames(reinf.time.coeffs), reinf.time.coeffs)
reinf.time.plot.coeffs <- reinf.time.coeffs

reinf.time.coeffs[, rateratio := or.ci(coef, `se(coef)`)]
reinf.time.coeffs[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
reinf.time.coeffs <- reinf.time.coeffs[, .(effect, rateratio, pvalue)]
reinf.time.coeffs <- pad.coeffs(reinf.time.coeffs, varnames=cox.varnames, dt=firstpos.split)
reinf.time.coeffs[, effect:= replace.names(effect)]

freqs.reinf <- univariate.tabulate(varnames=cox.varnames, outcome="event",
                                        data=firstpos.split, drop.reflevel=FALSE)
colnames(freqs.reinf)[1] <- gsub("^0", "Non-case intervals", colnames(freqs.reinf)[1])
colnames(freqs.reinf)[2] <- gsub("^1", "Cases", colnames(freqs.reinf)[2])
table.reinf <- data.table(reinf.time.coeffs[, 1], freqs.reinf,
                          reinf.time.coeffs[, -1], reinf.tests.coeffs[, -1])

####################################################################
reinf.coeffs <- rbind(reinf.time.plot.coeffs[, timescale := "calendar"],
                      reinf.tests.plot.coeffs[, timescale := "tests"])
reinf.coeffs[, yoffset := ifelse(timescale=="calendar", 0.1, -0.1)]
reinf.coeffs[, xmin := coef - 1.96 * `se(coef)`]
reinf.coeffs[, xmax := coef + 1.96 * `se(coef)`]
reinf.coeffs[, effect  := gsub("^occup", "Occupation: ", effect)]
reinf.coeffs[, effect  := gsub("^listedgr3", "", effect)]
reinf.coeffs[, effect  := gsub("^care\\.home", "", effect)]
reinf.coeffs[, effect  := gsub("^sexf", "", effect)]
reinf.coeffs[, effect  := gsub("^vax14.factor", "", effect)]
reinf.coeffs[, effect  := gsub("dose", "vaccine dose", effect)]
reinf.coeffs[, effect := replace.names(effect)]
reinf.coeffs[, effect := factor(effect, levels=effect[1:11])]
reinf.coeffs[, precision := `se(coef)`^-1]
breakpoints <- c(0.2, 0.3, 0.5, 1, 1.5, 2, 3, 5)
p.comparemodels <- ggplot(data=reinf.coeffs, aes(y=effect,
                                              x=coef,
                                              xmin=xmin,
                                              xmax= xmax, 
                                              color=timescale))  + 
    geom_point(shape=23, position = position_nudge(y = reinf.coeffs$yoffset)) + 
    geom_errorbarh(height=0.3,
                   position = position_nudge(y = reinf.coeffs$yoffset)) +
    geom_vline(xintercept=0, color="black", linetype="dashed", alpha=.5) +
    theme(legend.position = 'bottom') + 
    labs(x="Rate ratio for infection (log scale)", y="") +
    guides(color=guide_legend(title="Timescale")) + 
    scale_y_discrete(limits = rev(levels(reinf.coeffs$effect))) + 
    scale_x_continuous(breaks=log(breakpoints),
                       labels=breakpoints)
p.comparemodels

##############################################################
## Cox model for hospitalised / fatal reinfection
## we can re-use the firstpos.split, but with event indicator as hosporfatal

formulastring.hosp <- paste0("Surv(tstart, tstop, eventhosp) ~ ", paste(cox.varnames, collapse=" + "))
cox.formula.hosp <- with(firstpos.split, as.formula(formulastring.hosp))
reinf.coeffs.hosp <- summary(coxph(cox.formula.hosp,
                              data=firstpos.split))$coefficients
reinf.coeffs.hosp <- data.table(effect=rownames(reinf.coeffs.hosp), reinf.coeffs.hosp)
reinf.coeffs.hosp[, rateratio := or.ci(coef, `se(coef)`)]
reinf.coeffs.hosp[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
reinf.coeffs.hosp <- reinf.coeffs.hosp[, .(effect, rateratio, pvalue)]
reinf.coeffs.hosp <- pad.coeffs(reinf.coeffs.hosp, varnames=cox.varnames, dt=firstpos.split)
reinf.coeffs.hosp[, effect := replace.names(effect)]
print(reinf.coeffs.hosp)

table.reinf.hosp <- univariate.tabulate(varnames=cox.varnames, outcome="eventhosp",
                                        data=firstpos.split, drop.reflevel=FALSE)
table.reinf.hosp <- cbind(effect=reinf.coeffs.hosp$effect, table.reinf.hosp, reinf.coeffs.hosp[, -1])
colnames(table.reinf.hosp)[2] <- gsub("^0", "Non-case intervals", colnames(table.reinf.hosp)[2])
colnames(table.reinf.hosp)[3] <- gsub("^1", "Cases", colnames(table.reinf.hosp)[3])

###########################################################################

table.outcome <- with(cohort1,
                    table(outcome, agegr3, exclude=NULL))
table.outcome <- cbind(table.outcome, rowSums(table.outcome))
table.outcome <- paste.colpercent(table.outcome, digits=2)
colnames(table.outcome)[4] <- "All"
table.outcome

incidence1 <- cohort1[, list(personmonths=sum(tobs) * 12 / 365.25,
                             events=sum(event)), by=sexf]
incidence1[, rateper1000 := 1000 * events/personmonths]
averagefollowup1 <- sum(incidence1$personmonths) / nrow(cohort1)

incidence2 <- cohort2[, list(personmonths=sum(tobs) * 12 / 365.25,
                             events=sum(event)), by=sexf]
incidence2[, rateper1000 := 1000 * events/personmonths]
averagefollowup2 <- sum(incidence2$personmonths) / nrow(cohort2)

hosporfatal1 <- cohort1[, list(personmonths=sum(tobs) * 12 / 365.25,
                             events=sum(hosporfatal)), by=sexf]
hosporfatal1[, rateper1000 := 1000 * events/personmonths]

hosporfatal2 <- cohort2[, list(personmonths=sum(tobs) * 12 / 365.25,
                             events=sum(hosporfatal)), by=sexf]
hosporfatal2[, rateper1000 := 1000 * events/personmonths]

fatal1 <- cohort1[, list(personmonths=sum(tobs) * 12 / 365.25,
                             events=sum(as.integer(covid_ucod==1 & fatal28.reinf==1))), by=sexf]
fatal1[, ratepermillion := 1e6 * events/personmonths]

fatal2 <- cohort2[, list(personmonths=sum(tobs) * 12 / 365.25,
                             events=sum(covid_ucod)), by=sexf]
fatal2[, ratepermillion := 1e6 * events/personmonths]

cohort1.meanfollowup <- cohort1[, mean(tobs) * 12 /365.25]
cohort2.meanfollowup <- cohort2[, mean(tobs) * 12 /365.25]


################################################### 


## schoenfeld residuals returned as a matrix with one row per event
schoen <- as.data.table(residuals(reinf.tests.coxmodel, type="schoenfeld"))
schoen <- data.table(cohort1.astests[cumpos==1, .(testnum)], schoen)
schoen.means <- data.table(schoen[, lapply(.SD, mean), by=testnum],
                           N=schoen[, .N, by=testnum][["N"]])

schoen.means[, pointsize := 0.01 * sqrt( N / 1e4)]

schoen.plot <- function(varname, varlabel, legend=FALSE) {
    varname <- sym(varname)
    schoen.plot <- ggplot(data = schoen.means,
                          mapping = aes(x = testnum, y = !!varname,
                                        size=pointsize, weight=N)) +
        geom_point(aes(size=pointsize)) +
        geom_smooth() +
        xlab("Test number") +
        ylab(varlabel) + 
        theme_bw() + theme(legend.position = "none") +
     if(legend) {
        schoen.plot <- schoen.plot +
            scale_size_continuous(breaks=0.02 * sqrt(c(20, 50, 100, 1000) / 1e4),
                                  labels=c(20, 50, 100, 1000),
                                  guide=guide_legend(title="Number of events",
                                                     override.aes = list(linetype=0,
                                                                         color = "black"))) + 
            theme(legend.position = c(0.3, 0.2))
     }
    return(schoen.plot)
}

## plot mean schoenfeld residuals of each covariate against timescale
schoen.age <- schoen.plot("age_years", "Age")
schoen.sex <- schoen.plot("sexfFemale", "Female")
schoen.teacher <- schoen.plot("occupTeacher", "Teacher")
schoen.hcwnotpf <- schoen.plot("occupHealth care, not PF / undetermined", "HCW not PF")
schoen.hcwpf <- schoen.plot("occupHealth care PF", "HCW PF")
schoen.cev <- schoen.plot("listedgr3Eligible for shielding", "CEV status")

#######################################################

testscohort <- cohort1.astests[, .(anon_id, entrydate, exitdate, exitdate_followup, SpecimenDate, ecoss.result, cumpos, testnum, age_years, sex, care.home, occup, listedgr3, vaxdate_1, vaxdate_2)]

notestcohort <- cohort1[!(anon_id %in% testscohort$anon_id),
                        .(anon_id, entrydate, exitdate, vaxdate_1, vaxdate_2,
                          age_years, sex, care.home,
                          sector, occup,  listedgr3)]
notestcohort[, exitdate_followup := exitdate]
notestcohort[, cumpos := 0]
notestcohort[, testnum := 0]

testscohort <- rbind(testscohort, notestcohort, fill=TRUE) # 456121 records

setnames(testscohort, "cumpos", "testpos")
indivs <- unique(testscohort[, .(anon_id)]) # 143565 indivs
indivs[, id := .I] 
setkey(indivs, anon_id)
setkey(testscohort, anon_id)
testscohort <- indivs[testscohort]
testscohort[, anon_id := NULL]
testscohort[, exitdate_followup := as.integer(exitdate_followup)]
testscohort[, exitdate := as.integer(exitdate)]
testscohort[, entrydate := as.integer(entrydate)]
testscohort[, SpecimenDate := as.integer(SpecimenDate)]
testscohort[, vaxdate_1 := as.integer(vaxdate_1)]
testscohort[, vaxdate_2 := as.integer(vaxdate_2)]
testscohort[, anypos := max(testpos), by=id]
testscohort[anypos==0, exitdate := exitdate_followup]

table.defpos.2vax <- tobs.cohort.defpos[vaxstatus.atexit==2,
                                        .(vaxedfirst, tobs.months.vax2dose, events,
                                          eventshosp, rate.ci, hosprate.ci)]
table.defpos.2vax[, tobs.months.vax2dose := round(tobs.months.vax2dose, 1)]
table.defpos.2vax[, rate.ci := gsub("95% CI ", "", rate.ci)]
table.defpos.2vax[, hosprate.ci := gsub("95% CI ", "", hosprate.ci)]

colnames(table.defpos.2vax) <- c("Vaccination status at first infection", "Follow-up time (person-months) after second dose of vaccine", "Detected reinfections", "Hospitalised or fatal reinfections", "Detected reinfection rate", "Hospitalisation rate")
             
save(testscohort, file=paste0(datadir, "testscohort.RData")) 
rmarkdown::render("reinfections.Rmd")
