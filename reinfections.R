library(survival)
library(data.table)
library(ggplot2)

source("helperfunctions.R")

lastdate <- as.Date("2021-09-22")

datadir <- "./data/2021-09-22/"
if(!exists("cc.all")) {
    load(paste0(datadir, "cc.all.RData"))
}

load(paste0(datadir,"ecoss.alltests.RData"))
ecoss <- unique(ecoss)
setkey(ecoss, SpecimenDate)

## select positive records in ecoss based on > 1 test or Ct < 30
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

## alternatively restrict to ids.twicepos
ids.definite <- unique(c(ids.twicepos, ids.ctlt30)) 
#ids.definite <- unique(c(ids.twicepos)) 

## now restrict ecoss.firstpos to those in ids.definite

## create table of individuals with definite first positive result, not vaxed and survived 90 days
ecoss.firstpos.definite <- ecoss.firstpos[anon_id %in% ids.definite]
setkey(ecoss.firstpos.definite, anon_id, SpecimenDate)
setkey(cc.all, anon_id, specimen_date)
ecoss.firstpos.definite <-
    cc.all[, .(anon_id, specimen_date, date_of_death, vaxdate_1)][ecoss.firstpos.definite]
ecoss.firstpos.definite <- ecoss.firstpos.definite[is.na(vaxdate_1) |
                                                   vaxdate_1 > specimen_date]
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

#############################################################

setnames(reinfections, "specimen_date", "reinfection_date", skip_absent=TRUE)
reinfections <- reinfections[, .(anon_id, reinfection_date)]

## use rapid to derive time from specimen date to hospitalisation 
rapid.filename <- paste0(datadir, "CC_RAPID_ANON_2021-09-22.rds")
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

## left join all those who tested positive first time with reinfections to set up a cohort1 of all those at risk for reinfection

setkey(cc.all, anon_id)
setkey(reinfections, anon_id)
cohort1 <- 
    reinfections[cc.all[testpositive.case==TRUE,
                        .(anon_id, specimen_date, stratum, date_of_death, underlying_cause_of_death,
                          covid_ucod,
                          age_years, agegr3, agegr20, sex, care.home, hh.over18gr3,
                          hh.schoolage, occup, sector, listedgr3, shield.any,
                          adm.within14,
                          vaxdate_1, vaxdate_2)]]

# 508621 test-positive cases in cc.all
## cohort1 should be restricted to those with definite positive test first time who were unvaxed and survived at least 90 days
cohort1 <- cohort1[anon_id %in% ecoss.firstpos.definite$anon_id]
setnames(cohort1, "specimen_date", "firstinfection_date", skip_absent=TRUE)
cohort1[, case.reinf := as.integer(!is.na(reinfection_date))]
cohort1[, sexf := relevel(sex, ref="Male")]
table(cohort1$case.reinf)
numtestpos.definite <- nrow(cohort1)
# 231981 test-positive cases matched by anon_id in ecoss.firstpos.definite
###########################################################

## select matched controls of the cases who ended up in cohort1.
cohort2 <- cc.all[CASE==0 & stratum %in% cohort1$stratum,
                  .(anon_id, specimen_date, 
                    age_years, agegr3, agegr20, sex, care.home, hh.over18gr3,
                    hh.schoolage, occup, sector, listedgr3, shield.any,
                    vaxdate_1, vaxdate_2)]
## don't exclude those in cohort 1: on first testing positive they transfer after 90 days from cohort 2 to cohort 1 exit 
setnames(cohort2, "specimen_date", "controlsampling_date")
cohort2[, sexf := relevel(sex, ref="Male")]
## select unique anon_ids
setorder(cohort2, "controlsampling_date")
cohort2 <- unique(cohort2, by="anon_id")
setkey(cohort2, anon_id)
## now merge again on anon_id with cc.all[CASE==1] to get case, specimen_date, casegroup etc. 
setkey(cc.all, anon_id)
cohort2 <- cc.all[CASE==1, .(anon_id, specimen_date, fatalcase, casegroup, adm.within14, date_of_death, underlying_cause_of_death)][cohort2]
table(cohort2$casegroup, exclude=NULL)

###################################################################################

## restrict cohort1 to those who had not received first dose on or before specimen_date
cohort1 <- cohort1[is.na(vaxdate_1) | vaxdate_1 > firstinfection_date]
## restrict cohort2 similarly
cohort2 <- cohort2[is.na(vaxdate_1) | vaxdate_1 > controlsampling_date]

## assign variable for admission within 14 days of reinfection
cohort1[, adm.within14.reinf := as.integer(!is.na(daystoadmission) &
                                                         daystoadmission <= 28)]
table(cohort1$adm.within14.reinf, exclude=NULL)

## 28 fatal reinfections
cohort1[is.na(date_of_death) | case.reinf==0, fatal.reinf := 0]
cohort1[!is.na(date_of_death),
                             fatal.reinf := as.integer(as.integer(date_of_death -
                                                                  reinfection_date) <= 28)]
cohort1[covid_ucod==1, fatal.reinf := 1]
cohort1[is.na(reinfection_date), fatal.reinf := 0]
cohort1.ucod <- with(cohort1[fatal.reinf==1], table(covid_ucod))


with(cohort1, table(adm.within14.reinf, fatal.reinf, exclude=NULL))
with(cohort1, table(is.na(underlying_cause_of_death), fatal.reinf, exclude=NULL))
table(cohort1[fatal.reinf==1, underlying_cause_of_death])
# no deaths with underlying cause missing

## 71 hosp or fatal
cohort1[, hosporfatal := as.integer(adm.within14.reinf==1 | fatal.reinf==1)]
cohort2[, hosporfatal := as.integer(adm.within14==1 | fatalcase==1)]
cohort2[is.na(hosporfatal), hosporfatal := 0]
with(cohort2, table(is.na(underlying_cause_of_death), fatalcase, exclude=NULL))
cohort2[, covid_ucod := as.integer(fatalcase==1 &
                                   grepl("^U07[12]$", underlying_cause_of_death))]
cohort2.ucod <- with(cohort2[fatalcase==1], table(covid_ucod))
     
table(cohort1$case.reinf, exclude=NULL)
table(cohort1$adm.within14.reinf, exclude=NULL)
cohort1[, outcome := case.reinf + adm.within14.reinf]
cohort1[fatal.reinf==1, outcome := 3]
cohort1[fatal.reinf==1 & covid_ucod==0, outcome := 4]

cohort1[, outcome := car::recode(outcome,
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

table(cohort1$outcome, exclude=NULL)

table.outcome <- with(cohort1,
                    table(outcome, agegr3, exclude=NULL))
table.outcome <- cbind(table.outcome, rowSums(table.outcome))
table.outcome <- paste.colpercent(table.outcome, digits=2)
colnames(table.outcome)[4] <- "All"
table.outcome

########################################################################

## set up entry date, event, exitdate: enter at 90 days after first infection or control sampling
cohort1[, entrydate := firstinfection_date + 90]
cohort1 <- cohort1[is.na(date_of_death) | date_of_death > entrydate] 

cohort2[, entrydate := controlsampling_date + 90]
cohort2 <- cohort2[is.na(date_of_death) | date_of_death > entrydate] 

## event defined as detected reinfection
cohort1[, event := as.integer(!is.na(reinfection_date))]
cohort2[, event := as.integer(!is.na(casegroup))]

## 120 hospitalised / fatal reinfections
cohort1[, eventhosp := as.integer(!is.na(reinfection_date) & (adm.within14.reinf==1 | fatal.reinf==1))]
cohort2[, eventhosp := as.integer(casegroup=="A"| casegroup=="B")]
cohort2[is.na(casegroup), eventhosp := 0]

## assign exitdate
cohort1[, exitdate := lastdate]
cohort1[!is.na(date_of_death), exitdate := pmin(exitdate, date_of_death)]
cohort1[event==1, exitdate := pmin(exitdate, reinfection_date)]
cohort1 <- cohort1[entrydate < exitdate]
cohort1[, tstart := as.integer(entrydate)] ## this has to be an integer for survsplit 
cohort1[, tstop := as.integer(exitdate)]
cohort1[, tobs := tstop - tstart] 

cohort2[, exitdate := lastdate]
cohort2[!is.na(date_of_death), exitdate := pmin(exitdate, date_of_death)]
cohort2[event==1, exitdate := pmin(exitdate, specimen_date)]
cohort2 <- cohort2[entrydate < exitdate]
cohort2[, tstart := as.integer(entrydate)] ## this has to be an integer for survsplit 
cohort2[, tstop := as.integer(exitdate)]
cohort2[, tobs := tstop - tstart] 
setkey(cohort2, anon_id)


paste.rowpercent(with(cohort1, table(agegr3, event)), digits=1)
paste.rowpercent(with(cohort2, table(agegr3, event)), digits=1)
paste.rowpercent(with(cohort1, table(agegr3, eventhosp)), digits=2)
paste.rowpercent(with(cohort2, table(agegr3, eventhosp)), digits=2)
###################################################################

## merge cohort1 entry, exit, vax dates into ecoss 
setkey(ecoss, anon_id)
setkey(cohort1, anon_id)
## left join ecoss with cohort1 at risk of reinfection, then drop all records not matched in cohort1
cohort1.astests <- cohort1[, .(anon_id, firstinfection_date, entrydate, exitdate, vaxdate_1, vaxdate_2,
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
cohort1.astests[, months.sincefirstinfection := as.integer(SpecimenDate - firstinfection_date) *
            12 / 365.25]

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

## calculate observed time and testing rate for each individuals by vax status
cohort1[, tobs.unvax := as.integer(pmin(vaxdate_2, exitdate) - entrydate)]
## set tobs.unvax to tobs for those with no second dose
cohort1[is.na(tobs.unvax), tobs.unvax := tobs]
cohort1[tobs.unvax < 0, tobs.unvax := 0]
cohort1[, tobs.vax2dose := tobs - tobs.unvax]
cohort1[is.na(vaxdate_2), tobs.vax := 0]

## calculate observation time in months
cohort1[, tobs.months.unvax := tobs.unvax * 12 / 365.25]
cohort1[, tobs.months.vax2dose := tobs.vax2dose * 12 / 365.25]
cohort1[, Ntests := Ntests.unvax + Ntests.vax2dose]
cohort1[, testrate := Ntests / (tobs * 12 / 365.25)]
testrate1.permonth <- cohort1[, sum(Ntests) / sum(tobs)] * 365.25 / 12

## calculate observed time and testing rate for each individuals by vax status
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
testrate2.permonth <- cohort2[, sum(Ntests) / sum(tobs)] * 365.25 / 12

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
tests.listedgr3 <- tests.listedgr3[c(2, 1, 3), ]

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

## split cohort1 into person-time intervals by calendar time, censored at firstpos
## FIXME: eventhosp is being replicated over multiple person-time intervals

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
                             events=sum(as.integer(covid_ucod==1 & fatal.reinf==1))), by=sexf]
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

## FIXME: should include those who have not had a test 
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

save(testscohort, file=paste0(datadir, "testscohort.RData")) 
rmarkdown::render("reinfections.Rmd")
