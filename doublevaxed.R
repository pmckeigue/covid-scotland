library(survival)
library(data.table)
library(ggplot2)
library(survminer)
library(survMisc)

source("helperfunctions.R")

lastdate <- as.Date("2021-09-22")

datadir <- "./data/2021-09-22/"
if(!exists("cc.all")) {
    load(paste0(datadir, "cc.all.RData"))
}

load(paste0(datadir,"ecoss.alltests.RData"))
ecoss <- unique(ecoss)
setkey(ecoss, SpecimenDate)

## extract cohort of individuals sampled as controls and v2 on or after specimen_date
## entry date of this cohort is 14 days after vaxdate_2
v2.cohort <- cc.all[CASE==0 & !is.na(vaxdate_2) & specimen_date <= vaxdate_2 + 14,
                             .(anon_id, specimen_date, vaxdate_2, date_of_death,
                               age_years, agegr3, sex, care.home, occup5, hh.over18gr3, hh.schoolage,
                               hh.schoolagegr, listedgr3)]

setorder(v2.cohort, specimen_date)
## select first date sampled as control
v2.cohort <- unique(v2.cohort, by="anon_id") # 126804 individuals

setnames(v2.cohort, "specimen_date", "firstsampled.date")
setkey(cc.all, anon_id)

## merge with cases
v2.cohort.events <- cc.all[CASE==1 & anon_id %in% v2.cohort$anon_id,
                           .(anon_id, specimen_date, casegroup)]
## 83067 events in v2.cohort
setkey(v2.cohort, anon_id)
setkey(v2.cohort.events, anon_id)
v2.cohort <- v2.cohort.events[v2.cohort]
rm(v2.cohort.events)
v2.cohort[, event := as.integer(!is.na(specimen_date))] #still have 83067 events in 1262804 individuals

v2.cohort[, entrydate := pmax(firstsampled.date, vaxdate_2 + 14)]
v2.cohort[, exitdate := lastdate]
v2.cohort[!is.na(date_of_death), exitdate := pmin(exitdate, date_of_death)]
v2.cohort[!is.na(specimen_date), exitdate := pmin(exitdate, specimen_date)]
v2.cohort <- v2.cohort[entrydate < exitdate]
## this reduces v2.cohort to 7406 events in 1083844 individuals

v2.cohort[, sexf := relevel(sex, ref="Male")]
## FIXME: vaxdate2 on 11 dec 2020 
v2.cohort[, tstart := as.integer(entrydate)]
v2.cohort[, tstop := as.integer(exitdate)]
v2.cohort[, tobs := tstop - tstart]
## 7406 events
## 2803562 person-months
sum(v2.cohort$event) / (sum(v2.cohort$tobs) * 12 /365)

gc()

############################################################

## merge v2.cohort entry, exit, vax dates into ecoss 
setkey(ecoss, anon_id)
ecoss <- ecoss[!is.na(ecoss.result)]
ecoss <- unique(ecoss, by=c("anon_id", "SpecimenDate", "ecoss.result"))

ecoss.case.ids <- unique(ecoss[ecoss.result=="Positive", anon_id]) # 386894
cc.case.ids <- cc.all[CASE==1, anon_id] # 327178

rm(cc.all)
gc()


ecoss.extracase.ids <- ecoss.case.ids[!(ecoss.case.ids %in% cc.case.ids)]
ecoss.extracases <- ecoss[anon_id %in% ecoss.extracase.ids & ecoss.result=="Positive"]
setorder(ecoss.extracases, SpecimenDate)
ecoss.extracases <- unique(ecoss.extracases, by="anon_id") # 63981 test-positive cases not coded as cases in cc.all

## left join ecoss with v2.cohort at risk of reinfection, then drop all records not matched in v2.cohort
ecoss.cohort <- v2.cohort[, .(anon_id, entrydate, event, exitdate, vaxdate_2,
                       age_years, agegr3, sex, care.home, occup5, hh.over18gr3,
                       hh.schoolagegr, listedgr3)][ecoss]
setnames(ecoss.cohort, "event", "v2cohort.event")

ecoss.cohort[, calendar14day := floor(as.integer(SpecimenDate) / 14)]
ecoss.cohort <- ecoss.cohort[!is.na(entrydate)]
## drop tests before entrydate
ecoss.cohort <- ecoss.cohort[SpecimenDate >= entrydate]
## no need to drop tests before vaxdate2 -- already done in v2.cohort by dropping tests before entrydate
with(ecoss.cohort, table(SpecimenDate >= vaxdate_2 + 14))
ecoss.cohort <- unique(ecoss.cohort, by=c("anon_id", "SpecimenDate", "ecoss.result"))

## censor ecoss.cohort at first positive test
setorder(ecoss.cohort, SpecimenDate)
ecoss.cohort[, cumpos := cumsum(as.integer(ecoss.result=="Positive")), by=anon_id]

keep.tofirstnonzero <- function(x) {
    keep <- rep(TRUE, length(x))
    nonzero <- which(x > 0)
    if(length(nonzero) > 0) { 
        first.nonzero <- nonzero[1]
        if(first.nonzero < length(x)) {
            keep[(first.nonzero + 1):length(x)] <- FALSE
        }
    }
    return(keep)
}
    
ecoss.cohort[, keep := keep.tofirstnonzero(cumpos), by=anon_id]
ecoss.cohort <- ecoss.cohort[keep==TRUE]
table(ecoss.cohort[cumpos > 0, .N, by=anon_id][["N"]]) # 23577 with a positive test

setnames(ecoss.cohort, "cumpos", "event")  
with(ecoss.cohort, table(v2cohort.event, event, exclude=NULL))
setorder(ecoss.cohort, anon_id, SpecimenDate)

ecoss.cohort[v2cohort.event==0 & event==1, .(vaxdate_2, v2cohort.event, event, ecoss.result, SpecimenDate, exitdate)]

ecoss.cohort[, testnum := 1:.N, by=anon_id]
ecoss.cohort[, sexf := relevel(sex, ref="Male")]
setorder(ecoss.cohort, anon_id)

## should censor at say 50 or 100 tests
######################################################################################

cox.v2.varnames <- c("age_years", "sexf", "care.home", "occup5",
                     "hh.over18gr3", "hh.schoolagegr",
                     "listedgr3")

formulastring <- paste0("Surv(tstart, tstop, event) ~ ", paste(cox.v2.varnames, collapse=" + "))
cox.formula <- with(v2.cohort, as.formula(formulastring))
v2.time.coeffs <- summary(coxph(cox.formula,
                              data=v2.cohort))$coefficients
v2.time.coeffs[1, 1:2] <- 10 * v2.time.coeffs[1, 1:2] # multiply age coeffs by 10
v2.time.coeffs <- data.table(effect=rownames(v2.time.coeffs), v2.time.coeffs)
v2.time.coeffs[, rateratio := or.ci(coef, `se(coef)`)]
v2.time.coeffs[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
v2.time.coeffs <- v2.time.coeffs[, .(effect, rateratio, pvalue)]
v2.time.coeffs <- pad.coeffs(v2.time.coeffs, varnames=cox.v2.varnames, dt=v2.cohort)
v2.time.coeffs[, effect:= replace.names(effect)]

## Cox model for test-positive events with tests as timescale
formulastring.tests <- paste0("Surv(time=testnum, event=event) ~ ",
                        paste(c(cox.v2.varnames, "strata(calendar14day)"), collapse=" + "))
cox.formula.tests <- with(ecoss.cohort, as.formula(formulastring.tests))
tests.coxmodel <- coxph(formula=cox.formula.tests, data=ecoss.cohort)
v2.tests.coeffs <- summary(tests.coxmodel)$coefficients
v2.tests.coeffs[1, 1:2] <- 10 * v2.tests.coeffs[1, 1:2] # multiply age coeffs by 10
v2.tests.coeffs <- data.table(effect=rownames(v2.tests.coeffs), v2.tests.coeffs)
v2.tests.coeffs[, rateratio := or.ci(coef, `se(coef)`)]
v2.tests.coeffs[, pvalue := format.pvalue(z, `Pr(>|z|)`)]
v2.tests.coeffs <- v2.tests.coeffs[, .(effect, rateratio, pvalue)]
v2.tests.coeffs <- pad.coeffs(v2.tests.coeffs, varnames=cox.v2.varnames, dt=ecoss.cohort)

freqs.v2 <- univariate.tabulate(varnames=cox.v2.varnames, outcome="event",
                                        data=v2.cohort, drop.reflevel=FALSE)
colnames(freqs.v2)[1] <- gsub("^0", "Non-cases", colnames(freqs.v2)[1])
colnames(freqs.v2)[2] <- gsub("^1", "Cases", colnames(freqs.v2)[2])
table.v2 <- data.table(v2.time.coeffs[, 1], freqs.v2,
                          v2.time.coeffs[, -1], v2.tests.coeffs[, -1])
table.v2[1, effect := "Age (10-year scale)"]

###################################################

uncond.formulastring.tests <- paste0("event ~ ", paste(c(cox.v2.varnames,
                                                         "as.factor(calendar14day)"),
                                                   collapse=" + "))
glm.formula.tests <- with(ecoss.cohort, as.formula(uncond.formulastring.tests))
v2.tests.coeffs.uc <- summary(glm(formula=glm.formula.tests,
                                  #data=ecoss.cohort,
                                  data=ecoss.cohort[testnum <= 3],
                                  family="binomial"))$coefficients[2:13, ]
v2.tests.coeffs.uc[1, 1:2] <- 10 * v2.tests.coeffs.uc[1, 1:2] # multiply age coeffs by 10
v2.tests.coeffs.uc <- data.table(effect=rownames(v2.tests.coeffs.uc), v2.tests.coeffs.uc)
v2.tests.coeffs.uc[, oddsratio := or.ci(Estimate, `Std. Error`)]
v2.tests.coeffs.uc[, pvalue := format.pvalue(`z value`, `Pr(>|z|)`)]
v2.tests.coeffs.uc <- v2.tests.coeffs.uc[, .(effect, oddsratio, pvalue)]
v2.tests.coeffs.uc <- pad.coeffs(v2.tests.coeffs.uc, varnames=cox.v2.varnames, dt=ecoss.cohort)
v2.tests.coeffs.uc[1, effect := "Age (10-year scale)"]
####################################

# test if scaled Schoenfeld residuals are independent of time
test.ph <- survival::cox.zph(tests.coxmodel)
test.ph

ecoss.cohort[, resid_mart := residuals(tests.coxmodel, type = "martingale")]

## schoenfeld residuals returned as a matrix with one row per event
schoen <- as.data.table(residuals(tests.coxmodel, type="schoenfeld"))
schoen <- data.table(ecoss.cohort[event==1, .(testnum)], schoen)
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
schoen.age <- schoen.plot("age_years", "Age (years)")
schoen.sex <- schoen.plot("sexfFemale", "Female sex")
schoen.teacher <- schoen.plot("occup5Teacher, secondary/other", "Teacher status")
schoen.hcwnotpf <- schoen.plot("occup5Health care, not PF / undetermined", "Health care not PF")
schoen.hcwpf <- schoen.plot("occup5Health care PF", "Health care PF")
schoen.gt2adults <- schoen.plot("hh.over18gr33 or more", ">2 adults in household")
schoen.schoolage2 <- schoen.plot("hh.schoolagegr2 or more", ">1 school age children")
schoen.cev <- schoen.plot("listedgr3Eligible for shielding", "CEV status")


## plot martingale residuals against age
p.mart <- ggplot(data = ecoss.cohort, mapping = aes(x = age_years, y = resid_mart,
                                                    color=as.factor(event))) +
    geom_point(size=0.5) +
    geom_smooth(color="gray") +
    scale_y_continuous(limit = c(-1, 1)) +
    labs(x= "Age (years)",
         y="Martingale residual") +
    theme_bw() + theme(legend.key = element_blank())


## Cox-Snell residuals
## martingale residual m_i is cumulative observed minus expected events: y_i - lambda_i
## cox-snell "residual" is minus log survival prob i.e the fitted cumulative hazard rate
## survival prob is exp(-lambda_i) so this is just lambda_i 
## cox-snell resid is for each observation
ecoss.cohort[, resid_coxsnell := event - resid_mart]

## Fit model for observed event to Cox-Snell resids 
## resid_coxsnell is 
fit_coxsnell <- coxph(formula = Surv(resid_coxsnell, event) ~ 1,
                      data = ecoss.cohort,
                      ties = c("efron","breslow","exact")[1])

## Nelson-Aalen estimator for baseline hazard (all covariates zero)
## basehaz computes predicted survival curve for a Cox model
df_base_haz <- as.data.table(basehaz(fit_coxsnell, centered = FALSE))
df_base_haz[, fitminusexp := hazard - time]

p.coxsnell <- ggplot(data = df_base_haz, mapping = aes(x = hazard, y=time)) +
    geom_point(size=0.1) +
    labs(x = "Cox-Snell residual",
         y = "Cumulative hazard") +
    theme_bw() + theme(legend.key = element_blank())

p.coxsnell.obsminusexp <- ggplot(data = df_base_haz,
                                 mapping = aes(x = hazard, y=fitminusexp)) +
    geom_point(size=0.1) +
    scale_y_continuous(limit = c(-0.05, 0.05)) +
    labs(x = "Fitted cumulative hazard",
         y = "Fitted cumulative hazard minus standard exponential") +
    theme_bw() + theme(legend.key = element_blank())

if(FALSE) { # these functions hang
    library(survminer)
    ggcoxzph(test.ph)
    ggcoxdiagnostics(fit=tests.coxmodel, type="deviance", linear.predictions=FALSE)
    ggcoxdiagnostics(fit=tests.coxmodel, type = "martingale", linear.predictions = TRUE)
    ## plot martingale residuals
    ggcoxfunctional(formula=cox.formula.tests, data=ecoss.cohort)
}


## generate a table of number of tests and follow-up time for each individual
ecoss.cohort[, numtests := max(testnum), by=anon_id]
ecoss.cohort[, followupmonths := as.integer(exitdate - entrydate) * 12 / 365.25]
ecoss.indivs <- unique(ecoss.cohort, by="anon_id")

## loop over factor covariates to calculate testing rates by levels of each factor
## age group, sex, care home, occup5, listedgr3

tests.agegr3 <- ecoss.indivs[, list(followupmonths=sum(followupmonths),
                                   numtests=sum(numtests)), by=agegr3]
tests.agegr3[, testrate := numtests / followupmonths]
tests.agegr3 <- tests.agegr3[3:1, ]

tests.sexf <- ecoss.indivs[, list(followupmonths=sum(followupmonths),
                                   numtests=sum(numtests)), by=sexf]
tests.sexf[, testrate := numtests / followupmonths]
tests.sexf

tests.care.home <- ecoss.indivs[, list(followupmonths=sum(followupmonths),
                                   numtests=sum(numtests)), by=care.home]
tests.care.home[, testrate := numtests / followupmonths]
tests.care.home

tests.occup5 <- ecoss.indivs[, list(followupmonths=sum(followupmonths),
                                   numtests=sum(numtests)), by=occup5]
tests.occup5[, testrate := numtests / followupmonths]
tests.occup5

tests.listedgr3 <- ecoss.indivs[, list(followupmonths=sum(followupmonths),
                                   numtests=sum(numtests)), by=listedgr3]
tests.listedgr3[, testrate := numtests / followupmonths]
tests.listedgr3 <- tests.listedgr3[c(2, 1, 3), ]

## rbind into a table
table.testrates <- as.data.table(rbind(tests.agegr3, tests.sexf, tests.care.home,
                                       tests.occup5, tests.listedgr3,
                                       use.names=FALSE))

table.testrates[, testrate := round(testrate, 2)] 

rmarkdown::render("doublevaxed.Rmd")
