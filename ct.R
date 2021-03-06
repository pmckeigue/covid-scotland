library(ggplot2)
library(data.table)
library(survival)

source("helperfunctions.R")

tabulate.freqs.coxregressions <- function(varnames=covariates, stratumvar="week",
                                          data, split.data) {
    table.freqs <- univariate.tabulate(varnames=varnames, outcome="event",
                                   data=data,
                                   drop.reflevel=FALSE, drop.sparserows=FALSE,
                                   drop.emptyrows=TRUE, minrowsum=10)

    univariate.table <- univariate.clogit(varnames=varnames, outcome="event",
                                          data=split.data,
                                          add.reflevel=TRUE,
                                          model="cox", stratumvar="week")
    
    multivariate.table <- multivariate.clogit(varnames=varnames, outcome="event",
                                              data=split.data,
                                              add.reflevel=TRUE,
                                              model="cox", stratumvar="week")
    
    table.aug <- combine.tables3(table.freqs, univariate.table, multivariate.table)
    rownames(table.aug) <- replace.names(rownames(table.aug))
    return(table.aug)
}

## rollmean.dtN returns the variable date for the rolling mean 
rollmean.dtN <- function(dt, k) {
    ## add rows for missing dates by a left join
    N.dates <- dt[, .N, by=SPECIMENDATE]
    setkey(N.dates, SPECIMENDATE)
    all.dates <- data.table(date=seq(from=min(N.dates$SPECIMENDATE),
                                     to=max(N.dates$SPECIMENDATE),
                                     by=1))
    setkey(all.dates, date)
    all.dates <- N.dates[all.dates]
    all.dates[is.na(N), N := 0]
    return(data.table(date=dt$SPECIMENDATE,
                      avdaily=zoo::rollmean(dt$N, k=k, fill=NA)))
}

########################################################

datadir <- "./data/2021-02-18/"

load("ct.RData") ## FIXME - should have been saved to datadir

if(!exists("cc.all")) {
    load(paste0(datadir, "cc.all.RData"))
}

firstdate <- as.Date("2020-09-01")
lastdate <- as.Date("2021-02-01")
lastdate.tests <- lastdate.specimen - 7

##########################################################

ct.newnames <- grep("Ct$", names(ct), value=TRUE)

## diff2channels is N_Ct - ORF1ab Ct
## line of separation is -N_Ct = -ORF1abCT - 1 
## i.e. diff2channels = 1

## https://assets.publishing.service.gov.uk/government/uploads/system/uploads/attachment_data/file/949639/Technical_Briefing_VOC202012-2_Briefing_2_FINAL.pdf
## SGTF defined as S gene negative & ORF1ab Ct < 31  & N Ct <31
## try with separation line at diff2channels <= 1
ct[S_result=="POSITIVE", sgtf := 1]
ct[S_result=="NEGATIVE" & ORF1ab_result=="POSITIVE" &  N_result=="POSITIVE" &
   ORF1ab_Ct < 30 & N_Ct < 30 & diff2channels <= 1, sgtf := 3]
ct[S_result=="NEGATIVE" & ORF1ab_result=="POSITIVE" & N_result=="POSITIVE" &
   (ORF1ab_Ct >= 30 | N_Ct >= 30 | diff2channels > 1), sgtf := 2]

with(ct, tapply(diff2channels, Sgene.dropout, quantile, na.rm=TRUE))
table(ct$sgtf, ct$Sgene.dropout)

## create variables for max and min sgtf
ct[, max.sgtf := max(sgtf, na.rm=TRUE), by=ANON_ID]
ct[, min.sgtf := min(sgtf, na.rm=TRUE), by=ANON_ID]
ct[is.infinite(max.sgtf), max.sgtf := NA]
ct[is.infinite(-min.sgtf), min.sgtf := NA]
ct[max.sgtf==3, truepos := "Other test definite dropout"]
ct[min.sgtf==1, truepos := "Other test no dropout"]

## for tests with undetermined result but ground truth available from another test on same individual,
## plot S gene result against other two channel Ct values with colour-coding of true pos and true neg 
p.check <- ggplot(data=ct[sgtf.fix1=="Undetermined" & !is.na(truepos)], 
       aes(x=ORF1ab_Ct, y=N_Ct, color=truepos)) +
    geom_point() +
    scale_color_manual(values=c("red", "blue")) +
    scale_x_reverse(breaks=c(40, 30, 20, 10),
                    labels=c("40", "30", "20", "10")) +
    scale_y_reverse(breaks=c(40, 30, 20, 10),
                    labels=c("40", "30", "20", "10")) +
    # horiz line from 33 to 28 at y=31
    geom_segment(aes(x=33, y=31, xend=28, yend=31), linetype=3, color="black") +
    # vertical line from 31 to 10 at x=31
    geom_segment(aes(x=31, y=31, xend=31, yend=10), linetype=3, color="black") +
    geom_segment(aes(x=29.5, y=31, xend=8.5, yend=10), linetype=3, color="black") +
    # vertical line from 31 to 10 at x=31
    geom_segment(aes(x=33, y=31, xend=33, yend=10), color="green") +
    geom_segment(aes(x=28, y=31, xend=7, yend=10), color="green") +
    ylab("N gene cycle threshold (reverse scale)") +
    xlab("ORF1ab gene cycle threshold (reverse scale)") + 
    ggtitle("Repeat tests initially classified as Undetermined")
p.check

## plot S gene result against other two channel Ct values 
p.dropout <- ggplot(data=ct[ORF1ab_result=="POSITIVE" & N_result=="POSITIVE" &
                            !is.na(S_result)],
                    aes(x=ORF1ab_Ct, y=N_Ct, color=S_result)) +
    geom_point(size=0.5) +
    scale_color_manual(values=c("red", "blue")) +
    scale_x_reverse(breaks=c(40, 30, 20, 10),
                    labels=c("40", "30", "20", "10")) + 
    scale_y_reverse(breaks=c(40, 30, 20, 10),
                    labels=c("40", "30", "20", "10")) +
    ylab("N gene cycle threshold (reverse scale)") +
    xlab("ORF1ab gene cycle threshold (reverse scale)") + 
    theme(legend.position = c(0.2, 0.8)) +
    ggtitle("Relation of S gene detection to ORF and N gene Ct values")

###########################################

## plot median channel values by date
ct.mediansbydate <- ct[,
                       lapply(.SD, median, na.rm=TRUE),
                       keyby=SpecimenDate,
                       .SDcols=c(ct.newnames, "av2channels")]
ct.mediansbydate <- melt(ct.mediansbydate, id.vars="SpecimenDate",
                         measure.vars=ct.newnames)

p.mediansbydate <-
    ggplot(data=ct.mediansbydate[SpecimenDate >= as.Date("2020-10-01")],
           aes(x=SpecimenDate, y=value, group=variable,
               color=variable)) +
    scale_x_date(breaks = seq.Date(from = firstdate,
                                   to = lastdate.tests, by = "month"),
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = firstdate,
                                                  to = lastdate.tests, by = "month"),
                                         "%d %b")
                 ),
                 limits=c(firstdate, lastdate.tests)) +
    scale_color_manual(labels = c("ORF1ab", "N", "S", "MS2 (control)"),
                       values = c("red", "blue", "green", "black")) +
    scale_y_reverse(limits=c(25, 18), expand=c(0, 0)) + 
    labs(y="Median cycle threshold", 
         x="Specimen date",
         tag="(b)") + 
    theme(legend.position = c(0.1, 0.5))+
    theme(legend.title = element_blank()) + 
    geom_line()
p.mediansbydate

## area plot of time series by S gene dropout

setkey(ct, SpecimenDate)
ct.first <- ct[!duplicated(ANON_ID)]  ## restrict to first positive test
ct.Sgene <- ct.first[!is.na(Sgene.dropout), .N, by=c("SpecimenDate", "Sgene.dropout")]

p.variantbydate <- ggplot(data=ct.Sgene, aes(x=SpecimenDate, y=N,
                                              group=Sgene.dropout, fill=Sgene.dropout)) +
    geom_area() +
    scale_x_date(breaks = seq.Date(from = firstdate,
                                   to = lastdate.tests, by = "month"),
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = firstdate,
                                                  to = lastdate.tests, by = "month"),
                                         "%d %b")
                 ),
                 limits=c(firstdate, lastdate.tests)) +
    labs(x="Specimen date", 
         y="Number of test-positive cases",
         tag="(a)") + 
    scale_fill_manual(values=c("grey", "blue", "red")) +
    theme(legend.position = c(0.4, 0.7)) +
    theme(legend.title = element_blank())

p.variantbydate

########################################################
## Figure 1: severe cases by date of presentation
winsize.casedates <- 3
## use by= to get rolling means by category

casedates.bysource <- cc.severe[CASE==1, .N, by=.(SPECIMENDATE, exp.group)]
setkey(casedates.bysource, SPECIMENDATE)
casedates.bysource <- casedates.bysource[, rollmean.dtN(dt=.SD, k=winsize.casedates),
                                         by=exp.group,
                                         .SDcols=c("SPECIMENDATE", "N")]
                                 
p.bysource <- ggplot(data=casedates.bysource[date < as.Date("2021-01-20")],
                     aes(x=date,
                         y=avdaily, fill=exp.group)) +
    geom_area() +
    scale_x_date(breaks = seq.Date(from = firstdate,
                                   to = lastdate, by = "month"),
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = firstdate,
                                                  to = lastdate, by = "month"),
                                         "%d %b")
                 ),
                 limits=c(firstdate, lastdate), expand=c(0, 0)) +
    theme(legend.position = c(0.1, 0.8)) +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values=c("grey", "orange", "green")) +
    xlab(paste0("Specimen date: mid-point of ", winsize.casedates, "-day window")) +
         ylab("Daily severe cases")
p.bysource

############################################################################

cases.all <- cc.kept[!is.na(Sgene.dropout) & Sgene.dropout != "Undetermined" & SPECIMENDATE >= firstdate]
cases.all[, severe.case := as.integer(group == "A")]
cases.all[, hosp.case := as.integer(group == "A" | group == "B")]

## drop the redundant level (doesn't actually recode any values)
cases.all[, Sgene.dropout := car::recode(Sgene.dropout,
                                         "'Undetermined'='No dropout'", 
                                         as.factor=TRUE,
                                         levels=c("No dropout", "Definite dropout"))]

## about one quarter of those with nonmissing Sgene.dropout are not yet classified as cases
## but as long as they are linked we can include them in the analysis
with(cases.all, table(CASE, Sgene.dropout, exclude=NULL))

cases.all[, age.5yr := 5 * floor(AGE/5)]
cases.all[, week.specimen := as.factor(floor(as.integer(SPECIMENDATE) / 7))]
setkey(cases.all, week.specimen, sex, age.5yr)

cases.all[, stratum := .GRP, by=key(cases.all)]
cases.all[, stratum := as.factor(stratum)]
cat(length(unique(cases.all$stratum)), "strata based on week of specimen, sex and 5-year age group\n")  # 1213 strata

cases.all[, `:=`(N.severe = length(which(severe.case==1)),
              N.notsevere = length(which(severe.case==0))), by=stratum]
cases.all[, mutant := as.integer(Sgene.dropout=="Definite dropout")]

## tabulate cases by month
cases.casegr <-
    with(cc.all[CASE==1],
         table(casegr, lubridate::month(SpecimenDate, label=TRUE), exclude=NULL)[, c(9:12, 1)]
         )
colsums.cases.casegr <- colSums(cases.casegr)
casefatality.casegr <- round(100 * colSums(cases.casegr[3:6, ]) / colsums.cases.casegr, 1)
cases.casegr <- paste.colpercent(cases.casegr)
cases.casegr <- rbind(colsums.cases.casegr, cases.casegr, casefatality.casegr)
#cases.casegr[, ncol(cases.casegr)] <- "."

rownames(cases.casegr)[c(1, 8)] <- c("Total cases", "Cohort fatality rate (\\%)")

cases.testCt <-
    with(cc.all[CASE==1],   
         table(flag_lighthouse_labs_testing, lubridate::month(SpecimenDate, label=TRUE))[, c(9:12, 1)]
         )

colsums.cases.testCt <- colSums(cases.testCt)
cases.testCt <- paste.colpercent(cases.testCt)
cases.testCt <- rbind(colsums.cases.testCt, cases.testCt)
rownames(cases.testCt) <- c("Total test-positive cases",
                                 "Ct record, not Lighthouse",  "Lighthouse test",
                                 "NHS test")

cases.testsgtf <-
    with(cc.all[CASE==1 & !is.na(Sgene.dropout)],
         table(Sgene.dropout, lubridate::month(SpecimenDate, label=TRUE))[, c(9:12, 1)]
         )
colsums.cases.testsgtf <- colSums(cases.testsgtf)
cases.testsgtf <- paste.colpercent(cases.testsgtf)
cases.testsgtf <- rbind(colsums.cases.testsgtf, cases.testsgtf)
rownames(cases.testsgtf)[1] <- "Total cases with nonmissing S gene dropout status"

table.cases.bymonth <- rbind(cases.casegr, cases.testCt, cases.testsgtf)

cases.sgtf.fatality <-
        with(cc.all[CASE==1 & Sgene.dropout=="Definite dropout"],
         table(fatalcase, lubridate::month(SpecimenDate, label=TRUE))[, c(9:12, 1)]
         )
cases.sgtf.fatality <- paste.colpercent(cases.sgtf.fatality, digits=2)[2, ]  

cases.nosgtf.fatality <-
        with(cc.all[CASE==1 & Sgene.dropout=="No dropout"],
         table(fatalcase, lubridate::month(SpecimenDate, label=TRUE))[, c(9:12, 1)]
         )
cases.nosgtf.fatality <- paste.colpercent(cases.nosgtf.fatality, digits=2)[2, ]  

cases.sgtf.fatality <- rbind(cases.nosgtf.fatality, cases.sgtf.fatality)
#cases.sgtf.fatality[, ncol(cases.sgtf.fatality)] <- rep(".", 2)
rownames(cases.sgtf.fatality) <- c("No dropout", "Definite dropout")

table.cases.bymonth <- rbind(cases.casegr, cases.testCt, cases.testsgtf, cases.sgtf.fatality)

#####################################################################
## Surv object for deaths
## create event variable
cases.all[, death28 := ifelse(!is.na(daystodeath) & daystodeath <= 28, 1, 0)]
cases.all <- cases.all[!is.na(censoringdays.death) & censoringdays.death > 0]
## replace missing values and values > 28 with pmin (28, censoringdays.death)
cases.all[!is.na(censoringdays.death),
          daystodeath28 := pmin(28, censoringdays.death, daystodeath, na.rm=TRUE)]
death.Surv <- with(cases.all,
                 survival::Surv(time=daystodeath28, event=death28, type="right"))  
death.fit <- survival::survfit(death.Surv ~ Sgene.dropout, data=cases.all)

###############################################################
## Surv object for severe
## failure time is min of daystocritical, daystodeath
cases.all[, daystosevere := pmin(daystodeath, daystocritical, na.rm=TRUE)]
cases.all[, censoringdays.severe := pmin(censoringdays.death,
                                         censoringdays.critical, na.rm=TRUE)]
## create event variable
cases.all[, severe28 := ifelse(!is.na(daystosevere) & daystosevere <= 28, 1, 0)]
## replace missing values and values > 28
cases.all[, daystosevere28 := pmin(28, censoringdays.severe, daystosevere, na.rm=TRUE)]
severe.Surv <- with(cases.all,
                 survival::Surv(time=daystosevere28, event=severe28, type="right"))  
severe.fit <- survival::survfit(severe.Surv ~ Sgene.dropout, data=cases.all)

#####################################################################
## Surv object for hosp
## failure time is min of daystocritical, daystodeath, daystoadmission
cases.all[, daystohosp := pmin(daystodeath, 
                               daystoadmission, na.rm=TRUE)]
## FIXME -- some records shoudl have censoringdays.admission less than 9
cases.all[, censoringdays.hosp := pmin(censoringdays.death,
                                         censoringdays.admission, na.rm=TRUE)]
## create event variable
cases.all[, hosp28 := ifelse(!is.na(daystohosp) & daystohosp <= 28, 1, 0)]
## recode time values > 28 as 28
cases.all[, daystohosp28 := pmin(28, censoringdays.hosp, daystohosp, na.rm=TRUE)]
hosp.Surv <- with(cases.all,
                 survival::Surv(time=daystohosp28, event=hosp28, type="right"))  
hosp.fit <- survival::survfit(hosp.Surv ~ Sgene.dropout, data=cases.all)

#######################################################################
## survival curves
psurv.fatal <-
    survminer::ggsurvplot(death.fit, data = cases.all,
                          title="Death",
                          xlim=c(0, 28), ylim=c(0.99, 1), palette=c("blue", "red"),
                          xlab="", # "Days since testing positive",
                          break.time.by=7,
                          legend=c(0.2, 0.3), 
                          legend.title="S gene status:",
                          legend.labs=c("No dropout", "Definite dropout"),
                           axes.offset=FALSE,
                         risk.table.y.text = FALSE)
psurv.severe <-
    survminer::ggsurvplot(severe.fit, data = cases.all, 
                          title="Entry to critical care or death",
                          xlim=c(0, 28), ylim=c(0.95, 1), palette=c("blue", "red"),
                          xlab="", # "Days since testing positive",
                          break.time.by=7,
                          legend="none",
                          axes.offset=FALSE,
                          risk.table.y.text = FALSE)
psurv.hosp <-
    survminer::ggsurvplot(hosp.fit, data = cases.all, 
                          title="Hospitalisation or death",
                          xlim=c(0, 28), ylim=c(0.9, 1), palette=c("blue", "red"),
                          xlab="Days since testing positive",
                          break.time.by=7,
                          legend="none",
                          axes.offset=FALSE,
                          risk.table.y.text = FALSE)

###########################################################################
## tabulate data for survival curves

table.daystodeath <- cases.all[, .N, by=.(daystodeath28, death28)]
death28.table <- dcast(data=table.daystodeath, formula=daystodeath28 ~ death28, value.var="N")
colnames(death28.table) <- c("days.sincetest", "Censor.fordeath", "Deaths") 
setkey(death28.table, days.sincetest)

table.daystosevere <- cases.all[, .N, by=.(daystosevere28, severe28)]
severe28.table <- dcast(data=table.daystosevere, formula=daystosevere28 ~ severe28, value.var="N")
colnames(severe28.table) <- c("days.sincetest", "Censor.forsevere", "Severe") 
setkey(severe28.table, days.sincetest)

table.daystohosp <- cases.all[, .N, by=.(daystohosp28, hosp28)]
hosp28.table <- dcast(data=table.daystohosp, formula=daystohosp28 ~ hosp28, value.var="N")
colnames(hosp28.table) <- c("days.sincetest", "Censor.forhosp", "Hosp") 
setkey(hosp28.table, days.sincetest)

surv28.table <- severe28.table[hosp28.table]
surv28.table <- death28.table[surv28.table]

for (j in 2:7)
    set(surv28.table, which(is.na(surv28.table[[j]])), j, 0 )
totals <-  colSums(surv28.table)

surv28.table <- rbind(as.matrix(surv28.table), totals)

colnames(surv28.table) <- c("Days from testing positive", "Censored by last date of follow-up for deaths",  "Deaths",
                            "Censored by last date of follow-up for severe cases", "Severe cases",
                            "Censored by last date of follow-up for hospitalisation", "Hospitalised")

surv28.table[nrow(surv28.table), 1] <- NA

####################################################################

covariates=c("AGE", "sex", "region4", "care.home", "COVID.age",
             "qSIMD.integer", "listedgr3",  
             "hosp.recent", "av2channels", "Sgene.dropout")
covariates.split <- c("tstart", "tstop", "event", "week", covariates) 

cases.all[, startday := 1 + as.integer(SPECIMENDATE - firstdate)] # first startday is 1
cases.all[, startweek := ceiling(startday / 7)] # first startday is 1

##################################################################
## Cox regression for death within 28 days

covariates.death <- c("startday", "startweek", "daystodeath28", "death28", covariates)
cases.death <- cases.all[, ..covariates.death]
setnames(cases.death, "daystodeath28", "tstop")
setnames(cases.death, "death28", "event")
setkey(cases.death, startday)
cases.death.split <- split.strata(cases.death)
cases.death.split <- cases.death.split[, ..covariates.split]

## drop care home residents
covariates.death <- covariates[covariates != "care.home"]
cases.nocare.death <- cases.death[care.home=="Independent"][, care.home := NULL]
cases.nocare.death.split <- cases.death.split[care.home =="Independent"][, care.home := NULL]

table.death.coxph <-
    tabulate.freqs.coxregressions(varnames=covariates.death,
                                  stratumvar="week",
                                  data=cases.nocare.death,
                                  split.data=cases.nocare.death.split)
colnames(table.death.coxph)[1] <- gsub("X0", "Survived",
                                           colnames(table.death.coxph)[1])
colnames(table.death.coxph)[2] <- gsub("X1", "Death within 28 days", 
                                           colnames(table.death.coxph)[2])

######################################################################
## Cox regression for death or critical care within 28 days

covariates.severe <- c("startday", "startweek", "daystosevere28", "severe28", covariates)
cases.severe <- cases.all[, ..covariates.severe]
cases.severe <- cases.severe[daystosevere28 > 0] ## drop those with zero days to severe

                                        #cases.severe[, startweek := 1 + floor(startday / 7)]
setnames(cases.severe, "daystosevere28", "tstop")
setnames(cases.severe, "severe28", "event")
setkey(cases.severe, startday)

cases.severe.split <- split.strata(cases.severe)
cases.severe.split <- cases.severe.split[, ..covariates.split]

## drop care home residents
#covariates.severe <- covariates[covariates != "care.home"]
#cases.nocare.severe <- cases.severe[care.home=="Independent"][, care.home := NULL]
#cases.nocare.severe.split <- cases.severe.split[care.home =="Independent"][, care.home := NULL]

table.severe.coxph <-
    tabulate.freqs.coxregressions(varnames=covariates,
                                  stratumvar="week",
                                  data=cases.severe,
                                  split.data=cases.severe.split)
colnames(table.severe.coxph)[1] <- gsub("X0", "Survived, no critical care",
                                           colnames(table.severe.coxph)[1])
colnames(table.severe.coxph)[2] <- gsub("X1", "Critical care or death within 28 days", 
                                           colnames(table.severe.coxph)[2])

##############################################################################
## Cox regression for death or admission within 28 days

covariates.hosp <- c("startday", "startweek", "daystohosp28", "hosp28", covariates)
cases.hosp <- cases.all[, ..covariates.hosp]
cases.hosp <- cases.hosp[daystohosp28 > 0] ## drop those with zero days to hosp
setnames(cases.hosp, "daystohosp28", "tstop") 
setnames(cases.hosp, "hosp28", "event")
setkey(cases.hosp, startday)

cases.hosp.split <- split.strata(cases.hosp)
covariates.split <- c("tstart", "tstop", "event", "week", covariates) 
cases.hosp.split <- cases.hosp.split[, ..covariates.split]

## drop care home residents
#covariates.hosp <- covariates[covariates != "care.home"]
#cases.nocare.hosp <- cases.hosp[care.home=="Independent"][, care.home := NULL]
#cases.nocare.hosp.split <- cases.hosp.split[care.home =="Independent"][, care.home := NULL]

table.hosp.coxph <-
    tabulate.freqs.coxregressions(varnames=covariates,
                                  stratumvar="week",
                                  data=cases.hosp,
                                  split.data=cases.hosp.split)
colnames(table.hosp.coxph)[1] <- gsub("X0", "Survived, no hospital",
                                           colnames(table.hosp.coxph)[1])
colnames(table.hosp.coxph)[2] <- gsub("X1", "Hospitalisation or death within 28 days", 
                                           colnames(table.hosp.coxph)[2])
##################################################################################

## Cox regression for death or admission within 28 days by sliding time window

covariates.window=c("Sgene.dropout", "AGE", "sex", "region4", "care.home",
                    "qSIMD.integer", 
                    "hosp.recent", "av2channels")
covariates.window.hosp <- c("SPECIMENDATE", "startday", "startweek", "daystohosp28", "hosp28",
                            covariates.window)
cases.window.hosp <- cases.all[, ..covariates.window.hosp]
cases.window.hosp <- cases.window.hosp[daystohosp28 > 0] ## drop those with zero days to hosp
setnames(cases.window.hosp, "daystohosp28", "tstop") 
setnames(cases.window.hosp, "hosp28", "event")
setkey(cases.window.hosp, startday)

cases.window.hosp.split <- split.strata(cases.window.hosp)
covariates.window.split <- c("tstart", "tstop", "event", "week",
                             covariates.window.hosp) 

firstdate <- as.Date("2020-12-01")
lastdate <- as.Date("2021-02-04")
## for coeffs we have to specify explicitly the midpoint of the time window 
winsize <- 7 # should be an odd number for date.midpoint to be exactly centred
startdates <- with(cc.severe, firstdate:max(SPECIMENDATE) - winsize)
enddates <- startdates + winsize

## plot rate ratio associated with S gene dropout by time window

coeffs.timewindow <- NULL
for(timewin in 1:length(startdates)) {
    tdata <- cases.window.hosp.split[as.integer(SPECIMENDATE) >= startdates[timewin] &
                                            as.integer(SPECIMENDATE) <= enddates[timewin], ]
    if(length(with(tdata, table(event, Sgene.dropout))) > 2) {
        date.midpoint <- startdates[timewin] + floor(winsize / 2)
        coeffs <- c(date.midpoint,
                    tryCatch(summary(coxph(data=tdata,
                      formula=Surv(time=tstart, time2=tstop,
                                   event=event) ~ . -startday -startweek -SPECIMENDATE -week + strata(week)))$coefficients[1, ],
                             error=function(cond) return(matrix(rep(NA, 10), nrow=2))
                             )
                    )
        coeffs.timewindow <- rbind(coeffs.timewindow, coeffs)
    }
}

coeffs.timewindow <- as.data.table(coeffs.timewindow)
colnames(coeffs.timewindow) <- c("date.midpoint", "coeff", "rateratio", "se.coeff", "Z", "pvalue")
coeffs.timewindow[, date.midpoint := as.Date(date.midpoint, origin="1970-01-01")]

p.rateratio <-
    ggplot(data=coeffs.timewindow, aes(x=date.midpoint, y=coeff)) +
    geom_line(size=0.01 * coeffs.timewindow[, se.coeff]^-2) +
    scale_x_date(breaks = seq.Date(from = firstdate,
                                   to = lastdate, by = "week"),
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = firstdate,
                                                  to = lastdate, by = "week"),
                                         "%d %b")
                 ),
                 limits=c(firstdate, lastdate), expand=c(0, 0)) +
    scale_y_continuous(limits=log(c(0.8, 2.2)),
                       breaks=log(c(seq(0.8, 1, by=0.1), seq(1.2, 2.2, by=0.2))),
                           labels=c(seq(0.8, 1, by=0.1), seq(1.2, 2.2, by=0.2))) +
    labs(x=paste0("Presentation date: mid-point of ", winsize, "-day window"),
         y="Rate ratio (log scale)", 
         title="Rate ratio for hospitalisation associated with S gene dropout",
         tag="(a)")
p.rateratio


## tabulate last dates of outcome recording by health board
cases.all[, max(Date.Death, na.rm=TRUE), by=hb2019name]
cases.all[, max(Admission.Date, na.rm=TRUE), by=hb2019name]

##########################################################################################
table.symptoms <- tabulate.freqs.regressions(varnames=c("agegr20plus", 
                                                        covariates[covariates != "AGE"]),
                                             outcome="symptomatic",
                                             data=cases.all, model="logistic",
                                             stratumvar="week.specimen")
colnames(table.symptoms)[1] <- gsub("X0", "Asymptomatic",
                                    colnames(table.symptoms)[1])
colnames(table.symptoms)[2] <- gsub("X1", "Symptomatic",
                                    colnames(table.symptoms)[2])

covariates.lab <- c("agegr20plus", "sex", "care.home",
                    "qSIMD.integer", 
                    "hh.over18gr",
                    "preschool.any", "hh.schoolagegr",
                    "occup", "hosp.recent", "casegr")

cc.all[, withCtresult := as.integer(!is.na(Sgene.dropout))]

table.withCt <- univariate.tabulate(varnames=covariates.lab,
                                    outcome="withCtresult",
                                    data=cc.all[CASE==1 & testpositive.case==TRUE &
                                                SPECIMENDATE >= as.Date("2020-09-01")],
                                    drop.reflevel=FALSE)
table.withCt <- table.withCt[-nrow(table.withCt), ]
colnames(table.withCt) <- gsub("^0", "No Ct results", colnames(table.withCt))
colnames(table.withCt) <- gsub("^1", "Ct results", colnames(table.withCt))
rownames(table.withCt) <- replace.names(rownames(table.withCt))

## covid symptoms
with(cases.all, table(flag_covid_symptomatic, Sgene.dropout, exclude=NULL))

#print(paste.colpercent(with(cases.all[lastdate.rapid - SPECIMENDATE >=14],
#                            table(daystoadmission, Sgene.dropout, exclude=NULL)), digits=2)[1:21, ])

#print(paste.colpercent(with(cases.all[date_onset_of_symptoms <= SPECIMENDATE],
#                            table(SPECIMENDATE - date_onset_of_symptoms, Sgene.dropout, #exclude=NULL)), digits=2)[1:21, ])

######################################################################################

covariates.exposure <- c("agegr20plus", "sex", "region4", "care.home",
                         "COVID.age", 
                         "qSIMD.integer", "listedgr3",  
                         "hh.over18gr",
                         "preschool.any", "hh.schoolagegr",
                         "occup", "hosp.recent")

## unconditional logistic regression of variant status on covariates
table.coeffs.exposure <-
    tabulate.freqs.regressions(varnames=covariates.exposure,
                               data=cases.all[SPECIMENDATE >= as.Date("2020-12-01")], outcome="mutant",
                               model="logistic", stratumvar="week.specimen")
colnames(table.coeffs.exposure)[1] <- gsub("X0", "Original strain",
                                           colnames(table.coeffs.exposure)[1])
colnames(table.coeffs.exposure)[2] <- gsub("X1", "S gene dropout", 
                                           colnames(table.coeffs.exposure)[2])

covariates.exposure.glm <- c("mutant", covariates.exposure)

coeffs.exposure.timewindow <- NULL
for(timewin in 1:length(startdates)) {
    coeffs.exposure <-
         tryCatch(summary(glm(formula=mutant ~ .,
                    data=cases.all[as.integer(SPECIMENDATE) >= startdates[timewin] &
                                   as.integer(SPECIMENDATE) <= enddates[timewin],
                                   ..covariates.exposure.glm], 
                    family="binomial"))$coefficients[-1, , drop=FALSE],
                  error=function(cond) return(NULL))
    if(!is.null(coeffs.exposure)) {
        coeffs.exposure.timewindow <-
            rbind(coeffs.exposure.timewindow,
                  data.table(date.midpoint=rep(startdates[timewin] + floor(winsize / 2),
                                               nrow(table.coeffs.exposure)), 
                             effect=rownames(coeffs.exposure),
                             coeffs.exposure))
    }
}

colnames(coeffs.exposure.timewindow)[1:4] <- c("date.midpoint", "effect", "coeff", "se.coeff")
coeffs.exposure.timewindow[, date.midpoint := as.Date(date.midpoint, origin="1970-01-01")]

coeffs.region.timewindow <- coeffs.exposure.timewindow[substr(effect, 1, 7)=="region4"]
coeffs.region.timewindow[, effect := gsub("^region4", "", effect)]
p.region <-
    ggplot(data=coeffs.region.timewindow,
           aes(x=date.midpoint, y=coeff, color=effect)) +
    geom_line(size=0.005 * coeffs.region.timewindow[, se.coeff]^-2) +
    scale_x_date(breaks = seq.Date(from = firstdate,
                                   to = lastdate, by = "week"),
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = firstdate,
                                                  to = lastdate, by = "week"),
                                         "%d %b")
                 ),
                 limits=c(firstdate, lastdate), expand=c(0, 0)) +
    scale_y_continuous(limits=log(c(1, 6)),
                       breaks=log(c(seq(1, 4, by=0.5), 5., 6)),
                       labels=c(seq(1, 4, by=0.5), 5., 6)) +
    theme(legend.position = c(0.2, 0.25)) +
    theme(legend.title = element_blank()) + 
    labs(x=paste0("Presentation date: mid-point of ", winsize, "-day window"), 
         y="Rate ratio (log scale)", 
         title="Rate ratio for S gene dropout associated with region",
         tag="(b)")
p.region

##########################################################################

coeffs.hosp.timewindow <- coeffs.exposure.timewindow[effect=="hosp.recentTRUE"]

p.hosp <-
    ggplot(data=coeffs.hosp.timewindow,
           aes(x=date.midpoint, y=coeff)) + #, color=effect)) +
    geom_line(size=0.01 * coeffs.hosp.timewindow[, se.coeff]^-2) +
    scale_x_date(breaks = seq.Date(from = firstdate,
                                   to = lastdate, by = "week"),
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = firstdate,
                                                  to = lastdate, by = "week"),
                                         "%d %b")
                 ),
                 limits=c(firstdate, lastdate), expand=c(0, 0)) +
    scale_y_continuous(limits=log(c(0.5, 1.2)),
                       breaks=log(seq(0.5, 1.2, by=0.1)),
                       labels=c(seq(0.5, 1.2, by=0.1))) +
  #  theme(legend.position = c(0.45, 0.8)) +
 #   theme(legend.title = element_blank()) + 
    labs(x=paste0("Presentation date: mid-point of ", winsize, "-day window"), 
         y="Rate ratio (log scale)", 
         title="Rate ratio for S gene dropout associated with recent exposure to hospital",
         tag="(c)")
p.hosp


rmarkdown::render("ct.Rmd")
