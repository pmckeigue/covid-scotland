library(data.table)
library(ggplot2)
library(survival)

source("helperfunctions.R")

datadir <- "data/2021-07-28/"
if(!exists("cc.all")) {
    load(paste0(datadir, "cc.all.RData"))
}

firstdate <- as.Date("2020-12-01")
lastdate <- as.Date("2021-07-28")

## for coeffs we have to specify explicitly the midpoint of the time window 
winsize <- 42 
startdates <- with(cc.all, firstdate:lastdate - winsize)
enddates <- startdates + winsize

cc.sincedec2020 <- cc.all[specimen_date >= as.Date("2020-12-01") &
                                (casegroup=="A" | casegroup=="B"), 
                           .(CASE, casegroup, specimen_date, dayssincelastdose, care.home,
                             listedgr3, vacc_product_name_2, vaxgr, vax14.factor, vax14.dose, numdrugs.notcv, inpat.recent,
                             stratum)]
cc.sincedec2020[, atleast19may := specimen_date >= as.Date("2021-05-19")]

####################  time course of efficacy against severe disease ######
coeffs.timewindow <- NULL
for(timewin in 1:length(startdates)) {
    tdata <- cc.sincedec2020[casegroup=="A" & 
                    as.integer(specimen_date) >= startdates[timewin] &
                    as.integer(specimen_date) <= enddates[timewin],
                    .(CASE, care.home, listedgr3, numdrugs.notcv, vaxgr, vax14.dose, vax14.factor, inpat.recent, stratum)]
    date.midpoint <- startdates[timewin] + floor(winsize / 2)
    coeffs <- tryCatch(summary(clogit(formula=CASE ~ care.home + listedgr3 + 
                                          numdrugs.notcv + inpat.recent +
                                          vaxgr +
                                          strata(stratum), data=tdata))$coefficients,
                       error=function(cond) return(NULL)
                       )
    if(!is.null(coeffs)) {
        coeffs <- data.table(date.midpoint=rep(date.midpoint, nrow(coeffs)),
                             effect=rownames(coeffs), coeffs)
        coeffs.timewindow <- rbind(coeffs.timewindow, coeffs)
    }
}

coeffs.timewindow[, date.midpoint := as.Date(date.midpoint, origin="1970-01-01")]
#coeffs.timewindow <- coeffs.timewindow[grep("care.home", effect)]
coeffs.timewindow <- coeffs.timewindow[grep("vaxgr", effect)]
coeffs.timewindow[, effect := gsub("vaxgr", "", effect)]
coeffs.timewindow[, effect := replace.names(effect)]
coeffs.timewindow[, linesize := `se(coef)`^-2]
coeffs.timewindow[, linesize := 2 * linesize / max(linesize, na.rm=TRUE), by="effect"]
setorder(coeffs.timewindow, effect)

## Figure: rate ratios by date of presentation
breakpoints <- c(0.02, 0.05, 0.1, 0.2, 0.5, 1)
#breakpoints <- c(0.5, 1, 2, 5, 10)
p.rateratio <-
    ggplot(data=coeffs.timewindow, aes(x=date.midpoint, y=coef, color=effect)) +
    geom_line(size=coeffs.timewindow$linesize) +
    labs(x=paste0("Presentation date: mid-point of ", winsize, "-day window"),
         y="Rate ratio for severe COVID-19 (log scale)") + 
    scale_y_continuous(breaks=log(breakpoints),
                       labels=breakpoints,
                       limits=log(c(min(breakpoints), max(breakpoints))),
                       expand=c(0, 0)) +
    scale_x_date(breaks = seq.Date(from = firstdate,
                                   to = lastdate, by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = firstdate,
                                                  to = lastdate, by = "month"),
                                         "%d %b")
                             ),
                 limits=c(as.Date("2021-03-01"), lastdate)) +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.8, 0.85)) 
p.rateratio

##################### time course of efficacy against hospitalisation ###

coeffs.timewindow.hosp <- NULL
for(timewin in 1:length(startdates)) {
    tdata <- cc.sincedec2020[(casegroup=="A" | casegroup=="B") & 
                    as.integer(specimen_date) >= startdates[timewin] &
                    as.integer(specimen_date) <= enddates[timewin],
                    .(CASE, care.home, listedgr3, numdrugs.notcv, vaxgr, vax14.dose, vax14.factor, inpat.recent, stratum)]
    date.midpoint <- startdates[timewin] + floor(winsize / 2)
    coeffs <- tryCatch(summary(clogit(formula=CASE ~ care.home + listedgr3 + numdrugs.notcv + inpat.recent + vaxgr + strata(stratum), data=tdata))$coefficients,
                       error=function(cond) return(NULL)
                       )
    if(!is.null(coeffs)) {
        coeffs <- data.table(date.midpoint=rep(date.midpoint, nrow(coeffs)),
                             effect=rownames(coeffs), coeffs)
        coeffs.timewindow.hosp <- rbind(coeffs.timewindow.hosp, coeffs)
    }
}

coeffs.timewindow.hosp[, date.midpoint := as.Date(date.midpoint, origin="1970-01-01")]
coeffs.timewindow.hosp <- coeffs.timewindow.hosp[grep("vaxgr", effect)]
coeffs.timewindow.hosp[, effect := gsub("TRUE$", "", effect)]
coeffs.timewindow.hosp[, effect := gsub("vaxgr", "", effect)]
coeffs.timewindow.hosp[, effect := replace.names(effect)]
coeffs.timewindow.hosp[, linesize := `se(coef)`^-2]
coeffs.timewindow.hosp[, linesize := 2 * linesize / max(linesize, na.rm=TRUE), by="effect"]
setorder(coeffs.timewindow.hosp, effect)

## Figure: rate ratios by date of presentation
breakpoints <- c(0.05, 0.1, 0.2, 0.5, 1)
p.rateratio.hosp <-
    ggplot(data=coeffs.timewindow.hosp, aes(x=date.midpoint, y=coef, color=effect)) +
    geom_line(size=coeffs.timewindow.hosp$linesize) +
    labs(x=paste0("Presentation date: mid-point of ", winsize, "-day window"),
         y="Rate ratio for hospitalised/fatal COVID-19 (log scale)") + 
    scale_y_continuous(breaks=log(breakpoints),
                       labels=breakpoints,
                       limits=log(c(min(breakpoints), max(breakpoints))),
                       expand=c(0, 0)) +
    scale_x_date(breaks = seq.Date(from = firstdate,
                                   to = lastdate, by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = firstdate,
                                                  to = lastdate, by = "month"),
                                         "%d %b")
                             ),
                 limits=c(as.Date("2021-03-01"), lastdate)) +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.8, 0.85)) 
p.rateratio.hosp

###############################################################

## rollmean.dtN returns a data.table with columns date, avdaily (rolling mean) 
rollmean.dtN <- function(dt, k) {
    N.dates <- dt[, .N, by=specimen_date]
    setkey(N.dates, specimen_date)
    all.dates <- data.table(date=seq(from=min(N.dates$specimen_date),
                                     to=max(N.dates$specimen_date),
                                     by=1))
    setkey(all.dates, date)
    ## add rows for missing dates by a left join
    all.dates <- N.dates[all.dates]
    all.dates[is.na(N), N := 0]
    return(data.table(date=dt$specimen_date,
                      avdaily=zoo::rollmean(dt$N, k=k, fill=NA)))
}

## Figure: severe cases by date of presentation
winsize.casedates <- 3

casedates.allgr <-  cc.sincedec2020[CASE==1 & casegroup=="A", .N, by=specimen_date]
setkey(casedates.allgr, specimen_date)
casedates.allgr <- casedates.allgr[, rollmean.dtN(dt=.SD, k=winsize.casedates),
                                   .SDcols=c("specimen_date", "N")]

## use by= to get rolling means by category
casedates.byelig <- cc.sincedec2020[CASE==1 & casegroup=="A", .N, by=.(specimen_date, listedgr3)]
setkey(casedates.byelig, specimen_date)
casedates.byelig <- casedates.byelig[, rollmean.dtN(dt=.SD, k=winsize.casedates), by=listedgr3,
                       .SDcols=c("specimen_date", "N")]
setnafill(casedates.byelig, cols="avdaily", fill=0)
setorder(casedates.byelig, avdaily)

p.byelig <- ggplot(data=casedates.byelig,
                aes(x=date,
                    y=avdaily, fill=listedgr3)) +
    geom_area(position="stack") +
    scale_y_continuous(expand=c(0, 0)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-12-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = as.Date("2020-12-01"),
                                                  to = lastdate, by = "month"),
                                         "%d %b")
                 ),
                 limits=c(as.Date("2020-12-01"), lastdate)) +
    theme(legend.position = c(0.5, 0.7)) +
    scale_fill_manual(values=c("grey", "blue", "red")) +
    xlab(paste0("Presentation date: mid-point of ", winsize.casedates, "-day window")) +
         ylab("Daily severe cases")
p.byelig

##############  time since dose ################################################

## for coeffs we have to specify explicitly the midpoint of the time window 
winsize.sincedose <- 42

firstmidpoint <- floor(winsize.sincedose/2)
lastmidpoint <- max(cc.sincedec2020$dayssincelastdose, na.rm=TRUE) -
    floor(winsize.sincedose/2)

coeffs.timewindow.sincedose <- NULL
for(t in firstmidpoint:lastmidpoint) {
    tdata <- cc.sincedec2020[casegroup=="A" &
                             (vax14.dose == 0 | # include unvax as ref category
                              (dayssincelastdose >= t - floor(winsize.sincedose/2) & 
                               dayssincelastdose <= t + floor(winsize.sincedose/2)))]
    time.midpoint <- t
    coeffs <- tryCatch(summary(clogit(formula=CASE ~ care.home + listedgr3 + numdrugs.notcv +
                                          inpat.recent + vaxgr + strata(stratum),
                                      data=tdata))$coefficients,
                       error=function(cond) return(NULL)
                       )
    if(!is.null(coeffs)) {
        coeffs <- data.table(time.midpoint=rep(time.midpoint, nrow(coeffs)),
                             effect=rownames(coeffs), coeffs)
        coeffs.timewindow.sincedose <- rbind(coeffs.timewindow.sincedose, coeffs)
    }
}

coeffs.timewindow.sincedose <- na.omit(coeffs.timewindow.sincedose, cols="coef")
coeffs.timewindow.sincedose <- coeffs.timewindow.sincedose[coef > -5 & coef < 5]
coeffs.timewindow.sincedose <- coeffs.timewindow.sincedose[grep("vaxgr", effect)]
coeffs.timewindow.sincedose[, effect := gsub("vaxgr", "", effect)]
coeffs.timewindow.sincedose[, effect := replace.names(effect)]
coeffs.timewindow.sincedose[, time.midpoint := time.midpoint / 7]
summary(coeffs.timewindow.sincedose)

coeffs.timewindow.sincedose[, linesize := `se(coef)`^-2]
coeffs.timewindow.sincedose[, linesize := 2 * linesize / max(linesize, na.rm=TRUE), by="effect"]
setorder(coeffs.timewindow.sincedose, effect)
breakpoints <- c(0.02, 0.05, 0.1, 0.2, 0.5, 1)
p.rateratio.sincedose <-
    ggplot(data=coeffs.timewindow.sincedose, aes(x=time.midpoint, y=coef, color=effect)) +
    geom_line(size=coeffs.timewindow.sincedose$linesize) +
    labs(x=paste0("Weeks since last dose: mid-point of ", winsize.sincedose, "-day window"),
         y="Rate ratio for severe COVID-19 (log scale)") + 
    scale_y_continuous(breaks=log(breakpoints),
                       labels=breakpoints,
                       limits=log(c(min(breakpoints), max(breakpoints))),
                       expand=c(0, 0)) +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.7, 0.85)) 
p.rateratio.sincedose

################ weeks since dose for hospitalisation ###############
firstmidpoint <- floor(winsize.sincedose/2)
lastmidpoint <- max(cc.sincedec2020$dayssincelastdose, na.rm=TRUE) -
    floor(winsize.sincedose/2)

coeffs.timewindow.sincedose <- NULL
for(t in firstmidpoint:lastmidpoint) {
    tdata <- cc.sincedec2020[vax14.dose == 0 | # include unvax as ref category
                           (dayssincelastdose >= t - floor(winsize.sincedose/2) & 
                            dayssincelastdose <= t + floor(winsize.sincedose/2))]
    time.midpoint <- t
    coeffs <- tryCatch(summary(clogit(formula=CASE ~ care.home + listedgr3 + numdrugs.notcv +
                                          inpat.recent + vaxgr + strata(stratum),
                                      data=tdata))$coefficients,
                       error=function(cond) return(NULL)
                       )
    if(!is.null(coeffs)) {
        coeffs <- data.table(time.midpoint=rep(time.midpoint, nrow(coeffs)),
                             effect=rownames(coeffs), coeffs)
        coeffs.timewindow.sincedose <- rbind(coeffs.timewindow.sincedose, coeffs)
    }
}

coeffs.timewindow.sincedose <- na.omit(coeffs.timewindow.sincedose, cols="coef")
coeffs.timewindow.sincedose <- coeffs.timewindow.sincedose[coef > -5 & coef < 5]
coeffs.timewindow.sincedose <- coeffs.timewindow.sincedose[grep("vaxgr", effect)]
coeffs.timewindow.sincedose[, effect := gsub("vaxgr", "", effect)]
coeffs.timewindow.sincedose[, effect := replace.names(effect)]
coeffs.timewindow.sincedose[, time.midpoint := time.midpoint / 7]
summary(coeffs.timewindow.sincedose)

coeffs.timewindow.sincedose[, linesize := `se(coef)`^-2]
coeffs.timewindow.sincedose[, linesize := 2 * linesize / max(linesize, na.rm=TRUE), by="effect"]
setorder(coeffs.timewindow.sincedose, effect)
# Figure: rate ratios by date of presentations
breakpoints <- c(0.02, 0.05, 0.1, 0.2, 0.5, 1)
p.rateratio.sincedose.hosp <-
    ggplot(data=coeffs.timewindow.sincedose, aes(x=time.midpoint, y=coef, color=effect)) +
    geom_line(size=coeffs.timewindow.sincedose$linesize) +
    labs(x=paste0("Weeks since last dose: mid-point of ", winsize.sincedose, "-day window"),
         y="Rate ratio for hospitalised/fatal COVID-19 (log scale)") + 
    scale_y_continuous(breaks=log(breakpoints),
                       labels=breakpoints,
                       limits=log(c(min(breakpoints), max(breakpoints))),
                       expand=c(0, 0)) +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.7, 0.85)) 
p.rateratio.sincedose.hosp

table.vaxproduct2 <- with(cc.sincedec2020, table(vacc_product_name_2, vax14.dose))
table.vaxproduct2.colsums <- colSums(table.vaxproduct2)
table.vaxproduct2 <- paste.colpercent(table.vaxproduct2, digits=1)
table.vaxproduct2 <- rbind(table.vaxproduct2, table.vaxproduct2.colsums)

##### changepoint model for slope of rate ratio to time since 2nd dose  ########

## unvaccinated will be coded as vaxclass AZ
cc.sincedec2020[, vaxgr3 := car::recode(vaxgr,
                                        "c('1 dose mRNA vaccine', '1 dose AZ vaccine')=NA",
                                        as.factor=TRUE,
                                        levels=c("Not vaccinated",
                                                 "2 doses mRNA vaccine",
                                                 "2 doses AZ vaccine"))]

########### severe disease ###############################################
loglik.nocp <- clogit(data=cc.sincedec2020[casegroup=="A" & vax14.dose != 0.5],
                     formula=CASE ~ care.home + listedgr3 + numdrugs.notcv + inpat.recent +  
                         vaxgr3 + 
                         vax14.dose:dayssincelastdose +
                                     strata(stratum))$loglik[2]
loglik.changepoints <- NULL
week.cp <- 8:16
rows.coeffs.cp <- 8:9
for(week in week.cp) { 
    cc.sincedec2020[dayssincelastdose <= week * 7,
                    c("weeks.lt.changepoint", "weeks.ge.changepoint") :=
                        list(dayssincelastdose / 7, 0)]
    cc.sincedec2020[dayssincelastdose > week * 7,
                    c("weeks.lt.changepoint", "weeks.ge.changepoint") :=
                        list(week * 7, dayssincelastdose / 7 - week)]
    model.week <- clogit(data=cc.sincedec2020[casegroup=="A" & vax14.dose != 0.5],
                                 formula=CASE ~ care.home + listedgr3 + numdrugs.notcv + inpat.recent +  
                                     vaxgr3 + 
                                     vax14.dose:weeks.lt.changepoint + vax14.dose:weeks.ge.changepoint +
                                     strata(stratum))
    coeffs.week <- as.numeric(t(summary(model.week)$coefficients[rows.coeffs.cp, c(1, 3)]))
    
    loglik.changepoints <- rbind(loglik.changepoints, c(model.week$loglik[2], coeffs.week))
}
loglik.changepoints <- data.table(week.changepoint=week.cp, loglik.changepoints)
colnames(loglik.changepoints)[2:6] <- c("loglik", "coeff.before", "se.coeff.before", "coeff.after", "se.coeff.after") 
loglik.changepoints[, loglik := round(loglik - loglik.nocp, 2)]
loglik.changepoints[, effect.before := or.ci(coeff.before, se.coeff.before)]
loglik.changepoints[, z.before := coeff.before / se.coeff.before]
loglik.changepoints[, p.before := 2 * pnorm(-abs(z.before))]
loglik.changepoints[, p.before := format.pvalue(z.before, p.before)]

loglik.changepoints[, effect.after := or.ci(coeff.after, se.coeff.after)]
loglik.changepoints[, z.after := coeff.after / se.coeff.after]
loglik.changepoints[, p.after := 2 * pnorm(-abs(z.after))]
loglik.changepoints[, p.after := format.pvalue(z.after, p.after)]
loglik.changepoints <- loglik.changepoints[, .(week.changepoint, loglik,
                                               coeff.before, se.coeff.before, effect.before, p.before,
                                               coeff.after, se.coeff.after, effect.after, p.after)]
loglik.changepoints[, deviance := -2 * loglik]
print(loglik.changepoints[, .(week.changepoint, deviance, effect.before, p.before, effect.after, p.after)])

bestfit.cp <- loglik.changepoints[which.min(deviance), week.changepoint]
bestfit.effectbefore <- or.ci.text(loglik.changepoints[which.min(deviance), coeff.before],
                                   loglik.changepoints[which.min(deviance), se.coeff.before])
bestfit.effectafter <- or.ci.text(loglik.changepoints[which.min(deviance), coeff.after],
                                   loglik.changepoints[which.min(deviance), se.coeff.after])
                                   
########### hospitalisations ###############################################
loglik.nocp.hosp <- clogit(data=cc.sincedec2020[(casegroup=="A" | casegroup=="B")  &
                                                vax14.dose != 0.5],
                           formula=CASE ~ care.home + listedgr3 + numdrugs.notcv + inpat.recent +  
                               vaxgr3 + 
                               vax14.dose:dayssincelastdose +
                               strata(stratum))$loglik[2]

loglik.changepoints.hosp <- NULL
for(week in week.cp) { 
    cc.sincedec2020[dayssincelastdose <= week * 7,
                    c("weeks.lt.changepoint", "weeks.ge.changepoint") :=
                        list(dayssincelastdose / 7, 0)]
    cc.sincedec2020[dayssincelastdose > week * 7,
                    c("weeks.lt.changepoint", "weeks.ge.changepoint") :=
                        list(week * 7, dayssincelastdose / 7 - week)]
    model.week <- clogit(data=cc.sincedec2020[(casegroup=="A" | casegroup=="B") &
                                              vax14.dose != 0.5],
                                 formula=CASE ~ care.home + listedgr3 + numdrugs.notcv + inpat.recent +  
                                     vaxgr3 + 
                                     vax14.dose:weeks.lt.changepoint + vax14.dose:weeks.ge.changepoint +
                                     strata(stratum))
    coeffs.week <- as.numeric(t(summary(model.week)$coefficients[rows.coeffs.cp, c(1, 3)]))
    
    loglik.changepoints.hosp <- rbind(loglik.changepoints.hosp, c(model.week$loglik[2], coeffs.week))
}
loglik.changepoints.hosp <- data.table(week.changepoint=week.cp, loglik.changepoints.hosp)
colnames(loglik.changepoints.hosp)[2:6] <- c("loglik", "coeff.before", "se.coeff.before", "coeff.after", "se.coeff.after") 
loglik.changepoints.hosp[, loglik := round(loglik - loglik.nocp.hosp, 2)]
loglik.changepoints.hosp[, effect.before := or.ci(coeff.before, se.coeff.before)]
loglik.changepoints.hosp[, z.before := coeff.before / se.coeff.before]
loglik.changepoints.hosp[, p.before := 2 * pnorm(-abs(z.before))]
loglik.changepoints.hosp[, p.before := format.pvalue(z.before, p.before)]

loglik.changepoints.hosp[, effect.after := or.ci(coeff.after, se.coeff.after)]
loglik.changepoints.hosp[, z.after := coeff.after / se.coeff.after]
loglik.changepoints.hosp[, p.after := 2 * pnorm(-abs(z.after))]
loglik.changepoints.hosp[, p.after := format.pvalue(z.after, p.after)]
loglik.changepoints.hosp <- loglik.changepoints.hosp[, .(week.changepoint, loglik, coeff.before, se.coeff.before, effect.before, p.before, coeff.after, se.coeff.after, effect.after, p.after)]
loglik.changepoints.hosp[, deviance := -2 * loglik]
print(loglik.changepoints.hosp[, .(week.changepoint, deviance, effect.before, p.before, effect.after, p.after)])

bestfit.cp.hosp <- loglik.changepoints.hosp[which.min(deviance), week.changepoint]
bestfit.effectbefore.hosp <- or.ci.text(loglik.changepoints.hosp[which.min(deviance), coeff.before],
                                        loglik.changepoints.hosp[which.min(deviance), se.coeff.before])
bestfit.effectafter.hosp <- or.ci.text(loglik.changepoints.hosp[which.min(deviance), coeff.after],
                                       loglik.changepoints.hosp[which.min(deviance), se.coeff.after])

################### comparison of rate ratio for 2 doses versus 0 before and after 19 May Delta changepoint

coeffs.severe.19may <- summary(clogit(data=cc.sincedec2020[casegroup=="A" &
                                                                  vax14.dose != 0.5],
                                        formula=CASE ~ care.home + listedgr3 + numdrugs.notcv + inpat.recent +  
                                            atleast19may/vax14.dose + strata(stratum)))$coefficients

coeffs.hosp.19may <- summary(clogit(data=cc.sincedec2020[vax14.dose != 0.5],
                                    formula=CASE ~ care.home + listedgr3 + numdrugs.notcv + inpat.recent +  
                                        atleast19may/vax14.dose + strata(stratum)))$coefficients

estci.severe.before19may <- or.ci(coeffs.severe.19may[7, 1], coeffs.severe.19may[7, 3])
pval.severe.before19may <- format.pvalue(coeffs.severe.19may[7, 4], coeffs.severe.19may[7, 5])
estcipv.severe.before19may <- format.estcipv(estci.severe.before19may, pval.severe.before19may)

estci.severe.after19may <- or.ci(coeffs.severe.19may[8, 1], coeffs.severe.19may[8, 3])
pval.severe.after19may <- format.pvalue(coeffs.severe.19may[8, 4], coeffs.severe.19may[8, 5])
estcipv.severe.after19may <- format.estcipv(estci.severe.after19may, pval.severe.after19may)

estci.hosp.before19may <- or.ci(coeffs.hosp.19may[7, 1], coeffs.hosp.19may[7, 3])
pval.hosp.before19may <- format.pvalue(coeffs.hosp.19may[7, 4], coeffs.hosp.19may[7, 5])
estcipv.hosp.before19may <- format.estcipv(estci.hosp.before19may, pval.hosp.before19may)

estci.hosp.after19may <- or.ci(coeffs.hosp.19may[8, 1], coeffs.hosp.19may[8, 3])
pval.hosp.after19may <- format.pvalue(coeffs.hosp.19may[8, 4], coeffs.hosp.19may[8, 5])
estcipv.hosp.after19may <- format.estcipv(estci.hosp.after19may, pval.hosp.after19may)

table.severe <- tabulate.freqs.regressions(data=cc.sincedec2020[casegroup=="A"],
                                           varnames=c("care.home", "listedgr3", "numdrugs.notcv", "inpat.recent",
                                                      "vaxgr"))

table.hosp <- tabulate.freqs.regressions(data=cc.sincedec2020[specimen_date >= as.Date("2020-12-01")],,
                                         varnames=c("care.home", "listedgr3", "numdrugs.notcv", "inpat.recent",
                                                    "vaxgr"))

recent.coeffs <- coeffs.timewindow[date.midpoint==as.Date("2021-07-07") &
                                   grepl("2 doses", effect), .(effect, coef, `se(coef)`)]
az.efficacy.severe <- pctefficacy.ci.text(recent.coeffs[1, 2], recent.coeffs[1, 3])
mrna.efficacy.severe <- pctefficacy.ci.text(recent.coeffs[2, 2], recent.coeffs[2, 3])

recent.coeffs.hosp <- coeffs.timewindow.hosp[date.midpoint==as.Date("2021-07-07") &
                                   grepl("2 doses", effect), .(effect, coef, `se(coef)`)]
az.efficacy.hosp <- pctefficacy.ci.text(recent.coeffs.hosp[1, 2], recent.coeffs.hosp[1, 3])
mrna.efficacy.hosp <- pctefficacy.ci.text(recent.coeffs.hosp[2, 2], recent.coeffs.hosp[2, 3])

rmarkdown::render("vaxtrend.Rmd")
