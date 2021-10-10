library(data.table)
library(ggplot2)
library(survival)

source("helperfunctions.R")

datadir <- "data/2021-09-22/"

if(!exists("cc.all")) {
    load(paste0(datadir, "cc.all.RData"))
}

firstdate <- as.Date("2020-12-01")
lastdate <- as.Date("2021-09-08")

cc.sincedec2020 <- cc.all[specimen_date >= as.Date("2020-12-01") &
                                (casegroup=="A" | casegroup=="B"), 
                          .(CASE, casegroup, specimen_date, dayssincelastdose, care.home,
                            hh.over18, hh.schoolage, occup,
                            listedgr3, vacc_product_name_2, vaxgr, vax14.factor, vax14.dose,
                            numdrugs.notcv, inpat.recent,
                            stratum)]
cc.sincedec2020[, atleast19may := specimen_date >= as.Date("2021-05-19")]
cc.sincedec2020[, vaxgr4 := car::recode(as.character(vaxgr),
                                   "c('1 dose AZ vaccine', '1 dose mRNA vaccine')='1 dose any vaccine'",
                                   as.factor=TRUE,
                                   levels=c("Not vaccinated",
                                            "1 dose any vaccine",
                                            "2 doses AZ vaccine",
                                            "2 doses mRNA vaccine"))]

if(FALSE) {
## severe disease before Dec 2020
summary(clogit(formula=CASE ~ care.home +
                   hh.over18 + hh.schoolage + occup + 
                   listedgr3 + 
                   inpat.recent +
                   strata(stratum), data=cc.all[specimen_date < as.Date("2020-12-01") &
                                (casegroup=="A" )]))$coefficients

## severe disease after Dec 2020
summary(clogit(formula=CASE ~ care.home +
                   hh.over18 + hh.schoolage + occup +
                   listedgr3 + 
                   inpat.recent +
                   vaxgr +
                   strata(stratum), data=cc.sincedec2020[casegroup=="A"]))$coefficients
}

## for coeffs we have to specify explicitly the midpoint of the time window 
winsize <- 42 
# winsize <- 56  # for hh composition, occup, care home analysis
startdates <- with(cc.all, firstdate:lastdate - winsize)
enddates <- startdates + winsize

####################  time course of efficacy against severe disease ######

coeffs.timewindow <- NULL
for(timewin in 1:length(startdates)) {
    tdata <- cc.sincedec2020[casegroup=="A" & 
                    as.integer(specimen_date) >= startdates[timewin] &
                    as.integer(specimen_date) <= enddates[timewin],
                    .(CASE, care.home, occup,
                      hh.over18, hh.schoolage,
                      listedgr3, numdrugs.notcv, vaxgr4, vax14.dose, vax14.factor, inpat.recent, stratum)]
    date.midpoint <- startdates[timewin] + floor(winsize / 2)
    #print(as.Date(date.midpoint, origin="1970-01-01"))
    coeffs <- tryCatch(summary(clogit(formula=CASE ~ care.home +
                                          #hh.over18 + hh.schoolage +
                                          #occup + 
                                          listedgr3 +
                                          inpat.recent +
                                          vaxgr4 + 
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
coeffs.timewindow <- coeffs.timewindow[!is.na(coef) &
                                       date.midpoint >= firstdate &
                                       date.midpoint <= lastdate]

coeffs.timewindow <- coeffs.timewindow[grep("vaxgr4", effect)]
coeffs.timewindow[, effect := gsub("vaxgr4", "", effect)]
coeffs.timewindow[, effect := replace.names(effect)]
coeffs.timewindow[, linesize := `se(coef)`^-2]
coeffs.timewindow[, max(linesize, na.rm=TRUE), by="effect"]
coeffs.timewindow[, linesize := 2 * linesize / max(linesize, na.rm=TRUE), by="effect"]
setorder(coeffs.timewindow, effect)

calendardate.lowerlimit <- as.Date("2021-03-01")
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
                 limits=c(calendardate.lowerlimit, lastdate)) +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.8, 0.85)) 
##########################################################################

## plot for 1 dose any vaccine only, allowing an earlier time window

cc.sincedec2020[, vax14.1dose := as.integer(vax14.dose > 0)]

coeffs.1dose.timewindow <- NULL
for(timewin in 1:length(startdates)) {
    tdata <- cc.sincedec2020[casegroup=="A" & 
                    as.integer(specimen_date) >= startdates[timewin] &
                    as.integer(specimen_date) <= enddates[timewin] &
                           vax14.factor != "2",
                    .(CASE, care.home, occup,
                      hh.over18, hh.schoolage,
                      listedgr3, numdrugs.notcv, vaxgr, vax14.1dose, vax14.factor, inpat.recent, stratum)]
    date.midpoint <- startdates[timewin] + floor(winsize / 2)
    #print(as.Date(date.midpoint, origin="1970-01-01"))
    coeffs.1dose <- tryCatch(summary(clogit(formula=CASE ~ care.home +
                                          #hh.over18 + hh.schoolage +
                                          #occup + 
                                          listedgr3 +
                                          inpat.recent +
                                          vax14.1dose + 
                                          strata(stratum), data=tdata))$coefficients,
                       error=function(cond) return(NULL)
                       )
    if(!is.null(coeffs.1dose)) {
        coeffs.1dose <- data.table(date.midpoint=rep(date.midpoint, nrow(coeffs.1dose)),
                             effect=rownames(coeffs.1dose), coeffs.1dose)
        coeffs.1dose.timewindow <- rbind(coeffs.1dose.timewindow, coeffs.1dose)
    }
}
coeffs.1dose.timewindow[, date.midpoint := as.Date(date.midpoint, origin="1970-01-01")]

coeffs.1dose.timewindow <- coeffs.1dose.timewindow[grep("vax14", effect)]
#coeffs.1dose.timewindow[, effect := gsub("vax14", "", effect)]
coeffs.1dose.timewindow[, effect := replace.names(effect)]
coeffs.1dose.timewindow[, linesize := `se(coef)`^-2]
coeffs.1dose.timewindow[, linesize := 2 * linesize / max(linesize, na.rm=TRUE), by="effect"]
setorder(coeffs.1dose.timewindow, effect)

## Figure: rate ratios by date of presentation
breakpoints <- c(0.02, 0.05, 0.1, 0.2, 0.5, 1)
#breakpoints <- c(0.5, 1, 2, 5, 10)
p.rateratio.1dose <-
    ggplot(data=coeffs.1dose.timewindow, aes(x=date.midpoint, y=coef, color=effect)) +
    geom_line(size=coeffs.1dose.timewindow$linesize) +
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
                 limits=c(calendardate.lowerlimit, lastdate)) +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.8, 0.85)) 

########### plot household composition, occupation, care home  ################
winsize56 <- 56
lastdate56 <- as.Date("2021-08-31")
# winsize <- 56  # for hh composition, occup, care home analysis
startdates56 <- with(cc.all, firstdate:lastdate56 - winsize56)
enddates56 <- startdates56 + winsize56
coeffs.hh.timewindow <- NULL
for(timewin in 1:length(startdates56)) {
    tdata <- cc.sincedec2020[casegroup=="A" & 
                    as.integer(specimen_date) >= startdates[timewin] &
                    as.integer(specimen_date) <= enddates[timewin],
                    .(CASE, care.home, occup,
                      hh.over18, hh.schoolage,
                      listedgr3, numdrugs.notcv, vaxgr, vax14.dose, vax14.factor,
                      inpat.recent, stratum)]
    date.midpoint <- startdates[timewin] + floor(winsize56 / 2)
    #print(as.Date(date.midpoint, origin="1970-01-01"))
    coeffs.hh <- tryCatch(summary(clogit(formula=CASE ~ care.home +
                                          hh.over18 + hh.schoolage +
                                          occup + 
                                          listedgr3 +
                                          # inpat.recent +
                                          vax14.factor + # commented out for DS
                                          strata(stratum), data=tdata))$coefficients,
                       error=function(cond) return(NULL)
                       )
    if(!is.null(coeffs.hh)) {
        coeffs.hh <- data.table(date.midpoint=rep(date.midpoint, nrow(coeffs.hh)),
                             effect=rownames(coeffs.hh), coeffs.hh)
        coeffs.hh.timewindow <- rbind(coeffs.hh.timewindow, coeffs.hh)
    }
}
coeffs.hh.timewindow[, date.midpoint := as.Date(date.midpoint, origin="1970-01-01")]

coeffs.hh.tw <- coeffs.hh.timewindow[grep("care\\.home|hh|Teach|vax", effect)]
coeffs.hh.tw[, effect := gsub("care.home", "", effect)]
coeffs.hh.tw[, effect := gsub("TRUE", "", effect)]
coeffs.hh.tw[, effect := gsub("occup", "", effect)]
coeffs.hh.tw[, effect := gsub("listedgr3", "", effect)]
coeffs.hh.tw[, effect := gsub("vax14\\.factor", "vax dose ", effect)]
coeffs.hh.tw[, effect := replace.names(effect)]
coeffs.hh.tw[, linesize := `se(coef)`^-2]
coeffs.hh.tw[, linesize := 2 * linesize / max(linesize, na.rm=TRUE), by="effect"]
setorder(coeffs.hh.tw, effect)

## Figure: rate ratios by date of presentation
breakpoints <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)
#breakpoints <- c(0.5, 1, 2, 5, 10)
p.rateratio.hh <-
    ggplot(data=coeffs.hh.tw, aes(x=date.midpoint, y=coef, color=effect)) +
    geom_line(size=coeffs.hh.tw$linesize) +
    labs(x=paste0("Presentation date: mid-point of ", winsize56, "-day window"),
         y="Rate ratio for severe COVID-19 (log scale)") + 
    scale_y_continuous(breaks=log(breakpoints),
                       labels=breakpoints,
                       limits=log(c(min(breakpoints), max(breakpoints))),
                       expand=c(0, 0)) +
    scale_x_date(breaks = seq.Date(from = firstdate,
                                   to = lastdate56, by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = firstdate,
                                                  to = lastdate56, by = "month"),
                                         "%d %b")
                             ),
                 limits=c(as.Date("2020-12-01"), lastdate)) +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.7, 0.85)) 

table.coeffs.recent <- tabulate.freqs.regressions(varnames=c("care.home", "hh.over18",
                                                       "hh.schoolage", "occup",
                                                       "listedgr3", "inpat.recent",
                                                       "vax14.factor"),
                                            data=cc.sincedec2020[casegroup=="A" & 
                    as.integer(specimen_date) >= as.Date("2021-07-01") &
                    as.integer(specimen_date) <= as.Date("2021-08-31")])

################################################################################

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

coeffs.timewindow.hosp <- coeffs.timewindow.hosp[!is.na(coef)]
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
                 limits=c(calendardate.lowerlimit, lastdate)) +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.8, 0.85)) 

###############################################################
## Figure: severe cases by date of presentation
winsize.casedates <- 3

casedates.allgr <-  cc.sincedec2020[CASE==1 & casegroup=="A", .N, by=specimen_date]
setkey(casedates.allgr, specimen_date)
casedates.allgr <- casedates.allgr[, rollmean.dtN(dt=.SD, k=winsize.casedates,
                                                  datevar="specimen_date"),
                                   .SDcols=c("specimen_date", "N")]

## use by= to get rolling means by category
casedates.byelig <- cc.sincedec2020[CASE==1 & casegroup=="A", .N, by=.(specimen_date, listedgr3)]
setkey(casedates.byelig, specimen_date)
casedates.byelig <- casedates.byelig[, rollmean.dtN(dt=.SD, k=winsize.casedates, datevar="specimen_date"), by=listedgr3,
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
breakpoints <- c(0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
p.rateratio.sincedose <-
    ggplot(data=coeffs.timewindow.sincedose, aes(x=time.midpoint, y=coef, color=effect)) +
    geom_line(size=coeffs.timewindow.sincedose$linesize) +
    labs(x=paste0("Weeks since last dose: mid-point of ", winsize.sincedose, "-day window"),
         y="Rate ratio for severe COVID-19 (log scale)") +
    scale_x_continuous(limits=c(2, 30)) + 
    scale_y_continuous(breaks=log(breakpoints),
                       labels=breakpoints,
                       limits=log(c(min(breakpoints), max(breakpoints))),
                       expand=c(0, 0)) +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.7, 0.85)) 

################ weeks since dose for hospitalisation ###############
firstmidpoint <- floor(winsize.sincedose/2)
lastmidpoint <- max(cc.sincedec2020$dayssincelastdose, na.rm=TRUE) -
    floor(winsize.sincedose/2)

coeffs.timewindow.sincedose.hosp <- NULL
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
        coeffs.timewindow.sincedose.hosp <- rbind(coeffs.timewindow.sincedose.hosp, coeffs)
    }
}

coeffs.timewindow.sincedose.hosp <- na.omit(coeffs.timewindow.sincedose.hosp, cols="coef")
coeffs.timewindow.sincedose.hosp <- coeffs.timewindow.sincedose.hosp[coef > -5 & coef < 5]
coeffs.timewindow.sincedose.hosp <- coeffs.timewindow.sincedose.hosp[grep("vaxgr", effect)]
coeffs.timewindow.sincedose.hosp[, effect := gsub("vaxgr", "", effect)]
coeffs.timewindow.sincedose.hosp[, effect := replace.names(effect)]
coeffs.timewindow.sincedose.hosp[, time.midpoint := time.midpoint / 7]
summary(coeffs.timewindow.sincedose.hosp)

coeffs.timewindow.sincedose.hosp[, linesize := `se(coef)`^-2]
coeffs.timewindow.sincedose.hosp[, linesize := 2 * linesize / max(linesize, na.rm=TRUE), by="effect"]
setorder(coeffs.timewindow.sincedose.hosp, effect)
# Figure: rate ratios by date of presentations
breakpoints <- c(0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)
p.rateratio.sincedose.hosp <-
    ggplot(data=coeffs.timewindow.sincedose.hosp, aes(x=time.midpoint, y=coef, color=effect)) +
    geom_line(size=coeffs.timewindow.sincedose.hosp$linesize) +
    labs(x=paste0("Weeks since last dose: mid-point of ", winsize.sincedose, "-day window"),
         y="Rate ratio for hospitalised/fatal COVID-19 (log scale)") + 
    scale_x_continuous(limits=c(2, 30)) + 
    scale_y_continuous(breaks=log(breakpoints),
                       labels=breakpoints,
                       limits=log(c(min(breakpoints), max(breakpoints))),
                       expand=c(0, 0)) +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.7, 0.85)) 

table.vaxproduct2 <- with(cc.sincedec2020, table(vacc_product_name_2, vax14.dose))
table.vaxproduct2.colsums <- colSums(table.vaxproduct2)
table.vaxproduct2 <- paste.colpercent(table.vaxproduct2, digits=1)
table.vaxproduct2 <- rbind(table.vaxproduct2, table.vaxproduct2.colsums)

#### efficacy against hospitalisation by risk category  #########################

coeffs.byriskgr.timewindow.sincedose <- NULL
for(t in firstmidpoint:lastmidpoint) {
    tdata <- cc.sincedec2020[vax14.dose == 0 | # include unvax as ref category
                           (dayssincelastdose >= t - floor(winsize.sincedose/2) & 
                            dayssincelastdose <= t + floor(winsize.sincedose/2))]
    time.midpoint <- t
    coeffs.byriskgr <- tryCatch(summary(clogit(formula=CASE ~ care.home + numdrugs.notcv +
                                                   inpat.recent + listedgr3 / vaxgr + #vax14.factor +
                                                   strata(stratum),
                                      data=tdata))$coefficients,
                       error=function(cond) return(NULL)
                       )
    if(!is.null(coeffs.byriskgr)) {
        coeffs.byriskgr <- data.table(time.midpoint=rep(time.midpoint, nrow(coeffs.byriskgr)),
                             effect=rownames(coeffs.byriskgr), coeffs.byriskgr)
        coeffs.byriskgr.timewindow.sincedose <- rbind(coeffs.byriskgr.timewindow.sincedose, coeffs.byriskgr)
    }
}
coeffs.byriskgr.timewindow.sincedose <- na.omit(coeffs.byriskgr.timewindow.sincedose, cols="coef")

table(coeffs.byriskgr.timewindow.sincedose$effect)

summary(coeffs.byriskgr.timewindow.sincedose[grep("vaxgr2", effect)])

######################################################
coeffs.byriskgr.timewindow.sincedose <- coeffs.byriskgr.timewindow.sincedose[grep("vaxgr2", effect)]
coeffs.byriskgr.timewindow.sincedose[, effect := gsub("vaxgr", " ", effect)]
coeffs.byriskgr.timewindow.sincedose[, effect := gsub("listedgr3", "", effect)]


coeffs.byriskgr.timewindow.sincedose <- coeffs.byriskgr.timewindow.sincedose[coef > -5 & coef < 5]
coeffs.byriskgr.timewindow.sincedose[, time.midpoint := time.midpoint / 7]
coeffs.byriskgr.timewindow.sincedose[, linesize := `se(coef)`^-2]
coeffs.byriskgr.timewindow.sincedose[, linesize := 2 * linesize / max(linesize, na.rm=TRUE), by="effect"]

coeffs.byriskgr.timewindow.sincedose[, effect := replace.names(effect)]
summary(coeffs.byriskgr.timewindow.sincedose)

setorder(coeffs.byriskgr.timewindow.sincedose, effect)
# Figure: rate ratios by date of presentations
breakpoints <- c(0.02, 0.05, 0.1, 0.2, 0.5, 1)
p.byriskgr.rateratio.sincedose.hosp.AZ <-
    ggplot(data=coeffs.byriskgr.timewindow.sincedose[grep("AZ", effect)],
           aes(x=time.midpoint, y=coef, color=effect)) +
    geom_line(size=coeffs.byriskgr.timewindow.sincedose[grep("AZ", effect), linesize]) +
    labs(x=paste0("Weeks since second dose: mid-point of ", winsize.sincedose, "-day window"),
         y="Rate ratio for hospitalised/fatal COVID-19 (log scale)") + 
   scale_x_continuous(limits=c(2, 30)) + 
    scale_y_continuous(breaks=log(breakpoints),
                       labels=breakpoints,
                       limits=log(c(min(breakpoints), max(breakpoints))),
                       expand=c(0, 0)) +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.6, 0.2)) 
p.byriskgr.rateratio.sincedose.hosp.mRNA <-
    ggplot(data=coeffs.byriskgr.timewindow.sincedose[grep("mRNA", effect)],
           aes(x=time.midpoint, y=coef, color=effect)) +
    geom_line(size=coeffs.byriskgr.timewindow.sincedose[grep("mRNA", effect), linesize]) +
    labs(x=paste0("Weeks since second dose: mid-point of ", winsize.sincedose, "-day window"),
         y="Rate ratio for hospitalised/fatal COVID-19 (log scale)") + 
   scale_x_continuous(limits=c(2, 30)) + 
    scale_y_continuous(breaks=log(breakpoints),
                       labels=breakpoints,
                       limits=log(c(min(breakpoints), max(breakpoints))),
                       expand=c(0, 0)) +
    theme(legend.title = element_blank()) +
    theme(legend.position = c(0.6, 0.9)) 

##### model for decay of rate ratio to time since 2nd dose  ########

source("decay.R")

#####################################################################

#### comparison of rate ratio for 2 doses versus 0 before and after 19 May Delta changepoint

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

recent.coeffs <- coeffs.timewindow[date.midpoint==as.Date("2021-07-29") &
                                   grepl("2 doses", effect), .(effect, coef, `se(coef)`)]
az.efficacy.severe <- pctefficacy.ci.text(recent.coeffs[1, 2], recent.coeffs[1, 3])
mrna.efficacy.severe <- pctefficacy.ci.text(recent.coeffs[2, 2], recent.coeffs[2, 3])

recent.coeffs.hosp <- coeffs.timewindow.hosp[date.midpoint==as.Date("2021-07-29") &
                                   grepl("2 doses", effect), .(effect, coef, `se(coef)`)]
az.efficacy.hosp <- pctefficacy.ci.text(recent.coeffs.hosp[1, 2], recent.coeffs.hosp[1, 3])
mrna.efficacy.hosp <- pctefficacy.ci.text(recent.coeffs.hosp[2, 2], recent.coeffs.hosp[2, 3])


objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))


coeffs.timewindow[effect=="1 dose mRNA vaccine" & 
                  date.midpoint >= as.Date("2021-05-01") &
                  date.midpoint <= as.Date("2021-05-31"),
                  ]

rmarkdown::render("vaxtrend.Rmd")
