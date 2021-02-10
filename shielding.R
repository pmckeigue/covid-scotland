
## rollmean.dtN returns a data.table with columns date, avdaily (rolling mean) 
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

## tables and figures for shielding report

theme_set(theme_gray(base_size = 14))

lastdate <- as.Date("2020-12-31")

## Table 1 - shielding cohort by eligibility category x age

shielded.full[, shield.group6 := car::recode(shield.group,
                                             recode="'Pregnant with heart disease'=NA",
                                             levels=levels(shield.group)[-6])]
                                             
bygroup <- with(shielded.full, table(shield.group6))
bygroup <- as.integer(c(bygroup, sum(bygroup)))

shieldgroups.byage <- with(shielded.full,
                           as.matrix(as.data.frame.matrix(table(shield.group6, agegr20))))
shieldgroups.byage <- rbind(shieldgroups.byage, colSums(shieldgroups.byage))
rownames(shieldgroups.byage)[nrow(shieldgroups.byage)] <- "All shielding categories"

shieldgroups.byage.rowpct <- round(100 * shieldgroups.byage /
                                   rowSums(shieldgroups.byage))
shieldgroups.byage.rowpct <- paste0("(", shieldgroups.byage.rowpct, "\\%)")
table.shieldgroups.byage <- paste.vectortomatrix(shieldgroups.byage,
                                                 shieldgroups.byage.rowpct)
table.shieldgroups.byage <- data.frame(table.shieldgroups.byage,
                                       All=bygroup)
colnames(table.shieldgroups.byage) <- c(levels(shielded.full$agegr20), "All")


## Table 2 - shielding cohort by shielding category x outcome
table.shield.casegr <- univariate.tabulate(outcome="casegr3",
                                           varnames=c("shield.group6"),
                                           data=shielded.full, drop.reflevel=FALSE,
                                           digits=1, colpercent=FALSE)

table.shield.casegr.rowsums <- rowSums(with(shielded.full, table(shield.group6, casegr3)))

table.shield.casegr.rowsums <- c(sum(table.shield.casegr.rowsums),
                                     table.shield.casegr.rowsums)

## add row for all shielding categories with row percentages
table.casegr <- with(shielded.full, table(casegr3))
table.casegr <- paste0(table.casegr, " (",
                       round(100 * table.casegr / sum(table.casegr), 1),
                       "\\%)")
table.shield.casegr <- rbind(table.casegr, table.shield.casegr)
rownames(table.shield.casegr) <- replace.names(rownames(table.shield.casegr))
rownames(table.shield.casegr)[1] <- "All shielding categories"
colnames(table.shield.casegr) <- gsub("\\(.+\\)", "", colnames(table.shield.casegr))
table.shield.casegr <- data.frame(table.shield.casegr, All=table.shield.casegr.rowsums)
colnames(table.shield.casegr) <- gsub("\\.+", " ", colnames(table.shield.casegr))
colnames(table.shield.casegr) <- gsub(" non ", ", non-", colnames(table.shield.casegr))


###################################

## Table 3 - shielding cohort by date of letter
table.date.letter <- as.matrix(as.data.frame.matrix(
    with(shielded.full, table(shield.group, Date.Sent))))
col5name <- colnames(table.date.letter)[5]
table.date.letter <- rbind(table.date.letter, colSums(table.date.letter))
rownames(table.date.letter)[nrow(table.date.letter)] <- "All shielding categories"
                           
totals <-  rowSums(table.date.letter)
lastcol <- ncol(table.date.letter)

table.date.letter <- cbind(table.date.letter[, 1:4], rowSums(table.date.letter[, 5:lastcol]))
date.letter.rowpct <- round(100 * table.date.letter / totals)
date.letter.rowpct <- paste0("(", date.letter.rowpct, "\\%)")

table.date.letter <- paste.vectortomatrix(table.date.letter, date.letter.rowpct)
table.date.letter <- cbind(table.date.letter, totals)
colnames(table.date.letter)[5:6] <- c(col5name, "All")

table.date.letter <- as.data.frame(table.date.letter)
colnames(table.date.letter)[1:5] <- format(as.Date(colnames(table.date.letter)[1:5]),
                                           "%d %b")
colnames(table.date.letter) <- gsub("^0", "", colnames(table.date.letter))
colnames(table.date.letter)[5] <- paste(colnames(table.date.letter)[5], "or later")

#####################################
## Table 4: rate ratios asssociated with eligibility category by advice interval

table.shielding.byinterval <-
    cbind(
        tabulate.freqs.regressions(varnames=c("shield.group", "shield.any"),
                                   outcome="CASE",
                                   data=cc.severe[care.home=="Independent" &
                                               SPECIMENDATE <= as.Date("2020-04-17")])[, 1:4], 
        tabulate.freqs.regressions(varnames=c("shield.group", "shield.any"), 
                                   outcome="CASE",
                                   data=cc.severe[care.home=="Independent" &
                                               SPECIMENDATE > as.Date("2020-04-17")])[, 1:4]) 
rownames(table.shielding.byinterval)[nrow(table.shielding.byinterval)] <- "All shielding categories"

with(cc.all[shield.group != "No shielding"], table(hh.over18gr, agegr20))

## calculate estimated proportion of all cases infected after arrival of each batch of letters
severe.infected.date <- table(cc.severe[care.home=="Independent"][CASE==1, SPECIMENDATE - 7])
dateby.props <- cumsum(severe.infected.date) / sum(severe.infected.date)
dateby.props <- dateby.props[match(levels(as.factor(cc.all$Date.Sent)), names(dateby.props))]


#####################################

with(cc.all[shield.group != "No shielding"], table(hh.over18gr, agegr20))

## calculate estimated proportion of all cases infected after arrival of each batch of letters
severe.infected.date <- table(cc.severe[care.home=="Independent"][CASE==1, SPECIMENDATE - 7])
dateby.props <- cumsum(severe.infected.date) / sum(severe.infected.date)
dateby.props <- dateby.props[match(levels(as.factor(cc.all$Date.Sent)), names(dateby.props))]

## table by infection after letter
table.after.letter <- univariate.tabulate(varnames=c("shield.group", "adultsgt2",
                                                     "hosp.recent", "hb2019name"), 
                                            outcome="after.letter",
                                          data=cc.all[CASE==1 & care.home=="Independent" &
                                                      shield.group != "No shielding"],
                                          colpercent=FALSE)
rownames(table.after.letter) <- replace.names(rownames(table.after.letter))
colnames(table.after.letter) <- c("Letter sent at least days earlier",
                                  "Letter not sent at least 14 days earlier")

## same for severe cases
table.after.letter.severe <- univariate.tabulate(varnames=c("shield.group", "hh.over18gr",
                                                            "hosp.recent", "hb2019name"), 
                                            outcome="after.letter",
                                          data=cc.severe[CASE==1 & care.home=="Independent" &
                                                      shield.group != "No shielding"],
                                          colpercent=FALSE)
rownames(table.after.letter.severe) <- replace.names(rownames(table.after.letter.severe))
colnames(table.after.letter.severe) <-  c("Letter sent >= 14 days earlier",
                                          "Letter not sent at least 14 days earlier")

#########################################
## exclude care home residents

## FIXME: remove this later
#cc.severe[, occup := car::recode(var=occup,
#                                 recodes="'Health care PF'='Health care, not PF / undetermined'",
#                                 as.factor=TRUE,
#                                 levels=c("Other / undetermined",
#                                          "Teacher",
#                                          "Health care, not PF / undetermined"))]

## Table 5 - regression of severe case status on shielding group and covariates
table.shielded.severecases.nocare <-
    tabulate.freqs.regressions(varnames=c("shield.group", 
                                           "preschool.any", "hh.schoolagegr",
                                          "hh.over18gr", 
                                          "SIMD.quintile", "occup", "hosp.recent"),
                               outcome="CASE",
                               data=cc.severe[care.home=="Independent"])

## regression of severe case status on listedgr3 and covariates to get coeff for adultsgt1
table.listedgr3.severecases.nocare <-
    tabulate.freqs.regressions(varnames=c("listedgr3",
                                           "preschool.any", "hh.schoolagegr",
                                          "adultsgt1", 
                                          "SIMD.quintile", "occup", "hosp.recent"),
                               outcome="CASE",
                               data=cc.severe[care.home=="Independent"])

## regression of severe case status on listedgr3 and covariates to get coeff for hh.schoolage.any
table.listedgr3.severecases.nocare.anyschoolage <-
    tabulate.freqs.regressions(varnames=c("listedgr3",
                                           "preschool.any", "hh.schoolage.any",
                                          "hh.over18gr", 
                                          "SIMD.quintile", "occup", "hosp.recent"),
                               outcome="CASE",
                               data=cc.severe[care.home=="Independent"])


## Table 6 - regression in shielded eligible only
summary(cc.severe[listedgr3 == "Eligible for shielding" & care.home=="Independent",
                  .(shieldedonly.group, hh.over18gr, hosp.recent)])
table.shieldedonly.severecases.nocare <-
    tabulate.freqs.regressions(varnames=c("shieldedonly.group", 
                                           "preschool.any", "hh.schoolagegr",
                                          "hh.over18gr3", 
                                          "qSIMD.integer", "hosp.recent"),
                               outcome="CASE",
                               data=cc.severe[listedgr3 == "Eligible for shielding" &
                                              care.home=="Independent"])
####################################################


###########################################
## Figure 1: severe cases by date of presentation
winsize.casedates <- 3

casedates.allgr <-  cc.severe[CASE==1 & care.home=="Independent", .N, by=SPECIMENDATE]
setkey(casedates.allgr, SPECIMENDATE)
casedates.allgr <- casedates.allgr[, rollmean.dtN(dt=.SD, k=winsize.casedates),
                                   .SDcols=c("SPECIMENDATE", "N")]

## can use by= to get rolling means by category
casedates.byelig <- cc.severe[CASE==1 & care.home=="Independent", .N, by=.(SPECIMENDATE, listedgr3)]
setkey(casedates.byelig, SPECIMENDATE)
casedates.byelig <- casedates.byelig[, rollmean.dtN(dt=.SD, k=winsize.casedates), by=listedgr3,
                       .SDcols=c("SPECIMENDATE", "N")]
                                 
p.byelig <- ggplot(data=casedates.byelig,
                aes(x=date,
                    y=avdaily, fill=listedgr3)) +
    geom_area() +
    scale_y_continuous(expand=c(0, 0)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 0), 
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = as.Date("2020-03-01"),
                                                  to = lastdate, by = "month"),
                                         "%d %b")
                 ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
    theme(legend.position = c(0.5, 0.7)) +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values=c("grey", "blue", "red")) +
    xlab(paste0("Presentation date: mid-point of ", winsize.casedates, "-day window")) +
         ylab("Daily severe cases")
p.byelig

#######################################################################
## Figure 2: plot rate ratios associated with shielding eligibility, recent hospital exposure and number of children

firstdate <- as.Date("2020-03-01")

## for coeffs we have to specify explicitly the midpoint of the time window 
winsize <- 21 # should be an odd number for date.midpoint to be exactly centred
startdates <- with(cc.severe, firstdate:max(SPECIMENDATE) - winsize)
enddates <- startdates + winsize

## 2 (a) plot rate ratio associated with shielding eligibility by time window

coeffs.timewindow <- NULL
for(timewin in 1:length(startdates)) {
    tdata <- cc.severe[care.home=="Independent" &
                       as.integer(SPECIMENDATE) >= startdates[timewin] &
                       as.integer(SPECIMENDATE) <= enddates[timewin],
                       .(CASE, listedgr3, stratum)]
    if(length(with(tdata, table(CASE, listedgr3))) > 2) {
        date.midpoint <- startdates[timewin] + floor(winsize / 2)
        coeffs <- c(date.midpoint,
                    tryCatch(summary(clogit(formula=CASE ~ listedgr3 + strata(stratum), data=tdata))$coefficients,
                             error=function(cond) return(matrix(rep(NA, 10), nrow=2))
                             )
                    )
        coeffs.timewindow <- rbind(coeffs.timewindow, coeffs)
    }
}

coeffs.timewindow <- data.table(coeffs.timewindow)
coeffs.long <- melt(coeffs.timewindow, idvars=1, measure=list(2:3, 4:5, 6:7, 8:9, 10:11),
                    value.name=c("coeff", "rateratio", "se.coeff", "Z", "pvalue"))
colnames(coeffs.long)[1:2] <- c("date.midpoint", "riskgroup")
coeffs.long[, date.midpoint := as.Date(date.midpoint, origin="1970-01-01")]
coeffs.long[, riskgroup := car::recode(riskgroup,
                                       "1='Moderate risk condition'; 2='Eligible for shielding';",
                                       as.factor=TRUE)]
 

coeffs.long[date.midpoint >= as.Date("2020-06-01") &
                  date.midpoint < as.Date("2020-09-01"), coeff := NA]

num.letters <- with(shielded.full, table(Date.Sent))
num.letters <- num.letters / num.letters[1]
dates.letters <- as.Date(names(num.letters))
coeffs.long[date.midpoint >= as.Date("2020-06-01") & date.midpoint <= as.Date("2020-09-30"), coeff := NA]
                                  
p.rateratio <-
    ggplot(data=coeffs.long, aes(x=date.midpoint, y=coeff, color=riskgroup)) +
    geom_line(size=0.01 * coeffs.long[, se.coeff]^-2) +
    scale_color_manual(name="Risk category",
                       values=c("red", "blue")) + 
    theme(legend.position = c(0.5, 0.5)) +
    theme(legend.title = element_text(size=10), legend.text=element_text(size=10)) + 
    scale_y_continuous(breaks=log(c(4, 5, 6, 7, 8, 10, 12, 14)),
                       labels=c(4, 5, 6, 7, 8, 10, 12, 14),
                       limits=log(c(3, 14)), expand=c(0, 0)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 0), 
                 labels=gsub("^0", "", 
                     format.Date(seq.Date(from = as.Date("2020-03-01"),
                              to = lastdate, by = "month"),
                     "%d %b")
                 ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
    xlab(paste0("Presentation date: mid-point of ", winsize, "-day window")) + ylab("Rate ratio (log scale)") +
    ggtitle("Association with risk category") +
    annotate(geom="segment",
             x = dates.letters,
             xend = dates.letters,
             y=log(3),
             yend=log(3.5),
             size=as.numeric(num.letters), 
             arrow=arrow(ends="first", length=unit(0.1, "inches")))
p.rateratio

#####################################

## 2 (b) plot rate ratio associated with recent hospital exposure by time window

coeffs.hosp.timewindow <- NULL
for(timewin in 1:length(startdates)) {
    tdata <- cc.severe[care.home=="Independent" &
                       as.integer(SPECIMENDATE) >= startdates[timewin] &
                       as.integer(SPECIMENDATE) <= enddates[timewin], .(CASE, hosp.recent,
                                                                        stratum)]
    if(
        length(with(tdata, table(CASE))) > 1 & 
        length(with(tdata, which(hosp.recent==TRUE))) > 0
    ) {
        date.midpoint <- startdates[timewin] + floor(winsize / 2)
        coeffs.hosp <- c(date.midpoint,
                         summary(clogit(formula=CASE ~ hosp.recent + strata(stratum),
                                        data=tdata))$coefficients)
        coeffs.hosp.timewindow <- rbind(coeffs.hosp.timewindow, coeffs.hosp)
    }
}

colnames(coeffs.hosp.timewindow) <- c("date.midpoint", "coeff", "rateratio", "se.coeff", "Z", "pvalue")

coeffs.hosp.timewindow <- data.table(coeffs.hosp.timewindow)
coeffs.hosp.timewindow[, date.midpoint := as.Date(date.midpoint, origin="1970-01-01")]
coeffs.hosp.timewindow <- coeffs.hosp.timewindow[se.coeff < 0.5]
coeffs.hosp.timewindow[date.midpoint >= as.Date("2020-06-01") &
                  date.midpoint <= as.Date("2020-09-30"), coeff := NA]

p.hosp.rateratio <-
    ggplot(data=coeffs.hosp.timewindow, aes(x=date.midpoint, y=coeff)) +
    geom_line(size=0.01 * coeffs.hosp.timewindow[, se.coeff]^-2) +
    scale_y_continuous(breaks=log(c(5, 10, 20, 30, 50)),
                           labels=c(5, 10, 20, 30, 50)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 0), 
                 labels=gsub("^0", "", 
                     format.Date(seq.Date(from = as.Date("2020-03-01"),
                              to = lastdate, by = "month"),
                     "%d %b")
                 ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
     xlab(paste0("Presentation date: mid-point of ", winsize, "-day window")) + ylab("Rate ratio (log scale)") +
    ggtitle("Association with recent hospital exposure")

#################################################################################
## 2 (c) plot rate ratio associated with number of children by time window

coeffs.household.timewindow <- NULL
for(timewin in 1:length(startdates)) {
    tdata <- cc.severe[care.home=="Independent" &
                       as.integer(SPECIMENDATE) >= startdates[timewin] &
                       as.integer(SPECIMENDATE) <= enddates[timewin], .(CASE, hh.schoolage, hh.over18,
                                                                        qSIMD.integer,
                                                                        stratum)]
    if(
        length(with(tdata, table(CASE))) > 1 & 
        length(with(tdata, which(hh.schoolage > 0))) > 0
    ) {
        date.midpoint <- startdates[timewin] + floor(winsize / 2)
        coeffs.household <-
            data.frame(date.midpoint=rep(date.midpoint, 2),
                       hhvar=c("per child", "per adult"),
                       summary(clogit(formula=CASE ~ hh.schoolage + hh.over18 + qSIMD.integer 
                                      + strata(stratum),
                                      data=tdata))$coefficients[1:2, ])
        coeffs.household.timewindow <- rbind(coeffs.household.timewindow, coeffs.household)
    }
}

coeffs.household.timewindow <- data.table(coeffs.household.timewindow)

colnames(coeffs.household.timewindow) <- 
    c("date.midpoint", "RR", "coeff", "rateratio", "se.coeff", "Z", "pvalue")
coeffs.household.timewindow[, date.midpoint := as.Date(date.midpoint, origin="1970-01-01")]
coeffs.household.timewindow <- coeffs.household.timewindow[se.coeff < 0.5]
coeffs.household.timewindow[date.midpoint >= as.Date("2020-06-01") &
                  date.midpoint < as.Date("2020-09-01"), coeff := NA]

coeffs.adults.timewindow <- coeffs.household.timewindow[RR=="adults"]
coeffs.children.timewindow <- coeffs.household.timewindow[RR=="children"]
coeffs.household.timewindow[date.midpoint >= as.Date("2020-06-01") & date.midpoint <= as.Date("2020-09-30"), coeff := NA]

p.household.rateratio <-
    ggplot(data=coeffs.household.timewindow, aes(x=date.midpoint, y=coeff, color=RR)) +
    geom_line(size=0.05 * coeffs.household.timewindow[, se.coeff]^-1) +
    xlim(as.Date(c("2020-03-01", "2020-11-30"))) +   
    scale_y_continuous(breaks=log(c(0.6, 0.8, 1, 1.2, 1.4, 1.6)),
                       labels=c(0.6, 0.8, 1, 1.2, 1.4, 1.6), 
                       limits=log(c(0.6, 1.8)), expand=c(0, 0)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 0), 
                 labels=gsub("^0", "", 
                     format.Date(seq.Date(from = as.Date("2020-03-01"),
                              to = lastdate, by = "month"),
                     "%d %b")
                 ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
    xlab(paste0("Presentation date: mid-point of ", winsize, "-day window")) + ylab("Rate ratio (log scale)") +
    labs(color='Rate ratio') + 
    theme(legend.position = c(0.5, 0.6)) + 
    geom_hline(size=0.2, yintercept=0) +
    ggtitle("Association with adults and children in household")

p.children.rateratio <-
    ggplot(data=coeffs.children.timewindow, aes(x=date.midpoint, y=coeff)) +
    geom_line(size=0.01 * coeffs.children.timewindow[, se.coeff]^-2) +
    scale_y_continuous(breaks=log(c(0.6, 0.8, 1, 1.2)),
                       limits=log(c(0.6, 1.2)),
                       labels=c(0.6, 0.8, 1, 1.2)) + 
    xlim(as.Date(c("2020-03-01", "2020-11-30"))) +   
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 0), 
                 labels=gsub("^0", "", 
                     format.Date(seq.Date(from = as.Date("2020-03-01"),
                              to = lastdate, by = "month"),
                     "%d %b")
                 ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
    xlab(paste0("Presentation date: mid-point of ", winsize, "-day window")) + ylab("Rate ratio per child (log scale)") +
    geom_hline(size=0.2, yintercept=0) +
    ggtitle("Association with school-age children in household")


#################################################################################

## Figure 3: frequency of recent hospital exposure in controls(excluding care home residents, population attributable risk fraction excluding care home residents, and fraction of cases in care home residents

## 3 (a) frequency of recent hospital exposure
winsize.hosp <- 7
controls.hosp <- cc.severe[CASE==0 & care.home=="Independent", .N, by=.(SPECIMENDATE, hosp.recent, listedgr3)]
setkey(controls.hosp, SPECIMENDATE)

controls.hosp <- controls.hosp[, rollmean.dtN(dt=.SD, k=winsize.hosp), by=.(listedgr3, hosp.recent), 
                               .SDcols=c("SPECIMENDATE", "N")]

## compute sum of avdaily over hosp.recent TRUE and FALSE in each window by listedgr3
controls.N <- controls.hosp[, .(avdaily.total = sum(avdaily)), by=.(date, listedgr3)]
setkey(controls.N, date, listedgr3)

controls.hosp <- controls.hosp[hosp.recent==TRUE]
setnames(controls.hosp, "avdaily", "avdaily.exposed")
setkey(controls.hosp, date, listedgr3)
controls.hosp <- controls.N[controls.hosp]
controls.hosp[, expfreq := avdaily.exposed / avdaily.total]
#controls.hosp[avdaily.total < 5, expfreq := NA]
controls.hosp[avdaily.total < 5, avdaily.total := NA]
controls.hosp[, se.expfreq := sqrt(expfreq * (1 - expfreq) / avdaily.total)]
controls.hosp[date >= as.Date("2020-06-01") & date <= as.Date("2020-09-30"), expfreq := NA]


p.hosp <- ggplot(data=controls.hosp,
                 aes(x=date, y=expfreq,  color=listedgr3, group=listedgr3)) +
    geom_line(size=0.07 / controls.hosp$se.expfreq^0.4) +  
    guides(size=FALSE) + 
    scale_color_manual(name="Risk category", values=c("black", "blue", "red"),
                       guide = guide_legend(reverse = TRUE)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 0), 
                 labels=gsub("^0", "", 
                     format.Date(seq.Date(from = as.Date("2020-03-01"),
                              to = lastdate, by = "month"),
                     "%d %b")
                    ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
    scale_y_continuous(breaks=seq(0, 0.2, by=0.02), limits=c(0, 0.2), expand=c(0, 0)) +
    theme(legend.position = c(0.5, 0.7)) + 
    theme(legend.title = element_blank()) +
    xlab(paste0("Presentation date: midpoint of ", winsize.hosp, "-day window")) +
    ylab("Proportion exposed") +
    ggtitle("Recent hospital exposure in controls not resident in care homes")
p.hosp

###########################################################

## 3 (b) fraction attributable to recent hospital exposure excluding care home residents

## compute daily counts of cases and controls with recent hospital exposure
cc.hosp.recent <- cc.severe[CASE==0 & care.home=="Independent" & hosp.recent==TRUE,
                        .N, by=.(SPECIMENDATE, CASE)]
setkeyv(cc.hosp.recent, c("SPECIMENDATE", "CASE"))

## compute daily counts of all cases and controls
cc.counts.all <- cc.severe[care.home=="Independent", .N, by=.(SPECIMENDATE, CASE)]
setkeyv(cc.counts.all, c("SPECIMENDATE", "CASE"))

#########################

## compute rolling mean of avdaily in controls with recent hospital exposure
controls.hosp.rollmean <- cc.hosp.recent[CASE==0, rollmean.dtN(dt=.SD, k=winsize),
                                  .SDcols=c("SPECIMENDATE", "N")]
setnames(controls.hosp.rollmean, "avdaily", "avdaily.exposed")
setkey(controls.hosp.rollmean, date)
## compute rolling mean of avdaily in all controls
controls.all.rollmean <- cc.counts.all[CASE==0, rollmean.dtN(dt=.SD, k=winsize),
                               .SDcols=c("SPECIMENDATE", "N")]
setnames(controls.all.rollmean, "avdaily", "avdaily.total")
setkey(controls.all.rollmean, date)
## left join av counts in all controls with av counts in hosp.recent
controls.all.rollmean <- controls.hosp.rollmean[controls.all.rollmean]
## compute frequency of exposure in controls
controls.all.rollmean[, expfreq := avdaily.exposed / avdaily.total]
setkey(controls.all.rollmean, date)

######################################

## compute rolling mean of avdaily in cases with recent hospital exposure
cases.hosp.rollmean <- cc.hosp.recent[CASE==0, rollmean.dtN(dt=.SD, k=winsize),
                                  .SDcols=c("SPECIMENDATE", "N")]
setnames(cases.hosp.rollmean, "avdaily", "avdaily.exposed")
setkey(cases.hosp.rollmean, date)
## compute rolling mean of avdaily in all cases
cases.all.rollmean <- cc.counts.all[CASE==1, rollmean.dtN(dt=.SD, k=winsize),
                               .SDcols=c("SPECIMENDATE", "N")]
setnames(cases.all.rollmean, "avdaily", "avdaily.total")
setkey(cases.all.rollmean, date)
## left join av counts in all cases with av counts in hosp.recent
cases.all.rollmean <- cases.hosp.rollmean[cases.all.rollmean]
## compute frequency of exposure in cases
cases.all.rollmean[, case.expfreq := avdaily.exposed / avdaily.total]
setnames(cases.all.rollmean, "avdaily.exposed", "avdaily.exposed.case")
setnames(cases.all.rollmean, "avdaily.total", "avdaily.total.case")
setkey(cases.all.rollmean, date)

######################################

## merge with coeffs
setkey(coeffs.hosp.timewindow, date.midpoint)
coeffs.hosp.timewindow <- controls.all.rollmean[coeffs.hosp.timewindow]
coeffs.hosp.timewindow <- cases.all.rollmean[coeffs.hosp.timewindow]

## classic definition
#coeffs.hosp.timewindow[, popattr.frac := expfreq * (rateratio - 1) /
                        #                        (1 + expfreq * (rateratio - 1))]
## Miettinen's definition
coeffs.hosp.timewindow[, popattr.frac := case.expfreq * (rateratio - 1) / rateratio]

coeffs.hosp.timewindow[date.midpoint >= as.Date("2020-06-01") &
                      date.midpoint <= as.Date("2020-09-30"), popattr.frac := NA]

r.hcai <- summary(clogit(formula=CASE ~ prob.hcai + strata(stratum),
                                 data=cc.severe))$coefficients[2]

p.hcai <- with(cc.severe[CASE==0 & care.home=="Independent"], table(prob.hcai))
p.hcai <- p.hcai[2] / sum(p.hcai)
hcai.popattrfrac <- p.hcai * (r.hcai - 1) / (1 + p.hcai * (r.hcai - 1))
coeffs.hosp.timewindow[date >= as.Date("2020-06-01") & date <= as.Date("2020-09-30"), popattr.frac := NA]

p.fracattr <-
    ggplot(data=coeffs.hosp.timewindow, aes(x=date, y=popattr.frac)) +
    geom_line(size=0.01 * coeffs.hosp.timewindow[, se.coeff]^-2) +
    scale_y_continuous(limits=c(0, 0.55), expand=c(0, 0)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 0), 
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = as.Date("2020-03-01"),
                                                  to = lastdate, by = "month"),
                                         "%d %b")
                             ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
    xlab(paste0("Presentation date: mid-point of ", winsize, "-day window")) +
    ylab("Population attributable risk fraction") +
    ggtitle("Fraction of severe cases attributable to hospital exposure")

## 3 (c) fraction of severe cases in care home residents

cases.all <- cc.severe[CASE==1, .N, by=.(SPECIMENDATE, care.home)]
setkey(cases.all, SPECIMENDATE)
cases.all <- cases.all[, rollmean.dtN(dt=.SD, k=winsize),
                                       by=care.home, 
                                       .SDcols=c("SPECIMENDATE", "N")]

## compute sum of avdaily over care.home TRUE and FALSE in each window
cases.all.N <- cases.all[, .(avdaily.total = sum(avdaily)), by=date]
setkey(cases.all.N, date)
cases.all <- cases.all[care.home=="Care/nursing home"]
setnames(cases.all, "avdaily", "avdaily.care.home")
setkey(cases.all, date)
cases.all <- cases.all.N[cases.all]
cases.all[, care.frac := avdaily.care.home / avdaily.total]
setkey(cases.all, date)
cases.all[date.midpoint >= as.Date("2020-06-01") &
                  date.midpoint < as.Date("2020-09-01"), care.frac := NA]

p.carefrac <-
    ggplot(data=cases.all, aes(x=date, y=care.frac)) +
    geom_line() + # size=0.01 * coeffs.children.timewindow[, se.coeff]^-2) +
    scale_y_continuous(limits=c(0, 0.75), expand=c(0, 0)) + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 labels=gsub("^0", "", 
                     format.Date(seq.Date(from = as.Date("2020-03-01"),
                              to = lastdate, by = "month"),
                     "%d %b")
                 ),
                 limits=c(as.Date("2020-03-01"), lastdate)) +
    xlab(paste0("Presentation date: mid-point of sliding window")) +
    ylab("Proportion of severe cases") +
    ggtitle("Fraction of severe cases resident in care homes")

#############################################################################
## Supplementary Figure 1: case fatality rate

casedates.all <- cc.all[CASE==1, .N, by=.(SPECIMENDATE, fatalcase, shield.any)]
## cast to wide format
casedates.all <- dcast(casedates.all, SPECIMENDATE ~ fatalcase + shield.any, value.var="N")
## set NA to 0
for (j in 2:5) {
    set(casedates.all, which(is.na(casedates.all[[j]])), j, 0)
}
## compute sliding window
winsize.fatality <- 7
for(j in names(casedates.all)[2:5]) {
    casedates.all[, (j) := zoo::rollmean(casedates.all[[j]], k=winsize.fatality, fill=NA)]
}

colnames(casedates.all)[2:5] <- c("nonfatal.inelig", "nonfatal.elig",
                                  "fatal.inelig", "fatal.elig")

## compute N and  casefatality
casedates.all[, N.inelig := fatal.inelig + nonfatal.inelig]
casedates.all[, N.elig := fatal.elig + nonfatal.elig]

casedates.all[, casefatality.inelig := fatal.inelig / N.inelig]
casedates.all[, casefatality.elig := fatal.elig / N.elig]

## concatenate inelig and elig
casedates.inelig <- casedates.all[, .(SPECIMENDATE, N.inelig, casefatality.inelig)]
casedates.elig <- casedates.all[, .(SPECIMENDATE, N.elig, casefatality.elig)]
colnames(casedates.inelig) <- colnames(casedates.elig) <- c("SPECIMENDATE", "N", "casefatality")
casedates.inelig[, shield.elig := FALSE]
casedates.elig[, shield.elig := TRUE]
casedates.all <- rbind(casedates.inelig, casedates.elig)
setkey(casedates.all, SPECIMENDATE)

p.casefatality <- ggplot(data=casedates.all,
       aes(x=SPECIMENDATE, y=casefatality, group=shield.elig, color=shield.elig)) +
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = as.Date("2020-05-01"), by = "week"),
                 labels=gsub("^0", "", 
                     format.Date(seq.Date(from = as.Date("2020-03-01"),
                                          to = as.Date("2020-05-01"),
                                          by = "week"),
                                 "%d %b")
                     ),
                 limits=c(as.Date("2020-03-01"), as.Date("2020-05-01"))) +
    scale_y_continuous(expand=c(0, 0)) + 
    geom_line(size=ifelse(casedates.all[["shield.elig"]], 0.4, 0.02) *
                  casedates.all[["N"]]^0.5 
              ) +
    scale_color_manual(name="Shielding status",
                       labels=c("Ineligible","Eligible"),
                       values=c("blue", "red")) +
    theme(legend.position = c(0.6, 0.8)) + 
    xlab(paste0("Presentation date: midpoint of ", winsize.fatality, "-day window")) +
    ylab("Proportion of cases that were fatal") 

############################################################################
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
                 limits=c(firstdate, lastdate)) +
    theme(legend.position = c(0.1, 0.8)) +
    theme(legend.title = element_blank()) +
    scale_fill_manual(values=c("grey", "orange", "green")) +
    xlab(paste0("Specimen date: mid-point of ", winsize.casedates, "-day window")) +
         ylab("Daily severe cases")



#############################################################

rmarkdown::render("shielding.Rmd")
rmarkdown::render("shieldingslides.Rmd")
