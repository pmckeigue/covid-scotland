firstdate <- as.Date("2020-03-01")
## for coeffs we have to specify explicitly the midpoint of the time window 
winsize <- 42 
startdates <- with(cc.all, firstdate:max(specimen_date) - winsize)
enddates <- startdates + winsize
lastdate <- as.Date("2021-07-14")

coeffs.timewindow <- NULL
for(timewin in 1:length(startdates)) {
    tdata <- cc.all[casegroup=="A" & 
                       as.integer(specimen_date) >= startdates[timewin] &
                       as.integer(specimen_date) <= enddates[timewin],
                       .(CASE, care.home, listedgr3, vax14.dose, vax14.factor, occup, inpat.recent, stratum)]
    if(length(with(tdata, table(CASE, occup))) > 2) {
        date.midpoint <- startdates[timewin] + floor(winsize / 2)
        coeffs <- tryCatch(summary(clogit(formula=CASE ~ care.home + listedgr3 + inpat.recent + vax14.dose + strata(stratum), data=tdata))$coefficients,
                           error=function(cond) return(NULL)
                           )
        if(!is.null(coeffs)) {
            coeffs <- data.table(date.midpoint=rep(date.midpoint, nrow(coeffs)),
                                 effect=rownames(coeffs), coeffs)
            coeffs.timewindow <- rbind(coeffs.timewindow, coeffs)
        }
    }
}
coeffs.timewindow[, date.midpoint := as.Date(date.midpoint, origin="1970-01-01")]
coeffs.timewindow <- coeffs.timewindow[grep("care.home", effect, invert=TRUE)]
coeffs.timewindow[, effect := gsub("^care.home", "", effect)]
coeffs.timewindow[, effect := gsub("^listedgr3", "", effect)]
coeffs.timewindow[, effect := gsub("TRUE$", "", effect)]
coeffs.timewindow[, effect := replace.names(effect)]
coeffs.timewindow[, linesize := `se(coef)`^-1.5]
coeffs.timewindow[, linesize := 2 * linesize / max(linesize, na.rm=TRUE), by="effect"]
setorder(coeffs.timewindow, effect)

breakpoints <- c(0.2, 0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50)
p.rateratio <-
    ggplot(data=coeffs.timewindow, aes(x=date.midpoint, y=coef, color=effect)) +
        geom_line(size=coeffs.timewindow$linesize) +
        labs(x=paste0("Presentation date: mid-point of ", winsize, "-day window"),
             y="Rate ratio (log scale)") + 
        scale_y_continuous(breaks=log(breakpoints),
                           labels=breakpoints,
                           limits=log(c(min(breakpoints), max(breakpoints))),
                           expand=c(0, 0)) +
        scale_x_date(breaks = seq.Date(from = as.Date("2020-12-01"),
                                       to = lastdate, by = "month"),
                     expand=c(0, 10), 
                     labels=gsub("^0", "", 
                                 format.Date(seq.Date(from = as.Date("2020-12-01"),
                                                      to = lastdate, by = "month"),
                                             "%d %b")
                                 ),
                     limits=c(as.Date("2020-12-01"), lastdate)) +
        theme(legend.position = c(0.2, 0.2)) 
 p.rateratio

png("rateratiotimeplot.png")
p.rateratio
dev.off()

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

## Figure 1: severe cases by date of presentation
winsize.casedates <- 3

casedates.allgr <-  cc.all[CASE==1 & casegroup=="A", .N, by=specimen_date]
setkey(casedates.allgr, specimen_date)
casedates.allgr <- casedates.allgr[, rollmean.dtN(dt=.SD, k=winsize.casedates),
                                   .SDcols=c("specimen_date", "N")]

## use by= to get rolling means by category
casedates.byelig <- cc.all[CASE==1 & casegroup=="A", .N, by=.(specimen_date, listedgr3)]
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
    theme(legend.title = element_blank()) +
    scale_fill_manual(values=c("grey", "blue", "red")) +
    xlab(paste0("Presentation date: mid-point of ", winsize.casedates, "-day window")) +
         ylab("Daily severe cases")
p.byelig

png("dailysevere.png")
p.byelig
dev.off()
