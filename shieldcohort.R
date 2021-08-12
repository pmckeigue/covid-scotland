library(data.table)
library(ggplot2)
library(survival)

datadir <- "data/2021-07-28/"
lastdate <- as.Date("2021-07-14")
interval.length <- 28
    
load(paste0(datadir, "shielded.full.RData"))

shielded.full[, anycase := as.integer(!is.na(specimen_date))]
min.time <- min(as.integer(shielded.full$entrydate), na.rm=TRUE)
max.time <- max(as.integer(shielded.full$exitdate), na.rm=TRUE)

shielded.tsplit <-
    as.data.table(survSplit(formula=Surv(time=as.integer(entrydate),
                                         time2=as.integer(exitdate),	
                                         event=anycase) ~ age_years + sex + shield.group + casegr + casegr2,
                            data=shielded.full,
                            cut=seq(min.time, max.time, by=interval.length)))	
nrow(shielded.tsplit)

shielded.tsplit[, tobs := tstop - tstart]
shielded.tsplit[, severe := as.integer(anycase==1 & casegr2=="Severe")]

## any case ##################################
anycase.poissonmodel <- glm(formula = anycase ~
                                offset(log(tobs)) + splines::ns(tstart, 5) +
                                age_years + sex + shield.group,
                           family = poisson,
                           data = shielded.tsplit)
anycase.summary <- summary(anycase.poissonmodel)
## predict on scale of intensity of Poisson arrivals per interval
shielded.tsplit[, anycase.lambda := as.numeric(predict(anycase.poissonmodel,
                                                       newdata=shielded.tsplit,
                                                      type="response"))]

## severe case #################################################
severe.poissonmodel <- glm(formula = severe ~
                                offset(log(tobs)) + splines::ns(tstart, 5) +
                                age_years + sex + shield.group,
                           family = poisson,
                           data = shielded.tsplit)
severe.summary <- summary(severe.poissonmodel)
shielded.tsplit[, severe.lambda := as.numeric(predict(severe.poissonmodel,
                                                       newdata=shielded.tsplit,
                                                      type="response"))]

#######################################################################

shielded.tsplit[, anycase.prob := 1 - exp(-anycase.lambda * 365.25 / interval.length)]
shielded.tsplit[, severe.prob := 1 - exp(-severe.lambda * 365.25 / interval.length)]
shielded.tsplit[, date := as.Date(tstart, origin=as.Date("1970-01-01"))]

shield.anycase  <- shielded.tsplit[, .(riskyear=mean(anycase.prob, na.rm=TRUE)), by=date]
shield.severe <- shielded.tsplit[, .(riskyear=mean(severe.prob, na.rm=TRUE)), by=date]

save(anycase.summary, severe.summary, shield.anycase, shield.severe,
     file=paste0(datadir, "shieldcohort.models.RData"))
rm(shielded.tsplit)
gc()

########### plots ##################################
p.anycase <- ggplot(data=shield.anycase, aes(x=date, y=riskyear)) +
    geom_line() +
        labs(x=paste0("Presentation date: start of ", interval.length, "-day interval"),
             y="Rate per year") + 
    scale_x_date(breaks = seq.Date(from = as.Date("2020-03-01"),
                                   to = lastdate, by = "month"),
                 expand=c(0, 10), 
                 labels=gsub("^0", "", 
                             format.Date(seq.Date(from = as.Date("2020-03-01"),
                                                  to = lastdate, by = "month"),
                                         "%d %b")
                             ),
                 limits=c(as.Date("2020-12-01"), lastdate)) 
p.anycase

p.severe <- ggplot(data=shield.severe, aes(x=date, y=riskyear)) +
    geom_line()
p.severe

png("shieldedrisk.png")
gridExtra::grid.arrange(p.anycase, p.severe, nrow=2)
dev.off()
