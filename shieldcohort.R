library(data.table)
library(ggplot2)
library(survival)

source("helperfunctions.R")

datadir <- "data/2021-09-02/"
lastdate <- as.Date("2021-08-19")
interval.length <- 28
    
load(paste0(datadir, "shielded.linked.RData"))

shielded.full[, anycase := as.integer(!is.na(specimen_date))]
shielded.full <- shielded.full[!is.na(sex)] # 30 records with sex coded as empty string
shielded.full[, sex := factor(sex, levels=c("Male", "Female"))]

min.time <- min(as.integer(shielded.full$entrydate), na.rm=TRUE)
max.time <- max(as.integer(shielded.full$exitdate), na.rm=TRUE)

shielded.2vax <- shielded.full[!is.na(vaxdate_2)]
shielded.2vax[, entrydate := vaxdate_2 + 14]
shielded.2vax <- shielded.2vax[is.na(specimen_date) | specimen_date > entrydate]

## casegr has 5 categories
paste.rowpercent(with(shielded.2vax, table(shield.group, casegr)), digits=2)

## casegr2 has 2 categories
severe.2vax <- with(shielded.2vax, table(shield.group, casegr2))
severe.2vax <- data.table(shield.group=rownames(severe.2vax),
                               severe=severe.2vax[, 3])
setkey(severe.2vax, shield.group)

personmonths.2vax <- shielded.2vax[, as.numeric(sum(exitdate - entrydate)) * 12 / 365.25, by=shield.group]
names(personmonths.2vax)[2] <- "personmonths"
setkey(personmonths.2vax, shield.group)

severe.2vax <- severe.2vax[personmonths.2vax]

sumcols <- c("severe", "personmonths")
severe.2vax.sumcols <-
    severe.2vax[shield.group=="Additional conditions" |
                shield.group=="Rare diseases",
                lapply(.SD, sum),
                .SDcols = sumcols]
print(severe.2vax.sumcols)

severe.2vax[shield.group=="Additional conditions", severe := severe.2vax.sumcols[1, 1]]
severe.2vax[shield.group=="Additional conditions", personmonths := severe.2vax.sumcols[1, 2]]
severe.2vax <- severe.2vax[shield.group != "Rare diseases"]
severe.2vax[, shield.group := gsub("Additional conditions",
                                   "Additional / rare conditions", shield.group)]

severe.2vax[, rateper1000 := round(1000 * severe / personmonths, 2)]
setorder(severe.2vax, -rateper1000)
save(severe.2vax, file=paste0(datadir, "severe2vax.RData"))
     
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
                                offset(log(tobs)) + splines::ns(tstart, 6) +
                                age_years + sex + shield.group,
                           family = poisson,
                           data = shielded.tsplit)
anycase.summary <- summary(anycase.poissonmodel)
## predict on scale of log lambda and subtract offset term
shielded.tsplit[, anycase.loglambda := as.numeric(predict(anycase.poissonmodel,
                                                       newdata=shielded.tsplit,
                                                      type="link")) - log(tobs)]

## severe case #################################################
severe.poissonmodel <- glm(formula = severe ~
                                offset(log(tobs)) + splines::ns(tstart, 5) +
                                age_years + sex + shield.group,
                           family = poisson,
                           data = shielded.tsplit)
severe.summary <- summary(severe.poissonmodel)
shielded.tsplit[, severe.loglambda := as.numeric(predict(severe.poissonmodel,
                                                       newdata=shielded.tsplit,
                                                      type="link")) - log(tobs)]

#######################################################################

shielded.tsplit[, anycase.probmonth := 1 - exp(-exp(anycase.loglambda) * 365.25 / 12)]
shielded.tsplit[, severe.probmonth := 1 - exp(-exp(severe.loglambda) * 365.25 / 12)]
shielded.tsplit[, date := as.Date(tstart, origin=as.Date("1970-01-01"))]

shield.anycase  <- shielded.tsplit[, .(probmonth=mean(anycase.probmonth, na.rm=TRUE)), by=date]
shield.severe <- shielded.tsplit[, .(probmonth=mean(severe.probmonth, na.rm=TRUE)), by=date]

rm(shielded.tsplit)
gc()

shield.all <- rbind(shield.anycase[, casegr := "All cases"],
                    shield.severe[, casegr := "Severe cases"])

save(anycase.summary, severe.summary, shield.all,
     file=paste0(datadir, "shieldcohort.models.RData"))

