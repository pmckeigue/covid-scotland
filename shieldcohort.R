library(data.table)
library(ggplot2)
library(survival)

datadir <- "data/2021-07-28/"

load(paste0(datadir, "shielded.full.RData"))

min.time <- min(as.integer(shielded.full$entrydate), na.rm=TRUE)
max.time <- max(as.integer(shielded.full$exitdate), na.rm=TRUE)

interval.length <- 14

shielded.tsplit <-
			as.data.table(survSplit(formula=Surv(time=as.integer(entrydate),
					  time2=as.integer(exitdate),	
                			       	   event=severecase) ~ age_years + sex + shield.group,
			     	    data=shielded.full,
						   		cut=seq(min.time, max.time, by=interval.length)))	


shielded.tsplit[, tobs := tstop - tstart]
		  
severe.poissonmodel <- glm(formula = severecase ~
                                offset(log(tobs)) + splines::ns(tstart, 4) +
                                age_years + sex + shield.group,
                           family = poisson,
                           data = shielded.tsplit)

## predict on scale of intensity of Poisson arrivals per interval
shielded.tsplit[, severe.lambda := as.numeric(predict(severe.poissonmodel,
                                                       newdata=shielded.tsplit,
                                                       type="response"))]
shielded.tsplit <- shielded.tsplit[!is.na(severe.lambda)]
shielded.tsplit[, severe.prob := 1 - exp(-lambda * 365.25 / interval.length)]
shielded.tsplit[, date := as.Date(tstart, origin=as.Date("1970-01-01"))]

shield.meanrate  <- shielded.tsplit[, .(riskyear=mean(severe.prob)), by=date]
setorder(shield.meanrate, date)

ggplot(data=shield.meanrate, aes(x=date, y=riskyear)) +
    geom_line()
