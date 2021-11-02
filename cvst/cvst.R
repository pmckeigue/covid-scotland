library(data.table)
library(ggplot2)

source("../helperfunctions.R")

credint95.ratepermillion <- function(events, popatrisk) {
    N <- length(events)
    popatrisk <- as.numeric(unlist(popatrisk))
    ci.str <- NULL
    for(i in 1:N) {
        ci <- round(qgamma(p=c(0.025, 0.975),
                           shape=events[i], scale=1E6 / popatrisk[i]), 1)
        
        ci.str <- c(ci.str, paste0("(", ci[1], " to ", ci[2], ")"))
    }
    ci.str <- gsub("\\(0 to 0\\)", ".", ci.str)
    return(ci.str)
}

credint.rateratio <- function(r1, r0, N1, N0) {
    width=0.001
    theta <- seq(-5, 5, by=width)
    odds <- exp(theta) *  N1 / N0
    p <- odds / (1 + odds)
    lik <- p^r1 * (1 - p)^r0
    lik <- lik  / (width * sum(lik))
    posterior <- data.table(theta, density=lik, rateratio=exp(theta))
    posterior[, cumprob := cumsum(lik) / sum(lik)] 
    ci <- with(posterior,
               c(theta[which.min(abs(cumprob - 0.025))],
                 theta[which.min(abs(cumprob - 0.975))]))
    ci <- round(exp(ci), 1)
    ci <- paste(sprintf("%.1f", ci[1]), "to",
                      sprintf("%.1f", ci[2]))
}

mle.rateratio <- function(r1, r0, N1, N0) {
    width=0.001
    theta <- seq(-5, 5, by=width)
    odds <- exp(theta) *  N1 / N0
    p <- odds / (1 + odds)
    lik <- p^r1 * (1 - p)^r0
    lik <- lik  / (width * sum(lik))
    posterior <- data.table(theta, density=lik, rateratio=exp(theta))
    posterior[, cumprob := cumsum(lik) / sum(lik)] 
    mle <- posterior$theta[which.max(lik)]
    mle <- round(exp(mle), 1)           
}
    
rollsum.datewin <- function(dt, k, datevar, align="center") {
    ## returns a table of rolling sums of width k centred on date
    ## dt is a dataset of totals (N) by date, which may have no rows for some dates in the
    ## range of interest
    ## add rows for missing dates by a left join of all.dates with dt
    all.dates <- data.table(date=seq(from=min(dt[[datevar]]),
                                     to=max(dt[[datevar]]),
                                     by=1))
    setkey(all.dates, date)
    # setkeyv(dt, datevar) can't set physical key here
    all.dates <- dt[all.dates]
    all.dates[is.na(N), N := 0]

    return(data.table(date=all.dates[[datevar]],
                      rsum=zoo::rollsum(all.dates[, N], k=k, fill=0, align=align)))
}

getvaxpop.rsum <- function(vax.bydate) {
    ## calculate daily totals of those vaccinated in last 14 days
    ## for each vaccine
    vaxpop.rsum.az <- rollsum.datewin(vax.bydate[vacc_product_name=="AstraZeneca"],
                                      k=14, datevar="vaxdate",
                                      align="right")
    vaxpop.rsum.az[, product := "AstraZeneca"]
    vaxpop.rsum.pf <- rollsum.datewin(vax.bydate[vacc_product_name=="Pfizer"],
                                   k=14, datevar="vaxdate",
                                   align="right")
    vaxpop.rsum.pf[, product := "Pfizer"]
    vaxpop.rsum.mo <- rollsum.datewin(vax.bydate[vacc_product_name=="Moderna"],
                                      k=14, datevar="vaxdate",
                                   align="right")
    vaxpop.rsum.mo[, product := "Moderna"]
    
    vaxpop.rsum <- rbind(vaxpop.rsum.pf, vaxpop.rsum.az, vaxpop.rsum.mo)
    vaxpop.rsum[, rsum := as.numeric(rsum)]
    setnames(vaxpop.rsum, "rsum", "sum14")

    ## for daily totals of those vaccinated in last 15-28 days, increment date by 14
    vaxpop28 <- copy(vaxpop.rsum)
    vaxpop28[, date := date + 14]
    setnames(vaxpop28, "sum14", "sum28")
    setkey(vaxpop28, product, date)

    ## for daily totals of those vaccinated in last 29-42 days, increment date by 28
    vaxpop42 <- copy(vaxpop.rsum)
    vaxpop42[, date := date + 28]
    setnames(vaxpop42, "sum14", "sum42")
    setkey(vaxpop42, product, date)

    ## for daily totals of those vaccinated in last 43-56 days, increment date by 42
    vaxpop56 <- copy(vaxpop.rsum)
    vaxpop56[, date := date + 42]
    setnames(vaxpop56, "sum14", "sum56")
    setkey(vaxpop56, product, date)

    ## for daily totals of those vaccinated in last 57-70 days, increment date by 56
    vaxpop70 <- copy(vaxpop.rsum)
    vaxpop70[, date := date + 56]
    setnames(vaxpop70, "sum14", "sum70")
    setkey(vaxpop70, product, date)

    ## for daily totals of those vaccinated in last 71-84 days, increment date by 70
    vaxpop84 <- copy(vaxpop.rsum)
    vaxpop84[, date := date + 70]
    setnames(vaxpop84, "sum14", "sum84")
    setkey(vaxpop84, product, date)

    setkey(vaxpop.rsum, product, date)
    vaxpop.rsum <- vaxpop28[vaxpop.rsum]
    vaxpop.rsum <- vaxpop42[vaxpop.rsum]
    vaxpop.rsum <- vaxpop56[vaxpop.rsum]
    vaxpop.rsum <- vaxpop70[vaxpop.rsum]
    vaxpop.rsum <- vaxpop84[vaxpop.rsum]
    vaxpop.rsum[, sum28 := nafill(sum28, fill=0)]
    vaxpop.rsum[, sum42 := nafill(sum42, fill=0)]
    vaxpop.rsum[, sum56 := nafill(sum56, fill=0)]
    vaxpop.rsum[, sum70 := nafill(sum70, fill=0)]
    vaxpop.rsum[, sum84 := nafill(sum84, fill=0)]
    
    vaxpop.rsum[date > as.Date("2021-04-14") & date <= as.Date("2021-05-17"),
                lapply(.SD, function(x) 0.917 * x),
                .SDcols=c("sum14", "sum28", "sum42", "sum56", "sum70", "sum84")]
    vaxpop.rsum[date > as.Date("2021-05-17"),
                 lapply(.SD, function(x) 0.121 * x),
                .SDcols=c("sum14", "sum28", "sum42", "sum56", "sum70", "sum84")]
    return(vaxpop.rsum)
}    

## include code G08, exclude deaths or admissions with any mention of infection
cvst.matched <- function(x) as.integer(grepl("I636|I676|G08X?", x))
secondary.matched <- function(x) as.integer(grepl("^A|^B|^C7[012]|^D3[23]|G06[02]", x))

## read vaccine dataset
vaxbatch <- as.data.table(readRDS("data/vaccine_2021_07_19.rds")) 
vaxbatch <- unique(vaxbatch) # 6979166 records no CHI
with(vaxbatch, table(toupper(substr(vacc_batch_number, 1, 2)), vacc_product_name)[, -1])
with(vaxbatch[toupper(substr(vacc_batch_number, 1, 2))=="UN"],
     table(vacc_batch_number, vacc_product_name))
## Pfizer vaccines have batchnum beginning with E or F
## AZ vaccines have batchnum beginning with 4, A, K or P
## unknown batchnums coded as UNKNOWN or unknown
vaxbatch[toupper(substr(vacc_batch_number, 1, 2)) == "UN",
         vacc_batch_number := NA]
vaxbatch[toupper(substr(vacc_batch_number, 1, 1)) %in% c("E", "F"),
         productcorrect := "Covid-19 mRNA Vaccine Pfizer"]
vaxbatch[toupper(substr(vacc_batch_number, 1, 1)) %in% c("4", "A", "K", "P"),
         productcorrect := "Covid-19 Vaccine AstraZeneca"]
setkey(vaxbatch, vacc_product_name, vacc_occurence_time, vacc_dose_number, patient_sex, patient_date_of_birth)
vaxbatch[, count := .N, by=key(vaxbatch)]
vaxbatch <- vaxbatch[count==1] # 3204471 records uniquely identified by indirect matching
vaxbatch <- vaxbatch[!is.na(productcorrect) & productcorrect != vacc_product_name] # 1236 records

vaxpop <- as.data.table(readRDS("data/vaccine_data_2021_07_19.rds")) # 6930252 records with CHI
vaxpop[, source_system_patient_id := NULL]
vaxpop[, patient_derived_post_code := NULL]
vaxpop <- unique(vaxpop)
setkey(vaxpop, vacc_product_name, vacc_occurence_time, vacc_dose_number, patient_sex, patient_date_of_birth)
vaxpop <- vaxbatch[vaxpop]
vaxpop[!is.na(productcorrect), vacc_product_name := productcorrect]

setnames(vaxpop, "patient_derived_upi_number", "chi")
vaxpop[nchar(chi)==9, chi := paste0("0", chi)]
setnames(vaxpop, "patient_sex", "sex")
setnames(vaxpop, "patient_date_of_birth", "birth_date")
vaxpop[, sex := car::recode(sex, "'MALE'='Male'; 'FEMALE'='Female'")]
vaxpop[, birth_date := as.Date(birth_date)]
vaxpop[, vaxdate := as.Date(vacc_occurence_time)] # , format="%m/%d/%Y")]
vaxpop[vaxdate < as.Date("2020-12-01"), vaxdate := NA]
vaxpop[as.integer(vaxdate - birth_date) > 365.25 * 110,  birth_date := NA]
vaxpop <- vaxpop[vacc_product_name != "Not Applicable"] # code means dose not given
vaxpop[, vacc_product_name := gsub("Covid-19 ", "", vacc_product_name)]
vaxpop[, vacc_product_name := gsub("Vaccine ", "", vacc_product_name)]
vaxpop[, vacc_product_name := gsub("mRNA ", "", vacc_product_name)]
vaxpop[, vacc_product_name := factor(vacc_product_name,
                                     levels=c("Pfizer", "AstraZeneca", "Moderna"))]
vaxpop[, age.vax := floor(as.integer(vaxdate - birth_date) / 365.25)]
vaxpop[, vacc_dose_number := as.integer(vacc_dose_number)]
vaxpop[, vacc_occurence_time := NULL]
setkey(vaxpop, chi)

vax.bydate <- vaxpop[!is.na(vacc_product_name) &
                     !is.na(vaxdate) & vaxdate <= as.Date("2021-06-01"), .N,
                       by=list(vaxdate, vacc_product_name)]
setkey(vax.bydate, vaxdate)
vaxpop.rsum <- getvaxpop.rsum(vax.bydate)
pop.atrisk <- data.table(
    vaxpop.rsum[, sum(sum14) / 14, by=product],
    vaxpop.rsum[, sum(sum28) / 14, by=product][, 2],
    vaxpop.rsum[, sum(sum42 + sum56 + sum70 + sum84) / 14, by=product][, 2])
colnames(pop.atrisk)[2:4] <- c("1-14", "15-28", "29-84")
pop.atrisk <- pop.atrisk[c(3, 1, 2), ]
p28days <- as.integer(0.5 * sum(as.numeric(pop.atrisk[2, 2:3])))

vaxunder60.bydate <- vaxpop[!is.na(vacc_product_name) & !is.na(vaxdate) &
                            vaxdate <= as.Date("2021-06-01") & 
                            !is.na(age.vax) & age.vax < 60,
                            .N,
                            by=list(vaxdate, vacc_product_name)]
setkey(vaxunder60.bydate, vaxdate)

vaxpopunder60.rsum <- getvaxpop.rsum(vaxunder60.bydate)
#vaxpopunder60.rsum[, sum(sum14) / 14, by=product]

popunder60.atrisk <- data.table(
    vaxpopunder60.rsum[, sum(sum14) / 14, by=product],
    vaxpopunder60.rsum[, sum(sum28) / 14, by=product][, 2],
    vaxpopunder60.rsum[, sum(sum42 + sum56 + sum70 + sum84) / 14, by=product][, 2])
colnames(popunder60.atrisk)[2:4] <- c("1-14", "15-28", "29-84")
popunder60.atrisk <- popunder60.atrisk[c(3, 1, 2), ]
p28days.under60 <- as.integer(0.5 * sum(as.numeric(popunder60.atrisk[2, 2:3])))


p.vaxdates <- ggplot(data=vaxpop.rsum, aes(x=date, y=sum14, color=product)) + geom_line() +
    xlab("End date of 14-day window after vaccine dose") +
    ylab("Vaccinated in last 14 days") +
    scale_x_date(limits=as.Date(c("2020-12-07", "2021-05-01"))) + 
    scale_y_continuous(breaks=c(0, 2E5, 4E5), labels=c("0", "0.2", "0.4"))

## deaths
deaths.all <- RDStodt("./2021_05_19/all_deaths_from_Dec_2020.rds")
deaths.all[, DATE_OF_DEATH := as.Date(DATE_OF_DEATH, format="%d/%m/%Y")]
## code mention of CVST without infection
cols.cause <- grep("CAUSE", names(deaths.all), value=TRUE)
cvst.any <- rowSums(deaths.all[, lapply(.SD, cvst.matched), .SDcols=(cols.cause)]) > 1
cvst.secondary <- rowSums(deaths.all[, lapply(.SD, secondary.matched),
                                     .SDcols=(cols.cause)]) > 1
cvst.mention <- deaths.all[cvst.any > 0 & cvst.secondary==0]
## of 7 deaths with mention of CVST and no secondary cause, six had code G08 and one I676

deaths15 <- RDStodt("./2021_05_19/NRS_DEATHS_20210610.rds") # 15 obs
#setnames(deaths, "UPI_NUMBER", "CHI")
deaths15[, DATE_OF_DEATH := as.Date(DATE_OF_DEATH)]

## SMR01
smr01 <- RDStodt("2021_05_19/SMR01_2021-05-26.rds") # 6153 records since 2015
setnames(smr01, "UPI_NUMBER", "CHI")

smr01.cases <- RDStodt("./2021_05_19/SMR01_20210610.rds")
common_cols <- intersect(colnames(smr01), colnames(smr01.cases))
smr01 <- unique(rbind(smr01[, ..common_cols], smr01.cases[, ..common_cols]))
smr01[, SEX := car::recode(as.character(SEX), "'1'='Male'; '2'='Female'; '9'=NA",
                           as.factor=TRUE)]
smr01[, ADMISSION_DATE := as.Date(ADMISSION_DATE)]
smr01[, DISCHARGE_DATE := as.Date(DISCHARGE_DATE)]
setnames(smr01, "AGE_IN_YEARS", "age.seen")

## code mention of CVST without infection
cols.condition <- grep("CONDITION", names(smr01), value=TRUE) 
cols.cvst <- paste("cvst", cols.condition, sep=".")
cols.secondary <- paste("secondary", cols.condition, sep=".")

conditions <- smr01[, ..cols.condition]
sums.cvst.any <- rowSums(sapply(conditions, cvst.matched))
sums.secondary.any <- rowSums(sapply(conditions, secondary.matched))

disch <- data.table(smr01, cvst.any=sums.cvst.any, secondary.any=sums.secondary.any)
cvstany.disch <- disch[cvst.any > 0]
## select first record with a CVST diagnosis
setorder(cvstany.disch, ADMISSION_DATE)
cvstany.disch <- unique(cvstany.disch, by="CHI")

print(with(cvstany.disch[ADMISSION_DATE >= as.Date("2020-12-01")],
     table(CHI %in% disch[cvst.any > 0 & secondary.any > 0]$CHI)))

cvst.disch <- cvstany.disch[!(CHI %in% disch[cvst.any > 0 & secondary.any > 0]$CHI),
                            .(CHI, MAIN_CONDITION, SEX, age.seen, ADMISSION_DATE,
                              HBRES_CURRENTDATE)]
print(with(cvst.disch[ADMISSION_DATE >= as.Date("2020-12-01")],
     table(CHI %in% disch[cvst.any > 0 & secondary.any > 0]$CHI)))

setkey(cvst.disch, CHI)

## read unfiltered PACS scans
scans.pacs <- RDStodt("2021_05_19/scans3.rds")[STUDY_DATE >= as.Date("2020-12-01") &
                                               STUDY_DATE <= as.Date("2021-05-15")]
setkey(scans.pacs, STUDY_DATE)

scans.pacs.veno <- scans.pacs[venogram.code==TRUE, .N, by=STUDY_DATE]
setkey(scans.pacs.veno, STUDY_DATE)
scans.pacs.veno <- rollsum.datewin(dt=scans.pacs.veno[, list(N, STUDY_DATE)], k=7,
                                   datevar="STUDY_DATE")
scans.pacs.veno[, venogram := "Venogram"]

scans.pacs.nonveno <- scans.pacs[venogram.code==FALSE, .N, by=STUDY_DATE]
setkey(scans.pacs.nonveno, STUDY_DATE)
scans.pacs.nonveno <- rollsum.datewin(dt=scans.pacs.nonveno[, list(N, STUDY_DATE)], k=7,
                                   datevar="STUDY_DATE")
scans.pacs.nonveno[, venogram := "Other"]

scans.pacs.date <- rbind(scans.pacs.veno, scans.pacs.nonveno)

p.veno <- ggplot(data=scans.pacs.veno, aes(x=date, y=rsum)) +
    geom_line() +
    ylab("Number of venograms") + 
    scale_x_date(limits=as.Date(c("2020-12-07", "2021-05-01")))

p.nonveno <- ggplot(data=scans.pacs.nonveno, aes(x=date, y=rsum)) +
    geom_line() +
    ylab("Other head scans") + 
    scale_x_date(limits=as.Date(c("2020-12-07", "2021-05-01")))  +
    scale_y_continuous(limits=c(4000, 8000))
    
## read scans reviewed
scans.events <- RDStodt("2021_05_19/scans.all.events.rds") # 1665 events


## generate highest event code for each individual
scans.events[, assign.integer := as.integer(hfeventcode)]
scans.events[, min.assign := min(assign.integer), by=c("chi", "sex", "birth_date")]
setorder(scans.events, min.assign)

## retain distinct events only if they have different dates
scans.events <- unique(scans.events, by=c("hfeventonset", "study_date", "chi", "sex", "birth_date"))

source("cvst_withchi1.R")

scans.all <- unique(scans.events, by=c("chi", "sex", "birth_date"))
## some records still duplicated on chi
scans.all[anyDuplicated(scans.all[!is.na(chi), chi]), ] # one duplicated nonmissing chi
setkey(scans.all, chi)

scans.possible <- scans.all[study_date >= as.Date("2020-12-01") &
                            (hfeventcode=="primary acute" |
                             hfeventcode=="possible")]
scans.possible[, age.scan := floor(as.integer(study_date - birth_date) / 365.25)]
setorder(scans.possible, hfeventcode, na.last=TRUE)
setkey(scans.possible, chi)

## left join cvst.disch[ADMISSION_DATE >= as.Date("2020-12-01")] (29 rows) with scans.all on CHI

scans.join <- scans.all[!is.na(chi) & !duplicated(chi) & study_date >= as.Date("2020-12-01")]

cvst.leftjoin <-
    scans.join[, .(chi, site_name, study_date,
                hfeventonset,
                hfeventcode)][cvst.disch[ADMISSION_DATE >= as.Date("2020-12-01")]]
cvst.leftjoin <- unique(cvst.leftjoin) # 29 records
setorder(cvst.leftjoin, ADMISSION_DATE)

chi.scans.missing <- cvst.leftjoin[is.na(hfeventcode), chi]
scans.found <- scans.pacs[CHI %in% chi.scans.missing]
setorder(scans.found, CHI, STUDY_DATE)

source("cvst_withchi2.R")


cvst.leftjoin.keep <- cvst.leftjoin[is.na(hfeventcode) |
                                    !(hfeventcode=="chronic" | hfeventcode=="secondary" |
                                      hfeventcode=="negative")]
cvst.leftjoin.keep[, study_date := NULL]
cvst.leftjoin.keep <- unique(cvst.leftjoin.keep)

## outer join of cvst.leftjoin.keep with scans.possible
## cvst.leftjoin.keep has 22 records: chi, SEX, age.seen, hfeventonset, hfeventcode
## no missing chi
## scans.possible has 68 records (49 primary acute, 19 possible), one with missing chi

cvst.all <- merge(cvst.leftjoin.keep[ADMISSION_DATE >= as.Date("2020-12-01"),
                                     .(chi, ADMISSION_DATE, age.seen, SEX,
                                       MAIN_CONDITION, scancheck)],
                  scans.possible[!duplicated(chi)],
                  by.x="chi", by.y="chi", all=TRUE)
print(with(cvst.all, table(is.na(SEX), is.na(sex))))
cvst.all[is.na(SEX), SEX := sex] # 3 extra records with missing eventcode: now 71 records
cvst.all[, sex := NULL]
cvst.all[, dateseen := pmin(study_date, ADMISSION_DATE, hfeventonset, na.rm=TRUE)] 
cvst.all <- cvst.all[dateseen >= as.Date("2020-12-01")]
with(cvst.all, table(is.na(ADMISSION_DATE), is.na(study_date)))
cvst.all[, report_text := NULL]

table(scans.all[chi %in% deaths15$CHI, hfeventcode])
## 9 of 15 fatal cases in scans.all were scored as primary acute, 2 as possible and 4 as secondary

setorder(smr01, ADMISSION_DATE)
admissions.other <- smr01[CHI %in% cvst.all[is.na(ADMISSION_DATE) &
                                            hfeventcode=="primary acute", chi] &
                          ADMISSION_DATE >= as.Date("2020-12-01"), .(CHI, ADMISSION_DATE, MAIN_CONDITION)]
setkey(admissions.other, CHI)
setkey(scans.all, chi)
admissions.other.leftjoin <- scans.all[admissions.other]
admissions.other.leftjoin[, timediff := abs(as.integer(study_date - ADMISSION_DATE))]
setorder(admissions.other.leftjoin, timediff)
admissions.other.leftjoin <- unique(admissions.other.leftjoin, by="chi")
admissions.other.leftjoin <- admissions.other.leftjoin[timediff < 8,
                                                       .(chi, ADMISSION_DATE, MAIN_CONDITION)]
setnames(admissions.other.leftjoin, "ADMISSION_DATE", "admissiondate.other")
setnames(admissions.other.leftjoin, "MAIN_CONDITION", "diagnosis.other")
setkey(admissions.other.leftjoin, chi)

# left join cvst.all with vax
## 70 records in cvst.all: 51 primary acute, 18 possible 1 NA because no scan 
## one missing CHI reportdbid 37497977 male born 1968 scan date 2021-04-12 onset ? 12 days before
setkey(cvst.all, chi)
cvst.nochi <- cvst.all[is.na(chi)]
cvst.all <- vaxpop[cvst.all[!is.na(chi)]]
cvst.all <- rbind(cvst.all, cvst.nochi, fill=TRUE)

cvst.all[, daysfromvax := as.integer(dateseen - vaxdate)]
cvst.all[daysfromvax < 0 | daysfromvax > 1000, daysfromvax := NA]
cvst.all[is.na(daysfromvax), vacc_product_name := NA]
## select most recent vaccination date i.e. min(daysfromvax) excluding negative values
setorder(cvst.all, assign.integer, daysfromvax, na.last=TRUE)
cvst.all <- unique(cvst.all, by="chi")
table(cvst.all$hfeventcode, exclude=NULL) # 67 records: 48 primary acute, 18 possible, 1 NA

cvst.all[is.na(age.seen), age.seen := age.scan]
cvst.all[is.na(SEX), SEX := i.sex]
setorder(cvst.all, "dateseen")

rmarkdown::render("cvst_medrxiv.Rmd")
