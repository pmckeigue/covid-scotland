## save cases only table for diabetes paper

casesonly <- as.data.table(readRDS("./data/2020-09-09/all_cases.rds"), key="ANON_ID")
                                        # latest specimen date 29 Sep 20
                                        # latest death date 31 July 20
sicsag.new <- as.data.table(readRDS("./data/2020-09-09/CC_SICSAG_ANON.rds"), key="ANON_ID") # latest hospital date 6 Sep 2020
casesonly$DATE_OF_DEATH <- as.Date(casesonly$DATE_OF_DEATH)
casesonly <- casesonly[SPECIMENDATE <= as.Date("2020-08-25") &
                       (is.na(DATE_OF_DEATH) |
                        DATE_OF_DEATH <= as.Date("2020-07-31")), ]
sicsag.new <- sicsag.new[AdmitUnit <= as.Date("2020-08-23"), ]
casesonly <- sicsag.new[casesonly]
setnames(casesonly, "AgeYear", "AGE")
setnames(casesonly, "SEX", "sex")

severe.cases <- casesonly[(!is.na(covidICUorHDU) & (covidICUorHDU == 1 | covidICUorHDU ==3)) |
                          dead28 == 1 |
                          covid_cod == 1, ]
fatal.cases <- casesonly[dead28 == 1 | covid_cod == 1, ]

case.freqs <- with(severe.cases[, .(AGE, sex, SPECIMENDATE)],
                   table(AGE, sex, exclude=NULL))
death.freqs <- with(fatal.cases[, .(AGE, sex, SPECIMENDATE)], table(AGE, sex, exclude=NULL))
severe.cases <- severe.cases[,  .(AGE, sex, SPECIMENDATE)]
fatal.cases <- fatal.cases[,  .(AGE, sex, SPECIMENDATE)]

save(severe.cases, fatal.cases, case.freqs, death.freqs,
     file="casefreqs31July.RData") ## national cumulative cases and deaths by sex and one year age group
