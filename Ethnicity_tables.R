###--------------------------------------------------------------------------------------------------
### PART A - remove NRS deaths from severe tables and re-run
###--------------------------------------------------------------------------------------------------

cc.all$casegroup <- car::recode(cc.all$casegroup,
                                 "'A'='Critical care or fatal'; 'B'='Hospitalised, not severe'; 'C'='Test-positive, not hospitalised';
                                 'D'='Critical care or fatal'")

## tabulate ONOMAP ethnicity against SMR ethnicity
table.ethnic <- table(cc.all$ethnic5.onomap, cc.all$ethnic5.smr, exclude=NULL)

tn <- table(cc.all$ethnic5.onomap, cc.all$ethnic5.smr)
SouthAsian.sensitivity <- 100 * tn[5, 2] / sum(tn[, 2])
SouthAsian.specificity <- 100 * (sum(tn[, -2]) - sum(tn[5, ]) + tn[5, 2]) / sum(tn[, -2])
sum.xtabulate <- sum(tn)

table.ethnic <- paste.colpercent(table.ethnic)


testpositives.ethnic.smr <- paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                  table(ethnic5.smr, casegroup)), 1)
table.testpositives.demog.ethnicsmr <-
        tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home", "SIMD.quintile"),
                           data=cc.all[cc.all$nrs_covid_case!=1&!is.na(cc.all$ethnic5.smr), ])


table.hospitalized.demog.ethnicsmr <-
  tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home", "SIMD.quintile"),
                             data=cc.all[(cc.all$casegroup=="Hospitalised, not severe" |
                                            cc.all$casegroup=="Critical care or fatal") & cc.all$nrs_covid_case!=1 &
                                           !is.na(cc.all$ethnic5.smr), ])

testpositives.ethnic <- paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                              table(ethnic4.onomap, casegroup)), 1)


testpositives.carehome <- paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                table(ethnic4.onomap, care.home)), 0)


testpositives.healthboard <- t(paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                     table(ethnic4.onomap, hb2019name)), 0))

table.testpositives.demog <-
  tabulate.freqs.regressions(varnames=c("ethnic4.onomap", "care.home", "SIMD.quintile"),
                             data=cc.all[cc.all$nrs_covid_case!=1,])

table.hospitalized.demog <-
  tabulate.freqs.regressions(varnames=c("ethnic4.onomap", "care.home", "SIMD.quintile"),
                             data=cc.all[(cc.all$casegroup=="Hospitalised, not severe" |
                                            cc.all$casegroup=="Critical care or fatal") &cc.all$nrs_covid_case!=1, ])


table.ethnic5smr <-
  tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home",
                                        "SIMD.quintile"),
                             outcome="CASE",
                             data=cc.severe[cc.severe$nrs_covid_case!=1&!is.na(cc.severe$ethnic5.smr), ])
rownames(table.ethnic5smr) <- replace.names(rownames(table.ethnic5smr))


table.severe.demog <-
  tabulate.freqs.regressions(varnames=c("ethnic3.onomap", "care.home",
                                        "SIMD.quintile"),
                             data=cc.severe[cc.severe$nrs_covid_case!=1,])

##with diabetes
table.ethnic5smrDB <-
  tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home",
                                        "SIMD.quintile","diabetes.any"),
                             outcome="CASE",
                             data=cc.severe[cc.severe$nrs_covid_case!=1&!is.na(cc.severe$ethnic5.smr), ])
rownames(table.ethnic5smrDB) <- replace.names(rownames(table.ethnic5smrDB))

##with diabetes
table.severe.demogDB <-
  tabulate.freqs.regressions(varnames=c("ethnic3.onomap", "care.home",
                                        "SIMD.quintile","diabetes.any"),
                             data=cc.severe[cc.severe$nrs_covid_case!=1,])




table.ethnic5smrNRS <-
  tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home",
                                        "SIMD.quintile"),
                             outcome="CASE",
                             data=cc.severe[!is.na(cc.severe$ethnic5.smr), ])
rownames(table.ethnic5smrNRS) <- replace.names(rownames(table.ethnic5smrNRS))


table.severe.demogNRS <-
  tabulate.freqs.regressions(varnames=c("ethnic3.onomap", "care.home",
                                        "SIMD.quintile"),
                             data=cc.severe)

#with diabetes
table.ethnic5smrNRSDB <-
  tabulate.freqs.regressions(varnames=c("ethnic5.smr", "care.home",
                                        "SIMD.quintile","diabetes.any"),
                             outcome="CASE",
                             data=cc.severe[!is.na(cc.severe$ethnic5.smr), ])
rownames(table.ethnic5smrNRSDB) <- replace.names(rownames(table.ethnic5smrNRSDB))

#with diabetes
table.severe.demogNRSDB <-
  tabulate.freqs.regressions(varnames=c("ethnic3.onomap", "care.home",
                                        "SIMD.quintile","diabetes.any"),
                             data=cc.severe)

###--------------------------------------------------------------------------------------------------
###--------------------------------------------------------------------------------------------------
### PART B - remove NRS deaths from severe tables and re-run using disaggregated ethnicity
###--------------------------------------------------------------------------------------------------
###--------------------------------------------------------------------------------------------------

  ## tabulate ONOMAP ethnicity against SMR ethnicity
  table.ethnicB <- table(cc.all$ethnic8.onomap, cc.all$ethnic9.smr, exclude=NULL)
  
  #tn <- table(cc.all$ethnic5, cc.all$ethnic9.smr)
  #SouthAsian.sensitivity <- 100 * tn[5, 2] / sum(tn[, 2])
  #SouthAsian.specificity <- 100 * (sum(tn[, -2]) - sum(tn[5, ]) + tn[5, 2]) / sum(tn[, -2])
  #sum.xtabulate <- sum(tn)
  
  table.ethnicB <- paste.colpercent(table.ethnicB)
  

testpositives.ethnic.smrB <- paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                  table(ethnic9.smr, casegroup)), 1)

table.testpositives.demog.ethnicsmrB <-
  tabulate.freqs.regressions(varnames=c("ethnic9.smr", "care.home", "SIMD.quintile"),
                             data=cc.all[cc.all$nrs_covid_case!=1&!is.na(cc.all$ethnic9.smr), ])


table.hospitalized.demog.ethnicsmrB <-
  tabulate.freqs.regressions(varnames=c("ethnic9.smr", "care.home", "SIMD.quintile"),
                             data=cc.all[(cc.all$casegroup=="Hospitalised, not severe" |
                                            cc.all$casegroup=="Critical care or fatal") & cc.all$nrs_covid_case!=1 &
                                           !is.na(cc.all$ethnic9.smr), ])

  ## tabulate ethnicity by case group
  testpositives.ethnicB <- paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                table(ethnic8.onomap, casegroup)), 1)
  
  testpositives.carehomeB <- paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                  table(ethnic8.onomap, care.home)), 0)
  
  testpositives.healthboardB <- t(paste.colpercent(with(cc.all[cc.all$CASE==1&cc.all$nrs_covid_case!=1, ],
                                                       table(ethnic8.onomap, hb2019name)), 0))
  
  table.testpositives.demogB <-
    tabulate.freqs.regressions(varnames=c("ethnic8.onomap", "care.home", "SIMD.quintile"),
                               data=cc.all[cc.all$nrs_covid_case!=1,])
  
  table.hospitalized.demogB <-
    tabulate.freqs.regressions(varnames=c("ethnic8.onomap", "care.home", "SIMD.quintile"),
                               data=cc.all[(cc.all$casegroup=="Hospitalised, not severe" |
                                              cc.all$casegroup=="Critical care or fatal") &cc.all$nrs_covid_case!=1, ])

  table.hospitalized.demogB7 <-
    tabulate.freqs.regressions(varnames=c("ethnic7.onomap", "care.home", "SIMD.quintile"),
                               data=cc.all[(cc.all$casegroup=="Hospitalised, not severe" |
                                              cc.all$casegroup=="Critical care or fatal") &cc.all$nrs_covid_case!=1, ])
  
table.ethnic9smrB <-
  tabulate.freqs.regressions(varnames=c("ethnic9.smr", "care.home",
                                        "SIMD.quintile"),
                             outcome="CASE",
                             data=cc.severe[cc.severe$nrs_covid_case!=1&!is.na(cc.severe$ethnic9.smr), ])
rownames(table.ethnic9smrB) <- replace.names(rownames(table.ethnic9smrB))


table.severe.demogB7 <-
  tabulate.freqs.regressions(varnames=c("ethnic7.onomap", "care.home",
                                        "SIMD.quintile"),
                             data=cc.severe[cc.severe$nrs_covid_case!=1,])



table.ethnic9smrNRSB <-
  tabulate.freqs.regressions(varnames=c("ethnic9.smr", "care.home",
                                        "SIMD.quintile"),
                             outcome="CASE",
                             data=cc.severe[!is.na(cc.severe$ethnic9.smr), ])
rownames(table.ethnic9smrNRSB) <- replace.names(rownames(table.ethnic9smrNRSB))


table.severe.demogNRSB7 <-
  tabulate.freqs.regressions(varnames=c("ethnic7.onomap", "care.home",
                                        "SIMD.quintile"),
                             data=cc.severe)

###Before saving environment remove data don't need - list below ones to keep

rm(list=ls() [!(ls() %in% c("cc.all", "cc.severe","z","zz","zzz",
                            "table.ethnic","table.ethnicB",
                            "table.testpositives.demog.ethnicsmr","table.testpositives.demog.ethnicsmrB",
                            "table.hospitalized.demog.ethnicsmr","table.hospitalized.demog.ethnicsmr",
                            "table.ethnic5smr","table.ethnic5smrB",
                            "table.ethnic5.smrNRS","table.ethnic5.smrNRS",
                            "testpositives.ethnic.smr","testpositives.ethnic.smrB",
                            "table.testpositives.demog","table.testpositives.demogB",
                            "table.hospitalized.demog","table.hospitalized.demogB",
                            "table.severe.demog","table.severe.demogB7",
                            "table.severe.demogNRS","table.severe.demogNRSB7",
                            "testpositives.ethnic","testpositives.ethnicB",
                            "testpositives.carehome","testpositives.carehomeB",
                            "testpositives.healthboard","testpositives.healthboardB"))])


