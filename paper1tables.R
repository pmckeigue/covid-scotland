## tables for first paper reporting case-control study

########### variable lists for tabulating

if(old) {
    demog <- c("ethnic3.onomap", "SIMD.quintile", "care.home")
}

demog.smr <- c("ethnic4.smr", "SIMD.quintile", "care.home")

varnames.extended <- c("care.home", "scrip.any", "diag.any", "listed.any", "diag.other",
                     "scripordiag", listed.conditions)

varnames.listedpluscovs <- c("care.home", "scrip.any", "diag.any", listed.conditions)

bnfcols <- grep("^BNF", colnames(cc.severe))
bnf.chapternames <- colnames(cc.severe)[bnfcols]

icdcols <- grep("^Ch\\.", colnames(cc.severe))
icd.chapternames <- colnames(cc.severe)[icdcols]
## drop factors without at least 2 levels represented in cc.severe 
num.values <- unlist(lapply(cc.severe[, ..icd.chapternames], function(x) length(table(x))))
icd.chapternames <- icd.chapternames[num.values > 1]

full.varnames <- c(demog.smr, listed.conditions,
                   "diag.any", icd.chapternames, "scrip.any", bnf.chapternames)
## drop factors without at least 2 levels represented in cc.severe 
num.values <- unlist(lapply(cc.severe[, ..full.varnames], function(x) length(table(x))))
full.varnames <- full.varnames[num.values > 1]

## logical vectors for subsetting
nocare <- cc.severe$care.home=="Independent"
notlisted <- cc.severe$listed.any == 0
nocare.notlisted <- nocare & notlisted

##### tables for first paper  ######################
cat("Tables for first paper ...")

table.casegr <- univariate.tabulate(outcome="severe.casegr",
                                    varnames=c("AGE", "sex", "care.home",
                                               "scrip.any", "diag.any", listed.conditions),
                                    data=cc.severe[CASE==1])
rownames(table.casegr) <- replace.names(rownames(table.casegr))
rownames(table.casegr)[1] <- "Age [median (IQR)]"

##############################################
# table.severe.demog <- tabulate.freqs.regressions(varnames=demog, data=cc.severe)

table.scripordiag <- NULL
for(agegr in levels(cc.severe$agegr20)) {
    x <- with(cc.severe[cc.severe$agegr20==agegr, ],
              table(scripordiag, CASE))
    colnames(x) <- c("Controls", "Cases")
    x <- paste.colpercent(x)
    y <- with(cc.severe[cc.severe$agegr20==agegr, ],
              table(listed.any, CASE))
    y <- paste.colpercent(y)
    y <- y[2, , drop=FALSE]
    rownames(y) <- "Any listed condition"
    z <- with(cc.severe[cc.severe$agegr20==agegr, ],
              table(listed.any==0 & diag.any==1, CASE))
    z <- paste.colpercent(z)
    z <- z[2, , drop=FALSE]
    rownames(z) <- "No listed condition, but other diagnosis"
    
    x <- rbind(x[1, , drop=FALSE], y, z, x[2, , drop=FALSE])
    table.scripordiag <- cbind(table.scripordiag, x)
}

table.scripordiag.fatal <- NULL
for(agegr in levels(cc.severe$agegr3)) {
    x <- with(cc.severe[cc.severe$agegr3==agegr, ],
              table(scripordiag, fatalcase))
    colnames(x) <- c("Controls", "Fatal cases")
    x <- paste.colpercent(x)
    # x <- x[1, , drop=FALSE]
    rownames(x) <- c("No prescription or diagnosis", "Prescription or diagnosis")
    table.scripordiag.fatal <- rbind(table.scripordiag.fatal, x)
}


cat("by age group ... ")
## table.agegr shown in paper 1: frequencies by age group but no rate ratios
## varnames.listed
table.agegr <- NULL
for(agegr in levels(cc.severe$agegr20)) {
    x <- univariate.tabulate(varnames=varnames.extended, 
                             outcome="CASE",
                             data=cc.severe[cc.severe$agegr20==agegr, ],
                             drop.reflevel=FALSE)
    table.agegr <- cbind(table.agegr, x)
}
## all ages including 28-day fatality
freqs.all <- univariate.tabulate(varnames=c("deathwithin28", varnames.extended),
                             outcome="CASE",
                             data=cc.severe,
                             drop.reflevel=FALSE)

## list of tables by age group shown in paper 1
## varnames.listedpluscovs
tables.agegr <- vector("list", length(levels(cc.severe$agegr3)))
for(i in 1:length(levels(cc.severe$agegr3))) {
    agegr <- levels(cc.severe$agegr3)[i]
    tables.agegr[[i]] <-
        tabulate.freqs.regressions(varnames=varnames.listedpluscovs,
                                   data=cc.severe[cc.severe$agegr3==agegr, ])
}

## table.agegr.all shown in paper 1: rate ratios for all age groups combined
## varnames.listedpluscovs
table.agegr.all <- tabulate.freqs.regressions(varnames=varnames.listedpluscovs,
                                              data=cc.severe)
cat("demographic vars ... ")
## demographic vars

if(old) {
    table.demog.aug <- tabulate.freqs.regressions(varnames=demog, data=cc.severe)
}

## separate analysis using SMR ethnicity 
table.ethnicsmr <- univariate.tabulate(varnames="ethnic4.smr", outcome="CASE",
                                       data=cc.severe[!is.na(cc.severe$ethnic4.smr), ],
                                       drop.reflevel=FALSE)
univariate.ethnicsmr <-
    univariate.clogit(varnames="ethnic4.smr",
                      data=cc.severe[!is.na(cc.severe$ethnic4.smr), ],
                      add.reflevel=TRUE)
table.ethnicsmr.aug <- combine.tables2(table.ethnicsmr, univariate.ethnicsmr)
rownames(table.ethnicsmr.aug) <- replace.names(rownames(table.ethnicsmr.aug))

cat("listed conditions ... ")
## listed conditions
table.listed.conditions.lt60 <-
    tabulate.freqs.regressions(varnames=listed.conditions,
                               data=cc.severe[cc.severe$AGE < 60, ])
table.listed.conditions.ge60 <-
    tabulate.freqs.regressions(varnames=listed.conditions,
                               data=cc.severe[cc.severe$AGE >= 60 ])

## full multivariate model variables -- for this use dm.type rather than diabetes.any

summary(cc.severe[, ..full.varnames])

multivariate.all <-
    multivariate.clogit(varnames=full.varnames, data=cc.severe, add.reflevel=TRUE)

################# restrict to those without listed conditions #############
cat("restricting to those without listed conditions ...")

## conditions
table.conditions.aug <- tabulate.freqs.regressions(varnames=icd.chapternames, 
                                                   data=cc.severe[notlisted, ])
cat("done\n")

cat("Tabulating ICD subchapter diagnoses ...")
## tabulate subchapters in ICD chapters of interest
table.icdchapter2 <-
    tabulate.freqs.regressions(varnames=grep("^Ch_II:",
                                             colnames(cc.severe), value=TRUE),
                               data=cc.severe[notlisted, ])

table.icdchapter7 <-
    tabulate.freqs.regressions(varnames=grep("^Ch_VII:",
                                             colnames(cc.severe),
                                             value=TRUE),
                               data=cc.severe[notlisted, ])

table.icdchapter11 <-
    tabulate.freqs.regressions(varnames=grep("^Ch_XI:",
                                             colnames(cc.severe),
                                             value=TRUE),
                               data=cc.severe[notlisted, ])

icd.subchapternames <- grep("^Ch_", names(cc.severe), value=TRUE)
subchapters.freqs <- univariate.tabulate(varnames=icd.subchapternames, outcome="CASE",
                                         data=cc.severe[notlisted],
                                         drop.reflevel=FALSE, drop.sparserows=FALSE)
keep.subchapters <- rowSums(matrix(as.integer(gsub(" .+", "", subchapters.freqs)),
                                   ncol=2)) > 50
subchapters.freqs <- subchapters.freqs[keep.subchapters, ]
subchapters.clogit <- univariate.clogit(varnames=icd.subchapternames[keep.subchapters],
                                        outcome="CASE",
                                        data=cc.severe[notlisted], add.reflevel=TRUE)
table.icdsubchapters <- combine.tables2(subchapters.freqs, subchapters.clogit)

cat("done\n")

#########################################################################

## drugs 
table.drugs.aug <- tabulate.freqs.regressions(varnames=bnf.chapternames, 
                                              data=cc.severe[notlisted, ])

################################################################

## tabulate associations with drug chapters in those not in care homes and without listed conditions 
table.drugs.nocare.notlisted <- tabulate.freqs.regressions(varnames=bnf.chapternames, 
                                                           data=cc.severe[nocare.notlisted, ])
