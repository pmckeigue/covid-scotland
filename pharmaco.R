## pharmacoepi paper

####### derived variables  ##########################


##############################################################################

## tabulate rate ratios for each num drug group, by agegr 

table.numdrugsgr.agegr <- NULL
for(agegr in levels(cc.severe$agegr3)) {
    x <- tabulate.freqs.regressions(varnames="numdrugsgr",
                                    data=cc.severe[cc.severe$agegr3==agegr, ])[, 1:4]
    colnames(x)[1:2] <- c("controls", "cases")
    table.numdrugsgr.agegr <- rbind(table.numdrugsgr.agegr, x)
}
table.numdrugsgr.agegr <- data.frame(numdrugsgr=rep(levels(cc.severe$numdrugsgr), 3),
                               table.numdrugsgr.agegr)

## tabulate rate ratios for each num drug group by cats3 (carehome, listed conditions)

table.numdrugsgr.cats3 <- NULL
for(cat in levels(cc.severe$cats3)) {
    x <- tabulate.freqs.regressions(varnames="numdrugsgr",
                                    data=cc.severe[cc.severe$cats3==cat, ])[, 1:4]
    colnames(x)[1:2] <- c("controls", "cases")
    table.numdrugsgr.cats3 <- rbind(table.numdrugsgr.cats3, x)
}
table.numdrugsgr.cats3 <- data.frame(numdrugsgr=rep(levels(cc.severe$numdrugsgr), 3),
                                     table.numdrugsgr.cats3)

## fit joint model for num cardiovascular and num non-cardiovascular drugs

tj <- NULL
for(agegr in levels(cc.severe$agegr3)) {
    x <- summary(clogit(formula=CASE ~ numdrugs.cardiovasc + numdrugs.notcardiovasc + strata(stratum),
                        data=cc.severe[nocare.notlisted & cc.severe$agegr3==agegr, ]))$coefficients
    tj <- rbind(tj, x)
}

tj <- as.data.frame(tj)
tj$u.ci <- or.ci(tj[, 1],
                                    tj[, 3])
tj$pvalue <- pvalue.latex(tj[, 5])
tj$drugscat <- rep(c("Cardiovascular", "Non-cardiovascular"), 3)

means.cardio <- with(cc.severe[nocare.notlisted, ], tapply(numdrugs.cardiovasc, list(CASE, agegr3), mean)) 
means.notcardio <- with(cc.severe[nocare.notlisted, ], tapply(numdrugs.notcardiovasc, list(CASE, agegr3), mean))

means.numdrugs <- rbind(means.cardio, means.notcardio)
means.numdrugs <- t(matrix(means.numdrugs, nrow=2))
colnames(means.numdrugs) <- c("Controls", "Cases")

table.jointcardiovasc <- data.frame(drugscat=tj$drugscat,
                                    round(means.numdrugs, 1),
                                    tj[, c("u.ci", "pvalue")])

################################################################

subparacols <- grep("^subpara\\.", colnames(cc.severe))
y <- cc.severe[nocare.notlisted, "CASE"]
stratum <- cc.severe[nocare.notlisted, "stratum"]
subparacols.keep <- NULL
subpara.coeffs <- NULL
for(col in subparacols) {
    x <- cc.severe[nocare.notlisted, col]
    freqs <- table(as.integer(x))
    if(min(freqs) >= 10 & length(freqs) > 1) { # if at least 10 in each cell
        coeffs <- summary(clogit(formula=y ~ x + strata(stratum)))$coefficients
        if(coeffs[1, 5] < 0.001 & coeffs[1, 1] > 0) {
            subpara.coeffs <- rbind(subpara.coeffs,
                                    data.frame(subpara=colnames(cc.severe)[col],
                                               data.frame(coeffs)))
            subparacols.keep <- c(subparacols.keep, col)
        }
    }
}
print(subpara.coeffs)

subparacols.forstepwisedrop <- subparacols[match(subparacols.keep, subparacols)] ## reduces to 61 subpara codes

x <- cc.severe[nocare.notlisted, ][, subparacols.forstepwisedrop]
x <- matrix(as.integer(as.matrix(x)), nrow=nrow(x))
colnames(x) <- colnames(cc.severe)[subparacols.forstepwisedrop]
cat("Stepwise drop procedure over BNF subpara codes ...")
stepwise.drop.subparas <- stepwise.union.dropcols(x=x, y=y, stratum=stratum)
cat("done\n")

## tabulate associations with drug chapters in those not in care homes and without listed conditions 
table.drugs.notlisted <- tabulate.freqs.regressions(varnames=drugs, 
                                                           data=cc.severe[notlisted, ])

## tabulate proportion of effect of scrip.any that is explained by each chapter
## use cc.severe[notlisted, ]
table.anyscrip.chapter <-
    summary(clogit(formula=as.formula(paste("CASE ~ scrip.any + strata(stratum)")),
                                      data=cc.severe[notlisted, ]))$coefficients[1, 1:2, drop=FALSE]
for(bnfchapter in drugs) {
    ch.formula=as.formula(paste("CASE ~ scrip.any +", bnfchapter, "+ strata(stratum)"))
    table.anyscrip.chapter <-
        rbind(table.anyscrip.chapter,
              summary(clogit(ch.formula, data=cc.severe[notlisted, ]))$coefficients[1, 1:2])
}
rownames(table.anyscrip.chapter) <- c("Unadjusted", drugs)
table.anyscrip.chapter <- as.data.frame(table.anyscrip.chapter)
table.anyscrip.chapter$prop.explained <- round(with(table.anyscrip.chapter,
                                                    c(0, 1 - coef[-1] / coef[1])), 2)

## tabulate para or subpara codes in BNF chapters of interest
table.bnfchapter1 <- tabulate.bnfparas(chnum=1, data=cc.severe[notlisted, ])
table.bnfchapter2 <- tabulate.bnfsubparas(chnum=2, data=cc.severe[notlisted, ])
table.bnfchapter4 <- tabulate.bnfsubparas(chnum=4, data=cc.severe[notlisted, ])
table.bnfchapter9 <- tabulate.bnfsubparas(chnum=9, data=cc.severe[notlisted, ])
table.bnfchapter10 <- tabulate.bnfsubparas(chnum=10, data=cc.severe[notlisted, ])

## fix to tabulate BNF chemical substance

############# proton pump #########################

    ## calculate propensity score trained on cc.hosp
    
    cc.hosp <- readRDS("cchosp.rds")
    source("propensity.R")
    rm(cc.hosp)
    
    cc.severe$propensity <- propensity
    
    cc.severe$propensitygr <- ceiling(cc.severe$propensity)
    cc.severe$propensitygr <- as.factor(car::recode(cc.severe$propensitygr,
                                                    "5:hi='>4'"))
    cc.severe$propensitygr <- factor(cc.severe$propensitygr,
                                     levels=levels(cc.severe$propensitygr)[c(2:6, 1)])
}   


############# proton pump #########################

## calculate propensity score trained on cc.hosp
    
cc.hosp <- readRDS("cchosp.rds")
source("propensity.R")
rm(cc.hosp)

cc.severe$propensity <- propensity

cc.severe$propensitygr <- ceiling(cc.severe$propensity)
cc.severe$propensitygr <- as.factor(car::recode(cc.severe$propensitygr,
                                                "5:hi='>4'"))
cc.severe$propensitygr <- factor(cc.severe$propensitygr,
                                 levels=levels(cc.severe$propensitygr)[c(2:6, 1)])

####### effects of scrip and protonpump by care home / listed condition status ######################

table.protonpump.cats3 <- NULL
for(condition in levels(cc.severe$cats3)) {
    x <- tabulate.freqs.regressions(varnames="protonpump",
                                    data=cc.severe[cc.severe$cats3==condition, ])[, 1:4]
    rownames(x) <- condition
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.protonpump.cats3 <- rbind(table.protonpump.cats3, x)
}

######### tabulate ever-use effect

table.everuse.protonpump <- NULL
withcovariates.formula <- as.formula("CASE ~ protonpump + propensity + strata(stratum)")
coeff.row <- 1

for(agegr in levels(cc.severe$agegr3)) {
    x <- summary(clogit(formula=CASE ~ protonpump + strata(stratum), 
                        data=cc.severe[cc.severe$agegr3==agegr, ]))$coefficients[coeff.row, , drop=FALSE]
    x <- as.data.frame(x)
    x$u.ci <- or.ci(x[, 1], x[, 3])
    x$u.pvalue <- pvalue.latex(x[, 5])
    x <- x[, c("u.ci", "u.pvalue")]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.severe[cc.severe$agegr3==agegr, ]))$coefficients[coeff.row, , drop=FALSE]
    y <- as.data.frame(y)
    y$m.ci <- or.ci(y[, 1], y[, 3])
    y$m.pvalue <- pvalue.latex(y[, 5])
    y <- y[, c("m.ci", "m.pvalue")]
    x <- data.frame(x, y)               
    table.everuse.protonpump <- rbind(table.everuse.protonpump, x)
}

rownames(table.everuse.protonpump) <- NULL
freqs <- paste.colpercent(with(cc.severe, table(agegr3, CASE)))
colnames(freqs)[1:2] <- c("Controls", "Cases")
freqs <- as.data.frame(freqs)
table.everuse.protonpump  <- cbind(freqs, table.everuse.protonpump)

                                        # print(table.everuse.protonpump)

####### tabulate dose-response effect  ##########################################

## tabulate.freqs.regressions() can only be used with factor variables

table.dose.protonpump <- NULL
withcovariates.formula <- as.formula("CASE ~ DDDs.all + propensity + strata(stratum)")
coeff.row <- 1

for(agegr in levels(cc.severe$agegr3)) {
    x <- summary(clogit(formula=CASE ~ DDDs.all + strata(stratum), 
                        data=cc.severe[cc.severe$agegr3==agegr, ]))$coefficients[coeff.row, , drop=FALSE]
    x <- as.data.frame(x)
    x$u.ci <- or.ci(x[, 1], x[, 3])
    x$u.pvalue <- pvalue.latex(x[, 5])
    x <- x[, c("u.ci", "u.pvalue")]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.severe[cc.severe$agegr3==agegr, ]))$coefficients[coeff.row, , drop=FALSE]
    y <- as.data.frame(y)
    y$m.ci <- or.ci(y[, 1], y[, 3])
    y$m.pvalue <- pvalue.latex(y[, 5])
    y <- y[, c("m.ci", "m.pvalue")]
    x <- data.frame(x, y)               
    table.dose.protonpump <- rbind(table.dose.protonpump, x)
}

rownames(table.dose.protonpump) <- NULL
freqs <- paste.colpercent(with(cc.severe, table(agegr3, CASE)))
colnames(freqs)[1:2] <- c("Controls", "Cases")
freqs <- as.data.frame(freqs)
table.dose.protonpump  <- cbind(freqs, table.dose.protonpump)

                                        # print(table.dose.protonpump)

####### tabulate dose-response effect  ##########################################

## tabulate.freqs.regressions() can only be used with factor variables

table.dosegr.protonpump <- NULL
withcovariates.formula <- as.formula("CASE ~ DDDsgr + propensity + strata(stratum)")
coeff.rows <- 1:4

for(agegr in levels(cc.severe$agegr3)) {
    x <- summary(clogit(formula=CASE ~ DDDsgr + strata(stratum), 
                        data=cc.severe[cc.severe$agegr3==agegr, ]))$coefficients
    x <- as.data.frame(x)
    x$u.ci <- or.ci(x[, 1], x[, 3])
    x$u.pvalue <- pvalue.latex(x[, 5])
    x <- x[, c("u.ci", "u.pvalue")]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.severe[cc.severe$agegr3==agegr, ]))$coefficients[coeff.rows, ] # drop=FALSE]
    y <- as.data.frame(y)
    y$m.ci <- or.ci(y[, 1], y[, 3])
    y$m.pvalue <- pvalue.latex(y[, 5])
    y <- y[, c("m.ci", "m.pvalue")]
    x <- data.frame(x, y)               
   
    x <- rbind(rep(NA, ncol(x)), x)
    freqs <- paste.colpercent(with(cc.severe[cc.severe$agegr3==agegr, ], table(DDDsgr, CASE)))
    colnames(freqs)[1:2] <- c("Controls", "Cases")
    x <- data.frame(DDD.average=rownames(freqs), freqs, x)
    rownames(x) <- paste0(agegr, ": ", rownames(x))
    table.dosegr.protonpump <- rbind(table.dosegr.protonpump, x)
}
# print(table.dosegr.protonpump)

################# tabulate ever-use effect with adjustment for numdrugs

withcovariates.formula <- as.formula("CASE ~ protonpump + numdrugs.subpara + strata(stratum)")

coeff.row <- 1
table.everuse.protonpump.numdrugs <- NULL
for(agegr in levels(cc.severe$agegr20)) {
    x <- tabulate.freqs.regressions(varnames="protonpump",
                                    data=cc.severe[cc.severe$agegr20==agegr, ])[, 1:4]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.severe[cc.severe$agegr20==agegr, ]))$coefficients[coeff.row, , drop=FALSE]
    x$m.ci <- or.ci(y[, 1], y[, 3])
    x$m.pvalue <- pvalue.latex(y[, 5])
    rownames(x) <- agegr
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.everuse.protonpump.numdrugs <- rbind(table.everuse.protonpump.numdrugs, x)
}
x <- tabulate.freqs.regressions(varnames="protonpump", data=cc.severe)[, 1:4]
y <- summary(clogit(formula=withcovariates.formula,
                    data=cc.severe))$coefficients[coeff.row, , drop=FALSE]
x$m.ci <- or.ci(y[, 1], y[, 3])
x$m.pvalue <- pvalue.latex(y[, 5])
rownames(x) <- "All"
colnames(table.everuse.protonpump.numdrugs)[1:2] <- colnames(x)[1:2]
table.everuse.protonpump.numdrugs <- rbind(table.everuse.protonpump.numdrugs, x)

#### time window #######################

cc.severe$exposure.nonrecent <- as.integer(cc.severe$DDD.interval2 > 0 | cc.severe$DDD.interval3 > 0)
cc.severe$exposure.recent <- as.integer(cc.severe$DDD.interval1 > 0)

cc.severe$exposurecat <- integer(nrow(cc.severe))
cc.severe$exposurecat[cc.severe$exposure.nonrecent==0 & cc.severe$exposure.recent==0] <- 0
cc.severe$exposurecat[cc.severe$exposure.nonrecent==1 & cc.severe$exposure.recent==0] <- 1
cc.severe$exposurecat[cc.severe$exposure.nonrecent==0 & cc.severe$exposure.recent==1] <- 2
cc.severe$exposurecat[cc.severe$exposure.nonrecent==1 & cc.severe$exposure.recent==1] <- 3

cc.severe$exposurecat <- as.factor(car::recode(cc.severe$exposurecat,
                                "0='No prescriptions'; 1='Non-recent only'; 2='Recent only';
                                 3='Prescriptions in both time windows'"))
cc.severe$exposurecat <- factor(cc.severe$exposurecat,
                                levels=levels(cc.severe$exposurecat)[c(4, 2, 3, 1)])

summary(clogit(formula=CASE ~ exposurecat + strata(stratum), 
               data=cc.severe[cc.severe$AGE < 75, ]))$coefficients

 
###########################################

## tabulate effect of each drug
ppinames <- grep("PRAZOLE$", names(cc.severe), value=TRUE)
ppinames.formula <- as.formula(paste("CASE ~", paste(ppinames, collapse="+"), "+ strata(stratum)"))

ppi.means <- cbind(
    round(apply(cc.severe[cc.severe$CASE==0, ppinames], 2, mean), 4),
    round(apply(cc.severe[cc.severe$CASE==1, ppinames], 2, mean), 4)
)

table.ppinames <- NULL
for(ppi in ppinames) {
    x  <- summary(clogit(formula=as.formula(paste("CASE ~", ppi, "+ strata(stratum)")),
                         data=cc.severe))$coefficients
    table.ppinames <- rbind(table.ppinames, x)
}
y <- summary(clogit(formula=ppinames.formula, 
                    data=cc.severe))$coefficients
x <- as.data.frame(table.ppinames)
x$u.ci <- or.ci(x[, 1], x[, 3])
x$u.pvalue <- pvalue.latex(x[, 5])
x$m.ci <- or.ci(y[, 1], y[, 3])
x$m.pvalue <- pvalue.latex(y[, 5])
x <- x[, c("u.ci", "u.pvalue", "m.ci", "m.pvalue")]
table.ppinames <- data.frame(ppi.means, x)

################# tabulate effect of any proton pump exposure by age group #############################

    withcovariates.formula <- as.formula("CASE ~ protonpump + propensity + strata(stratum)")
coeff.row <- 1

table.protonpump <- NULL
for(agegr in levels(cc.severe$agegr20)) {
    x <- tabulate.freqs.regressions(varnames="protonpump",
                                    data=cc.severe[cc.severe$agegr20==agegr, ])[, 1:4]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.severe[cc.severe$agegr20==agegr, ]))$coefficients[coeff.row, , drop=FALSE]
    x$m.ci <- or.ci(y[, 1], y[, 3])
    x$m.pvalue <- pvalue.latex(y[, 5])
    rownames(x) <- agegr
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.protonpump <- rbind(table.protonpump, x)
}
    x <- tabulate.freqs.regressions(varnames="protonpump", data=cc.severe)[, 1:4]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.severe))$coefficients[coeff.row, , drop=FALSE]
x$m.ci <- or.ci(y[, 1], y[, 3])
    x$m.pvalue <- pvalue.latex(y[, 5])
    rownames(x) <- "All"
colnames(table.protonpump)[1:2] <- colnames(x)[1:2]
    table.protonpump <- rbind(table.protonpump, x)

############# effect of protonpump by propensity ##############################

table.propensitygr <- tabulate.freqs.regressions(varnames="propensitygr",
                                                 data=cc.severe)[, 1:4]
colnames.casecontrol <- colnames(table.propensitygr)[1:2]
    
    table.protonpump.bypropensity <- NULL
    for(pr in levels(cc.severe$propensitygr)) {
        x <- tabulate.freqs.regressions(varnames="protonpump",
                                        data=cc.severe[cc.severe$propensitygr==pr, ])
        rownames(x) <- pr
        x <- as.data.frame(x)
        colnames(x)[1:2] <- c("controls", "cases")
        table.protonpump.bypropensity <- rbind(table.protonpump.bypropensity, x)
    }

    table.protonpump.bypropensity <- table.protonpump.bypropensity[, 1:4]
    colnames(table.protonpump.bypropensity)[1:2] <- colnames.casecontrol

############ effect of protonpump by numdrugs ####################

    table.numdrugs.protonpump <- paste.colpercent(with(cc.severe[cc.severe$CASE==0, ], 
                                                       table(numdrugs.notppi.gr, protonpump)))

    colnames(table.numdrugs.protonpump) <- paste0("(",
                                                  table(cc.severe$CASE, cc.severe$protonpump)[1, ],
                                                  ")")
    
    colnames(table.numdrugs.protonpump) <- paste(c("Non-users", "Ever-users"),
                                                 colnames(table.numdrugs.protonpump)
                                                 )

    
    freqs <- NULL
    for(numdrugs in levels(cc.severe$numdrugs.notppi.gr)) {
        x <-  with(cc.severe[cc.severe$numdrugs.notppi.gr==numdrugs, ],
                   table(protonpump, CASE))
        freqs <- rbind(freqs, x)
    }
    
    
    table.protonpump.bynumdrugs <- NULL
    for(numdrugs in levels(cc.severe$numdrugs.notppi.gr)) {
    x <- tabulate.freqs.regressions(varnames="protonpump",
                                    data=cc.severe[cc.severe$numdrugs.notppi.gr==numdrugs, ])
    rownames(x) <- numdrugs
        colnames(x)[1:2] <- c("Controls", "Cases")
    x <- as.data.frame(x)
    table.protonpump.bynumdrugs <- rbind(table.protonpump.bynumdrugs, x)
    }
    table.protonpump.bynumdrugs <- table.protonpump.bynumdrugs[, 1:4]


#############################################################################

    
table.dose.propensity <- NULL
for(pr in levels(cc.severe$propensitygr)) {
    x <- summary(clogit(formula=CASE ~ DDDs.all + strata(stratum),
                data=cc.severe[cc.severe$propensitygr==pr, ]))$coefficients
    rownames(x) <- pr
    table.dose.propensity <- rbind(table.dose.propensity, x)
}

    cc.severe$exposurecat <- factor(cc.severe$exposurecat,
                                    levels=levels(cc.severe$exposurecat)[c(4, 1, 2, 3)])

    table.timewindow.protonpump <- tabulate.freqs.regressions(varnames="exposurecat",
                                                          data=cc.severe)[, 1:4]

    table.timewindow.protonpump.propensity0 <-
        tabulate.freqs.regressions(varnames="exposurecat",
                                   data=cc.severe[cc.severe$propensity==0, ])[, 1:4]
    
#####################################
    
## tabulate fatal cases  by age group
table.fatal.protonpump <- NULL
for(agegr in levels(cc.severe$agegr3)) {
    x <- tabulate.freqs.regressions(varnames="protonpump", outcome="fatalcase",
                                    data=cc.severe[cc.severe$agegr3==agegr, ])[, 1:4]
    rownames(x) <- agegr
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.fatal.protonpump <- rbind(table.fatal.protonpump, x)
}
x <- tabulate.freqs.regressions(varnames="protonpump", data=cc.severe)[, 1:4]
rownames(x) <- "All"
colnames(x)[1:2] <- c("Controls", "Cases")
table.fatal.protonpump <- rbind(table.fatal.protonpump, x)

