## pharmacoepi paper

########################################################
## tabulate rate ratios for each num drug group by listed.any
## this shows that the polypharmacy association is only in those without listed conditions

table.numdrugsgr.listed.any <- NULL
for(cat in levels(cc.severe$listed.any)) {
    x <- tabulate.freqs.regressions(varnames="numdrugsgr",
                                    data=cc.severe[cc.severe$listed.any==cat, ])[, 1:4]
    colnames(x)[1:2] <- c("controls", "cases")
    table.numdrugsgr.listed.any <- rbind(table.numdrugsgr.listed.any, x)
}
table.numdrugsgr.listed.any <- data.frame(numdrugsgr=rep(levels(cc.severe$numdrugsgr), 2),
                                     table.numdrugsgr.listed.any)
colnames(table.numdrugsgr.listed.any)[2:3] <- paste0(c("Controls (", "Cases ("), as.integer(table(cc.severe$CASE)), rep(")", 2))


########################################################################
## in those without listed conditions, by age group

table.numdrugsgr.agegr2 <- NULL
for(cat in levels(cc.severe$agegr2)) {
    x <- tabulate.freqs.regressions(varnames="numdrugsgr",
                                    data=subset(cc.severe, notlisted & agegr2==cat))[, 1:4]
    colnames(x)[1:2] <- c("controls", "cases")
    table.numdrugsgr.agegr2 <- rbind(table.numdrugsgr.agegr2, x)
}
table.numdrugsgr.agegr2 <- data.frame(numdrugsgr=rep(levels(cc.severe$numdrugsgr), 2),
                                     table.numdrugsgr.agegr2)
colnames(table.numdrugsgr.agegr2)[2:3] <- paste0(c("Controls (", "Cases ("), as.integer(table(cc.severe$CASE[notlisted])), rep(")", 2))

########################################################################

## in those without listed conditions, fit joint model for num cardiovascular and num non-cardiovascular drugs

tj <- summary(clogit(formula=CASE ~ numdrugs.cardiovasc + numdrugs.notcardiovasc + strata(stratum),
                     data=cc.severe[notlisted, ]))$coefficients

tj <- as.data.frame(tj)
tj$u.ci <- or.ci(tj[, 1], tj[, 3])
tj$pvalue <- pvalue.latex(tj[, 5])
means.cardio <- with(cc.severe[notlisted, ], tapply(numdrugs.cardiovasc, CASE, mean)) 
means.notcardio <- with(cc.severe[notlisted, ], tapply(numdrugs.notcardiovasc, CASE, mean))

means.numdrugs <- rbind(means.cardio, means.notcardio)
drugcat <- c("Cardiovascular", "Other")

table.jointcardiovasc <- data.frame(drugcat, round(means.numdrugs, 1),
                                    tj[, c("u.ci", "pvalue")])
colnames(table.jointcardiovasc)[2:3] <- paste0(c("Controls (", "Cases ("),
                                               table(cc.severe[notlisted, ]$CASE), rep(")", 2))

################################################################
## generate subpara.coeffs and subparacols.keep
## consider restricting to age < 75
## FIXME - rewrite code below to use subset
restrict <- with(cc.severe, listed.any==0)

## select subpara codes excluding cardiovascular (chapter 2)
subparacols <- grep("^subpara\\.[^2]", colnames(cc.severe))
#notcardiovasc.subparacols <- !grepl("^subpara\\.2",  colnames(cc.severe)[subparacols.forstepwisedrop]

y <- subset(cc.severe, restrict)[["CASE"]]
stratum <- subset(cc.severe, restrict)[["stratum"]]
subparacols.keep <- NULL
subpara.coeffs <- NULL
for(col in subparacols) {
    x <- subset(cc.severe, restrict)[[col]]
    freqs.cc <- table(as.integer(x), y)
    if(nrow(freqs.cc) > 1 & ncol(freqs.cc) > 1) {
        if(sum(freqs.cc[2, ]) >= 20) { # if at least 20 in row
            coeffs <- summary(clogit(formula=y ~ x + strata(stratum)))$coefficients
            if(coeffs[1, 5] < 0.001) {
                subpara.coeffs <- rbind(subpara.coeffs,
                                        data.frame(subpara=colnames(cc.severe)[col],
                                                   controls=freqs.cc[2, 1],
                                                   cases=freqs.cc[2, 2], 
                                                   rate.ratio=or.ci(coeffs[, 1], coeffs[, 3]),
                                                   pvalue=pvalue.latex(coeffs[, 5])))
                subparacols.keep <- c(subparacols.keep, col)
            }
        }
    }
}
subpara.coeffs$subpara <- gsub("^subpara\\.", "", subpara.coeffs$subpara)
freqs <- as.integer(with(subset(cc.severe, restrict), table(CASE)))
colnames(subpara.coeffs)[2:3] <- paste0(c("Controls (", "Cases ("),
                                        freqs, rep(")", 2))
print(subpara.coeffs)

#############################################################################

colnames.drugclasses.forstepwise <- colnames(cc.severe)[subparacols.keep]

table.subpara.coeffs <- tabulate.freqs.regressions(colnames.drugclasses.forstepwise,
                                                   data=subset(cc.severe, restrict))

## stepwise drop procedure over subparacols.keep

## replace original subparas with nonopiate and anyopiate (scrips for compound preparations assigned to opiate)
#colnames.drugclasses.forstepwise <-
#    grep("\\.407020\\.", colnames.drugclasses.forstepwise, invert=TRUE, value=TRUE)

subparacols.forstepwisedrop <- match(colnames.drugclasses.forstepwise,
                                     colnames(cc.severe)) ## reduces to 61 columns

x <- subset(cc.severe, subset=restrict, select=subparacols.forstepwisedrop)
x <- matrix(as.integer(as.matrix(x)), nrow=nrow(x))
colnames(x) <- colnames(cc.severe)[subparacols.forstepwisedrop]
y <- subset(cc.severe, restrict)[["CASE"]]
stratum <- subset(cc.severe, restrict)[["stratum"]]
covariates.subparas <- subset(cc.severe, restrict)[["care.home"]] 
cat("Stepwise drop procedure over BNF subpara codes adjusted for",
    names(covariates.subparas), "...")
stepwise.drop.subparas <- stepwise.union.dropcols(x=x, y=y, covariates=covariates.subparas, stratum=stratum)
cat("done\n")

print(stepwise.drop.subparas)

################################################################

## tabulate dose-response effects by age group


#### opiates  #####################################
####### tabulate dose-response effect by opiate dose group ###############################

opiate.covariates <- c("care.home", "SIMD.quintile", "neoplasm.any")
opiate.covariates.string <- paste(opiate.covariates, collapse="+")
coeff.rows <- 1:3
withcovariates.formula <- as.formula(paste("CASE ~ dosegr.opiate +",
                                           opiate.covariates.string, "+ strata(stratum)"))

table.dosegr.opiate <- NULL
for(agegr in levels(cc.severe$agegr2)) {
    x <- summary(clogit(formula=CASE ~ dosegr.opiate + strata(stratum), 
                        data=cc.severe[cc.severe$agegr2==agegr, ]))$coefficients
    x <- as.data.frame(x)
    x$u.ci <- or.ci(x[, 1], x[, 3])
    x$u.pvalue <- pvalue.latex(x[, 5])
    x <- x[, c("u.ci", "u.pvalue")]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.severe[cc.severe$agegr2==agegr, ]))$coefficients[coeff.rows, ] # drop=FALSE]
    y <- as.data.frame(y)
    y$m.ci <- or.ci(y[, 1], y[, 3])
    y$m.pvalue <- pvalue.latex(y[, 5])
    y <- y[, c("m.ci", "m.pvalue")]
    x <- data.frame(x, y)               
   
    x <- rbind(rep(NA, ncol(x)), x)
    freqs <- paste.colpercent(with(cc.severe[cc.severe$agegr2==agegr, ], table(dosegr.opiate, CASE)))
    x <- data.frame(DDD.average=rownames(freqs), freqs, x)
    rownames(x) <- paste0(agegr, ": ", rownames(x))
    table.dosegr.opiate <- rbind(table.dosegr.opiate, x)
}

colnames(table.dosegr.opiate)[2:3] <- paste0(c("Controls (", "Cases ("),
                                                 as.integer(table(cc.severe$CASE)),
                                                 rep(")", 2))

##############################################################################
## tabulate.freqs.regressions() can only be used with factor variables

compound.opiates.varnames <- c("dosegr.compound.opiates", "care.home", "neoplasm.any")
table.dosegr.compound.opiates <- tabulate.freqs.regressions(varnames=compound.opiates.varnames,
                                                            data=cc.severe)

####### tabulate proton pump dose-response effect by age group ###############################

ppi.covariates <- c("care.home", "SIMD.quintile", "esoph.stomach.duod", "nsaid",
                    "antiplatelet", "anticoagulant.any")

ppi.covariates.string <- paste(ppi.covariates, collapse="+")
withcovariates.formula <- as.formula(paste("CASE ~ DDDsgr +",
                                           ppi.covariates.string, "+ strata(stratum)"))

table.dosegr.protonpump <- NULL
coeff.rows <- 1:4
for(agegr in levels(cc.severe$agegr2)) {
    x <- summary(clogit(formula=CASE ~ DDDsgr + strata(stratum), 
                        data=cc.severe[cc.severe$agegr2==agegr, ]))$coefficients
    x <- as.data.frame(x)
    x$u.ci <- or.ci(x[, 1], x[, 3])
    x$u.pvalue <- pvalue.latex(x[, 5])
    x <- x[, c("u.ci", "u.pvalue")]
    y <- summary(clogit(formula=withcovariates.formula,
                        data=cc.severe[cc.severe$agegr2==agegr, ]))$coefficients[coeff.rows, ] # drop=FALSE]
    y <- as.data.frame(y)
    y$m.ci <- or.ci(y[, 1], y[, 3])
    y$m.pvalue <- pvalue.latex(y[, 5])
    y <- y[, c("m.ci", "m.pvalue")]
    x <- data.frame(x, y)               
   
    x <- rbind(rep(NA, ncol(x)), x)
    freqs <- paste.colpercent(with(cc.severe[cc.severe$agegr2==agegr, ], table(DDDsgr, CASE)))
    x <- data.frame(DDD.average=rownames(freqs), freqs, x)
    rownames(x) <- paste0(agegr, ": ", rownames(x))
    table.dosegr.protonpump <- rbind(table.dosegr.protonpump, x)
}
colnames(table.dosegr.protonpump)[2:3] <- paste0(c("Controls (", "Cases ("),
                                                 as.integer(table(cc.severe$CASE)),
                                                 rep(")", 2))

## print(table.dosegr.protonpump)

#### opiate effect by time window ##################################

## excluding neoplasms

table.timewindow.opiate <-
    tabulate.freqs.regressions(varnames="opiate.exposurecat",
                               data=cc.severe[cc.severe$neoplasm.any==0, ])[, 1:4]
 
############# proton pump effect by time window #####################

table.timewindow.protonpump <-
    tabulate.freqs.regressions(varnames="protonpump.exposurecat",
                               data=cc.severe)[, 1:4]

#######################################################################

######### tabulate ever-use effect with adjustment for prespecified covariates


table.everuse.protonpump <- NULL
withcovariates.formula <- as.formula(paste("CASE ~ protonpump +",
                                           ppi.covariates.string, "+ strata(stratum)"))
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


####### tabulate dose-response effect with DDDs as continuous, adj for covariates ##########

## tabulate.freqs.regressions() can only be used with factor variables

table.dose.protonpump <- NULL
withcovariates.formula <- as.formula(paste("CASE ~ DDDs.all +",
                                           ppi.covariates.string, "+ strata(stratum)"))
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
colnames(freqs)[1:2] <- paste0(c("Controls (", "Cases ("),
                               as.integer(table(cc.severe$CASE)), rep(")", 2))
freqs <- as.data.frame(freqs)
table.dose.protonpump  <- cbind(freqs, table.dose.protonpump)
colnames(table.dose.protonpump)[1:2] <- paste0(c("Controls (", "Cases ("),
                                               as.integer(table(cc.severe$CASE)), rep(")", 2))



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


## tabulate effect of each drug
ppinames <- grep("PRAZOLE$", names(cc.severe), value=TRUE)
ppinames.formula <- as.formula(paste("CASE ~", paste(ppinames, collapse="+"), "+ strata(stratum)"))

ppi.means <- cbind(
    round(apply(subset(cc.severe, subset=CASE==0, select=ppinames), 2, mean), 4),
    round(apply(subset(cc.severe, subset=CASE==1, select=ppinames), 2, mean), 4)
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

############## tabulate effect of any proton pump exposure by age group ###############

withcovariates.formula <- as.formula(paste("CASE ~ protonpump +",
                                           ppi.covariates.string,
                                           "+ strata(stratum)"))
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

summary(clogit(formula=CASE ~ protonpump + numdrugs.notppi + strata(stratum),
               data=cc.severe[nocare.notlisted & cc.severe$AGE < 60, ]))$coefficients

#############################################################################

## tabulate fatal cases  by age group
table.fatal.protonpump <- NULL
for(agegr in levels(cc.severe$agegr3)) {
    x <- tabulate.freqs.regressions(varnames="protonpump", outcome="fatalcase",
                                    data=subset(cc.severe, subset=agegr3==agegr))[, 1:4]
    rownames(x) <- agegr
    colnames(x)[1:2] <- c("Controls", "Cases")
    table.fatal.protonpump <- rbind(table.fatal.protonpump, x)
}
x <- tabulate.freqs.regressions(varnames="protonpump", data=cc.severe)[, 1:4]
rownames(x) <- "All"
colnames(x)[1:2] <- c("Controls", "Cases")
table.fatal.protonpump <- rbind(table.fatal.protonpump, x)

#################################################################
## tabulate para or subpara codes in BNF chapters of interest
table.bnfchapter1 <- tabulate.bnfparas(chnum=1, data=cc.severe, minrowsum=50)
table.bnfchapter2 <- tabulate.bnfsubparas(chnum=2, data=cc.severe, minrowsum=50)
table.bnfchapter4 <- tabulate.bnfsubparas(chnum=4, data=cc.severe, minrowsum=50)
## chapter 10 has to be disaggregated to drug names as grouping of DMARDs is too broad
table.bnfchapter10 <- tabulate.bnfchemicals(chnum=10, data=cc.severe, minrowsum=50)

nocare.drugfreqs <-
    sapply(subset(cc.severe,
                  subset=care.home=="Independent" & cc.severe$CASE==0,
                  select=subparacols),
           function(x) table(x)[2] / sum(table(x))) %>% sort(decreasing=TRUE) %>% head(20)

care.drugfreqs <-
    sapply(subset(cc.severe,
                  subset=care.home=="Care/nursing home" & cc.severe$CASE==0,
                  select=subparacols),
           function(x) table(x)[2] / sum(table(x))) %>% sort(decreasing=TRUE) %>% head(20)

#############################

coeff.anticoagulant.univariate <- summary(clogit(formula=CASE ~ anticoagulant.any + strata(stratum),
               data=cc.severe))$coefficients[1, ]

coeff.anticoagulant.adjusted <- summary(clogit(formula=CASE ~ anticoagulant.any +  heart.other.any + protonpump + strata(stratum),
               data=cc.severe))$coefficients[1, ]

############################

coeff.hydroxychlor.univariate <- summary(clogit(formula=CASE ~ hydroxychloroquine + strata(stratum),
               data=cc.severe))$coefficients

coeff.hydroxychlor.adjusted <- summary(clogit(formula=CASE ~ hydroxychloroquine + protonpump + strata(stratum),
               data=cc.severe))$coefficients[1, ]

##############################################################################
