## pharmacoepi paper


## list of logical vectors to be used as argument to compare.timewindows()
subsets.laporte.scrips <- list(
    ppi = with(scrips, bnf_paragraph_code == "0103050"),
    antispasm.gi = with(scrips, bnf_paragraph_code == "0102000"),
    antihist = with(scrips, bnf_paragraph_code == "0304010"),
    
    hypno = with(scrips, bnf_paragraph_code == "0401010"),
    anxio = with(scrips, bnf_paragraph_code == "0401020"),
    antipsych = with(scrips, bnf_paragraph_code == "0402010"),
    manic = with(scrips, bnf_paragraph_code == "0402030"),
    
    tricyclic = with(scrips, bnf_paragraph_code == "0403010"),
    ssri = with(scrips, bnf_paragraph_code == "0403030"),
    antidepr.other = with(scrips, bnf_paragraph_code == "0403040"),
    
    nausea = with(scrips, bnf_paragraph_code == "0406000"),
    opioid = with(scrips, bnf_paragraph_code == "0407020"),
    gaba = with(scrips, approved_name=="GABAPENTIN" | approved_name=="PREGABALIN"),
    antiepilep.other = with(scrips, bnf_paragraph_code=="0408010" &
                                    !(approved_name=="GABAPENTIN" | approved_name=="PREGABALIN")),
    antispasm.ur = with(scrips, bnf_paragraph_code == "0406000"),
    nsaid = with(scrips, bnf_paragraph_code == "0406000")
)


########################################################
## tabulate rate ratios for each num drug group by listed.any
## this shows that the polypharmacy association is only in those without listed conditions

table.notcv.numdrugsgr.listed.any <- NULL
for(cat in levels(cc.severe$listed.any)) {
    x <- tabulate.freqs.regressions(varnames="notcv.numdrugsgr",
                                    data=cc.severe[listed.any==cat])[, 1:4]
    colnames(x)[1:2] <- c("controls", "cases")
    table.notcv.numdrugsgr.listed.any <- rbind(table.notcv.numdrugsgr.listed.any, x)
}
table.notcv.numdrugsgr.listed.any <- data.frame(notcv.numdrugsgr=rep(levels(cc.severe$notcv.numdrugsgr), 2),
                                     table.notcv.numdrugsgr.listed.any)
colnames(table.notcv.numdrugsgr.listed.any)[2:3] <- paste0(c("Controls (", "Cases ("), as.integer(table(cc.severe$CASE)), rep(")", 2))

## tabulate by care home status

table.notcv.numdrugsgr.care.home <- NULL
for(cat in levels(cc.severe$care.home)) {
    x <- tabulate.freqs.regressions(varnames="notcv.numdrugsgr",
                                    data=cc.severe[care.home==cat])[, 1:4]
    colnames(x)[1:2] <- c("controls", "cases")
    table.notcv.numdrugsgr.care.home <- rbind(table.notcv.numdrugsgr.care.home, x)
}
table.notcv.numdrugsgr.care.home <- data.frame(notcv.numdrugsgr=rep(levels(cc.severe$notcv.numdrugsgr), 2),
                                     table.notcv.numdrugsgr.care.home)
colnames(table.notcv.numdrugsgr.care.home)[2:3] <- paste0(c("Controls (", "Cases ("), as.integer(table(cc.severe$CASE)), rep(")", 2))


## in those notresident in care homes, tabulate by listed.any

table.notcv.numdrugsgr.listed.any <- NULL
for(cat in levels(cc.severe[care.home=="Independent"]$listed.any)) {
    x <- tabulate.freqs.regressions(varnames="notcv.numdrugsgr",   data=cc.severe[care.home=="Independent" & listed.any==cat])[, 1:4]
    colnames(x)[1:2] <- c("controls", "cases")
    table.notcv.numdrugsgr.listed.any <- rbind(table.notcv.numdrugsgr.listed.any, x)
}
table.notcv.numdrugsgr.listed.any <- data.frame(notcv.numdrugsgr=rep(levels(cc.severe$notcv.numdrugsgr), 2),
                                                table.notcv.numdrugsgr.listed.any)
colnames(table.notcv.numdrugsgr.listed.any)[2:3] <- paste0(c("Controls (", "Cases ("), as.integer(table(cc.severe[nocare]$CASE)), rep(")", 2))


########################################################################
## in those not resident in care homes, tabulate by age group

table.numdrugsgr.agegr2 <- NULL
for(cat in levels(cc.severe$agegr2)) {
    x <- tabulate.freqs.regressions(varnames="numdrugsgr",
                       data=cc.severe[care.home=="Independent" & agegr2==cat])[, 1:4]
    colnames(x)[1:2] <- c("controls", "cases")
    table.numdrugsgr.agegr2 <- rbind(table.numdrugsgr.agegr2, x)
}
table.numdrugsgr.agegr2 <- data.frame(numdrugsgr=rep(levels(cc.severe$numdrugsgr), 2),
                                      table.numdrugsgr.agegr2)
colnames(table.numdrugsgr.agegr2)[2:3] <- paste0(c("Controls (", "Cases ("), as.integer(table(cc.severe[nocare, CASE])), rep(")", 2))

########################################################################

## in those without listed conditions, fit joint model for num cardiovascular and num non-cardiovascular drugs

tj <- summary(clogit(formula=CASE ~ numdrugs.cardiovasc + numdrugs.notcardiovasc + strata(stratum),
                     data=cc.severe[notlisted]))$coefficients

tj <- as.data.frame(tj)
tj$u.ci <- or.ci(tj[, 1], tj[, 3])
tj$pvalue <- pvalue.latex(tj[, 5])
means.cardio <- with(cc.severe[notlisted], tapply(numdrugs.cardiovasc, CASE, mean)) 
means.notcardio <- with(cc.severe[notlisted], tapply(numdrugs.notcardiovasc, CASE, mean))

means.numdrugs <- rbind(means.cardio, means.notcardio)
drugcat <- c("Cardiovascular", "Other")

table.jointcardiovasc <- data.frame(drugcat, round(means.numdrugs, 1),
                                    tj[, c("u.ci", "pvalue")])
colnames(table.jointcardiovasc)[2:3] <- paste0(c("Controls (", "Cases ("),
                                               table(cc.severe[notlisted, ]$CASE), rep(")", 2))

## joint model with num cardiovascular and noncardiovascular drugs grouped 

coeffs.cvnotcv <- summary(clogit(CASE ~ cv.numdrugsgr + notcv.numdrugsgr + strata(stratum), data=cc.severe))$coefficients

coeffs.data <- data.frame(numdrugsgr=c(0:3, 0:5),
                              Indication=c(rep("Cardiovascular", 4), rep("Other", 6)),
                              coeff=c(0, coeffs.cvnotcv[1:3, 1], 0, coeffs.cvnotcv[4:8, 1]),
                              se=c(0, coeffs.cvnotcv[1:3, 3], 0, coeffs.cvnotcv[4:8, 3]))

coeffs.cvnotcv.notlisted <- summary(clogit(CASE ~ cv.numdrugsgr + notcv.numdrugsgr + strata(stratum), data=cc.severe[notlisted]))$coefficients

coeffs.cvnotcv.listed <- summary(clogit(CASE ~ cv.numdrugsgr + notcv.numdrugsgr + strata(stratum), data=cc.severe[!notlisted]))$coefficients

################################################################
## generate subpara.coeffs 
restrict <- nocare.notlisted

## select subpara codes excluding cardiovascular (chapter 2)
subparacols <- grep("^subpara\\.[^2]", colnames(cc.severe))

y <- cc.severe[restrict, CASE]
stratum <- cc.severe[restrict, stratum]
subparacols.keep <- NULL
subpara.coeffs <- NULL
for(col in subparacols) {
    x <- cc.severe[restrict][[col]]
    freqs.cc <- table(as.integer(x), y)
    if(nrow(freqs.cc) > 1 & ncol(freqs.cc) > 1) {
        if(sum(freqs.cc[2, ]) >= 50) { # if at least 20 in row
            coeffs <- summary(clogit(formula=y ~ x + strata(stratum)))$coefficients
            print(coeffs)
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
freqs <- as.integer(with(cc.severe[restrict], table(CASE)))
colnames(subpara.coeffs)[2:3] <- paste0(c("Controls (", "Cases ("),
                                        freqs, rep(")", 2))
print(subpara.coeffs)

#############################################################################

################################################################

## tabulate dose-response effects by age group

#### opiates  #####################################
####### tabulate dose-response effect by opiate dose group ###############################

opiate.covariates <- c("care.home", "SIMD.quintile", "neoplasm.any")
opiate.covariates.string <- paste(opiate.covariates, collapse="+")
coeff.rows <- 1:3
withcovariates.formula <- as.formula(paste("CASE ~ dosegr.opiate +",
                                           opiate.covariates.string, "+ strata(stratum)"))

table.dosegr.opioid <- NULL
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
    x <- data.frame(MME.average=rownames(freqs), freqs, x)
    rownames(x) <- paste0(agegr, ": ", rownames(x))
    table.dosegr.opioid <- rbind(table.dosegr.opioid, x)
}

colnames(table.dosegr.opioid)[2:3] <- paste0(c("Controls (", "Cases ("),
                                             as.integer(table(cc.severe$CASE)),
                                             rep(")", 2))

x <- summary(clogit(formula=CASE ~ dosegr.opiate + strata(stratum), 
                    data=subset(cc.severe, notlisted)))$coefficients
x <- as.data.frame(x)
x$u.ci <- or.ci(x[, 1], x[, 3])
x$u.pvalue <- pvalue.latex(x[, 5])
x <- x[, c("u.ci", "u.pvalue")]
y <- summary(clogit(formula=withcovariates.formula,
                    data=subset(cc.severe, notlisted)))$coefficients[coeff.rows, ]
y <- as.data.frame(y)
y$m.ci <- or.ci(y[, 1], y[, 3])
y$m.pvalue <- pvalue.latex(y[, 5])
y <- y[, c("m.ci", "m.pvalue")]
x <- data.frame(x, y)               
x <- rbind(rep(NA, ncol(x)), x)
freqs <- paste.colpercent(with(subset(cc.severe, notlisted), table(dosegr.opiate, CASE)))
x <- data.frame(MME.average=rownames(freqs), freqs, x)
table.dosegr.opioid.notlisted <- x 

colnames(table.dosegr.opioid.notlisted)[2:3] <-
    paste0(c("Controls (", "Cases ("),
           as.integer(table(subset(cc.severe, notlisted)[["CASE"]])),
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

cc.severe$DDDsgr4 <- car::recode(cc.severe$DDDsgr,
                                 "c('1.5', '2 or more')='1.5 or more'")

withcovariates.formula <- as.formula(paste("CASE ~ DDDsgr4 +",
                                           ppi.covariates.string, "+ strata(stratum)"))
coeff.rows <- 1:3

x <- summary(clogit(formula=CASE ~ DDDsgr4 + strata(stratum), 
                    data=cc.severe[notlisted & AGE < 75]))$coefficients
x <- as.data.frame(x)
x$u.ci <- or.ci(x[, 1], x[, 3])
x$u.pvalue <- pvalue.latex(x[, 5])
x <- x[, c("u.ci", "u.pvalue")]
y <- summary(clogit(formula=withcovariates.formula,
                    data=cc.severe[notlisted & AGE < 75]))$coefficients[coeff.rows, ]
y <- as.data.frame(y)
y$m.ci <- or.ci(y[, 1], y[, 3])
y$m.pvalue <- pvalue.latex(y[, 5])
y <- y[, c("m.ci", "m.pvalue")]
x <- data.frame(x, y)               
x <- rbind(rep(NA, ncol(x)), x)
freqs <- paste.colpercent(with(cc.severe[notlisted & AGE < 75], table(DDDsgr4, CASE)))
x <- data.frame(DDDs.average=rownames(freqs), freqs, x)
table.dosegr.protonpump.notlisted.agelt75 <- x 

colnames(table.dosegr.protonpump.notlisted.agelt75)[2:3] <-
    paste0(c("Controls (", "Cases ("),
           as.integer(table(cc.severe[notlisted & AGE < 75][["CASE"]])),
           rep(")", 2))

#### opiate effect by time window ##################################

## excluding neoplasms

table.timewindow.opioid <-
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
## exclude care home residents otherwise associations dominated by care home prescribing

table.bnfchapter1 <- tabulate.bnfparas(chnum=1, data=cc.severe[nocare], minrowsum=50)
table.bnfchapter2 <- tabulate.bnfsubparas(chnum=2, data=cc.severe[nocare], minrowsum=50)
table.bnfchapter4 <- tabulate.bnfsubparas(chnum=4, data=cc.severe[nocare], minrowsum=50)
## chapter 10 has to be disaggregated to drug names as grouping of DMARDs is too broad
table.bnfchapter10 <- tabulate.bnfchemicals(chnum=10, subpara.exclude="1003020",
                                            data=cc.severe[nocare],
                                            minrowsum=50)

## antipsychotics
table.bnfchapter4.antipsych <- tabulate.bnfchemicals(chnum=4, subpara="0402010", data=cc.severe[nocare & neoplasm.lastyear=="0"], minrowsum=20)
##  anti-nauseants
table.bnfchapter4.nauseavertigo <- tabulate.bnfchemicals(chnum=4, subpara="0406000", data=cc.severe[nocare & neoplasm.lastyear=="0"], minrowsum=20)
## opioids
table.bnfchapter4.opioid <- tabulate.bnfchemicals(chnum=4, subpara="0407020", data=cc.severe[nocare & neoplasm.lastyear=="0"], minrowsum=20)
## anti-epileptics
table.bnfchapter4.antiepilep <- tabulate.bnfchemicals(chnum=4, subpara="0408010",
                     data=cc.severe[nocare & neoplasm.lastyear=="0"], minrowsum=20)

table.bnfchapter4.chem <- rbind(table.bnfchapter4.antipsych,
                                table.bnfchapter4.nauseavertigo,
                                table.bnfchapter4.opioid,
                                table.bnfchapter4.antiepilep)[, 1:4]

############################################################################

coeff.anticoagulant.univariate <- summary(clogit(formula=CASE ~ anticoagulant.any + strata(stratum),
               data=cc.severe[nocare]))$coefficients[1, ]

coeff.anticoagulant.adjusted <- summary(clogit(formula=CASE ~ anticoagulant.any + IHD.any +  heart.other.any + protonpump + strata(stratum),
               data=cc.severe[nocare]))$coefficients[1, ]

############################

coeff.hydroxychlor.univariate <- summary(clogit(formula=CASE ~ hydroxychloroquine + strata(stratum),
               data=cc.severe[nocare]))$coefficients

coeff.hydroxychlor.adjusted <- summary(clogit(formula=CASE ~ hydroxychloroquine + protonpump + strata(stratum),
               data=cc.severe[nocare]))$coefficients[1, ]

##############################################################################
## tabulate LaPorte-Healy list of drug classes

table.laporte <- tabulate.freqs.regressions(varnames=c("SIMD.quintile", subparas.laporte),
                                            data=cc.severe[care.home=="Independent"])
rownames(table.laporte) <- gsub("^subpara\\.", "", rownames(table.laporte))
varnames=c("SIMD.quintile", subparas.laporte)
cc.laporte <- cc.severe[care.home=="Independent"] %>%
        select(CASE, stratum, all_of(varnames)) %>% na.omit()

laporte.formula <- as.formula(CASE ~ . + strata(stratum))
laporte.model <- clogit(formula=laporte.formula, data=cc.laporte)

recode.as0 <- function(x) {
    x[x=="1"] <- "0"
    return(x)
}

cc.nodrugs <- copy(cc.laporte)
cc.nodrugs[, (subparas.laporte) := lapply(.SD, recode.as0),
           .SDcols=subparas.laporte]
laporte.matrix <-  model.matrix(laporte.formula, cc.laporte)
nodrugs.matrix <- model.matrix(laporte.formula, cc.nodrugs)
keep.cols <- grep("stratum", colnames(laporte.matrix), invert=TRUE)[-1]
laporte.matrix <- laporte.matrix[, keep.cols]
nodrugs.matrix <- nodrugs.matrix[, keep.cols]
coeffs.matrix <- matrix(laporte.model$coefficients[-1], ncol=1)
lp.laporte  <- laporte.matrix %*% coeffs.matrix
lp.nodrugs  <- nodrugs.matrix %*% coeffs.matrix
mean.oddsratio <- mean(exp(lp.laporte - lp.nodrugs))
aetiologic.fraction <- (mean.oddsratio - 1) / mean.oddsratio

###############################################################
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

coeffs.withlaporte <- summary(clogit(CASE ~ numdrugs.cardiovasc + numdrugs.notcardiovasc + numdrugs.laporte + strata(stratum), data=cc.severe))$coefficients
print(coeffs.withlaporte)

summary(clogit(CASE ~ numdrugs.cardiovasc + numdrugs.notcardiovasc + numdrugs.laporte + strata(stratum), data=cc.severe))$coefficients

table.timewindows <- NULL

## subsets.laporte.scrips is a list of 16 logical vectors each of length 3763801
keep <- rep(FALSE, length(subsets.laporte.scrips))
for(i in 1:length(subsets.laporte.scrips)) {
    compare <- compare.timewindows(scrips[subsets.laporte.scrips[[i]]], data=cc.severe, recent.cutoff=120)
    num.nonrecentonly <- sum(as.integer(gsub(" .+", "", as.matrix(compare["Non-recent only", 1:2]))))
    if(num.nonrecentonly >= 500) {
        table.timewindows <- rbind(table.timewindows, compare[1:2, ])
        keep[i] <- TRUE
    }
}

table.timewindows <- 
data.frame(twin=rep(c("Non-recent only", "Recent only"), nrow(table.timewindows) / 2), 
 table.timewindows)
drugclasses.timewindows <- subparas.laporte[keep]

##########################################
## table with care home as response variable, drugs as covariates

cc.severe[, care.home.int := as.integer(care.home) -1]

new.carehometable <- FALSE

if(new.carehometable) {
    sfreqs <- lapply(cc.severe[CASE==0 & care.home.int==1, ..subparacols], table)

    ## keep.carehome is logical vector of rows to keep
    keep.carehome <- unlist(lapply(sfreqs,
                                    function(x) ifelse(length(x)==2 & x[2] > 50, TRUE, FALSE)))
    subparacols.carehome <- subparacols[keep.carehome]
    table.drugs.carehome <-
        tabulate.freqs.regressions(varnames=names(cc.severe)[subparacols.carehome],
                                   outcome="care.home.int",
                                   data=cc.severe[CASE==0 & AGE >= 75])
    
    save(table.drugs.carehome, file="table.drugs.carehome.RData")
} else {
    load(file="table.drugs.carehome.RData")
}

table.drugs.carehome <- table.drugs.carehome[grep("ensuremath",
                                                  table.drugs.carehome$u.pval), 1:4]
rownames(table.drugs.carehome) <- gsub("subpara\\.", "", rownames(table.drugs.carehome))
table.drugs.carehome <- table.drugs.carehome[grep("^1[123][0-9]+{5}",
                                                  rownames(table.drugs.carehome),
                                                  invert=TRUE), ]
                                                                  
ssri <- tabulate.freqs.regressions(varnames="subpara.403030.Selective serotonin re-uptake inhibitors",
                           data=cc.severe)[, 1:4]

    if(FALSE) {
        scrips.colchicine <- 
            subset(scrips, approved_name=="COLCHICINE",
                   select=c("ANON_ID", "SPECIMENDATE",
                            "dispensed_date", "item_strength",
                            "quantity"))            
        scrips.last15.colchicine <- 
            subset(scrips.last15, approved_name=="COLCHICINE",
                   select=c("ANON_ID", "SPECIMENDATE",
                            "dispensed_date", "item_strength",
                            "quantity"))
        scrips.colchicine <- rbind(scrips.colchicine, scrips.last15.colchicine) %>%
            as.data.table(key="ANON_ID")
        scrips.colchicine[, end.date := dispensed_date + floor(quantity / 2)]
        scrips.colchicine[, on.colchicine := dispensed_date < SPECIMENDATE &
                                SPECIMENDATE <= end.date]
        scrips.colchicine <- scrips.colchicine[on.colchicine==TRUE]
        cc.all[, colchicine.current := as.factor(as.integer(ANON_ID %in%
                                                            scrips.colchicine$ANON_ID))]
    }
}

