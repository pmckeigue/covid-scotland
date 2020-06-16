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

## fit conditional logistic regression model for CASE ~ care.home + numdrugs
cc.severe$numdrugs.sq <- cc.severe$numdrugs.subpara^2
cc.severe.nonmissing <- cc.severe[nonmissing.obs(cc.severe, "SIMD.quintile"), ]

multi.formula <- as.formula(paste("CASE ~ care.home + SIMD.quintile + diag.any +",
                                  paste(listed.conditions, collapse="+"),
                                  "+ numdrugsgr + numdrugs.cardiovasc + strata(stratum)"))

model.numdrugs <- clogit(formula=multi.formula, data=cc.severe.nonmissing)
print(summary(model.numdrugs)$coefficients)

theta <- model.numdrugs$linear.predictors
unnorm.p <- 1 / (1 + exp(-theta))
norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                        stratum=cc.severe.nonmissing$stratum,
                                       y=cc.severe.nonmissing$CASE)
rm(cc.severe.nonmissing)

norm.predicted <- norm.predicted[norm.predicted$prior.p < 1, ]

## FIXME -- fitted densities are not compatible with Turing identity
numdrugs.densities <- with(norm.predicted,
                           Wdensities(y, posterior.p, prior.p,
                                      recalibrate=FALSE, adjust.bw=0.5))
pander(summary(numdrugs.densities), table.style="multiline",
       split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
       caption="Prediction of severe COVID-19 from number of drugs dispensed")

# plotWdists(numdrugs.densities)

## tabulate rate ratios for each numdrugs group, by num.icdchapters.gr.  

table.numdrugsgr.num.icdch <- NULL
for(numicd in levels(cc.severe$num.icdchapters.gr)) {
    x <- tabulate.freqs.regressions(varnames="numdrugsgr",
                                    data=cc.severe[notlisted & cc.severe$num.icdchapters.gr==numicd, ])[, 1:4]
    colnames(x)[1:2] <- c("controls", "cases")
    table.numdrugsgr.num.icdch <- rbind(table.numdrugsgr.num.icdch, x)
}
table.numdrugsgr.num.icdch <- data.frame(numdrugsgr=rep(levels(cc.severe$numdrugsgr), 3),
                                         table.numdrugsgr.num.icdch)
colnames(table.numdrugsgr.num.icdch)[2:3] <- paste0(c("Controls (", "Cases ("), as.integer(table(cc.severe[notlisted, ]$CASE)), rep(")", 2))

####################################################################################
## tabulate rate ratios for each numdrugs group, by agegr 

table.numdrugsgr.agegr <- NULL
for(agegr in levels(cc.severe$agegr3)) {
    x <- tabulate.freqs.regressions(varnames="numdrugsgr",
                                    data=cc.severe[cc.severe$agegr3==agegr, ])[, 1:4]
    colnames(x)[1:2] <- c("controls", "cases")
    table.numdrugsgr.agegr <- rbind(table.numdrugsgr.agegr, x)
}
table.numdrugsgr.agegr <- data.frame(numdrugsgr=rep(levels(cc.severe$numdrugsgr), 3),
                               table.numdrugsgr.agegr)
colnames(table.numdrugsgr.agegr)[2:3] <- paste0(c("Controls (", "Cases ("), as.integer(table(cc.severe$CASE)), rep(")", 2))

## tabulate rate ratios for each num drug group by care home in those without listed conditions)

table.numdrugsgr.carehome <- NULL
for(cat in levels(cc.severe$care.home)) {
    x <- tabulate.freqs.regressions(varnames="numdrugsgr",
                                    data=cc.severe[cc.severe$listed.any==0 &
                                                   cc.severe$care.home==cat, ])[, 1:4]
    colnames(x)[1:2] <- c("controls", "cases")
    table.numdrugsgr.carehome <- rbind(table.numdrugsgr.carehome, x)
}
table.numdrugsgr.carehome <- data.frame(numdrugsgr=rep(levels(cc.severe$numdrugsgr), 2),
                                     table.numdrugsgr.carehome)


################################################################

subparacols <- grep("^subpara\\.", colnames(cc.severe))

y <- cc.severe[nocare.notlisted, ][["CASE"]]
stratum <- cc.severe[nocare.notlisted, "stratum"]
subparacols.keep <- NULL
subpara.coeffs <- NULL
for(col in subparacols) {
    x <- cc.severe[nocare.notlisted, ][[col]]
            freqs.cc <- table(as.integer(x),
                              cc.severe[nocare.notlisted, ][["CASE"]])
    
    if(nrow(freqs.cc) > 1 & ncol(freqs.cc) > 1) {
        if(sum(freqs.cc[2, ]) >= 20) { # if at least 10 in each cell
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
freqs <- as.integer(with(cc.severe[nocare.notlisted, ], table(CASE)))
colnames(subpara.coeffs)[2:3] <- paste0(c("Controls (", "Cases ("),
                                        freqs, rep(")", 2))
print(subpara.coeffs)

#############################################################################
restrict <- notlisted

subparacols.forstepwisedrop <- subparacols[match(subparacols.keep, subparacols)] ## reduces to 61 subpara codes
notcardiovasc.subparacols <- !grepl("^subpara\\.2",
                                   colnames(cc.severe)[subparacols.forstepwisedrop])
subparacols.forstepwisedrop <- subparacols.forstepwisedrop[notcardiovasc.subparacols]

x <- subset(cc.severe, subset=restrict, select=subparacols.forstepwisedrop)
x <- matrix(as.integer(as.matrix(x)), nrow=nrow(x))
colnames(x) <- colnames(cc.severe)[subparacols.forstepwisedrop]
y <- cc.severe[restrict, ][["CASE"]]
stratum <- cc.severe[restrict, "stratum"]
covariates.subparas <- cc.severe[restrict, c("care.home")] 
cat("Stepwise drop procedure over BNF subpara codes adjusted for",
    names(covariates.subparas), "...")
stepwise.drop.subparas <- stepwise.union.dropcols(x=x, y=y, covariates=covariates.subparas, stratum=stratum)
cat("done\n")

print(stepwise.drop.subparas)

## tabulate associations with drug chapters in those not in care homes and without listed conditions 
table.drugs.nocare.notlisted <- tabulate.freqs.regressions(varnames=drugs, 
                                                           data=cc.severe[nocare.notlisted, ])

#### opiates  #####################################
####### tabulate dose-response effect by opiate dose group ###############################

opiates.varnames <- c("dosegr.opiate", "care.home", "neoplasm.any")

table.dosegr.opiate <- tabulate.freqs.regressions(varnames=opiates.varnames,
                                                  data=cc.severe)

compound.opiates.varnames <- c("dosegr.compound.opiates", "care.home", "neoplasm.any")

table.dosegr.compound.opiates <- tabulate.freqs.regressions(varnames=compound.opiates.varnames,
                                                            data=cc.severe)


#### opiate effect by time window #######################

cc.severe$opiate.exposure.nonrecent <- as.integer(cc.severe$opiateMME.interval2 > 0 | cc.severe$opiateMME.interval3 > 0)
cc.severe$opiate.exposure.recent <- as.integer(cc.severe$opiateMME.interval1 > 0)

cc.severe$opiate.exposurecat <- integer(nrow(cc.severe))
cc.severe$opiate.exposurecat[cc.severe$opiate.exposure.nonrecent==0 & cc.severe$opiate.exposure.recent==0] <- 0
cc.severe$opiate.exposurecat[cc.severe$opiate.exposure.nonrecent==1 & cc.severe$opiate.exposure.recent==0] <- 1
cc.severe$opiate.exposurecat[cc.severe$opiate.exposure.nonrecent==0 & cc.severe$opiate.exposure.recent==1] <- 2
cc.severe$opiate.exposurecat[cc.severe$opiate.exposure.nonrecent==1 & cc.severe$opiate.exposure.recent==1] <- 3

cc.severe$opiate.exposurecat <- as.factor(car::recode(cc.severe$opiate.exposurecat,
                                "0='No prescriptions'; 1='Non-recent only'; 2='Recent only';
                                 3='Prescriptions in both time windows'"))
cc.severe$opiate.exposurecat <- factor(cc.severe$opiate.exposurecat,
                                levels=levels(cc.severe$opiate.exposurecat)[c(1, 2, 4, 3)])

table.timewindow.opiate <-
    tabulate.freqs.regressions(varnames="opiate.exposurecat",
                               data=cc.severe[cc.severe$neoplasm.any==0, ])[, 1:4]
 
############# proton pump #########################

## calculate propensity score trained on cc.hosp

cc.hosp <- readRDS("cchosp.rds")
source("propensity.R")
numrows.cchosp <- nrow(cc.hosp)
rm(cc.hosp)

cc.severe$propensity <- propensity
cc.severe$propensitygr <- ceiling(cc.severe$propensity)
cc.severe$propensitygr <- as.factor(car::recode(cc.severe$propensitygr,
                                                "1='>0 to 1'; 2='>1 to 2';
                                                 3='>2 to 3'; 4=;'>3 to 4'; 5:hi='>4'"))

cc.severe$propensitygr <- factor(cc.severe$propensitygr,
                                 levels=levels(cc.severe$propensitygr)[c(5, 1:4)])

######### tabulate ever-use effect with adjustment for prespecified covariates

ppi.covariates <- c("care.home + SIMD.quintile + esoph.stomach.duod + nsaid + antiplatelet + anticoagulant.any")

ppi.covariates.string <- paste(ppi.covariates, collapse="+")

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


####### tabulate dose-response effect  by DDDs group ###############################

## tabulate.freqs.regressions() can only be used with factor variables

table.dosegr.protonpump <- NULL
withcovariates.formula <- as.formula(paste("CASE ~ DDDsgr +",
                                           ppi.covariates.string, "+ strata(stratum)"))
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
    x <- data.frame(DDD.average=rownames(freqs), freqs, x)
    rownames(x) <- paste0(agegr, ": ", rownames(x))
    table.dosegr.protonpump <- rbind(table.dosegr.protonpump, x)
}
colnames(table.dosegr.protonpump)[2:3] <- paste0(c("Controls (", "Cases ("),
                                                 as.integer(table(cc.severe$CASE)), rep(")", 2))
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

cc.severe$ppi.exposure.nonrecent <- as.integer(cc.severe$DDD.interval2 > 0 | cc.severe$DDD.interval3 > 0)
cc.severe$ppi.exposure.recent <- as.integer(cc.severe$DDD.interval1 > 0)

cc.severe$ppi.exposurecat <- integer(nrow(cc.severe))
cc.severe$ppi.exposurecat[cc.severe$ppi.exposure.nonrecent==0 & cc.severe$ppi.exposure.recent==0] <- 0
cc.severe$ppi.exposurecat[cc.severe$ppi.exposure.nonrecent==1 & cc.severe$ppi.exposure.recent==0] <- 1
cc.severe$ppi.exposurecat[cc.severe$ppi.exposure.nonrecent==0 & cc.severe$ppi.exposure.recent==1] <- 2
cc.severe$ppi.exposurecat[cc.severe$ppi.exposure.nonrecent==1 & cc.severe$ppi.exposure.recent==1] <- 3

cc.severe$ppi.exposurecat <- as.factor(car::recode(cc.severe$ppi.exposurecat,
                                "0='No prescriptions'; 1='Non-recent only'; 2='Recent only';
                                 3='Prescriptions in both time windows'"))
cc.severe$ppi.exposurecat <- factor(cc.severe$ppi.exposurecat,
                                levels=levels(cc.severe$ppi.exposurecat)[c(4, 2, 3, 1)])

############## time window analysis ######################

table.timewindow.protonpump <- tabulate.freqs.regressions(varnames="ppi.exposurecat",
                                                          data=cc.severe)[, 1:4]

table.timewindow.protonpump.nocare.notlisted <-
    tabulate.freqs.regressions(varnames="ppi.exposurecat",
                               data=cc.severe[nocare.notlisted, ])[, 1:4]

###########################################

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

summary(clogit(formula=CASE ~ protonpump + numdrugs.notppi + strata(stratum),
               data=cc.severe[nocare.notlisted & cc.severe$AGE < 60, ]))$coefficients

#############################################################################
    
table.dose.propensity <- NULL
for(pr in levels(cc.severe$propensitygr)) {
    x <- summary(clogit(formula=CASE ~ DDDs.all + strata(stratum),
                data=cc.severe[cc.severe$propensitygr==pr, ]))$coefficients
    rownames(x) <- pr
    table.dose.propensity <- rbind(table.dose.propensity, x)
}


###################################################
    
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

