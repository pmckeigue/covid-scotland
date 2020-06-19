## stepwise regressions for casecontrol.R

lower.varnames <- "care.home"
demog.varnames <- c("care.home", "SIMD.quintile")
listed.varnames <- c(demog.varnames, "dm.type", listed.conditions)
full.varnames <- c(listed.varnames, "diag.any", "scrip.any", drugs)

lower.formula <- as.formula(paste0("CASE ~ ",
                            paste(lower.varnames, collapse="+"), 
                            " + strata(stratum)"))
demog.formula <- as.formula(paste0("CASE ~ ",
                            paste(demog.varnames, collapse="+"), 
                            " + strata(stratum)"))
listed.formula <- as.formula(paste0("CASE ~ ",
                            paste(listed.varnames, collapse="+"), 
                            " + strata(stratum)"))
full.formula <- as.formula(paste0("CASE ~ ",
                            paste(full.varnames, collapse="+"), 
                            " + strata(stratum)"))

## fit models
nonmissing <- nonmissing.obs(cc.severe, full.varnames)
cc.nonmissing <- subset(cc.severe, subset=nonmissing)

demog.model <- clogit(formula=demog.formula, data=cc.nonmissing)
listed.model <- clogit(formula=listed.formula, data=cc.nonmissing)
full.model <- clogit(formula=full.formula, data=cc.nonmissing)

if(stepwise) { # stepwise for full variable set
    stepwise.full <- step(full.model,
                          scope=list(lower=lower.formula, upper=full.formula),
                          direction="both", method="approximate", trace=-1)
    stepwise.full <- summary(stepwise.full)$coefficients
    rownames(stepwise.full) <- replace.names(rownames(stepwise.full))

########## cross-validation of stepwise regression ######
    ## each stratum is assigned to one of nfold disjoint test folds
    test.folds <- testfolds.bystratum(stratum=cc.nonmissing[["stratum"]],
                                      y=cc.nonmissing[["CASE"]], nfold=nfold)
    ## test.folds has one row per stratum, one column per fold containing indicator vars
    ## merge by stratum adds single column test.fold
    cv.data <- merge(cc.nonmissing, test.folds, by="stratum")
    cv.data <- subset(cv.data, select=c("test.fold", "CASE", "stratum", full.varnames))
    cat("cv.data has", nrow(cv.data), "rows\n")
   
###########################################
    demog.predicted <- cv.predict(nfold, cv.data, lower.formula, demog.formula) 
    demog.predicted <- demog.predicted[demog.predicted$prior.p < 1, ]
    demog.densities <- with(demog.predicted,
                            Wdensities(y, posterior.p, prior.p,
                                       recalibrate=FALSE))
    pander(summary(demog.densities), table.style="multiline",
           split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
           caption="Prediction of severe COVID-19 from demographic variables only")
    
########################################
    listed.predicted <- cv.predict(nfold, cv.data, lower.formula, listed.formula) 
    listed.predicted <- listed.predicted[listed.predicted$prior.p < 1, ]
    listed.densities <- with(listed.predicted,
                             Wdensities(y, posterior.p, prior.p,
                                        recalibrate=FALSE))
    pander(summary(listed.densities), table.style="multiline",
           split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
           caption="Prediction of severe COVID-19 from demographic + listed variables")
    
######################################
    full.predicted <- cv.predict(nfold, cv.data, lower.formula, full.formula) 
    full.predicted <- full.predicted[full.predicted$prior.p < 1, ]
    full.densities <- with(full.predicted,
                           Wdensities(y, posterior.p, prior.p,
                                      recalibrate=FALSE))
    pander(summary(full.densities), table.style="multiline",
           split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
           caption="Prediction of severe COVID-19 from full variable set")

    if(old) {
    save(stepwise.full,
         demog.densities, listed.densities, full.densities, 
         file="./data/stepwise_15May.RData")
    } else {
            save(stepwise.full,
         demog.densities, listed.densities, full.densities, 
         file="./data/stepwise.RData")

    }

} else {
    if(old) {
        load("./data/stepwise_15May.RData")
    } else {
         load("./data/stepwise.RData")
    }
}

