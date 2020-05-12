## stepwise regressions for casecontrol.R

lower.varnames <- "care.home"
listed.varnames <- c("care.home", "dm.type", listed.conditions)
full.varnames <- c(listed.varnames, "diag.any", drugs)

lower.formula <- as.formula(paste0("CASE ~ ",
                            paste(lower.varnames, collapse="+"), 
                            " + strata(stratum)"))
listed.formula <- as.formula(paste0("CASE ~ ",
                            paste(listed.varnames, collapse="+"), 
                            " + strata(stratum)"))
full.formula <- as.formula(paste0("CASE ~ ",
                            paste(full.varnames, collapse="+"), 
                            " + strata(stratum)"))

## fit models
care.model <- clogit(formula=lower.formula, data=cc.severe)
listed.model <- clogit(formula=listed.formula, data=cc.severe)
full.model <- clogit(formula=full.formula, data=cc.severe)

if(stepwise) {
## stepwise for listed conditions
    stepwise.listed <- step(listed.model,
                            scope=list(lower=lower.formula, upper=listed.formula),
                            direction="both", method="approximate", trace=-1)
    stepwise.listed <- summary(stepwise.listed)$coefficients
    rownames(stepwise.listed) <- replace.names(rownames(stepwise.listed))
    print(stepwise.listed)
    
    ## stepwise for full variable set
    stepwise.full <- step(full.model,
                          scope=list(lower=lower.formula, upper=full.formula),
                          direction="both", method="approximate", trace=-1)
    stepwise.full <- summary(stepwise.full)$coefficients
    rownames(stepwise.full) <- replace.names(rownames(stepwise.full))
    print(stepwise.full)

################ cross-validation of stepwise regression
    
    cc.nonmissing <- nonmissing.obs(cc.severe, full.varnames)
    
    test.folds <- testfolds.bystratum(stratum=cc.nonmissing$stratum,
                                      y=cc.nonmissing$CASE, nfold=nfold)
    cv.data <- merge(cc.nonmissing, test.folds, by="stratum")
    
    ##demog.predicted <- foreach(i=1:nfold, .combine=rbind) %do% {
    care.predicted <- NULL
    for(i in 1:nfold) {
        test.data  <- cv.data[cv.data$test.fold == i, ]
        train.data <- cv.data[cv.data$test.fold != i, ]
        xcare.model <- clogit(formula=lower.formula, data=train.data)
        ## normalize within each stratum
        unnorm.p <- predict(object=xcare.model, newdata=test.data,
                            na.action="na.pass", 
                            type="risk", reference="sample")
        norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                                stratum=test.data$stratum,
                                                y=test.data$CASE)
        listed.predicted <- rbind(care.predicted, norm.predicted)
    }
    listed.densities <- with(listed.predicted,
                             Wdensities(y, posterior.p, prior.p,
                                        recalibrate=FALSE))
    pander(summary(listed.densities), table.style="multiline",
           split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
           caption="Prediction of severe COVID-19 from care home only")
    
    ##demog.predicted <- foreach(i=1:nfold, .combine=rbind) %do% {
    listed.predicted <- NULL
    for(i in 1:nfold) {
        test.data  <- cv.data[cv.data$test.fold == i, ]
        train.data <- cv.data[cv.data$test.fold != i, ]
        start.model <- clogit(formula=listed.formula, data=train.data)
        stepwise.model <- step(start.model,
                               scope=list(lower=lower.formula, upper=listed.formula),
                               direction="both", trace=-1, method="approximate")
        ## normalize within each stratum
        unnorm.p <- predict(object=stepwise.model, newdata=test.data,
                            na.action="na.pass", 
                            type="risk", reference="sample")
        norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                                stratum=test.data$stratum,
                                                y=test.data$CASE)
        listed.predicted <- rbind(listed.predicted, norm.predicted)
    }
    listed.densities <- with(listed.predicted,
                             Wdensities(y, posterior.p, prior.p,
                                        recalibrate=FALSE))
    pander(summary(listed.densities), table.style="multiline",
           split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
           caption="Prediction of severe COVID-19 from care home + listed conditions only")
    
    
    full.predicted <- NULL
    for(i in 1:nfold) {
        test.data  <- cv.data[cv.data$test.fold == i, ]
        train.data <- cv.data[cv.data$test.fold != i, ]
        start.model <- clogit(formula=full.formula, data=train.data)
        stepwise.model <- step(start.model,
                               scope=list(lower=lower.formula, upper=full.formula),
                               direction="both", trace=-1, method="approximate")
        ## normalize within each stratum
        unnorm.p <- predict(object=stepwise.model, newdata=test.data,
                            na.action="na.pass", 
                            type="risk", reference="sample")
        norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                                stratum=test.data$stratum,
                                                y=test.data$CASE)
        full.predicted <- rbind(full.predicted, norm.predicted)
    }
    full.densities <- with(full.predicted,
                           Wdensities(y, posterior.p, prior.p,
                                      recalibrate=FALSE))
    
    pander(summary(full.densities), table.style="multiline",
           split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
           caption="Prediction of severe COVID-19 from full variable set")
    
    save(stepwise.listed, stepwise.full,
         listed.densities, full.densities, 
         file="./data/stepwise.RData")
} else {
    load("./data/stepwise.RData")
}

