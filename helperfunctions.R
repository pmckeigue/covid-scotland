
## helperfunctions for COVID analyses


pvalue.latex <- function(x, n=1, nexp=1) {
  sapply(x, function(z) {
    if (is.na(z)) {
        return(NA)
    } else if(z > 0.001) {
        return(signif(z, 1))
    } else {
    z <- sprintf("%.*E", 0, signif(z, n)) # default is 1 sig fig
#    z <- sprintf("%.*E", n, signif(z, n + 1))
    z <- as.numeric(unlist(strsplit(z, "E"))) # split z at E
    sprintf("\\ensuremath{%.*f\\times 10^{%0*d}}", 0, z[1], nexp, z[2])
                                        #    sprintf("\\ensuremath{%.*f\\times 10^{%0*d}}", n, z[1], nexp + 1, z[2])
    }
  })
}
## create and format a 95% confidence interval around an odds/hazard ratio
or.ci <- function(coeff, se, ndigits=2) {
  ci.lo <- exp(coeff - 1.96 * se)
  ci.up <- exp(coeff + 1.96 * se)
  sprintf("%.*f (%.*f, %.*f)", ndigits, exp(coeff), ndigits, ci.lo, ndigits, ci.up)
}

univariate.tabulate <- function(varnames, outcomevar, data) {
    outcome <- data[, match(outcomevar, names(data))]
    table.varnames <- NULL
    for(i in 1:length(varnames)) {
        ## test whether variable is factor or numeric
        z <- data[, match(varnames[i], names(data))] 
        if(is.numeric(z)) { # median (IQR) for numeric variables 
            x <- tapply(z, outcome,
                        function(x) return(paste0(median(x, na.rm=TRUE),
                                                  " (", quantile(x, probs=0.25, na.rm=TRUE),
                                                  "-", quantile(x, probs=0.75, na.rm=TRUE), ")")))
            x <- matrix(x, nrow=1)
            rownames(x) <- varnames[i]
        } else { # freqs for factor variables
            x <- table(z, outcome)
            x <- paste.colpercent(x)
            x <- x[-1, , drop=FALSE] # drop reference category
            ## FIXME: should keep this line if length(levels(z)) > 2, instead add extra line to univariate
        }
        if(nrow(x)==1) {
            rownames(x) <- varnames[i] # use variable name if original variable had only 2 levels 
        }
        table.varnames <- rbind(table.varnames, x)
    }
    return(table.varnames)
}

replace.names <- function(varnames) {
    names <- varnames
    found <- names %in% lookup.names$varname
    names[found] <- as.character(lookup.names$longname[match(names[found],
                                                             lookup.names$varname)])
    return(names)
}

traintest.fold <- function(i, cv.data,
                           lower.formula, upper.formula) { # stepwise clogit on training fold, predict on test fold
    test.data  <- cv.data[cv.data$test.fold == i, ]
    train.data <- cv.data[cv.data$test.fold != i, ]
    start.model <- clogit(formula=upper.formula, data=train.data)
    cat("train.data dimensions", dim(train.data), "\n")
    stepwise.model <- step(start.model,
                           scope=list(lower=lower.formula, upper=upper.formula),
                           direction="both", trace=-1, method="approximate")
    ## normalize within each stratum
    unnorm.p <- predict(object=stepwise.model, newdata=test.data,
                        na.action="na.pass", 
                        type="risk", reference="sample")
    norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                            stratum=test.data$stratum,
                                            y=test.data$CASE)
    return(norm.predicted)
}

nonmissing.obs <- function(x, varnames) { ## subset rows nonmissing for varnames
    keep <- rep(TRUE, nrow(x))
    for(j in 1:length(varnames)) {
        keep[is.na(x[, match(varnames[j], colnames(x))])] <- FALSE
    }
    return(x[keep, ])
}

normalize.predictions <- function(unnorm.p, stratum, y) { # format a pvalue in latex
    ## normalize probs so that they sum to 1 within each stratum
    ## returns data frame with 4 columns: stratum, normconst, prior.p, posterior.p

    unnorm.p <- as.numeric(unnorm.p)
    stratum <- as.integer(as.character(stratum))
    nonmissing <- !is.na(unnorm.p)

    ## keep only strata that contain a single case 
    table.strata <- tapply(y, stratum, sum) == 1
    strata.onecase <- as.integer(names(table.strata)[as.integer(table.strata)==1])
    keep <- !is.na(unnorm.p) & stratum %in% strata.onecase

    cat(length(which(!keep)), "predictions dropped because stratum does not contain one case\n")

    unnorm.p <- unnorm.p[keep]
    stratum <- as.factor(stratum[keep])
    y <- y[keep]

    ## compute normalizing constant and prior for each stratum, then merge
    norm.constant <- as.numeric(tapply(unnorm.p, stratum, sum))
    stratum.size <- as.numeric(tapply(unnorm.p, stratum, length))
    print(table(stratum.size))
    prior.p <- 1/stratum.size
    norm.constant <- data.frame(stratum=levels(stratum),
                                normconst=norm.constant,
                                prior.p=as.numeric(prior.p))
    
    norm.predicted <- data.frame(unnorm.p, y, stratum)
    norm.predicted <- merge(norm.predicted, norm.constant, by="stratum")
    # normalize probs
    norm.predicted$posterior.p <- norm.predicted$unnorm.p / norm.predicted$normconst
    return(norm.predicted)
}

append.emptyrow <- function(x) {
    empty=matrix(c(rep.int(NA, length(x))), nrow=1, ncol=ncol(x))  
    colnames(empty) = colnames(x)
    return(rbind(x, empty))
}

paste.colpercent <- function(x, digits=0) { # paste column percentages into freq table
    x.colpct <- paste0("(", round(100 * prop.table(x, 2), digits), "%)")
    z <- matrix(paste(x, x.colpct), nrow=nrow(x),
                dimnames=dimnames(x))
    return(z)
}

testfolds.bystratum <- function(stratum, y, nfold) {
    ## keep only strata that contain a single case 
    table.strata <- tapply(y, stratum, sum) == 1
    strata.onecase <- as.integer(names(table.strata)[as.integer(table.strata)==1])
    keep <- stratum %in% strata.onecase

    stratum.unique <- unique(stratum[keep])
    N <- length(stratum.unique)
    test.folds <- data.frame(stratum=stratum.unique,
                             test.fold=1 + sample(1:N, size=N) %% nfold)
    return(test.folds)
}
