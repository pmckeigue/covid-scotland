
## helperfunctions for COVID analyses

tabulate.bnfsection <- function(chnum) {
    ## get sectioncode in wide format, one col per subchapter
    scrips.ch.wide <- reshape2::dcast(scrips[scrips$chapternum==chnum, ],
                                         ANON_ID ~ sectioncode, fun.aggregate=length,
                                         value.var="sectioncode")
    colnames(scrips.ch.wide)[-1] <- paste0("section.",
                                              as.integer(colnames(scrips.ch.wide)[-1]))
    ## drop rare sectioncodes
    scrips.ch.wide <- scrips.ch.wide[, colSums(scrips.ch.wide) > 20]
    
    cc.bnf.ch <- merge(cc.severe[, c("ANON_ID", "CASE", "stratum")],
                       scrips.ch.wide,
                       by="ANON_ID", all.x=TRUE)
    ## now fix colnames and set missing to 0
    bnfcols <- grep("^section.", colnames(cc.bnf.ch))
    for(j in bnfcols) {
        cc.bnf.ch[, j][is.na(cc.bnf.ch[, j])] <- 0
        cc.bnf.ch[, j][cc.bnf.ch[, j] > 1] <- 1
        cc.bnf.ch[, j] <- as.factor(cc.bnf.ch[, j])
    }
    
    bnf.ch <- colnames(cc.bnf.ch)[-(1:3)]
    ## FIXME: regressions should be fixed to use only rows retained by univariate.tabulate
    table.bnf.ch <- univariate.tabulate(varnames=bnf.ch, outcomevar="CASE", data=cc.bnf.ch,
                                        drop.sparserows=TRUE)
    if(nrow(table.bnf.ch) >= 1) { 
        bnf.ch <- rownames(table.bnf.ch)  
        subsectioncodes <- as.integer(gsub("section.", "", rownames(table.bnf.ch)))
        rownames(table.bnf.ch) <-
            bnfcodes$sectionname[match(subsectioncodes, bnfcodes$sectioncode)]
        univariate.bnf.ch <- NULL
        for(i in 1:length(bnf.ch)) {
            univariate.formula <- as.formula(paste("CASE ~ ", bnf.ch[i], " + strata(stratum)"))
            x <- summary(clogit(formula=univariate.formula, data=cc.bnf.ch))$coefficients
            univariate.bnf.ch <- rbind(univariate.bnf.ch, x)
        }
        
        table.bnf.ch.aug <- combine.tables2(table.bnf.ch, univariate.bnf.ch)
        return(table.bnf.ch.aug) 
    } else {
        return(NULL)
    }
}

tabulate.icdsubchapter <- function(chapternum) {
    ## get chapternum in wide format, one col per subchapter
    diagnoses.ch.wide <- reshape2::dcast(diagnoses[diagnoses$chapter==chapternum, ],
                                         ANON_ID ~ subchapter, fun.aggregate=length,
                                         value.var="subchapter")
    colnames(diagnoses.ch.wide)[-1] <- paste0("subCh.",
                                              as.integer(colnames(diagnoses.ch.wide)[-1]))
    ## drop rare subchapters
    diagnoses.ch.wide <- diagnoses.ch.wide[, colSums(diagnoses.ch.wide) > 20]
    
    cc.icd.ch <- merge(cc.severe[, c("ANON_ID", "CASE", "stratum")],
                       diagnoses.ch.wide,
                       by="ANON_ID", all.x=TRUE)
    ## now fix colnames and set missing to 0
    icdcols <- grep("^subCh.", colnames(cc.icd.ch))
    for(j in icdcols) {
        cc.icd.ch[, j][is.na(cc.icd.ch[, j])] <- 0
        cc.icd.ch[, j][cc.icd.ch[, j] > 1] <- 1
        cc.icd.ch[, j] <- as.factor(cc.icd.ch[, j])
    }
    
    icd.ch <- colnames(cc.icd.ch)[-(1:3)]
    ## FIXME: this function should be fixed to drop rows with small numbers
    ## regressions should be fixed to use only rows retained by univariate.tabulate
    table.icd.ch <- univariate.tabulate(varnames=icd.ch, outcomevar="CASE", data=cc.icd.ch,
                                        drop.sparserows=TRUE)
    #colnames(table.icd.ch) <- c("Controls", "Cases")
    #colnames(table.icd.ch) <- paste(colnames(table.icd.ch),
    #                                gsub("([0-9]+)", "\\(N = \\1\\)",
    #                                     as.integer(table(cc.icd.ch$CASE))))
    ## drop rows with small numbers of controls or cases
    # keep <- table.icd.ch[, 1] + table.icd.ch[, 2] > 10
    # table.icd.ch <- table.icd.ch[keep, , drop=FALSE]
    if(nrow(table.icd.ch) >= 1) { 
        icd.ch <- rownames(table.icd.ch)  
        subchapternums <- as.integer(gsub("subCh.", "", rownames(table.icd.ch)))
        rownames(table.icd.ch) <- icdsubchapters$name[subchapternums]
        univariate.icd.ch <- NULL
        for(i in 1:length(icd.ch)) {
            univariate.formula <- as.formula(paste("CASE ~ ", icd.ch[i], " + strata(stratum)"))
            x <- summary(clogit(formula=univariate.formula, data=cc.icd.ch))$coefficients
            univariate.icd.ch <- rbind(univariate.icd.ch, x)
        }
        
        table.icd.ch.aug <- combine.tables2(table.icd.ch, univariate.icd.ch)
        return(table.icd.ch.aug) 
    } else {
        return(NULL)
    }
}

combine.tables3 <- function(ftable, utable, mtable)  {# returns single table from freqs, univariate, multivariate 
    
    u.ci <- or.ci(utable[, 1], utable[, 3]) 
    u.pvalue <- signif(utable[, 5], 1)
    # pvalue <- pvalue.latex(pvalue)
    
    mult.ci <- or.ci(mtable[, 1], mtable[, 3])
    mult.pvalue <- signif(mtable[, 5], 1)
    
    table.aug <- data.frame(ftable,
                            u.ci, u.pvalue,
                            mult.ci, mult.pvalue)
    return(table.aug)
}    

combine.tables2 <- function(ftable, utable, mtable)  {# returns single table from freqs, univariate 
    
    u.ci <- or.ci(utable[, 1], utable[, 3]) 
    u.pvalue <- signif(utable[, 5], 1)
    # pvalue <- pvalue.latex(pvalue)
    
   table.aug <- data.frame(ftable,
                            u.ci, u.pvalue)
    return(table.aug)
}    

pvalue.latex <- function(x, n=1, nexp=1) {
  sapply(x, function(z) {
    if (is.na(z)) {
        return(NA)
    } else if(z > 0.001) {
        return(signif(z, 1))
    } else {
    z <- sprintf("%.*E", 0, signif(z, n)) # default is 1 sig fig
    z <- as.numeric(unlist(strsplit(z, "E"))) # split z at E
    sprintf("\\ensuremath{%.*f\\times 10^{%0*d}}", 0, z[1], nexp, z[2])
    }
  })
}
## create and format a 95% confidence interval around an odds/hazard ratio
or.ci <- function(coeff, se, ndigits=2) {
  ci.lo <- exp(coeff - 1.96 * se)
  ci.up <- exp(coeff + 1.96 * se)
  sprintf("%.*f (%.*f, %.*f)", ndigits, exp(coeff), ndigits, ci.lo, ndigits, ci.up)
}

univariate.tabulate <- function(varnames, outcomevar, data, drop.referencelevel=TRUE,
                                drop.sparserows=FALSE) {
    outcome <- data[, match(outcomevar, names(data))]
    table.varnames <- NULL
    for(i in 1:length(varnames)) {
        keep.x <- TRUE
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
            keep.x <- !any(rowSums(x) < 10)
            x <- paste.colpercent(x)
            ## rownames are labelled with levels(varname)
            ## if two levels OR drop.reference level, drop reference level
            ## if single row left, label rows with varname
            ## else keep reference level
            if(nrow(x) ==2  | drop.referencelevel) {
                x <- x[-1, , drop=FALSE] # drop reference category
                if(nrow(x) == 1)
                    rownames(x) <- varnames[i]
            } 
       
            ## FIXME: should fix univariate and multivariate regression outputs to
            ## work with drop.referencelevel=FALSE
        }        
        if(!drop.sparserows | keep.x) { # exclude this table
            table.varnames <- rbind(table.varnames, x)
        }
    }
    if(ncol(table.varnames) ==2) {
        colnames(table.varnames) <- c("Controls", "Cases")
        colnames(table.varnames) <- paste0(colnames(table.varnames),
                                           " (", as.integer(table(outcome)), ")")
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
