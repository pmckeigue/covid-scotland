## analysis script for case-control study

replace.names <- function(varnames) {
    names <- varnames
    found <- names %in% lookup.names$varname
    names[found] <- as.character(lookup.names$longname[match(names[found],
                                                             lookup.names$varname)])
    return(names)
}

traintest.fold <- function(train.data, test.data, start.model,
                           lower.formula, upper.formula) { # stepwise clogit on training fold, predicted probs on test fold
    start.model <- clogit(formula=upper.formula, data=train.data)
    stepwise.model <- step(start.model,
                           scope=list(lower=lower.formula, upper=upper.formula),
                           direction="both", trace=-1, method="approximate")
    ## use reference=sample, then normalize separately within each stratum
    unnorm.p <- predict(object=stepwise.model, newdata=test.data,
                        na.action="na.pass", 
                        type="risk", reference="sample")
    norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                            stratum=test.data$stratum,
                                            y=test.data$CASE)
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
    ## also compute prior probs
    ## this is appropriate if there is only 1 case per stratum
    unnorm.p <- as.numeric(unnorm.p)
    stratum <- as.factor(as.integer(as.character(stratum)))
    nonmissing <- !is.na(unnorm.p)

    unnorm.p <- unnorm.p[nonmissing]
    stratum <- stratum[nonmissing]
    y <- y[nonmissing]

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

library(car)
library(survival)
library(MASS)
library(wevid)
library(rmarkdown)
library(pander)
library(ggplot2)
        
cc.all <- readRDS("./cc_linked_file_20200417_onomap.rds")
names(cc.all) <- gsub("CASE_NO", "stratum", names(cc.all))
names(cc.all) <- gsub("^SEX$", "sex", names(cc.all))
names(cc.all) <- gsub("imumune", "immune", names(cc.all))
names(cc.all) <- gsub("^ethnic$", "ethnic.old", names(cc.all))
names(cc.all) <- gsub("^CAREHOME$", "care.home", names(cc.all))

cc.all$stratum <- as.integer(cc.all$stratum)
simd.missing <- is.na(cc.all$simd2020_sc_quintile)
cc.all$SIMD_quintile <- as.factor(cc.all$simd2020_sc_quintile)
#cc.all$simd2020_sc_quintile[simd.missing] <- NA

cc.all$PIS <- as.factor(cc.all$PIS)
cc.all$SMR01 <- as.factor(cc.all$SMR01)

######## coding ethnicity ##############################

source("ethnic_assign.R")

## variable ETHNIC_smr is based on SMR, collapsed to 5 categories
## variable eth5 is based on Onomap only
## classifies all those with non-local Muslim names as South Asian, and misclassifies most  Black Caribbean as White

cc.all$ethnic <- as.factor(cc.all$eth5) 

####################################################################

cc.all$sex <- recode(as.factor(cc.all$sex), "1='Male'; 2='Female'")
cc.all <- within(cc.all, sex <- relevel(sex, ref="Female"))

cc.all$care.home <- recode(as.factor(cc.all$care.home), "0='Independent'; 1='Care home resident'")
cc.all <- within(cc.all, care.home <- relevel(care.home, ref="Independent"))

## assign controls to same group as cases i.e. A, B, C
casegroups <- cc.all[cc.all$CASE==1, ][, c("stratum", "group")]
colnames(casegroups)[2]<-"casegroup"
cc.all <- merge(cc.all, casegroups, by=c("stratum"), all.x=T)
table(cc.all$CASE, cc.all$casegroup)

with(cc.all[cc.all$CASE==0, ], table(casegroup, dead28, exclude=NULL))
with(cc.all[cc.all$CASE==1, ], table(casegroup, dead28, exclude=NULL))

with(cc.all[cc.all$CASE==1 & cc.all$icu.hdu.ccu==1, ], table(ethnic, dead28, exclude=NULL))

                                   
cc.all$casegroup <- recode(cc.all$casegroup,
                           "'A'='Critical care or fatal'; 'B'='Hospitalised, not severe'; 'C'='Test-positive, not hospitalised'")
cc.all$casegroup <- as.factor(cc.all$casegroup)
cc.all <- within(cc.all, casegroup <- relevel(casegroup, ref="Critical care or fatal"))

## recode diabetes type
cc.all$dm.type <- as.integer(cc.all$dm.type)
## missing recoded as zero
cc.all$dm.type[is.na(cc.all$dm.type)] <- 0
cc.all$dm.type <- 
  car::recode(cc.all$dm.type, 
              "0='Not diabetic'; 1='Type 1'; 2='Type 2'; 3:hi='Other DM'")

## anticoagulants6 same as anticoagulants_protamine6

cc.all$kidney.any <- as.integer(cc.all$kidney.advanced==1 | cc.all$kidney.other==1)

cc.all$asthma.any <- as.integer(cc.all$asthma12.bnf==1 | cc.all$asthma==1)

cc.all$cysticfibrosis.any <- as.integer(cc.all$cystic_fibrosis12.bnf==1 |
                                        cc.all$cystic.fibrosis==1)

cc.all$tuberculosis.any <- as.integer(cc.all$tuberculosis==1 |
                                      cc.all$tuberculosis12.bnf==1)

cc.all$otherresp.any <- as.integer(cc.all$resp.other==1 |
                                  cc.all$other_chronic_lrd12.bnf)

cc.all$respinf.orTB <- as.integer(cc.all$resp.inf==1 |
                                  cc.all$tuberculosis.any)

cc.all$HIV.any <- as.integer(cc.all$hiv12.bnf==1 |
                             cc.all$HIV==1)

cc.all$immune.any <- as.integer(cc.all$immune==1 |
                                cc.all$transplant.not.kidney==1 |
                                cc.all$cytotoxic_immune6.bnf==1 |
                                cc.all$HIV.any)

cc.all$IHD.any <- as.integer(cc.all$ihd12.bnf==1 | cc.all$IHD==1)

cc.all$connectivetissue.any <- as.integer(cc.all$connective_tissue12.bnf==1 | cc.all$connective==1)

cc.all$otherneuro.any <- as.integer(cc.all$ms12.bnf==1 | cc.all$neuro.other==1)

cc.all$antihypertensive.any <- with(cc.all,
                                    as.integer(vasodilator_ah6.bnf==1 |
                                               centrally_acting_ah6.bnf==1 |
                                               adrenergic_neurone_block6.bnf==1 |
                                               alpha_adrenoceptor6.bnf==1 |
                                               ace6.bnf==1 |
                                               angio6.bnf==1 |
                                               renin_angiotensin6.bnf==1 |
                                               thiazides6.bnf==1 |
                                               calcium_channel6.bnf==1))
cc.all$antihypertensive.other <- with(cc.all,
                                      as.integer(vasodilator_ah6.bnf==1 |
                                                 centrally_acting_ah6.bnf==1))

antihypertensive.classes <- c("vasodilator_ah6.bnf", 
                              "centrally_acting_ah6.bnf",  
                              "adrenergic_neurone_block6.bnf",  
                              "alpha_adrenoceptor6.bnf", 
                              "ace6.bnf",  
                              "angio6.bnf",  
                              "renin_angiotensin6.bnf",  
                              "thiazides6.bnf",  
                              "calcium_channel6.bnf")

antihypertensives <- c(antihypertensive.classes[4:9], "antihypertensive.other")

demog <- c("ethnic", "simd2020_sc_decile", "care.home", "PIS", "SMR01")
demog.smr <- c("ETHNIC_smr", "simd2020_sc_decile", "care.home")

conditions <- c("dm.type", "antihypertensive.any", "IHD.any", "CVD", "heart.other", "circulatory.other",
                "asthma.any", "respinf.orTB", "resp.other",
                "cysticfibrosis.any",
                "kidney.any", "connectivetissue.any", 
                "epilepsy", "mono.poly.neuro", "otherneuro.any",
                "blood.cancer", "lung.cancer", "other.cancer",
                "immune.any") # , "HIV.any")

drugs <- eval(c("antihypertensive.any",
                antihypertensives, # "antihypertensive.other",
                "anticoagulants6.bnf", "nsaids6.bnf",
                "lipid_regulating6.bnf", "statins6.bnf",
                "hydroxychloroquine6.bnf"))

admissions <- "SMR01"

lookup.names <- data.frame(varname=c("PIS", "SMR01", "simd2020_sc_decile",
                                     "ace6.bnf", "angio6.bnf"),
                           longname=c("Any prescription", "Any admission", "SIMD decile",
                                      "ACE inhibitor", "Angiotensin-II receptor blocker"))

########################################################################################

## tabulate test positives

table.ethnic <- table(cc.all$eth5, cc.all$ETHNIC_smr, exclude=NULL)
table.ethnic <- paste.colpercent(table.ethnic, 1)

testpositives.ethnic <- paste.colpercent(with(cc.all[cc.all$CASE==1, ], table(ethnic, casegroup)), 1)
testpositives.ethnic.smr <- paste.colpercent(with(cc.all[cc.all$CASE==1, ], table(ETHNIC_smr, casegroup)), 1)

#####################################################################################
cc.testpositives <- cc.all[cc.all$CASE==1 & !is.na(cc.all$ethnic), ]
sex.demog <- c("sex", demog) 
table.testpositives <-
    with(cc.testpositives,
     tapply(AGE,
            casegroup,
            function(x) return(paste0(median(x),
                                      " (", quantile(x, probs=0.25),
                                      "-", quantile(x, probs=0.75), ")"))))
for(i in 1:length(sex.demog)) {
    x <- table(cc.testpositives[, match(sex.demog[i], names(cc.testpositives))],
               cc.testpositives$casegroup)
    x <- paste.colpercent(x, 1)
    x <- x[-1, , drop=FALSE]
    if(nrow(x)==1) {
        rownames(x) <- sex.demog[i]
    }
    table.testpositives <- rbind(table.testpositives, x)
}
colnames(table.testpositives) <- paste(colnames(table.testpositives),
                                       gsub("([0-9]+)", "\\(N = \\1\\)",
                                            as.integer(table(cc.testpositives$casegroup))))
rownames(table.testpositives)[1] <- "Age"
rownames(table.testpositives)[2] <- "Male"

sex.demog.smr <- c("sex", demog.smr) 
cc.testpositives <- cc.all[cc.all$CASE==1 & !is.na(cc.all$ETHNIC_smr), ]
table.testpositives.smr <-
    with(cc.testpositives, 
     tapply(AGE,
            casegroup,
            function(x) return(paste0(median(x),
                                      " (", quantile(x, probs=0.25),
                                      "-", quantile(x, probs=0.75), ")"))))
for(i in 1:length(sex.demog.smr)) {
    x <- table(cc.testpositives[, match(sex.demog.smr[i], names(cc.testpositives))],
               cc.testpositives$casegroup)
    x <- paste.colpercent(x, 1)
    x <- x[-1, , drop=FALSE]
    if(nrow(x)==1) {
        rownames(x) <- sex.demog[i]
    }
    table.testpositives.smr <- rbind(table.testpositives.smr, x)
}

colnames(table.testpositives.smr) <- paste(colnames(table.testpositives.smr),
                                       gsub("([0-9]+)", "\\(N = \\1\\)",
                                            as.integer(table(cc.testpositives$casegroup))))
                                                
rownames(table.testpositives.smr)[1] <- "Age"
rownames(table.testpositives.smr)[2] <- "Male"

########### restrict to severe cases and matched controls ###################### 
 
cc.severe <- cc.all[cc.all$casegroup=="Critical care or fatal", ]
                                        #cc.severe <- cc.all[cc.all$casegroup=="Severe" & cc.all$AGE < 75, ]

## case-control analysis for ethnicity 

cc.ethnic <- paste.colpercent(table(cc.severe$ethnic, cc.severe$CASE), 1)
univariate.ethnic <- summary(clogit(formula=CASE ~ ethnic + strata(stratum),
                                    data=cc.severe))$coefficients
univariate.ethnic <- rbind(rep(NA, 2), univariate.ethnic[, c(2, 5)])
cc.ethnic <- data.frame(cc.ethnic, univariate.ethnic)
colnames(cc.ethnic)[1:2] <- paste(c("Controls", "Cases"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(cc.severe$CASE[!is.na(cc.severe$ethnic)]))))
colnames(cc.ethnic)[3:4] <- c("Rate ratio", "p-value")

## case-control analysis for SMR ethnicity
cc.ethnic.smr <- paste.colpercent(table(cc.severe$ETHNIC_smr, cc.severe$CASE), 1)
univariate.ethnic.smr <- summary(clogit(formula=CASE ~ ETHNIC_smr + strata(stratum),
                                        data=cc.severe))$coefficients
univariate.ethnic.smr <- rbind(rep(NA, 2), univariate.ethnic.smr[, c(2, 5)])
cc.ethnic.smr <- data.frame(cc.ethnic.smr, univariate.ethnic.smr)
colnames(cc.ethnic.smr)[1:2] <- paste(c("Controls", "Cases"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(cc.severe$CASE[!is.na(cc.severe$ETHNIC_smr)]))))
colnames(cc.ethnic.smr)[3:4] <- c("Rate ratio", "p-value")

cc.severe$ethnic <- recode(cc.severe$ethnic, "'Black'='Other'; 'Chinese'='Other'")
cc.severe <- within(cc.severe, ethnic <- relevel(as.factor(ethnic), ref="White"))

## univariate associations
table.demog <- NULL
for(i in 1:length(demog)) {
    ## test whether variable is factor or numeric
    z <- cc.severe[, match(demog[i], names(cc.severe))] 
    if(is.numeric(z)) {
        x <- tapply(z, cc.severe$CASE,
                    function(x) return(paste0(median(x, na.rm=TRUE),
                                              " (", quantile(x, probs=0.25, na.rm=TRUE),
                                              "-", quantile(x, probs=0.75, na.rm=TRUE), ")")))
        x <- matrix(x, nrow=1)
        rownames(x) <- demog[i]
    } else {
        x <- table(z, cc.severe$CASE)
        x <- paste.colpercent(x)
        x <- x[-1, , drop=FALSE] # drop reference category
        ## FIXME: should keep this line and pad univariate
    }
    if(nrow(x)==1) {
        rownames(x) <- demog[i]
    }
    table.demog <- rbind(table.demog, x)
}
colnames(table.demog) <- c("Controls", "Cases")
colnames(table.demog) <- paste(colnames(table.demog),
                                       gsub("([0-9]+)", "\\(N = \\1\\)",
                                            as.integer(table(cc.severe$CASE))))

univariate.demog <- NULL
for(i in 1:length(demog)) {
    univariate.formula <- as.formula(paste("CASE ~ ", demog[i], "+ strata(stratum)"))
    x <- summary(clogit(formula=univariate.formula, data=cc.severe))$coefficients
    univariate.demog <- rbind(univariate.demog, x)
}

############################ conditions #################

table.conditions <- NULL
for(i in 1:length(conditions)) {
    x <- table(cc.severe[, match(conditions[i], names(cc.severe))],
               cc.severe$CASE)
    x <- paste.colpercent(x)
    x <- x[-1, , drop=FALSE]
    if(nrow(x)==1) {
        rownames(x) <- conditions[i]
    }
    table.conditions <- rbind(table.conditions, x)
}
colnames(table.conditions) <- c("Controls", "Cases")
colnames(table.conditions) <- paste(colnames(table.conditions),
                                       gsub("([0-9]+)", "\\(N = \\1\\)",
                                            as.integer(table(cc.severe$CASE))))

univariate.conditions <- NULL
for(i in 1:length(conditions)) {
    univariate.formula <- as.formula(paste("CASE ~ ", conditions[i], "+ strata(stratum)"))
    x <- summary(clogit(formula=univariate.formula, data=cc.severe))$coefficients
    univariate.conditions <- rbind(univariate.conditions, x)
}
    
################# drugs ######################################
table.drugs <- NULL
for(i in 1:length(drugs)) {
    x <- table(cc.severe[, match(drugs[i], names(cc.severe))],
               cc.severe$CASE)
    x <- paste.colpercent(x)
    x <- x[-1, , drop=FALSE]
    if(nrow(x)==1) {
        rownames(x) <- drugs[i]
    }
    table.drugs <- rbind(table.drugs, x)
}
colnames(table.drugs) <- c("Controls", "Cases")
colnames(table.drugs) <- paste(colnames(table.drugs),
                                       gsub("([0-9]+)", "\\(N = \\1\\)",
                                            as.integer(table(cc.severe$CASE))))

univariate.drugs <- NULL
for(i in 1:length(drugs)) {
    univariate.formula <- as.formula(paste("CASE ~ ", drugs[i], "+ strata(stratum)"))
    x <- summary(clogit(formula=univariate.formula, data=cc.severe))$coefficients
    univariate.drugs <- rbind(univariate.drugs, x)
}

## fit model and evaluate prediction: demog only
demog.formula <- as.formula(paste("CASE ~",
                                  paste(demog, collapse=" + "),
                                  " + strata(stratum)"))
cc.severe.demog.nonmissing <- nonmissing.obs(cc.severe, (demog))
demog.model <- clogit(formula=demog.formula, data=cc.severe.demog.nonmissing,
                                 method="approximate")
stepwise.demog <- step(demog.model, direction="both",
                       trace=2)

unnorm.p <- predict(object=demog.model, newdata=cc.severe.demog.nonmissing,
                    na.action="na.pass", 
                    type="risk", reference="strata", )
demog.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                         stratum=cc.severe.demog.nonmissing$stratum,
                                        y=cc.severe.demog.nonmissing$CASE)
demog.densities <- with(demog.predicted, Wdensities(y, posterior.p, prior.p))
pander(summary(demog.densities), table.style="multiline",
       split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
       caption="Prediction of severe COVID-19 from demographic variables")

## fit model with demog + conditions
conditions.formula <- as.formula(paste("CASE ~",
                                   paste(demog, collapse=" + "), "+",
                                   paste(conditions, collapse=" + "), 
                                   "+ strata(stratum)"))
cc.severe.conditions.nonmissing <- nonmissing.obs(cc.severe, c(demog, conditions))
conditions.model <- clogit(formula=conditions.formula,
                                      data=cc.severe.conditions.nonmissing,
                                      method="approximate")
stepwise.conditions <- step(conditions.model, direction="both", trace=2)

unnorm.p <- predict(object=conditions.model, newdata=cc.severe.conditions.nonmissing,
                    na.action="na.pass", 
                    type="risk", reference="strata", )
norm.predicted <- normalize.predictions(unnorm.p=unnorm.p,
                                        stratum=cc.severe.conditions.nonmissing$stratum,
                                        y=cc.severe.conditions.nonmissing$CASE)
conditions.densities <- with(norm.predicted, Wdensities(y, posterior.p, prior.p))
pander(summary(conditions.densities), table.style="multiline",
       split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
       caption="Prediction of severe COVID-19 from demographic + diagnoses")

############# fit full model #####################################
lower.formula <- as.formula("CASE ~ care.home + strata(stratum)")
all.formula <- as.formula(paste("CASE ~",
                                paste(demog, collapse=" + "), "+",
                                paste(conditions, collapse=" + "), "+", 
                                paste(drugs, collapse=" + "), 
                                "+ strata(stratum)"))
cc.severe.all.nonmissing <- nonmissing.obs(cc.severe, c(demog, conditions, drugs,
                                                        "PIS", "SMR01"))
all.model <- clogit(formula=all.formula, data=cc.severe.all.nonmissing)
stepwise.all <- step(all.model,
                     scope=list(lower=lower.formula, upper=all.formula),
                     direction="both", method="approximate", trace=2)
stepwise.all <- summary(stepwise.all)$coefficients

rownames(stepwise.all) <- replace.names(rownames(stepwise.all))
print(stepwise.all)

## cross-validation of stepwise procedure
cv.data <- cc.severe.all.nonmissing
stratum.unique <- unique(cv.data$stratum)
N <- length(stratum.unique)
nfold <- 4
test.folds <- data.frame(stratum=stratum.unique,
                         test.fold=1 + sample(1:N, size=N) %% nfold)
cv.data <- merge(cv.data, test.folds, by="stratum")

## loop over test folds to obtain predictions from corresponding training folds
norm.predicted <- NULL
for(i in 1:nfold) {
    test.data <- cv.data[cv.data$test.fold==i, ]
    train.data <- cv.data[cv.data$test.fold != i, ]
    norm.predicted <- rbind(norm.predicted,
                            traintest.fold(train.data=train.data,
                                           test.data=test.data,
                                           start.model=start.model,
                                           lower.formula=demog.formula,
                                           upper.formula=all.formula))
}

stepwise.densities <- with(norm.predicted, Wdensities(y, posterior.p, prior.p))

pander(summary(stepwise.densities), table.style="multiline",
       split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
       caption="Prediction of severe COVID-19 from demographic + diagnoses")

plotWdists(stepwise.densities) + theme_grey(base_size = 18) + theme(legend.title=element_blank())
