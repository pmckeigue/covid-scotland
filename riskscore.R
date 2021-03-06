library(data.table)
library(readxl)
library(mgcv) # gam
library(ggplot2)

source("helperfunctions.R")

datadir <- "./data/2021-02-18/"

load(paste0(datadir, "casefreqs.4cats.agesex.RData"))

scotpop <- as.data.table(read_excel("./Scotland_midyearpop_est2019.xlsx"))

## fit gam model for age and sex to narrow.case.freqs
narrow.case.freqs <- data.table(Age=as.integer(rownames(narrow.case.freqs)),
                         Females=as.integer(narrow.case.freqs[, 1]),
                         Males=as.integer(narrow.case.freqs[, 2]))
case.long <- melt(narrow.case.freqs, id="Age")
colnames(case.long) <- c("Age", "Sex", "Cases")
case.long <- case.long[Age < 90]

scotpop <- scotpop[Age < 90]
scotpop.long <- melt(scotpop[, -2], id="Age")
colnames(scotpop.long) <- c("Age", "Sex", "Population")
scotpop.long[, Population := as.integer(Population)]

setkey(case.long, Age, Sex)
setkey(scotpop.long, Age, Sex)
discrim <- scotpop.long[case.long]

discrim[, NonCases := Population - Cases]
discrim[, Sex := car::recode(Sex, recode="'Females'='Female'; 'Males'='Male'",
                            as.factor=TRUE)]
y.cases <- as.matrix(discrim[, .(Cases, NonCases)])

gam.model <- mgcv::gam(formula=y.cases ~ s(Age) + s(Age, by=Sex) + Sex,
                      family=binomial("logit"),
                     data=discrim)

gam.logit <- logit(gam.model$fitted.values)

discrim[, gam.logit := ..gam.logit]

dt.lookup <- discrim[, .(Age, Sex, gam.logit)]
setnames(dt.lookup, "Age", "COVID.age")

##################################################################
## plot logistic gam model for severe cases by age and sex

breaks.gam <- c("1/20,000", "1/10,000", "1/5,000",
                "1/2000", "1/1000", "1/500",
                "1/200", "1/100")
theme_set(
    theme_gray(base_size = 12)
)

p.incidence <- ggplot(data=discrim, aes(x=Age, y=gam.logit, colour=Sex)) +
    geom_line() +
    scale_x_continuous(name="Age", limits=c(20, 90),
                       breaks=seq(20, 80, by=20),
                       labels=seq(20, 80, by=20)) + 
    scale_x_continuous(limits=c(20, 75), expand=c(0, 0)) +
    scale_y_continuous(name="Risk (logit scale)",
                       breaks=logit(c(5E-5, 1E-4, 2E-4,
                                           5E-4, 1E-3, 2E-3,
                                           5E-3, 1E-2)),
                       labels=breaks.gam,
                       limits=logit(c(5E-5, 5E-3))) +
    theme(axis.title=element_text(size=12)) +
    theme(legend.position=c(0.85, 0.3))
p.incidence

###########################################

## calibrate risk score to equate observed and predicted cases using age-sex distribution of risk score in controls

load(paste0(datadir, "betadraws.riskscore.RData"))
load(paste0(datadir, "cc.severe.lt75.RData"))

covariate.names.agesex <- c("AGE", "sex", covariate.names) # no need to retain ID
covariate.names.caseagesex <- c("CASE", "AGE", "sex", covariate.names) # no need to retain ID
cc.lt75 <- cc.lt75[, ..covariate.names.caseagesex]

controls <- cc.lt75[CASE==0, ..covariate.names.agesex]
controls <- na.omit(controls)
rm(cc.lt75)
setnames(controls, "AGE", "Age")
setnames(controls, "sex", "Sex")

setkey(controls, Age, Sex)
setkey(discrim, Age, Sex)
controls <- discrim[controls]
controls.predict <- controls[, .(Age, Sex, gam.logit)]

## deselect specimendate, case, stratum
covariateonly.names <- covariate.names[-(1:3)]
X.controls <- model.matrix(object=~ ., data=controls[, ..covariateonly.names])[, -1]
rm(controls)

cat("calculating unnormalized logits ...")
## FIXME -- should calculate means with logsumexp()
#beta.means <- colMeans(beta.draws)
#logit.unnorm <- X.controls %*% matrix(beta.means, ncol=1)
unnorm.lograte.draws <- X.controls %*% t(beta.draws)
rm(X.controls)
unnorm.lograte <- matrixStats::rowLogSumExps(unnorm.lograte.draws) -
    log(ncol(unnorm.lograte.draws))
rm(unnorm.lograte.draws)
controls.predict[, logit.unnorm := ..unnorm.lograte]
cat("done\n")

intercept.gam.unnorm <-  optimize(f=findintercept, interval=c(-10, 5),
                   logit1=controls.predict$gam.logit,
                   logit2=controls.predict$logit.unnorm)$minimum

cat("Intercept set to", intercept.gam.unnorm, "\n")

if(FALSE) {
    controls.predict[, logit.norm := ..alpha + gam.logit + logit.unnorm]
    controls.predict <- controls.predict[, id := .I]
    controls.predict <- controls.predict[, .(id, Age, Sex, vaxnow, logit.norm)]
    controls.covidage <- getCOVIDage(x=controls.predict, y=dt.lookup)
}
