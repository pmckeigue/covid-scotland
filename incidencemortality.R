
## https://www.nrscotland.gov.uk/statistics-and-data/statistics/statistics-by-theme/population/population-estimates/mid-year-population-estimates/mid-2019

## https://www.nrscotland.gov.uk/files//statistics/population-estimates/mid-19/mid-year-pop-est-2019-figures.xlsx

scotpop <- read_excel("./Scotland_midyearpop_est2019.xlsx")

# load("casefreqs.RData") ## national cumulative cases and deaths by sex and one year age group

####### incidence and mortality using national population estimates ######################

narrow.case.freqs <- data.frame(Age=as.integer(rownames(narrow.case.freqs)),
                         Females=as.integer(narrow.case.freqs[, 1]),
                         Males=as.integer(narrow.case.freqs[, 2]))
case.long <- reshape2::melt(narrow.case.freqs, id="Age")
colnames(case.long) <- c("Age", "Sex", "Cases")

narrow.death.freqs <- data.frame(Age=as.integer(rownames(narrow.death.freqs)),
                         Females=as.integer(narrow.death.freqs[, 1]),
                         Males=as.integer(narrow.death.freqs[, 2]))
death.long <- reshape2::melt(narrow.death.freqs, id="Age")
colnames(death.long) <- c("Age", "Sex", "Deaths")

scotpop.long <- reshape2::melt(scotpop[, -2], id="Age")
colnames(scotpop.long) <- c("Age", "Sex", "Population")

discrim <- merge(scotpop.long, case.long, by=c("Age", "Sex"), all.x=TRUE)
discrim <- merge(discrim, death.long, by=c("Age", "Sex"), all.x=TRUE)
discrim$Cases[is.na(discrim$Cases)] <- 0
discrim$Deaths[is.na(discrim$Deaths)] <- 0

discrim$Sex <- as.factor(discrim$Sex)

discrim$Noncases <- discrim$Population - discrim$Cases
y.cases <- cbind(as.integer(discrim$Cases), as.integer(discrim$Noncases))

discrim$Survivors <- discrim$Population - discrim$Deaths
y.deaths <- cbind(as.integer(discrim$Deaths), as.integer(discrim$Survivors))

cases.model <- glm(formula=y.cases ~ Sex + Age, family="binomial", data=discrim)
deaths.model <- glm(formula=y.deaths ~ Sex + Age, family="binomial", data=discrim)

cases.model.coeffs <- summary(cases.model)$coefficients
deaths.model.coeffs <- summary(deaths.model)$coefficients
logistic.coeffs <- data.frame(severecase=cases.model.coeffs[, 1],
                              death=deaths.model.coeffs[, 1])

male <- discrim$Sex=="Males"
female <- discrim$Sex=="Females"
gam.model.MaleDeaths <- gam::gam(formula=y.deaths[male, ] ~ s(Age), family=binomial("logit"),
                                 data=discrim[male, ])
gam.model.FemaleDeaths <- gam::gam(formula=y.deaths[female, ] ~ s(Age), family=binomial("logit"),
                                   data=discrim[female, ])
gam.model.MaleCases<- gam::gam(formula=y.cases[male, ] ~ s(Age), family=binomial("logit"),
                               data=discrim[male, ])
gam.model.FemaleCases <- gam::gam(formula=y.cases[female, ] ~ s(Age), family=binomial("logit"),
                                  data=discrim[female, ])

gam.male <- data.frame(Cases=logit(gam.model.MaleCases$fitted.values),
                       Deaths=logit(gam.model.MaleDeaths$fitted.values),
                       Age=discrim$Age[male])
gam.male.long <- reshape2::melt(data=gam.male, id="Age")
colnames(gam.male.long)[2] <- "Status"
gam.male.long$Sex <- "Males"

gam.female <- data.frame(Cases=logit(gam.model.FemaleCases$fitted.values),
                       Deaths=logit(gam.model.FemaleDeaths$fitted.values),
                       Age=discrim$Age[female])
gam.female.long <- reshape2::melt(data=gam.female, id="Age")
colnames(gam.female.long)[2] <- "Status"
gam.female.long$Sex <- "Females"
gam <- rbind(gam.male.long, gam.female.long)

###############################################################

logodds.posterior <- predict(object=cases.model, newdata=discrim, type="link")
logodds.prior <- log(sum(discrim$Cases) / sum(discrim$Noncases))
log.likratio <- logodds.posterior - logodds.prior
discrim$W <- log.likratio / log(2)
lambda1 <- sum(discrim$W * discrim$Cases) / sum(discrim$Cases)
lambda0 <- sum(-discrim$W * discrim$Noncases) / sum(discrim$Noncases)
cases.Lambda.agesex <- 0.5 * (lambda0 +  lambda1)


logodds.posterior <- predict(object=deaths.model, newdata=discrim, type="link")
logodds.prior <- log(sum(discrim$Deaths) / sum(discrim$Survivors))
log.likratio <- logodds.posterior - logodds.prior
discrim$W <- log.likratio / log(2)
lambda1 <- sum(discrim$W * discrim$Deaths) / sum(discrim$Deaths)
lambda0 <- sum(-discrim$W * discrim$Survivors) / sum(discrim$Survivors)
deaths.Lambda.agesex <- 0.5 * (lambda0 +  lambda1)


get.covidage <- function(popmodel, sex, logitrisk) {
    ## get interpolated age given logitrisk and sex
    ## popmodel has columns sex, age, logit
    N <- length(sex)
    covidage <- numeric(N)
    model.m <- subset(popmodel, sex=="Males")
    model.f <- subset(popmodel, sex=="Females")
    covidage[sex=="Males"] <- with(model.m, approx(x=logit,  y=age, xout=logitrisk)$y)
    covidage[sex=="Females"] <- with(model.f, approx(x=logit,  y=age, xout=logitrisk)$y)
    covidage[logitrisk < min(popmodel$logit)] <- 0
    covidage[logitrisk > max(popmodel$logit)] <- max(popmodel$age)
    return(covidage)
}

getlogitrisk <- function(popmodel, sex, age) {
    ## get logit risk given sex and age, from population model
    age.xout <- age
    N <- length(age)
    logitrisk <- numeric(N)
    model.m <- subset(popmodel, sex=="Males")
    model.f <- subset(popmodel, sex=="Females")
    logitrisk[sex=="Males"] <- with(model.m, approx(x=model.m$age,  y=logit, xout=age.xout)$y)
    logitrisk[sex=="Females"] <- with(model.f, approx(x=model.f$age,  y=logit,
                                                      xout=age.xout)$y)
    return(logitrisk)
}

calibrate.condlogit <- function(popmodel, sex, age, c.logit) {
    ## find value of intercept that equates observed and expected risk
    ## where the linear predictor is the sum of the logit from popmodel and
    ## conditional logit c.logit from case-control study

    ## left join data.table(sex, age, clogit) with popmodel

    ## sum the logits

    ## equate observed and predicted risk over all individuals
    
}

popmodel <- subset(gam, Status=="Deaths", select=c("Sex", "Age", "value"))
colnames(popmodel) <- c("sex", "age", "logit")
      
covidage <- get.covidage(popmodel=popmodel, sex="Males", logitrisk=-10)

print(covidage)

logitrisk <- getlogitrisk(popmodel=popmodel, sex="Males", age=covidage)

print(logitrisk)

## function to calibrate the sum of the logit from the population model and the logit from the conditional logistic regression model by equating observed and expected for each sex 


