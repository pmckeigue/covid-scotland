get.sexcoding <- function(sexvar) {
    sexcodes <- names(table(sexvar))
    malepos <- grep("^M", sexcodes, ignore.case=TRUE)
    if(malepos==2) sexcodes <- rev(sexcodes)
    return(sexcodes)
}

get.covidage <- function(popmodel, sex, logitrisk) {
    ## get interpolated age given logitrisk and sex
    ## popmodel has columns sex, age, logit
    N <- length(sex)
    covidage <- numeric(N)
    sexcodes <- get.sexcoding(sex)
    m <- sex==sexcodes[1]
    f <- sex==sexcodes[2]
    mpop <- popmodel
    p.sexcodes <- get.sexcoding(mpop$sex)
    colnames(mpop) <- paste0("popmodel.", colnames(mpop))
    covidage[m] <- with(subset(mpop, popmodel.sex==p.sexcodes[1]),
                        approx(x=popmodel.logit,  y=popmodel.age, xout=logitrisk[m],
                               yleft=0, yright=90))$y
    covidage[f] <- with(subset(mpop, popmodel.sex==p.sexcodes[2]),
                        approx(x=popmodel.logit,  y=popmodel.age, xout=logitrisk[f],
                               yleft=0, yright=90)$y)
    return(covidage)
}

get.logitrisk <- function(popmodel, sex, age) {
    ## get u.logit given sex and age, from population model
    N <- length(sex)
    sexcodes <- get.sexcoding(sex)
    logitrisk <- numeric(N)
    m <- sex==sexcodes[1]
    f <- sex==sexcodes[2]
    mpop <- popmodel
    p.sexcodes <- get.sexcoding(mpop$sex)
    colnames(mpop) <- paste0("popmodel.", colnames(mpop))
    logitrisk[m] <- with(subset(mpop, popmodel.sex==p.sexcodes[1]),
                         approx(x=popmodel.age,  y=popmodel.logit, xout=age[m],
                                yleft=0, yright=90))$y
    logitrisk[f] <- with(subset(mpop, popmodel.sex==p.sexcodes[2]),
                         approx(x=popmodel.age,  y=popmodel.logit, xout=age[f],
                                yleft=0, yright=90)$y)
    return(logitrisk)
}

popmodel <- subset(gam, Status=="Deaths", select=c("Sex", "Age", "value"))
colnames(popmodel) <- c("sex", "age", "logit")
      
ominuse.sq <- function(alpha, logit.p, observed.p) { # objective function to be minimized
    logit.p <- alpha + logit.p
    expected.p <- mean(invlogit(logit.p)) # inverts logit function log p / (1-p)
    obj  <-  (observed.p - expected.p)^2  # square it so that the minimum will be zero.
    return(obj)
}

calibrate.clogit <- function(popmodel,observed.p, c.logits) {
    ## find value of intercept that equates observed and predicted risk
    ## where the linear predictor uc.logit is the sum of u.logit from population model and
    ## c.logit from case-control study
    ## c.logits has columns sex, age, logit

    u.logit <- getlogitrisk(popmodel, c.logits$sex, c.logits$age)
    uc.logit <- u.logit + c.logits$logit

    ## equate observed and predicted risk over all individuals
    optim.result <-  optim(par=0, fn=ominuse.sq,
                           logit.p=uc.logit, observed.p=observed.p, method="L-BFGS-B")
    alpha <- optim.result$par # value to be added to logit.p
    cat("observed risk", observed.p,
        "expected risk", mean(invlogit(alpha + uc.logit)), "\n")
    return(alpha)
}

## controls are representative of the population at risk
c.logits <- cc.nonmissing[, .(CASE, age=AGE, sex, logit=demog.model$linear.predictors)]
c.logits <- c.logits[CASE==0 & age < 90]

observed.numdeaths <- sum(discrim$Deaths) ## 3055
observed.p <- observed.numdeaths / nrow(c.logits)

alpha <- calibrate.clogit(popmodel=popmodel, observed.p=observed.p, c.logits=c.logits)

print(alpha)

get.covidage.uc <- function(popmodel, alpha, c.logits) {
    ## get covid age given sex, age and linear predictor conditional on sex and age
    ## covid age is alpha + get.logitrisk(popmodel, sex, age) + get.covidac.logit)
    ## c.logits has columns sex, age, logit
    u.logit <- get.logitrisk(popmodel, c.logits$sex, c.logits$age) 
    uc.logit <- alpha  + u.logit + c.logits$logit
    covidage <- get.covidage(popmodel, c.logits$sex, uc.logit)
    return(covidage)
}

covidage.all <- get.covidage.uc(popmodel=popmodel, alpha=alpha,
                                c.logits=c.logits)

print(summary(covidage.all))
print(summary(c.logits))
      
