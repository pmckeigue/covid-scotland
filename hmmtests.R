library(data.table)
library(survival)
library(rstan)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

datadir <- "./data/2021-09-22/"
load(paste0(datadir, "testscohort.RData"))

mu.recover <- 0.5 # arrival rate of recoveries per week
p_recov <- 1 - exp(-mu.recover) # prob of recovery in 1 week

## replace anon_id with integer id
testscohort <- testscohort[1:20000, ]

## generate vector of cutpoints for calendar time
cutpoints <- seq(min(testscohort$entrydate - 7),
                 max(testscohort$exitdate + 7), by=7)

# 1.2 million records
tsplit <- as.data.table(survSplit(Surv(entrydate, exitdate, testpos) ~ .,
                                  data=testscohort,
                                  cut = cutpoints,
                                  episode="cutpointnum"))

## derive observed 3-state test variable
tsplit[, y := 1]
tsplit[SpecimenDate > entrydate & SpecimenDate <= exitdate & ecoss.result=="Negative", y := 2]
tsplit[SpecimenDate > entrydate & SpecimenDate <= exitdate & ecoss.result=="Positive", y := 3]

## keep highest observed test result in each person-time interval
setorder(tsplit, -y)
tsplit <- unique(tsplit, by=c("id", "cutpointnum"))
setorder(tsplit, id, cutpointnum)
N.intervals <- nrow(tsplit)

## derive time-updated vax status at start of each interval
tsplit[, vax14.1dose := as.integer(!is.na(vaxdate_1) &
                                        entrydate - vaxdate_1 >= 14)]
tsplit[is.na(vax14.1dose), vax14.1dose := 0]
tsplit[, vax14.2dose := as.integer(!is.na(vaxdate_2) &
                                        entrydate - vaxdate_2 >= 14)]
tsplit[is.na(vax14.2dose), vax14.2dose := 0]

tsplit[, vax14.dose := vax14.1dose + vax14.2dose]
tsplit[, vax14.factor := car::recode(vax14.dose,
                                             "0='Unvaccinated'; 1='1 dose'; 2='2 doses'",
                                             as.factor=TRUE,
                                             levels=c("Unvaccinated", "1 dose", "2 doses"))]
tsplit <- na.omit(tsplit[, .(id, cutpointnum, y, age_years, sex, care.home, occup,
                             listedgr3, vax14.factor)])
summary(tsplit)


X.intervals <- model.matrix(object=as.formula(~ age_years + sex + care.home + occup +
                                                  listedgr3 + vax14.factor),
                            data=tsplit)[, -1]

## number the first and last rows for each individual in the data table
setorder(tsplit, id, cutpointnum)
tsplit[, rownum := .I]
tsplitrows <- setorder(tsplit[, .(id, rownum)], rownum)
firstrows <- unique(tsplitrows[, .(id, rownum)], by="id")
setnames(firstrows, "rownum", "firstrow")

tsplitrows <- setorder(tsplit[, .(id, rownum)], -rownum)
lastrows <- unique(tsplitrows[, .(id, rownum)], by="id")
setnames(lastrows, "rownum", "lastrow")
setkey(lastrows, id)
setkey(firstrows, id)
indivs <- firstrows[lastrows]
indivs[, followuptime := 1 + lastrow - firstrow] 
N.indivs <- nrow(indivs)
N.intervals <- nrow(X.intervals)


# X.intervals <- X.intervals[, c(1:2, 9:10)]
#X.intervals <- X.intervals[, c(1:5)]

P <- ncol(X.intervals) # number of covariates
K <- 3 # number of hidden states
V <- 3 # number of observed states

data.stan <- list(N=N.intervals, N_indivs=N.indivs, K=K, V=V, P=P,
                  p_recov=p_recov,
                  coolness=0.1,
                  params_alpha=c(10, 0.2, 2),
                  followuptime=indivs$followuptime,
                  firstrow=indivs$firstrow, lastrow=indivs$lastrow,
                  test=tsplit$y,
                  X=X.intervals)

#summary(glm(data=tsplit, formula=y>1 ~ age_years + sex + care.home + occup +
#                             listedgr3 + vax14.factor, family=poisson))
#summary(glm(data=tsplit, formula=y>2 ~ age_years + sex + care.home + occup +
#                             listedgr3 + vax14.factor, family=poisson))

##########################################################

vaxeff.stanmodel <- stan_model(file="vaxeff.stan")

####################################################################
    
init.vb <- list(theta0 = -10, theta0_test = -2, prior_alpha=c(0.76, 0.02, 0.22), gamma=3)
vaxeff.samples.vb <- rstan::vb(vaxeff.stanmodel,
              #algorithm="fullrank",
              algorithm="meanfield",
              data=data.stan, # seed=12345,
              init=init.vb,
              pars=c("theta0", "theta", "theta0_test", "theta_test", "prior_alpha", "gamma",
                     "beta", "beta_test")
              )
vaxeff.coeffs.vb <- summary(vaxeff.samples.vb,
                                probs=c(0.1, 0.9))$summary 
round(vaxeff.coeffs.vb, 3)

gc()
objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

############################################################
rnormv <- function(x)
    rnorm(length(x), mean=x)

# function to generate initial values for a single chain
initf2 <- function(chain_id = 1) {
    list(theta0 = rnormv(-8),
         theta0_test = rnormv(-2),
         gamma=rnormv(3.3),
         prior_alpha=c(0.74, 0.04, 0.22),
         chainid = chain_id)
}

# generate list of lists to specify initial values
n_chains <- 4
init_ll <- lapply(1:n_chains, function(id) initf2(chain_id = id))

vaxeff.samples.mc <- rstan::sampling(vaxeff.stanmodel,
              data=data.stan, seed=12345,
              pars=c("theta0", "theta", "theta0_test", "theta_test", "beta", "beta_test",
                     "prior_alpha", "gamma"), 
              iter=200, warmup=100,
              chains=n_chains, open_progress=FALSE,
              init=init_ll,
              control=list(adapt_delta=0.8)
              )

vaxeff.coeffs <- summary(vaxeff.samples.mc,
                                probs=c(0.1, 0.9))$summary 
round(vaxeff.coeffs, 2)

pairs(vaxeff.samples.mc, pars=c("alpha", "theta", "alpha_test", "theta_test"))

traceplot(vaxeff.samples.mc, pars=c("alpha", "beta", "alpha_test", "beta_test", "lp__"), inc_warmup=TRUE)
