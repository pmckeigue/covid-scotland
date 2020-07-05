## train model to predict fatal disease

library(rstan)

rstan.sampling <- function(model=rstan.model, data=stan.data, pars=pars,
                           iter=1000, warmup=500,
                           open_progress=FALSE,
                           adapt_delta=0.8, 
                           chains=4,
                           seed=12345) {
    
    cat("sampling via rstan ... ")
    samples.mc <- rstan::sampling(object=model, data=data,
                                  pars=pars, 
                           iter=iter, warmup=warmup,
                           control=list(adapt_delta=adapt_delta),
                           chains=chains, seed=seed,
                           open_progress=open_progress)
    cat("done\n")
    return(samples.mc)
}

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
qr <- TRUE

clogit.model <- stan_model(file="~/hsstan/inst/stan/hs_clogit.stan")

## restrict to fatal cases + matched controls, and complete data on covariates
select.cols <- c("CASE", "stratum", "care.home", "diabetes.any",
                 grep("^Ch\\.", colnames(cc.severe), value=TRUE))


cc.nonmissing <- na.omit(cc.severe[fatal.casegroup==1], cols=select.cols)
cc.nonmissing <- cc.nonmissing[, ..select.cols]
covariate.names <- select.cols[-(1:2)]
str(cc.nonmissing)
cc.nonmissing <- cc.nonmissing[, lapply(.SD, as.numeric), by=row.names(cc.nonmissing)][, -1]
str(cc.nonmissing)
              
## drop strata with num cases not equal to 1, or < 2 observations 
cc.drop <- cc.nonmissing[CASE==1, .(.N), by = .(stratum)][N != 1, ]
cc.nonmissing <- cc.nonmissing[!(stratum %in% cc.drop$stratum)]
cc.drop <- cc.nonmissing[, .(.N), by = .(stratum)][N < 2, ]
cc.nonmissing <- cc.nonmissing[!(stratum %in% cc.drop$stratum)]
              
table(cc.nonmissing[, .(.N), by = .(stratum)][, N])

## renumber levels of stratum   
cc.nonmissing[, stratum := as.integer(as.factor(as.integer(stratum)))]
setkey(cc.nonmissing, stratum)

y <- cc.nonmissing[, CASE]
N <- length(y)
X <- as.matrix(cc.nonmissing[, ..covariate.names])
P <- ncol(X)
U <- 1
starts <- c(1, which(diff(as.integer(cc.nonmissing[, stratum])) > 0))
stops <- c(starts[-1] - 1, N)
S <- length(starts)
table(stops - starts)

ycat <- integer(S)
for(s in 1:S) {
    yind <- y[starts[s]:stops[s]]
    y.gt0 <- which(yind > 0)
    if(length(y.gt0) != 1)
        stop("numcases in stratum != 1")
    ycat[s] <- which(yind > 0)
}

   ## thin QR decomposition
    if (P > N) qr <- FALSE
    if (qr) {
        qr.dec <- qr(X)
        Q.qr <- qr.Q(qr.dec)
        R.inv <- qr.solve(qr.dec, Q.qr) * sqrt(N - 1)
        Q.qr <- Q.qr * sqrt(N - 1)
    }

train.data <- list(N=N, P=P, U=U, S=S, 
                   scale_u=10, ycat=ycat,
                   X=if (qr) Q.qr else X,
                   starts=starts, stops=stops,
                   global_scale=0.5 / sqrt(S), 
                   regularized=1, nu=1, 
                   par.ratio=0.05, global_df=1, slab_scale=2, slab_df=4) 

cat("Sampling clogit model ... ")
clogit.samples.mc <-
    rstan.sampling(model=clogit.model,
                   data=train.data,
                   pars=c("tau", "c2", "beta_u", "beta_p"), 
                   adapt_delta=0.95,
                   iter=1000, warmup=500,
                   chains=4, seed=12345)
cat("done\n")

## assign covariate names
par.idx <- grep("^beta_[up]", names(clogit.samples))
stopifnot(length(par.idx) == ncol(X))
names(samples)[par.idx] <- covariate.names

if (qr) { # invert QR transformation
    pars <- grep("beta_", clogit.samples.mc@sim$pars_oi, value=TRUE)
    stopifnot(pars[1] == "beta_u")
    beta.tilde <- rstan::extract(clogit.samples.mc, pars=pars,
                                 inc_warmup=TRUE, permuted=FALSE)
    B <- apply(beta.tilde, 1:2, FUN=function(z) R.inv %*% z)
    chains <- ncol(beta.tilde)
    for (chain in 1:chains) {
        for (p in 1:P)
            clogit.samples.mc@sim$samples[[chain]][[par.idx[p]]] <- B[p, , chain]
    }
}

print(summary(do.call(rbind, get_sampler_params(clogit.samples.mc)), digits = 2))
mc.summary <- round(summary(clogit.samples.mc,
                            probs=c(0.1, 0.9))$summary, 3)
print(mc.summary)


## within each stratum, normalize predicted values are a simplex that sums to 1.

## we want to fit a model to these predicted values

## in other words to fit a multinomial logistic model.

## 

## for an online app, we want easy to use variables
## age, sex, care home+++
## diagnoses of each listed condition++++
## up to 10 drug classes
## up to 10 extra hospital diagnoses
## household size

## for a PHS algo, we can use any variables accessible to PHS

## care home, SIMD_quintile
## all drug classes
## all hospital diagnoses
## maybe don't include listed conditions except from registries
