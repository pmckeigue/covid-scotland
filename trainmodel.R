rm(list=ls())
gc()

## train model to predict fatal disease

library(rstan)
library(cmdstanr)
library(data.table)
library(wevid)
library(pander)

cmdstan.sampling <- function(model=cmdstan.model,
                             data=stan.data,
                             pars=NULL,
                             open_progress=FALSE,
                             iter=500, warmup=200,
                             adapt_delta=0.95, 
                             chains=4,
                             threads_per_chain=2, 
                             seed=12345) { # returns a CmdStanMCMC object
    ## can't use the pars argument with cmdstanr
    cat("sampling via cmdstanr ... ")
    samples.cmdstan <- model$sample(data=data,
                                    iter_sampling=iter - warmup, iter_warmup=warmup,
                                    adapt_delta=adapt_delta,
                                    chains=chains,
                                    threads_per_chain=threads_per_chain,
                                    seed=seed)
    cat("done\n")
    samples.mc <- rstan::read_stan_csv(samples.cmdstan$output_files())
    return(samples.mc)
}

rstan.sampling <- function(model=NULL,
                           data=NULL,
                           pars=pars,
                           open_progress=FALSE,
                           iter=500, warmup=200,
                           adapt_delta=0.95, 
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

hsstan <- function(clogit.model,
                   X, N, P, U, S, starts, stops, ycat,
                   clogit=TRUE, 
                   penalized=TRUE, family=NULL,
                   iter=500, warmup=floor(iter / 2),
                   scale.u=2,  regularized=TRUE,
                   nu=ifelse(regularized, 1, 3),
                   global.scale=NULL, 
                   global.df=1, slab.scale=2, slab.df=4,
                   par.ratio=0.1,
                   qr=TRUE, seed=12345, adapt.delta=0.95,
                   keep.hs.pars=FALSE,
                   grainsize=1) {

    if(is.null(global.scale)) {

    }
    
    ## parameter names not to include by default in the stanfit object
    hs.pars <- c("lambda", "tau", "z", "c2")
    if (keep.hs.pars)
        hs.pars <- NA
    
    ## thin QR decomposition
    if (P > N) qr <- FALSE
    if (qr) {
        qr.dec <- qr(X)
        Q.qr <- qr.Q(qr.dec)
        R.inv <- qr.solve(qr.dec, Q.qr) * sqrt(N - 1)
        Q.qr <- Q.qr * sqrt(N - 1)
      }
    
    cat("writing data.input as list ...")
    data.input <- list(N=N, P=P, U=U, S=S, 
                       scale_u=scale.u, ycat=ycat,
                       X=X,
                       starts=starts, stops=stops,
                       global_scale=global.scale, 
                       regularized=1, nu=1, 
                       global_df=global.df,
                       slab_scale=slab.scale, slab_df=slab.df,
                       grainsize=grainsize)
    if(qr) {
        data.input$X <- Q.qr
    }
    print(summary(data.input$X))
    cat("done\n")
     
    cat("Sampling clogit model ... ")
    samples <-
       cmdstan.sampling(model=clogit.model,
       # rstan.sampling(model=clogit.model,
                       data=data.input,
                       #pars=c("tau", "c2", "beta_u", "beta_p"), 
                       adapt_delta=adapt.delta,
                       iter=iter, warmup=floor(0.4 * iter),
                       chains=4, threads_per_chain=2, seed=12345)
    cat("done\n")
    ## label beta parameters with covariate names / levels
    par.idx <- grep("^beta", names(samples))
    stopifnot(length(par.idx) == ncol(X))
    names(samples)[par.idx] <- covariatelevels

    if (qr) { # invert QR transformation
        pars <- grep("beta", samples@sim$pars_oi, value=TRUE)
        stopifnot(pars[1] == "beta")
        beta.tilde <- rstan::extract(samples, pars=pars,
                                     inc_warmup=TRUE, permuted=FALSE)
        B <- apply(beta.tilde, 1:2, FUN=function(z) R.inv %*% z)
        numchains <- ncol(beta.tilde)
        
        for (chain in 1:numchains) {
            for (p in 1:P)
                samples@sim$samples[[chain]][[par.idx[p]]] <- B[p, , chain]
        }
    }
    return(samples)
}

normalize.predictions <- function(unnorm.p, stratum, y) { ## normalize probs within each stratum
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


############################################

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Stan settings
iter=500
adapt.delta <- 0.9
seed=12345
qr <- TRUE

## shrinkage settings
scale.u <- 2
regularized <- TRUE
nu <- 1
par.ratio <- 0.1
global.df <- 1
slab.scale <- 2
slab.df <- 4

clogit.model.rstan <- stan_model(file="~/hsstan/inst/stan/hs_clogit.stan")
clogit.model.cmdstan <- cmdstanr::cmdstan_model(
                                      stan_file="~/hsstan/inst/stan/hs_clogit_parallel.stan",
                                      #force_recompile = TRUE,
                                      cpp_options = list(stan_threads = TRUE))

## load case-control dataset
## exclude care home residents and restrict to first wave
cc.severe.lt75 <- as.data.table(readRDS("data/cc.severe.lt75.rds"), key="ANON_ID")
cc.severe.lt75 <- cc.severe.lt75[care.home=="Independent" &
                                 SPECIMENDATE < as.Date("2020-06-01")]

## for an online app, we want easy to use variables
## age, sex, care home+++
## diagnoses of each listed condition++++
## up to 10 drug classes
## up to 10 extra hospital diagnoses
## household size

## for a PHS algo, we can use any variables accessible to PHS

## care home, SIMD_quintile, children in household
## all drug classes ## could maybe omit BNF chapters 11-13 (eye, ear, miscellaneous)
## all hospital diagnoses

## select covariates
select.cols <- c("CASE", "stratum", "hh.over18", "qSIMD.integer", "dm.type3", "shield.group",
                 grep("^Ch_", colnames(cc.severe.lt75), value=TRUE),
                 grep("subpara\\.", colnames(cc.severe.lt75), value=TRUE)
                 )

## FIXME: calculate number of unpenalized columns as length(levels()) - 1 for factors
U <- 10 # number of unpenalized columns

## other possible variables to add in: number of children, health-care worker, teacher

## select columns and exclude rows with missing values
cc.nonmissing <- na.omit(cc.severe.lt75[, ..select.cols])

## drop strata without at least one case and one control
cc.drop <- cc.nonmissing[, drop := length(unique(CASE)) == 1,
                         by = .(stratum)][drop==TRUE, .(CASE, stratum)]     
cc.nonmissing <- cc.nonmissing[!(stratum %in% cc.drop$stratum)]
cc.nonmissing[, drop := NULL]

## drop strata with numcases not equal to 1
## FIXME: should drop only the extra cases
cc.numcases.not1 <- cc.nonmissing[CASE==1, .(.N), by = .(stratum)][N != 1, ]
cc.nonmissing <- cc.nonmissing[!(stratum %in% cc.numcases.not1$stratum)]
## renumber levels of stratum   
cc.nonmissing[, stratum := as.integer(as.factor(as.integer(stratum)))]
setkey(cc.nonmissing, stratum)

## check we have one case per stratum
table(cc.nonmissing[CASE==1, .(.N), by = .(stratum)][, N])

## exclude columns without at least 2 levels
keep.cols <- unlist(lapply(cc.nonmissing, function(col) length(unique(col)) > 1))
cc.nonmissing <- cc.nonmissing[, ..keep.cols]

## exclude covariates with < 10 in one level
keep.cols <- unlist(lapply(cc.nonmissing, function(x) min(table(x)) >= 10))
keep.cols[1:2] <- TRUE
cc.nonmissing <- cc.nonmissing[, ..keep.cols]

covariate.names <- names(cc.nonmissing)[-(1:2)]
y <- cc.nonmissing[, CASE]
N <- length(y)
## X matrix without intercept column
X <- model.matrix(object=~ ., data=cc.nonmissing[, ..covariate.names])[, -1]
covariatelevels <- colnames(X)

P <- ncol(X)
K <- P - U
starts <- c(1, which(diff(as.integer(cc.nonmissing[, stratum])) > 0))
stops <- c(starts[-1] - 1, N)
S <- length(starts)
global.scale <- par.ratio / sqrt(S)
table(stops - starts)

ycat <- integer(S)
for(s in 1:S) { # set ycat to the index of the case within each stratum
    yind <- y[starts[s]:stops[s]]
    #print(yind)
    y.gt0 <- which(yind > 0)
    if(length(y.gt0) != 1) {
        stop("numcases in stratum != 1")
    }
    ycat[s] <- which(yind > 0)
}

############################################################################


grainsize <- 1 # ceiling(S / 16)
    
samples <- hsstan(X=X,
                  clogit.model=clogit.model.cmdstan,
                  # clogit.model=clogit.model.rstan,
                  ycat=ycat,
                  S=S,
                  N=N, P=P, U=U, starts=starts, stops=stops,
                  global.scale=global.scale, global.df=global.df,
                  iter=iter, seed=12345,
                  scale.u=scale.u, par.ratio=par.ratio,
                  grainsize=grainsize)

## store the hierarchical shrinkage settings
opts <- list(adapt.delta=adapt.delta, qr=qr, seed=seed, scale.u=scale.u)
if (K > 0)
    opts <- c(opts, regularized=regularized, nu=nu, par.ratio=par.ratio,
              global.scale=global.scale, global_df=global.df,
              slab.scale=slab.scale, slab.df=slab.df)

save(samples, opts, file="samplesopts.RData")

## compute the posterior means of the regression coefficients
betas <- colMeans(as.matrix(samples, pars="beta"))

print(summary(do.call(rbind, get_sampler_params(samples)), digits = 2))
mc.summary <- round(summary(samples,
                            probs=c(0.1, 0.9))$summary, 3)
print(mc.summary[, 1])

beta.draws <- as.matrix(samples)
beta.draws <- beta.draws[, match(covariatelevels, colnames(beta.draws))] # matrix with one row per variable
linearpredictor.draws <- X %*% t(beta.draws)
unnorm.p <- rowMeans(exp(linearpredictor.draws))

predictions <- as.data.table(normalize.predictions(unnorm.p, cc.nonmissing[, stratum], y))

## flatten extreme predictions
predictions[posterior.p < 0.005, posterior.p := 0.005]
predictions[posterior.p > 0.995, posterior.p := 0.995]

W.densities <- with(predictions,
                    Wdensities(y, posterior.p, prior.p,
                               recalibrate=FALSE))

pander(summary(W.densities), table.style="multiline",
       split.cells=c(5, 5, 5, 5, 5, 5, 5), split.table="Inf",
       caption="Prediction of severe COVID-19 from conditional logistic regression")

p.wdists <- plotWdists(densities=W.densities, distlabels=c("Crude", "Adjusted")) +
    theme_grey(base_size = 16) +
    theme(axis.title=element_text(size=14)) +
    theme(legend.title=element_blank()) +
    theme(legend.position=c(0.85, 0.75))

png("wdists.png")
p.wdists
dev.off()

## to generate risk scores in a new dataset, we need
## (1) vector covariate.names
## (2) matrix beta.draws
## then generate X matrix with no intercept
## X <- model.matrix(object=~ ., data=newdata[, ..covariate.names])[, -1]
##
## linearpredictor.draws <- X %*% t(beta.draws)
## unnorm.p <- rowMeans(exp(linearpredictor.draws))

save(covariate.names, beta.draws, file="betadraws.riskscore.RData")
