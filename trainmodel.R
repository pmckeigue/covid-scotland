
rm(list=ls())
gc()

## train model on first wave to predict severe disease, test on second wave 

library(rstan)
library(cmdstanr)
library(data.table)
library(wevid)
library(pander)

source("helperfunctions.R")

glmm.report <- function(samples.mc, report.title, exclude.pars=c("lambda", "z")) {
    pdf(paste0("diagnostics.",
               report.title, ".pdf"), title=report.title)
    etime <- round(sum(colMeans(get_elapsed_time(samples.mc))) / 3600, 2)
    cat("Sampler params for", report.title, "with adapt_delta =",
        samples.mc@stan_args[[1]]$control$adapt_delta, "at", samples.mc@date,
        "with elapsed time", etime, "hours\n")
    print(summary(do.call(rbind, get_sampler_params(samples.mc)), digits = 2))
    pairs(samples.mc, pars=c("aux1_global", "aux2_global", "caux", "lp__")) 
    stan_trace(samples.mc, pars=c("aux1_global", "aux2_global", "caux", "lp__"))
    samples.mc.densplot <- stan_dens(samples.mc,
                                     pars=c("aux1_global", "aux2_global", "caux", "lp__"))
    print(samples.mc.densplot + labs(title=""))
  
    logp <- get_logposterior(samples.mc)
    N.iterations <- length(logp[[1]])
    logp.chains <- data.frame(iteration=1:N.iterations,
                              chain1=logp[[1]], 
                              chain2=logp[[2]], 
                              chain3=logp[[3]], 
                              chain4=logp[[4]]) 
    logp.chains <- reshape2::melt(logp.chains, id="iteration")  # convert to long format
    p.conv <- ggplot(data=logp.chains,
                     aes(x=iteration, y=value, colour=variable)) + 
        geom_line() +
        theme(text = element_text(size = 20)) +
        ggtitle("Trace plot of log posterior excluding warmup")
    print(p.conv)

    mc.summary <- round(summary(samples.mc,
                                probs=c(0.1, 0.9))$summary, 3)
    keep.rows <- grep("lambda|z\\[|_local\\[", rownames(mc.summary), invert=TRUE)
    mc.summary <- mc.summary[keep.rows, ]
    print(head(mc.summary, 20))
    dev.off()
}

cmdstan.sampling <- function(model=cmdstan.model,
                             data=stan.data,
                             pars=NULL,
                             open_progress=FALSE,
                             iter=500, warmup=200,
                             adapt_delta=0.95, 
                             chains=4,
                             threads_per_chain=2, 
                             seed=12345) { # returns a stanfit object
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
                   global.scale=0.01, 
                   global.df=1, slab.scale=2, slab.df=4,
                   par.ratio=0.1,
                   qr=TRUE, seed=12345, adapt.delta=0.95,
                   keep.hs.pars=FALSE,
                   grainsize=1) {

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
    cat("done\n")
    if(qr) {
        data.input$X <- Q.qr
        cat("QR transformed variables: summary statistics for means and SDs:\n")
        print(summary(apply(Q.qr, 2, mean)))
        print(summary(apply(Q.qr, 2, sd)))
    }
    
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
    names(samples)[par.idx] <- colnames(X)
    
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

############################################
datadir <- "./data/2021-02-18/"

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Stan settings
iter=1200
warmup=floor(iter/4)
adapt.delta <- 0.8
seed=12345
qr <- TRUE
regularized <- TRUE
scale.u <- 2  # scale param for prior on unpenalized, ignored if QR transform used 
minfreq <- 1000 # minimum number of observations in each level of a binary covariate 

## shrinkage settings

# prior on beta_i is gaussian with scale param lambda.tilde_i * tau 
## lambda.tilde_i = sqrt( c2 * lambda^2 / (c2 + tau * lambda^2) )

## half-t prior on global scale parameter tau with df global.df and scale param global.scale
## Piironen and Vehtari argue that for a given number of nonzero coeffs, tau must scale inversely with the number of observations
## they recommend setting prior on tau so as to give a number of nonzero coeffs that is consistent with prior beliefs
par.ratio <- 0.1 ## can set global.scale = par.ratio / sqrt(num obs)
global.df <- 1 ## df for half-t prior on global scale param tau 
global.scale <- 0.1  ## fixed for now

## half-t prior on local shrinkage params lambda_i with df nu and scale param 1
nu <- 2 # (nu=1 gives half-Cauchy prior, but we might want less heavy tails) 

## scale param for prior on c:  c2 = slab.scale * caux
## large value of slab.scale gives prior approaching the original horseshoe (sparsity)
## small value of slab.scale gives a gaussian prior that shrinks large effects
slab.scale <- 0.1 

## controls shape of inverse gamma prior on caux
## large value will give less skewed prior on regularizing param c2
slab.df <- 4

clogit.model.rstan <- stan_model(file="~/hsstan/inst/stan/hs_clogit.stan")
clogit.model.cmdstan <- cmdstanr::cmdstan_model(
                                      stan_file="~/hsstan/inst/stan/hs_clogit_parallel.stan",
                                      #force_recompile = TRUE,
                                      cpp_options = list(stan_threads = TRUE))

## load case-control dataset
## exclude care home residents and restrict to first wave
load(paste0(datadir, "cc.severe.lt75.RData"))

## for an online app, we want easy to use variables
## age, sex
## diagnoses of each listed condition++++
## up to 10 drug classes
## up to 10 extra hospital diagnoses
## household size

## for a PHS algo, we can use any variables accessible to PHS

## care home, SIMD_quintile, children in household
## all drug classes ## could maybe omit BNF chapters 11-13 (eye, ear, miscellaneous)
## all hospital diagnoses

## select covariates
select.cols <- c("SPECIMENDATE", "CASE", "stratum",
                 "hh.upto4", "hh.5to11", "hh.12to17", "hh.over18", "qSIMD.integer",
                 "dm.type3", "shield.group",
                 grep("^Ch_", colnames(cc.lt75), value=TRUE),
                 grep("subpara\\.", colnames(cc.lt75), value=TRUE)
                 )
## other possible variables to add in: number of children, health-care worker, teacher

cc.select <- cc.lt75[, ..select.cols]
## exclude columns without at least 2 levels
keep.col <- cc.select[, sapply(cc.select, function(x) length(unique(x)) > 1)]
# keep.cols <- unlist(lapply(cc.lt75, function(col) length(unique(col)) > 1))
cc.select <- cc.select[, ..keep.col]

## exclude binary covariates from col4 onwards with < minfreq observations in one level
keep.covariate <- cc.select[, sapply(cc.select,
                                     function(x) length(unique(x)) > 2 | min(table(x)) >= minfreq)]
keep.covariate[1:3] <- TRUE

## TRIAL RUN
## keep.covariate <- keep.covariate[1:15]

cc.select <- cc.select[, ..keep.covariate]
covariate.names <- names(cc.select)

## training dataset
cc.train <- cc.select[SPECIMENDATE < as.Date("2020-09-01")]

## exclude rows with missing values
cc.train <- na.omit(cc.train)
## drop strata without at least one case and one control
cc.drop <- cc.train[, drop := length(unique(CASE)) == 1,
                         by = .(stratum)][drop==TRUE, .(CASE, stratum)]     
cc.train <- cc.train[!(stratum %in% cc.drop$stratum)]
cc.train[, drop := NULL]
## drop strata with numcases not equal to 1
## FIXME: should drop only the extra cases
cc.numcases.not1 <- cc.train[CASE==1, .(.N), by = .(stratum)][N != 1, ]
if(nrow(cc.numcases.not1) > 0) {
    cc.train <- cc.train[!(stratum %in% cc.numcases.not1$stratum)]
}
## renumber levels of stratum   
cc.train[, stratum := as.integer(as.factor(as.integer(stratum)))]
setkey(cc.train, stratum)

y <- cc.train[, CASE]
N <- length(y)
## deselect specimendate, case, stratum
covariateonly.names <- covariate.names[-(1:3)]
## X matrix without intercept column
X.train <- model.matrix(object=~ ., data=cc.train[, ..covariateonly.names])[, -1]
covariatelevels <- colnames(X.train)

P <- ncol(X.train)
## FIXME: calculate number of unpenalized columns as length(levels()) - 1 for factors
U <- 10 # number of unpenalized columns, ignored if QR transform used
K <- P - U
## elements of diff begin at first diff, so to get index of increments we must add one
starts <- c(1, 1 + which(diff(as.integer(cc.train[, stratum])) > 0))
stops <- c(starts[-1] - 1, N)
S <- length(starts)
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
samples <- hsstan(X=X.train,
                  clogit.model=clogit.model.cmdstan,
                  # clogit.model=clogit.model.rstan,
                  ycat=ycat,
                  S=S,
                  N=N, P=P, U=U, starts=starts, stops=stops,
                  global.scale=global.scale, global.df=global.df,
                  iter=iter, warmup=warmup, seed=12345,
                  scale.u=scale.u, par.ratio=par.ratio,
                  grainsize=grainsize)

## store the hierarchical shrinkage settings
opts <- list(adapt.delta=adapt.delta, qr=qr, seed=seed, scale.u=scale.u)
if (K > 0)
    opts <- c(opts, regularized=regularized, nu=nu, par.ratio=par.ratio,
              global.scale=global.scale, global_df=global.df,
              slab.scale=slab.scale, slab.df=slab.df)

save(samples, opts, file="samplesopts.RData")
rm(X.train)
glmm.report(samples, report.title="trainriskscore")
print(summary(do.call(rbind, get_sampler_params(samples)), digits = 2))

## compute the posterior means of the regression coefficients
betas <- colMeans(as.matrix(samples, pars="beta"))

mc.summary <- round(summary(samples, probs=c(0.1, 0.9))$summary, 3)
varname <- gsub('`', '', rownames(mc.summary))
mc.summary <- data.table(varname=varname, mc.summary)
mc.summary <- mc.summary[!grepl("z\\[|local|lambda|lp__", varname)]

## generate short varnames for bar chart
beta.summary <- mc.summary[varname %in% gsub('`', '', colnames(X.train))]
varname.short <- beta.summary$varname
wordsgt5 <- sapply(strsplit(varname.short, "\\s+"), length) > 5
varname.short[wordsgt5] <- stringr::word(varname.short[wordsgt5], 1, 5)
beta.summary[, varname.short := varname.short]
beta.summary[, varname.short := factor(varname.short, levels=rev(varname.short))]
p.betas <- ggplot(data=beta.summary,
                  aes(x=varname.short, y=mean)) +
    geom_col() + coord_flip() 
p.betas

beta.draws <- as.matrix(samples)
beta.draws <- beta.draws[, match(covariatelevels, colnames(beta.draws))] # matrix with one row per variable


## generate test dataset
cc.test <- cc.lt75[SPECIMENDATE >= as.Date("2020-09-01"), ..covariate.names]
## exclude rows with missing values
cc.test <- na.omit(cc.test)
## drop strata without at least one case and one control
cc.drop <- cc.test[, drop := length(unique(CASE)) == 1,
                         by = .(stratum)][drop==TRUE, .(CASE, stratum)]     
cc.test <- cc.test[!(stratum %in% cc.drop$stratum)]
cc.test[, drop := NULL]
## drop strata with numcases not equal to 1
## FIXME: should drop only the extra cases
cc.numcases.not1 <- cc.test[CASE==1, .(.N), by = .(stratum)][N != 1, ]
if(nrow(cc.numcases.not1) > 0) {
    cc.test <- cc.test[!(stratum %in% cc.numcases.not1$stratum)]
}
## renumber levels of stratum   
cc.test[, stratum := as.integer(as.factor(as.integer(stratum)))]
setkey(cc.test, stratum)

## calculate normalized predictions in test dataset
X.test <- model.matrix(object=~ ., data=cc.test[, ..covariateonly.names])[, -1]
unnorm.lograte.draws <- X.test %*% t(beta.draws)
rm(X.test)
unnorm.lograte <- matrixStats::rowLogSumExps(unnorm.lograte.draws) -
    log(ncol(unnorm.lograte.draws))
cc.predictions <- cc.test[, .(CASE, stratum, unnorm.lograte)]        
predictions <- with(cc.predictions, normalize.logrates(unnorm.lograte, stratum, CASE))

Cstat <- mean(predictions[, pROC::auc(response=y,
                                      predictor=posterior.p),
                          by=stratum][["V1"]])
cat("C-statistic computed as mean AUC over strata:", Cstat, "\n")

## flatten extreme predictions
predictions[posterior.p < 0.002, posterior.p := 0.002]
predictions[posterior.p > 0.998, posterior.p := 0.998]

W.densities <- with(predictions,
                    Wdensities(y, posterior.p, prior.p, adjust.bw=.1, recalibrate=FALSE))

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

save(covariate.names, covariateonly.names, beta.draws,
     file=paste0(datadir, "betadraws.riskscore.RData"))
