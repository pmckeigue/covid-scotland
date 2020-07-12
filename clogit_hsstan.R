## train model to predict fatal disease

library(hsstan)
library(data.table)
library(parallel)

cc.severe <- as.data.table(readRDS("data/cc.severe.rds"))
                           
## restrict to fatal cases + matched controls, and complete data on covariates
select.cols <- c("CASE", "stratum", "care.home", "diabetes.any",
                 grep("^Ch\\.", colnames(cc.severe), value=TRUE))
covariates <- select.cols[-(1:2)][1:8]

cc.nonmissing <- na.omit(cc.severe[fatal.casegroup==1], cols=select.cols)
cc.nonmissing <- cc.nonmissing[, ..select.cols]
covariate.names <- select.cols[-(1:2)]
## convert all columns to numeric
cc.nonmissing <- cc.nonmissing[, lapply(.SD, as.numeric), by=row.names(cc.nonmissing)][, -1]
cc.nonmissing[, (covariate.names) := lapply(.SD, scale), .SDcols=covariate.names] 
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
cc.nonmissing <- as.data.frame(cc.nonmissing)

covs.model <- as.formula(paste("CASE ~", paste(covariates, collapse="+")))  

penalized <- covariates[-1:2]
iter <- 1000
warmup <- 500
regularized <- TRUE
nu <- 1

############################################################################

hs.clogit <- hsstan(x=cc.nonmissing, covs.model=covs.model,
                   penalized=penalized, family="cox",
                   iter=iter, warmup=warmup,
                   scale.u=2, regularized=TRUE, nu=ifelse(regularized, 1, 3),
                   par.ratio=0.05, global.df=1, slab.scale=2, slab.df=4,
                   qr=TRUE, seed=123, adapt.delta=0.9,
                   keep.hs.pars=FALSE)

save(hs.clogit, file=paste0("hs.clogit_", length(covariates), "_covariates.RData"))
print(hs.clogit)
sampler.stats(hs.clogit)

invlogit <- function(x) 1/(1+exp(-x))
hs.clogit$family <- list(family="cox", linkinv=invlogit)

sel.vars <- projsel(hs.clogit, start.from="diabetes.any")
print(sel.biom, digits=4)

                                        #  
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
