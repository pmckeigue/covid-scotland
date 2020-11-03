## train model to predict fatal disease

library(hsstan)
library(data.table)
library(parallel)

devtools::load_all("../hsstan")

newrun <- TRUE

cc.severe <- as.data.table(readRDS("data/cc.severe.rds"))

drugclasses <-  grep("^subpara\\.", colnames(cc.severe), value=TRUE)
## restrict to BNF subparas with at least 50 exposed
keep.drugclasses <- logical(length(drugclasses))
for(j in 1:length(keep.drugclasses)) {
    x <- length(which(cc.severe[[drugclasses[j]]]=="1"))
    if(x >= 50) {
        keep.drugclasses[j] <- TRUE
    }
}
drugclasses <- drugclasses[keep.drugclasses]

subchapters <- grep("^Ch_", colnames(cc.severe), value=TRUE)
## restrict to ICD subchapters with at least 50 exposed
keep.subchapters <- logical(length(subchapters))
for(j in 1:length(keep.subchapters)) {
    x <- length(which(cc.severe[[subchapters[j]]]=="1"))
    if(x >= 50) {
        keep.subchapters[j] <- TRUE
    }
}
subchapters <- subchapters[keep.subchapters]

## restrict to fatal cases + matched controls, and complete data on covariates
select.cols <- c("CASE", "stratum", "care.home", subchapters, drugclasses)

cc.nonmissing <- na.omit(cc.severe[fatal.casegroup==1, ..select.cols], cols=select.cols)
rm(cc.severe)

## drop strata with num cases not equal to 1, or < 2 observations 
cc.drop <- cc.nonmissing[CASE==1, .(.N), by = .(stratum)][N != 1, ]
cc.nonmissing <- cc.nonmissing[!(stratum %in% cc.drop$stratum)]
cc.drop <- cc.nonmissing[, .(.N), by = .(stratum)][N < 2, ]
cc.nonmissing <- cc.nonmissing[!(stratum %in% cc.drop$stratum)]
table(cc.nonmissing[, .(.N), by = .(stratum)][, N])
## renumber levels of stratum   
cc.nonmissing[, stratum := as.integer(as.factor(as.integer(stratum)))]
setkey(cc.nonmissing, stratum)

## clean up column names
## surrounding names with backticks returns error later when backticks are stripped out
# colnames(cc.nonmissing) <- paste0("`", colnames(cc.nonmissing), "`")

colnames(cc.nonmissing) <- gsub("[ ,':()/]|\\-|\\[|\\]", "_", colnames(cc.nonmissing))
#colnames(cc.nonmissing) <- gsub("-", "_", colnames(cc.nonmissing))

print(colnames(cc.nonmissing))

## rename covariates
colnames(cc.nonmissing)[1:3] <- c("y", "stratum", "care.home")
covariate.names <- colnames(cc.nonmissing)[-(1:2)]

## convert all columns to numeric
cc.nonmissing <- cc.nonmissing[, lapply(.SD, as.numeric), by=row.names(cc.nonmissing)][, -1]

## scale and centre all covariates
cc.nonmissing[, (covariate.names) := lapply(.SD, scale), .SDcols=covariate.names] 

covs.model <- as.formula(paste("y ~",
                               paste(covariate.names[1], collapse=" + ")))
## removes backticks

penalized <- covariate.names[-1]
iter <- 1200
warmup <- 400
regularized <- TRUE
nu <- 1

############################################################################
if(newrun) {
    hs.clogit <- hsstan(x=cc.nonmissing, covs.model=covs.model,
                        penalized=penalized, 
                        family="clogit",
                        iter=iter, warmup=warmup,
                        scale.u=2, regularized=TRUE, nu=ifelse(regularized, 1, 3),
                        par.ratio=0.05, global.df=1, slab.scale=2, slab.df=4,
                        qr=TRUE, seed=123, adapt.delta=0.9,
                        keep.hs.pars=FALSE)
    
    save(hs.clogit, file=paste0("hs.clogit_", length(covariate.names), "_covariates.RData"))
} else {
    load(file=paste0("hs.clogit_", length(covariate.names), "_covariates.RData"))
} 
rm(cc.nonmissing)

print(hs.clogit)
sampler.stats(hs.clogit)

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))
gc()

sel.vars <- projsel(hs.clogit, start.from=c("care.home"))
gc()

print(sel.vars, digits=4)
save(sel.vars, file=paste0("projsel_", length(covariate.names), "_covariates.RData"))

pdf("plot_projsel.pdf")
 plot.projsel(sel.vars, covariates.from=FALSE)
dev.off()

                                     #  
## for an online app, we want easy to use variables
## age, sex, care home+++
## diagnoses of each listed condition++++
## up to 10 drug classes
## up to 10 extra hospital diagnoses

## for a PHS algo, we can use any variables accessible to PHS

## care home, SIMD_quintile
## all drug classes
## all hospital diagnoses
## maybe don't include listed conditions except from registries
