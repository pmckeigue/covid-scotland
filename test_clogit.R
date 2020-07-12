library("extraDistr")

data(diabetes, package="hsstan")

map.varnames <- grep("map", names(diabetes), value=TRUE)
P <- length(map.varnames)
N <- nrow(diabetes)
beta <- matrix(rnorm(P, sd=0.5), ncol=1)
Xbeta <- as.matrix(diabetes[, map.varnames]) %*% beta
S <- ceiling(N / 5)
diabetes$stratum <- rep(1:S, each=5)[1:N]

seed <- 12345

diabetes$y.clogit <- numeric(N)
for(s in 1:S) {
    xbeta <- Xbeta[diabetes$stratum==s, ]
    ycat <- extraDistr::rcat(n=1, prob=exp(xbeta) / sum(exp(xbeta)) )

    ## for s = 1 this will assign y=1 to element ycat
    diabetes$y.clogit[match(s, diabetes$stratum) - 1 + ycat] <- 1
}

penalized <- map.varnames[-(1:2)]
regularized <- TRUE

covs.model <- as.formula(paste("y.clogit ~", paste(map.varnames, collapse="+")))  

hs.clogit <- hsstan(x=diabetes, covs.model=covs.model,
                   penalized=penalized, family="cox",
                   iter=1000, warmup=500,
                   scale.u=2, regularized=TRUE, nu=ifelse(regularized, 1, 3),
                   par.ratio=0.05, global.df=1, slab.scale=2, slab.df=4,
                   qr=TRUE, seed=123, adapt.delta=0.9,
                   keep.hs.pars=FALSE)

print(hs.clogit)

sampler.stats(hs.clogit)

invlogit <- function(x) 1/(1+exp(-x))

hs.clogit$family <- list(family="cox", linkinv=invlogit)

sel.vars <- projsel(hs.clogit, start.from="map.2")

print(sel.vars, digits=4)

