## model for slope of rate ratio to time since 2nd dose  ########

profile.ci <- function(loglik, x) {
    N <- length(x)
    row.max <- which.max(loglik)
    if(row.max==1 | row.max==N)
       warning("loglik maximal at first or last value of param: use wider range\n")
    loglik.diff <- loglik - max(loglik)
    rows.inner <- which(loglik.diff >= -1.92)
    N.inner <- length(rows.inner)
    ## corner cases
    if(row.max==1) {
        warning("log lik maximal at lowest param value\n")
        rows.inner <- c(1, rows.inner[N.inner])
        rows.outer <- rows.inner + 0:1
    } else if(row.max==N) {
        warning("log lik maximal at highest param value\n")
        rows.inner <- c(rows.inner[1], N)
        rows.outer <- rows.inner + 1:0
    }  else {
        rows.inner <- c(rows.inner[1], rows.inner[N.inner])
        rows.outer <- rows.inner + c(-1, 1) 
        if(rows.inner[1]==1 ) {
            warning("lowest param value is above profile lik threshold\n")
            rows.outer[1] <- 1
        }
        if(rows.inner[2]==N) {
            warning("highest param value is above profile lik threshold\n")
            rows.outer[2] <- N.inner
        }
    }

    x.inner <- x[rows.inner] 
    x.outer <- x[rows.outer]
    ll.inner <- loglik.diff[rows.inner] 
    ll.outer <- loglik.diff[rows.outer]
    gradient <- (x.inner - x.outer) / ifelse(abs(ll.inner - ll.outer) > 0, ll.inner - ll.outer, 1)
    delta.ll <- -1.92 - ll.inner
    delta.x <- gradient * delta.ll
    x.interpolated <- round(x.inner + delta.x)
    x.interpolated[x.interpolated > exp(100)] <- Inf
    return(x.interpolated)
}
#loglik <- decay.nonzero.severe$loglik.decay.mRNA
#x <- decay.nonzero.severe$halflife
#x[which.max(loglik)]
#profile.ci(loglik, x)

## unvaccinated will be coded as vaxclass AZ
cc.sincedec2020[, vaxgr3 := car::recode(vaxgr,
                                        "c('1 dose mRNA vaccine', '1 dose AZ vaccine')=NA",
                                        as.factor=TRUE,
                                        levels=c("Not vaccinated",
                                                 "2 doses mRNA vaccine",
                                                 "2 doses AZ vaccine"))]

halflives <-  round(exp(c(seq(log(3), log(500), length=50), Inf))) # lambda = log 2 / halflife
halflives <- unique(halflives)

########### severe disease, waning to zero ##############

loglik.decay <- numeric(length(halflives))
coeff.decay.AZ <- numeric(length(halflives))
coeff.decay.mRNA <- numeric(length(halflives))
loglik.decay.AZ <- numeric(length(halflives))
loglik.decay.mRNA <- numeric(length(halflives))
for(i in 1:length(halflives)) {
    cc.sincedec2020[, decay.sincelastdose := exp(-(log(2)/ halflives[i]) *
                                                 dayssincelastdose)]
    cc.sincedec2020[, mRNA.decay := as.integer(vaxgr3=="2 doses mRNA vaccine") * decay.sincelastdose]
    cc.sincedec2020[, AZ.decay := as.integer(vaxgr3=="2 doses AZ vaccine") * decay.sincelastdose]
    cc.sincedec2020[vaxgr3=="Not vaccinated", AZ.decay := 0]
    cc.sincedec2020[vaxgr3=="Not vaccinated", mRNA.decay := 0]
    decay.model.AZ <-  clogit(data=cc.sincedec2020[casegroup=="A" &
                                                vaxgr3 != "2 doses mRNA vaccine"],
                           formula=CASE ~ care.home + listedgr3 + numdrugs.notcv +
                               inpat.recent +  
                              AZ.decay + # decay term
                               strata(stratum))
    decay.model.mRNA <-  clogit(data=cc.sincedec2020[casegroup=="A" &
                                              vaxgr3 != "2 doses AZ vaccine"],
                             formula=CASE ~ care.home + listedgr3 + numdrugs.notcv +
                               inpat.recent +  
                               mRNA.decay + # decay term
                               strata(stratum))
    coeff.decay.AZ[i] <- summary(decay.model.AZ)$coefficients[6, 1]
    coeff.decay.mRNA[i] <- summary(decay.model.mRNA)$coefficients[6, 1]
    loglik.decay.AZ[i] <- decay.model.AZ$loglik[2]
    loglik.decay.mRNA[i] <- decay.model.mRNA$loglik[2]
}
decaytozero.severe <- data.table(halflife=halflives, loglik.decay.AZ, loglik.decay.mRNA, 
                                 coeff.decay.AZ, coeff.decay.mRNA)

halflife.ci.AZ <- with(decaytozero.severe, profile.ci(loglik=loglik.decay.AZ, x=halflife))
halflife.mle.AZ <- decaytozero.severe[loglik.decay.AZ==max(loglik.decay.AZ), halflife]
halflife.severe.waningtozero.estci.AZ <- paste0(halflife.mle.AZ, " (95% CI ", halflife.ci.AZ[1], " to ", halflife.ci.AZ[2], ")")
loglik.max.severe.waningtozero.AZ <- decaytozero.severe[, max(loglik.decay.AZ)]

halflife.ci.mRNA <- with(decaytozero.severe, profile.ci(loglik=loglik.decay.mRNA, x=halflife))
halflife.mle.mRNA <- decaytozero.severe[loglik.decay.mRNA==max(loglik.decay.mRNA), halflife]
halflife.severe.waningtozero.estci.mRNA <- paste0(halflife.mle.mRNA, " (95% CI ", halflife.ci.mRNA[1], " to ", halflife.ci.mRNA[2], ")")
loglik.max.severe.waningtozero.mRNA <- decaytozero.severe[, max(loglik.decay.mRNA)]

halflife.severe.waningtozero.estci.AZ
loglik.max.severe.waningtozero.AZ
halflife.severe.waningtozero.estci.mRNA
loglik.max.severe.waningtozero.mRNA

########## severe disease, waning to a nonzero value ###############

loglik.decay <- numeric(length(halflives))
coeff.const.AZ <- numeric(length(halflives))
coeff.const.mRNA <- numeric(length(halflives))
coeff.decay.AZ <- numeric(length(halflives))
coeff.decay.mRNA <- numeric(length(halflives))
for(i in 1:length(halflives)) {
    cc.sincedec2020[, decay.sincelastdose := exp(-(log(2)/ halflives[i]) *
                                                 dayssincelastdose)]
    cc.sincedec2020[, mRNA.decay := as.integer(vaxgr3=="2 doses mRNA vaccine") * decay.sincelastdose]
    cc.sincedec2020[, AZ.decay := as.integer(vaxgr3=="2 doses AZ vaccine") * decay.sincelastdose]
    cc.sincedec2020[vaxgr3=="Not vaccinated", AZ.decay := 0]
    cc.sincedec2020[vaxgr3=="Not vaccinated", mRNA.decay := 0]
    decay.model.AZ <-  clogit(data=cc.sincedec2020[casegroup=="A" &
                                               vaxgr3 != "2 doses mRNA vaccine"],
                           formula=CASE ~ care.home + listedgr3 + numdrugs.notcv +
                               inpat.recent + vaxgr3 + 
                              AZ.decay + # decay term
                               strata(stratum))
    decay.model.mRNA <-  clogit(data=cc.sincedec2020[casegroup=="A" &
                                              vaxgr3 != "2 doses AZ vaccine"],
                           formula=CASE ~ care.home + listedgr3 + numdrugs.notcv +
                               inpat.recent + vaxgr3 + 
                               mRNA.decay + # decay term
                               strata(stratum))
    coeff.const.AZ[i] <- summary(decay.model.AZ)$coefficients[7, 1]
    coeff.const.mRNA[i] <- summary(decay.model.mRNA)$coefficients[6, 1]
    coeff.decay.AZ[i] <- summary(decay.model.AZ)$coefficients[8, 1]
    coeff.decay.mRNA[i] <- summary(decay.model.mRNA)$coefficients[8, 1]
    loglik.decay.AZ[i] <- decay.model.AZ$loglik[2]
    loglik.decay.mRNA[i] <- decay.model.mRNA$loglik[2]
}
decay.nonzero.severe <- data.table(halflife=halflives, loglik.decay.AZ, loglik.decay.mRNA,
                                   coeff.const.mRNA, coeff.const.AZ,
                                   coeff.decay.mRNA, coeff.decay.AZ)

halflife.ci.AZ <- as.integer(with(decay.nonzero.severe, profile.ci(loglik=loglik.decay.AZ, x=halflife)))
halflife.mle.AZ <- as.integer(decay.nonzero.severe[loglik.decay.AZ==max(loglik.decay.AZ), halflife])
halflife.severe.waningnonzero.estci.AZ <- paste0(halflife.mle.AZ, " (95% CI ", halflife.ci.AZ[1], " to ", halflife.ci.AZ[2], ")")
loglik.max.severe.waningnonzero.AZ <- decay.nonzero.severe[, max(loglik.decay.AZ)]

halflife.ci.mRNA <- as.integer(with(decay.nonzero.severe, profile.ci(loglik=loglik.decay.mRNA, x=halflife)))
halflife.mle.mRNA <- as.integer(decay.nonzero.severe[loglik.decay.mRNA==max(loglik.decay.mRNA), halflife])
halflife.severe.waningnonzero.estci.mRNA <- paste0(halflife.mle.mRNA, " (95% CI ", halflife.ci.mRNA[1], " to ", halflife.ci.mRNA[2], ")")
loglik.max.severe.waningnonzero.mRNA <- decay.nonzero.severe[, max(loglik.decay.mRNA)]

#################  repeat for hosp  #######################################

loglik.decay <- numeric(length(halflives))
coeff.decay.AZ <- numeric(length(halflives))
coeff.decay.mRNA <- numeric(length(halflives))
loglik.decay.AZ <- numeric(length(halflives))
loglik.decay.mRNA <- numeric(length(halflives))
for(i in 1:length(halflives)) {
    cc.sincedec2020[, decay.sincelastdose := exp(-(log(2)/ halflives[i]) *
                                                 dayssincelastdose)]
    cc.sincedec2020[, mRNA.decay := as.integer(vaxgr=="2 doses mRNA vaccine") * decay.sincelastdose]
    cc.sincedec2020[, AZ.decay := as.integer(vaxgr=="2 doses AZ vaccine") * decay.sincelastdose]
    cc.sincedec2020[vaxgr=="Not vaccinated", AZ.decay := 0]
    cc.sincedec2020[vaxgr=="Not vaccinated", mRNA.decay := 0]
    decay.model.AZ <-  clogit(data=cc.sincedec2020[(casegroup=="A" | casegroup=="B") &
                                                vax14.factor != "1"],
                           formula=CASE ~ care.home + listedgr3 + numdrugs.notcv +
                               inpat.recent +  
                              AZ.decay + # decay term
                               strata(stratum))
    decay.model.mRNA <-  clogit(data=cc.sincedec2020[casegroup=="A" &
                                                vax14.factor != "1"],
                           formula=CASE ~ care.home + listedgr3 + numdrugs.notcv +
                               inpat.recent +  
                               mRNA.decay + # decay term
                               strata(stratum))
    coeff.decay.AZ[i] <- summary(decay.model.AZ)$coefficients[6, 1]
    coeff.decay.mRNA[i] <- summary(decay.model.mRNA)$coefficients[6, 1]
    loglik.decay.AZ[i] <- decay.model.AZ$loglik[2]
    loglik.decay.mRNA[i] <- decay.model.mRNA$loglik[2]
}
decaytozero.hosp <- data.table(halflife=halflives, loglik.decay.AZ, loglik.decay.mRNA, 
                                 coeff.decay.AZ, coeff.decay.mRNA)
decaytozero.hosp[, loglik.decay.AZ.diff := loglik.decay.AZ - max(loglik.decay.AZ)]
decaytozero.hosp[, loglik.decay.mRNA.diff := loglik.decay.mRNA - max(loglik.decay.mRNA)]

halflife.ci.AZ <- with(decaytozero.hosp, profile.ci(loglik=loglik.decay.AZ, x=halflife))
halflife.mle.AZ <- decaytozero.hosp[loglik.decay.AZ==max(loglik.decay.AZ), halflife]
halflife.hosp.waningtozero.estci.AZ <- paste0(halflife.mle.AZ, " (95% CI ", halflife.ci.AZ[1], " to ", halflife.ci.AZ[2], ")")
loglik.max.hosp.waningtozero.AZ <- decaytozero.hosp[, max(loglik.decay.AZ)]

halflife.ci.mRNA <- with(decaytozero.hosp, profile.ci(loglik=loglik.decay.mRNA.diff, x=halflife))
halflife.mle.mRNA <- decaytozero.hosp[loglik.decay.mRNA==max(loglik.decay.mRNA), halflife]
halflife.hosp.waningtozero.estci.mRNA <- paste0(halflife.mle.mRNA, " (95% CI ", halflife.ci.mRNA[1], " to ", halflife.ci.mRNA[2], ")")
loglik.max.hosp.waningtozero.mRNA <- decaytozero.hosp[, max(loglik.decay.mRNA)]

########## hosp disease, waning to a nonzero value ###############

loglik.decay <- numeric(length(halflives))
coeff.const.AZ <- numeric(length(halflives))
coeff.const.mRNA <- numeric(length(halflives))
coeff.decay.AZ <- numeric(length(halflives))
coeff.decay.mRNA <- numeric(length(halflives))
for(i in 1:length(halflives)) {
    cc.sincedec2020[, decay.sincelastdose := exp(-(log(2)/ halflives[i]) *
                                                 dayssincelastdose)]
    cc.sincedec2020[, mRNA.decay := as.integer(vaxgr=="2 doses mRNA vaccine") * decay.sincelastdose]
    cc.sincedec2020[, AZ.decay := as.integer(vaxgr=="2 doses AZ vaccine") * decay.sincelastdose]
    cc.sincedec2020[vaxgr=="Not vaccinated", AZ.decay := 0]
    cc.sincedec2020[vaxgr=="Not vaccinated", mRNA.decay := 0]
    decay.model.AZ <-  clogit(data=cc.sincedec2020[(casegroup=="A" | casegroup=="B") &
                                                vax14.factor != "1"],
                           formula=CASE ~ care.home + listedgr3 + numdrugs.notcv +
                               inpat.recent + vaxgr3 + 
                              AZ.decay + # decay term
                               strata(stratum))
    decay.model.mRNA <-  clogit(data=cc.sincedec2020[casegroup=="A" &
                                                vax14.factor != "1"],
                           formula=CASE ~ care.home + listedgr3 + numdrugs.notcv +
                               inpat.recent + vaxgr3 + 
                               mRNA.decay + # decay term
                               strata(stratum))
    coeff.const.AZ[i] <- summary(decay.model.AZ)$coefficients[7, 1]
    coeff.const.mRNA[i] <- summary(decay.model.mRNA)$coefficients[6, 1]
    coeff.decay.AZ[i] <- summary(decay.model.AZ)$coefficients[8, 1]
    coeff.decay.mRNA[i] <- summary(decay.model.mRNA)$coefficients[8, 1]
    loglik.decay.AZ[i] <- decay.model.AZ$loglik[2]
    loglik.decay.mRNA[i] <- decay.model.mRNA$loglik[2]
}
decay.nonzero.hosp <- data.table(halflife=halflives, loglik.decay.AZ, loglik.decay.mRNA,
                                   coeff.const.mRNA, coeff.const.AZ,
                                   coeff.decay.mRNA, coeff.decay.AZ)
decay.nonzero.hosp[, loglik.decay.AZ.diff := loglik.decay.AZ - max(loglik.decay.AZ)]
decay.nonzero.hosp[, loglik.decay.mRNA.diff := loglik.decay.mRNA - max(loglik.decay.mRNA)]

halflife.ci.AZ <- as.integer(with(decay.nonzero.hosp, profile.ci(loglik=loglik.decay.AZ, x=halflife)))
halflife.mle.AZ <- as.integer(decay.nonzero.hosp[loglik.decay.AZ==max(loglik.decay.AZ), halflife])
halflife.hosp.waningnonzero.estci.AZ <- paste0(halflife.mle.AZ, " (95% CI ", halflife.ci.AZ[1], " to ", halflife.ci.AZ[2], ")")
loglik.max.hosp.waningnonzero.AZ <- decay.nonzero.hosp[, max(loglik.decay.AZ)]

halflife.ci.mRNA <- as.integer(with(decay.nonzero.hosp, profile.ci(loglik=loglik.decay.mRNA, x=halflife)))
halflife.mle.mRNA <- as.integer(decay.nonzero.hosp[loglik.decay.mRNA==max(loglik.decay.mRNA), halflife])
halflife.hosp.waningnonzero.estci.mRNA <- paste0(halflife.mle.mRNA, " (95% CI ", halflife.ci.mRNA[1], " to ", halflife.ci.mRNA[2], ")")
loglik.max.hosp.waningnonzero.mRNA <- decay.nonzero.hosp[, max(loglik.decay.mRNA)]

######################################################

halflife.severe.waningtozero.estci.AZ
halflife.severe.waningnonzero.estci.AZ
halflife.severe.waningtozero.estci.mRNA
halflife.severe.waningnonzero.estci.mRNA
halflife.hosp.waningtozero.estci.AZ
halflife.hosp.waningtozero.estci.mRNA
halflife.hosp.waningnonzero.estci.AZ
halflife.hosp.waningnonzero.estci.mRNA

loglik.max.severe.waningtozero.AZ
loglik.max.severe.waningnonzero.AZ
loglik.max.severe.waningtozero.mRNA
loglik.max.severe.waningnonzero.mRNA
loglik.max.hosp.waningtozero.AZ
loglik.max.hosp.waningnonzero.AZ
loglik.max.hosp.waningtozero.mRNA
loglik.max.hosp.waningnonzero.mRNA

table.models <- data.table(model=rep(c("To zero", "To constant efficacy"), 4),
                           product=rep(c(rep("AZ", 2), rep("mRNA", 2)), 2),
                           outcome=c(rep("Severe", 4), rep("Hospitalised", 4)),
                           loglik=c(loglik.max.severe.waningtozero.AZ, 
                                    loglik.max.severe.waningnonzero.AZ, 
                                    loglik.max.severe.waningtozero.mRNA, 
                                    loglik.max.severe.waningnonzero.mRNA, 
                                    loglik.max.hosp.waningtozero.AZ, 
                                    loglik.max.hosp.waningnonzero.AZ, 
                                    loglik.max.hosp.waningtozero.mRNA, 
                                    loglik.max.hosp.waningnonzero.mRNA),
                           halflife=c(halflife.severe.waningtozero.estci.AZ, 
                                      halflife.severe.waningnonzero.estci.AZ, 
                                      halflife.severe.waningtozero.estci.mRNA, 
                                      halflife.severe.waningnonzero.estci.mRNA, 
                                      halflife.hosp.waningtozero.estci.AZ, 
                                      halflife.hosp.waningnonzero.estci.AZ, 
                                      halflife.hosp.waningtozero.estci.mRNA, 
                                      halflife.hosp.waningnonzero.estci.mRNA)
                           )

subtract <- table.models[model=="To zero", loglik]
table.models[model=="To zero", loglik := loglik - subtract]
table.models[model=="To constant efficacy", loglik := loglik - subtract]
table.models[, loglik := round(loglik, 1)]
table.models[, halflife := gsub("to N.*$", "to Inf)", halflife)]

########################################################################
days <- 1:240

coeffs.max.severe.waning.AZ <- decay.nonzero.severe[which.max(loglik.decay.AZ),
                                                        .(halflife, coeff.const.AZ, coeff.decay.AZ)]
fitted.severe.waning.AZ <- with(coeffs.max.severe.waning.AZ,
                                   data.table(days, product="AZ",
                                              coeff=coeff.const.AZ +
                                                  coeff.decay.AZ * exp( -(log(2) / halflife) * days)))
coeffs.max.severe.waning.mRNA <- decay.nonzero.severe[which.max(loglik.decay.mRNA),
                                                        .(halflife,
                           coeff.const.mRNA, coeff.decay.mRNA)]
fitted.severe.waning.mRNA <- with(coeffs.max.severe.waning.mRNA,
                                   data.table(days, product="mRNA",
                                              coeff=coeff.const.mRNA + coeff.decay.mRNA *
                                                  exp( -(log(2) / halflife) * days)))
fitted.severe.waning <- rbind(fitted.severe.waning.AZ,
                              fitted.severe.waning.mRNA)
fitted.severe.waning[, outcome := "Severe COVID-19"]
fitted.severe.waning[, weeks := days / 7]

coeffs.max.hosp.waning.AZ <- decay.nonzero.hosp[which.max(loglik.decay.AZ),
                                                        .(halflife,
                           coeff.const.AZ, coeff.decay.AZ)]
fitted.hosp.waning.AZ <- with(coeffs.max.hosp.waning.AZ,
                                   data.table(days, product="AZ",
                                              coeff=coeff.const.AZ + coeff.decay.AZ *
                                                  exp( -(log(2) / halflife) * days)))
coeffs.max.hosp.waning.mRNA <- decay.nonzero.hosp[which.max(loglik.decay.mRNA),
                                                        .(halflife,
                           coeff.const.mRNA, coeff.decay.mRNA)]
fitted.hosp.waning.mRNA <- with(coeffs.max.hosp.waning.mRNA,
                                   data.table(days, product="mRNA",
                                              coeff=coeff.const.mRNA + coeff.decay.mRNA *
                                                  exp( -(log(2) / halflife) * days)))
fitted.hosp.waning <- rbind(fitted.hosp.waning.AZ,
                                    fitted.hosp.waning.mRNA)
fitted.hosp.waning[, outcome := "Hosp COVID-19"]
fitted.hosp.waning[, weeks := days / 7]

####################################################################
fitted.all <- rbind(fitted.severe.waning,
                    fitted.hosp.waning)

breakpoints <- c(0.02, 0.05, 0.1, 0.2, 0.5, 1)
p.decay <- ggplot(data=fitted.all, aes(x=weeks, y=coeff, color=product, linetype=outcome)) + geom_line() +
    labs(x=paste0("Weeks since second dose"),
         y="Rate ratio, unvaccinated as reference category (log scale)") +
    scale_x_continuous(limits=c(2, 30)) + 
    scale_y_continuous(breaks=log(breakpoints),
                       labels=breakpoints,
                       limits=log(c(min(breakpoints), max(breakpoints))),
                       expand=c(0, 0)) +
    theme(legend.position = c(0.2, 0.8)) +
    guides(color = guide_legend(reverse=TRUE)) +
    guides(linetype = guide_legend(reverse=TRUE)) +
    labs(color="Vaccine product", linetype="Outcome")
p.decay

