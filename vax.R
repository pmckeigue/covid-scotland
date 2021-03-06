## vaccination analysis

vax.firstdoses <- vacc[, .N, by=c("weekdose1", "vax_type_1")]
ggplot(data=vax.firstdoses, aes(x=weekdose1, y=N, color=vax_type_1)) +
    geom_line() 

vax.varnames <- c("care.home", "COVID.age", "listedgr3", "hosp.recent", 
                  "SIMD.quintile", "hh.over18gr", "hh.schoolagegr", 
                  "occup", "vaxgr")
tabulate.freqs.regressions(data=cc.severe, outcome="CASE", varnames=vax.varnames)

controls.vax <- cc.kept[CASE==0 & care.home=="Independent", .(AGE, COVID.age, vax_dose_1)]
controls.vax[, vaxnow := !is.na(vax_dose_1)]
controls.vax[, vaxnow := car::recode(vaxnow,
                                   recodes="FALSE='Not vaccinated';
                                            TRUE='At least one dose'",
                                   as.factor=TRUE,
                                   levels=c("Not vaccinated", "At least one dose"))]

theme_set(theme_gray(base_size = 16))
p.covidage <- ggplot(data=controls.vax, aes(x=AGE, y=COVID.age, color=vaxnow)) +
    geom_point(size=0.5, position="jitter") +
    scale_x_continuous(breaks=seq(20, 80, by=10), limits=c(15, 90)) + 
    scale_y_continuous(breaks=seq(20, 80, by=10), limits=c(15, 90)) +
    theme(legend.title = element_blank()) +
    theme(legend.position=c(0.15, 0.85)) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(title=paste("COVID age versus calendar age by vaccination status on",
                        format.Date(max(vacc$vax_dose_1), '%d %b %Y')),
         caption="Controls matched to severe cases, not resident in care homes. Limits of COVID age are 16 and 90 years",
         x="Calendar age (years)",
         y="COVID age (years)")
p.covidage

png("COVIDage_vaxstatus.png")
p.covidage
dev.off()

rm(p.covidage)
