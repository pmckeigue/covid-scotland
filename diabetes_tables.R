##### Cond log regression - all cases A+B+C/ hospitalised cases A+B / severe cases A /
## Any COVID death inc death within 28 days  /severe A + all COVID Deaths 
## Restrict data to Age bands - <50, 50-70, >70 to get age specific OR
# Using ONO ethnicity 

##NRS cases are FALSE testpositive.cases - nrs_covid_case =1 for NRS deaths and matched controls
#cc.all$nrs_covid_case<-rep(0,length(cc.all[,1]))
#cc.all[cc.all$CASE==1&cc.all$testpositive.case==FALSE,]$nrs_covid_case<-1
#cc.all[cc.all$stratum%in%cc.all[cc.all$nrs_covid_case==1,]$stratum,]$nrs_covid_case<-1
#table(cc.all$CASE==1,cc.all$nrs_covid_case)


### Amend Paul's helper functions
### take restrictions off the pvalue, remove / , and more DP
combine.tables2 <- function(ftable, utable)  {# returns single table from freqs, univariate 
  
  u.ci <- or.ci(utable[, 1], utable[, 3]) 
  u.pvalue <- signif(utable[, 5], 3)
  #pvalue <- pvalue.latex(u.pvalue)
  
  table.aug <- data.frame(ftable,
                          u.ci, u.pvalue)
  return(table.aug)
}


paste.colpercent <- function(x, digits=1, escape.pct=TRUE) { # paste column percentages into freq table
  x.colpct <- paste0("(", round(100 * prop.table(x, 2), digits))
  if(escape.pct) {
    x.colpct <- paste0(x.colpct, "%)")
  } else {
    x.colpct <- paste0(x.colpct, "%)")
  }
  z <- matrix(paste(x, x.colpct), nrow=nrow(x),
              dimnames=dimnames(x))
  return(z)
}

paste.rowpercent <- function(x, digits=1, escape.pct=TRUE) { # paste column percentages into freq table
  x.colpct <- paste0("(", round(100 * prop.table(x, 1), digits))
  if(escape.pct) {
    x.colpct <- paste0(x.colpct, "%)")
  } else {
    x.colpct <- paste0(x.colpct, "%)")
  }
  z <- matrix(paste(x, x.colpct), nrow=nrow(x),
              dimnames=dimnames(x))
  return(z)
}





## Diabetes plot ---------------------------------------------------------------------------------
## In people with diabetes - plot the number of COVID-19 cases per week for each outcome
# 1. All cases - All Test Positives + NRS deaths
# 2. Severe cases (remove NRS deaths)
# 3. NRS deaths
# 4. Hospital cases not severe (remove NRS deaths)
# 5. Test Positives NOT Hospital NOT severe (remove NRS deaths)

#use death date for nrs deaths, not spec date
library(isoweek)
library(lubridate)
cc.all$tpnrswk<-cc.all$SPECDATE
cc.all[cc.all$CASE==1&cc.all$nrs_covid_case==1,]$tpnrswk<-cc.all[cc.all$CASE==1&cc.all$nrs_covid_case==1,]$Date.Death
cc.all$tpnrs.wk<-isoweek(ymd(cc.all$tpnrswk))

data1<-cc.all[cc.all$diabetes.any=="Diabetic"&cc.all$CASE==1,]
data2<-cc.all[cc.all$diabetes.any=="Diabetic"&cc.all$CASE==1&cc.all$casegroup=="A"&cc.all$nrs_covid_case!=1,]
data3<-cc.all[cc.all$diabetes.any=="Diabetic"&cc.all$CASE==1&cc.all$nrs_covid_case==1,]
data4<-cc.all[cc.all$diabetes.any=="Diabetic"&cc.all$CASE==1&cc.all$casegroup=="B"&cc.all$nrs_covid_case!=1,]
data5<-cc.all[cc.all$diabetes.any=="Diabetic"&cc.all$CASE==1&cc.all$casegroup=="C"&cc.all$nrs_covid_case!=1,]

plot(table(data1$tpnrs.wk), ylim=c(0,500), type="l",col="red"
     ,xlab=c("Day/Month 2020"), ylab=c("Number of COVID-19 cases among those with diabetes"),xaxt = "n")
lines(table(data2$tpnrs.wk), type="l",col="blue")
lines(table(data3$tpnrs.wk), type="l",col="green")
lines(table(data4$tpnrs.wk), type="l",col="orange")
lines(table(data5$tpnrs.wk), type="l",col="magenta")
legend(19, 500, legend=c("All COVID", "Severe COVID","NRS COVID Deaths","Hospitalised, Not Severe","Test Positive, Not Hospitlalised"),
       col=c("red", "blue","green","orange","magenta"), lty=1, cex=0.7)
axis(1, at=10:24, labels=c("02/03","09/03","16/03","23/03","30/03",
                           "06/04","13/04","20/04","27/04","04/05","11/05",
                           "18/05","25/05","01/06","0806"))



### Diabetes tables ----------------------------------------------------------------------------------------------------------


# Recode AGE into Stuarts age bands
cc.all$agegr3 <- as.factor(car::recode(as.integer(cc.all$AGE),
                                       "0:49='0-49 <50'; 50:70='50-70';
                                       71:hi='71 or more >70'"))



vars1ono<-c("ethnic3.onomap","SIMD.quintile","care.home","dm.type") #"HAI"

cc.all$casegroup2 <- car::recode(cc.all$casegroup,
                                 "'A'='Critical care or fatal'; 'B'='Hospitalised, not severe'; 'C'='Test-positive, not hospitalised';
                                 'D'='Critical care or fatal'")


#ALL TP cases #--------------------------------------------------------------------------------------

data.dm<-cc.all[cc.all$nrs_covid_case!=1,]           # 1. all test positives (A+B+C, inclues all deaths) #remove NRS deaths/controls

table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoALL<-table.ono

# uni/multi reg##-------------

# All ages

table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm,   #            data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
AllAge<-temp4; colnames(AllAge)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")


#Age <50-----------------------------------------------------------------------------------------------
data.dm<-cc.all[cc.all$nrs_covid_case!=1&cc.all$agegr3=="0-49 <50",]           # 1. all test positives (A+B+C, inclues all deaths) #remove NRS deaths/controls

table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoALL50<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm,   #data=data.dm[!is.na(data.dm$ethnic3.onomap), ]
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
AllAge50<-temp4; colnames(AllAge50)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")

#Age 50-70 -----------------------------------------------------------------------------------------------
data.dm<-cc.all[cc.all$nrs_covid_case!=1&cc.all$agegr3=="50-70",]           # 1. all test positives (A+B+C, inclues all deaths) #remove NRS deaths/controls
table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoALL5070<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm,   #data=data.dm[!is.na(data.dm$ethnic3.onomap), ]
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
AllAge5070<-temp4; colnames(AllAge5070)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")


#Age >70-----------------------------------------------------------------------------------------------
data.dm<-cc.all[cc.all$nrs_covid_case!=1&cc.all$agegr3=="71 or more >70",]           # 1. all test positives (A+B+C, inclues all deaths) #remove NRS deaths/controls
table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoALL70<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm,  # data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
AllAge70<-temp4; colnames(AllAge70)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")

Allunimulti<-rbind(AllAge,AllAge50,AllAge5070,AllAge70)
rownames(Allunimulti)<-c("Not Diabetic","Type1", "Type 2","Other/unknown type",
                         "Not Diabetic <50","Type1 <50", "Type 2 <50","Other/unknown type <50",
                         "Not Diabetic 50-70","Type1 50-70", "Type 2 50-70","Other/unknown type 50-70",
                         "Not Diabetic >70","Type1 >70", "Type 2 >70","Other/unknown type >70")

colnames(Allunimulti)<-c("Controls","Cases","Uni OR (95% CI)","P-value","Multi1 OR (95% CI)","P-value","Multi2 OR (95% CI)","P-value")

Allunimulti[,4]<-ifelse(round(Allunimulti[,4],3)<0.001, paste("<0.001"), round(Allunimulti[,4],3))
Allunimulti[,6]<-ifelse(round(Allunimulti[,6],3)<0.001, paste("<0.001"), round(Allunimulti[,6],3))
Allunimulti[,8]<-ifelse(round(Allunimulti[,8],3)<0.001, paste("<0.001"), round(Allunimulti[,8],3))

write.csv(Allunimulti, "./Amanda files/DiabetesUniMultiAgeALL_180620extract_260620.csv")


# Hospitalised or severe cases #----------------------------------------------------------------------------------------------------

data.dm<-cc.all[(cc.all$casegroup2=="Hospitalised, not severe" | cc.all$casegroup2=="Critical care or fatal")&cc.all$nrs_covid_case!=1,]           # 2.hospitalised(A+B, inclues deaths>28d), remove nrs deaths

table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoALLHOSP<-table.ono

# uni/multi reg##-------------

# All ages

table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm, #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
AllAgeHOSP<-temp4; colnames(AllAgeHOSP)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")



#Age <50-----------------------------------------------------------------------------------------------
data.dm<-cc.all[(cc.all$casegroup2=="Hospitalised, not severe" | cc.all$casegroup2=="Critical care or fatal")&cc.all$nrs_covid_case!=1&cc.all$agegr3=="0-49 <50",]# 2.hospitalised(A+B, inclues deaths>28d), remove nrs deaths

table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoHOSP50<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm, #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
HospAge50<-temp4; colnames(HospAge50)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")

#Age 50-70-----------------------------------------------------------------------------------------------
data.dm<-cc.all[(cc.all$casegroup2=="Hospitalised, not severe" | cc.all$casegroup2=="Critical care or fatal")&cc.all$nrs_covid_case!=1&cc.all$agegr3=="50-70",]# 2.hospitalised(A+B, inclues deaths>28d), remove nrs deaths

table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoHOSP5070<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm, #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
HospAge5070<-temp4; colnames(HospAge5070)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")


#Age >70-----------------------------------------------------------------------------------------------
data.dm<-cc.all[(cc.all$casegroup2=="Hospitalised, not severe" | cc.all$casegroup2=="Critical care or fatal")&cc.all$nrs_covid_case!=1&cc.all$agegr3=="71 or more >70",]           # 2.hospitalised(A+B, inclues deaths>28d), remove nrs deaths
table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoHOSP70<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm, #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
HospAge70<-temp4; colnames(HospAge70)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")

Hospunimulti<-rbind(AllAgeHOSP, HospAge50,HospAge5070,HospAge70)
rownames(Hospunimulti)<-c("Not Diabetic","Type1", "Type 2","Other/unknown type",
                          "Not Diabetic <50","Type1 <50", "Type 2 <50","Other/unknown type <50",
                          "Not Diabetic 50-70","Type1 50-70", "Type 2 50-70","Other/unknown type 50-70",
                          "Not Diabetic >70","Type1 >70", "Type 2 >70","Other/unknown type >70")

colnames(Hospunimulti)<-c("Controls","Cases","Uni OR (95% CI)","P-value","Multi1 OR (95% CI)","P-value","Multi2 OR (95% CI)","P-value")

Hospunimulti[,4]<-ifelse(round(Hospunimulti[,4],3)<0.001, paste("<0.001"), round(Hospunimulti[,4],3))
Hospunimulti[,6]<-ifelse(round(Hospunimulti[,6],3)<0.001, paste("<0.001"), round(Hospunimulti[,6],3))
Hospunimulti[,8]<-ifelse(round(Hospunimulti[,8],3)<0.001, paste("<0.001"), round(Hospunimulti[,8],3))

write.csv(Hospunimulti, "./Amanda files/DiabetesUniMultiAgeHOSP_180620extract_260620.csv")


# Severe cases #------------------------------------------------------------------------------------------------

data.dm<-cc.all[cc.all$casegroup2=="Critical care or fatal"&cc.all$nrs_covid_case!=1,]           #3. severe and no nRS deaths

table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoALLSEV<-table.ono

# uni/multi reg##-------------

# All ages

table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm, #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
AllAgeSEV<-temp4; colnames(AllAgeSEV)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")


#Age <50-----------------------------------------------------------------------------------------------
data.dm<-cc.all[cc.all$casegroup2=="Critical care or fatal"&cc.all$nrs_covid_case!=1&cc.all$agegr3=="0-49 <50",]# 3. severe and no nRS deaths

table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoSEV50<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm, #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
SevAge50<-temp4; colnames(SevAge50)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")

#Age 50-70-----------------------------------------------------------------------------------------------
data.dm<-cc.all[cc.all$casegroup2=="Critical care or fatal"&cc.all$nrs_covid_case!=1&cc.all$agegr3=="50-70",]# 3. severe no NRS deaths

table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoSEV5070<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm,  #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
SevAge5070<-temp4; colnames(SevAge5070)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")


#Age >70-----------------------------------------------------------------------------------------------
data.dm<-cc.all[cc.all$casegroup2=="Critical care or fatal"&cc.all$nrs_covid_case!=1&cc.all$agegr3=="71 or more >70",]           # 3. severe no NRS deaths
table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoSEV70<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm,  #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
SevAge70<-temp4; colnames(SevAge70)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")

Sevunimulti<-rbind(AllAgeSEV, SevAge50,SevAge5070,SevAge70)
rownames(Sevunimulti)<-c("Not Diabetic","Type1", "Type 2","Other/unknown type",
                         "Not Diabetic <50","Type1 <50", "Type 2 <50","Other/unknown type <50",
                         "Not Diabetic 50-70","Type1 50-70", "Type 2 50-70","Other/unknown type 50-70",
                         "Not Diabetic >70","Type1 >70", "Type 2 >70","Other/unknown type >70")

colnames(Sevunimulti)<-c("Controls","Cases","Uni OR (95% CI)","P-value","Multi1 OR (95% CI)","P-value","Multi2 OR (95% CI)","P-value")

Sevunimulti[,4]<-ifelse(round(Sevunimulti[,4],3)<0.001, paste("<0.001"), round(Sevunimulti[,4],3))
Sevunimulti[,6]<-ifelse(round(Sevunimulti[,6],3)<0.001, paste("<0.001"), round(Sevunimulti[,6],3))
Sevunimulti[,8]<-ifelse(round(Sevunimulti[,8],3)<0.001, paste("<0.001"), round(Sevunimulti[,8],3))

write.csv(Sevunimulti, "./Amanda files/DiabetesUniMultiAgeSEV_180620extract_260620.csv")


# Any case with a COVID death #------------------------------------------------------------------------------------------------

data.dm<-cc.all[cc.all$stratum%in%cc.all[cc.all$CASE==1&!(is.na(cc.all$Date.Death))&!(cc.all$deathwithin28==FALSE&cc.all$covid_cod==0),]$stratum,]# 4. any COVID death - all <28 days counted as COVID and >28 days with COVID as COD

table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoALLDEATH<-table.ono

# uni/multi reg##-------------

# All ages

table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm,  #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
AllAgeDEATH<-temp4; colnames(AllAgeDEATH)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")



#Age <50-----------------------------------------------------------------------------------------------
data.dm<-cc.all[cc.all$stratum%in%cc.all[cc.all$CASE==1&!(is.na(cc.all$Date.Death))&!(cc.all$deathwithin28==FALSE&cc.all$covid_cod==0),]$stratum&cc.all$agegr3=="0-49 <50",]# 4. any COVID death

table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoDEATH50<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm,  #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
DeathAge50<-temp4; colnames(DeathAge50)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")

#Age 50-70-----------------------------------------------------------------------------------------------
data.dm<-cc.all[cc.all$stratum%in%cc.all[cc.all$CASE==1&!(is.na(cc.all$Date.Death))&!(cc.all$deathwithin28==FALSE&cc.all$covid_cod==0),]$stratum&cc.all$agegr3=="50-70",]# 3. severe no NRS deaths

table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoDEATH5070<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm,  #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
DeathAge5070<-temp4; colnames(DeathAge5070)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")


#Age >70-----------------------------------------------------------------------------------------------
data.dm<-cc.all[cc.all$stratum%in%cc.all[cc.all$CASE==1&!(is.na(cc.all$Date.Death))&!(cc.all$deathwithin28==FALSE&cc.all$covid_cod==0),]$stratum&cc.all$agegr3=="71 or more >70",]           # 3. severe no NRS deaths
table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoDEATH70<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm, #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
DeathAge70<-temp4; colnames(DeathAge70)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")

Deathunimulti<-rbind(AllAgeDEATH,DeathAge50,DeathAge5070,DeathAge70)
rownames(Deathunimulti)<-c("Not Diabetic","Type1", "Type 2","Other/unknown type",
                           "Not Diabetic <50","Type1 <50", "Type 2 <50","Other/unknown type <50",
                           "Not Diabetic 50-70","Type1 50-70", "Type 2 50-70","Other/unknown type 50-70",
                           "Not Diabetic >70","Type1 >70", "Type 2 >70","Other/unknown type >70")

colnames(Deathunimulti)<-c("Controls","Cases","Uni OR (95% CI)","P-value","Multi1 OR (95% CI)","P-value","Multi2 OR (95% CI)","P-value")

Deathunimulti[,4]<-ifelse(round(Deathunimulti[,4],3)<0.001, paste("<0.001"), round(Deathunimulti[,4],3))
Deathunimulti[,6]<-ifelse(round(Deathunimulti[,6],3)<0.001, paste("<0.001"), round(Deathunimulti[,6],3))
Deathunimulti[,8]<-ifelse(round(Deathunimulti[,8],3)<0.001, paste("<0.001"), round(Deathunimulti[,8],3))

write.csv(Deathunimulti, "./Amanda files/DiabetesUniMultiAgeDEATH_180620extract_260620.csv")



# Severe case plus ANY COVID deaths #------------------------------------------------------------------------------------------------

data.dm<-cc.all[(cc.all$casegroup2=="Critical care or fatal" | cc.all$stratum%in%cc.all[cc.all$CASE==1&!(is.na(cc.all$Date.Death))&!(cc.all$deathwithin28==FALSE&cc.all$covid_cod==0),]$stratum),]# 4. any COVID death


table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoALLSEVDEATH<-table.ono

# uni/multi reg##-------------

# All ages

table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm, #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic", "Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
AllAgeSEVDEATH<-temp4; colnames(AllAgeSEVDEATH)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")



#Age <50-----------------------------------------------------------------------------------------------
data.dm<-cc.all[(cc.all$casegroup2=="Critical care or fatal" | cc.all$stratum%in%cc.all[cc.all$CASE==1&!(is.na(cc.all$Date.Death))&!(cc.all$deathwithin28==FALSE&cc.all$covid_cod==0),]$stratum)&cc.all$agegr3=="0-49 <50",]# 4. any COVID death

table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoSEVDEATH50<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm,  #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
SevDEATHAge50<-temp4; colnames(SevDEATHAge50)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")

#Age 50-70-----------------------------------------------------------------------------------------------
data.dm<-cc.all[(cc.all$casegroup2=="Critical care or fatal" | cc.all$stratum%in%cc.all[cc.all$CASE==1&!(is.na(cc.all$Date.Death))&!(cc.all$deathwithin28==FALSE&cc.all$covid_cod==0),]$stratum)&cc.all$agegr3=="50-70",]# 4. any COVID death

table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoSEVDEATH5070<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm,  #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
SevDEATHAge5070<-temp4; colnames(SevDEATHAge5070)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")


#Age >70-----------------------------------------------------------------------------------------------
data.dm<-cc.all[(cc.all$casegroup2=="Critical care or fatal" | cc.all$stratum%in%cc.all[cc.all$CASE==1&!(is.na(cc.all$Date.Death))&!(cc.all$deathwithin28==FALSE&cc.all$covid_cod==0),]$stratum)&cc.all$agegr3=="71 or more >70",]# 4. any COVID death

table.ono <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[1]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
for(i in 2:length(vars1ono)){
  x <- paste.colpercent(table(data.dm[!is.na(data.dm$ethnic3.onomap), ][[vars1ono[i]]], data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE), 1)
  table.ono <- rbind(table.ono, x)
}
colnames(table.ono)[1:2] <- paste(c("Control", "Case"),
                                  gsub("([0-9]+)", "\\(N = \\1\\)",
                                       as.integer(table(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE[!is.na(data.dm[!is.na(data.dm$ethnic3.onomap), ]$CASE)]))))
table.onoSEVDEATH70<-table.ono
# uni/multi reg##-------------
table.dm.ono <- univariate.tabulate(varnames=c("dm.type"), outcome="CASE",
                                    data=data.dm,  #data.dm[!is.na(data.dm$ethnic3.onomap), ],
                                    drop.reflevel=FALSE)
univariate.dm.ono <-
  univariate.clogit(varnames=c("dm.type"),
                    data=data.dm[!is.na(data.dm$ethnic3.onomap), ],
                    add.reflevel=TRUE)
table.dm.aug.ono <- combine.tables2(table.dm.ono, univariate.dm.ono)
rownames(table.dm.aug.ono) <- replace.names(rownames(table.dm.aug.ono))
multivariate.dm.ono1 <-
  multivariate.clogit(varnames=c(vars1ono),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
multivariate.dm.ono2 <-
  multivariate.clogit(varnames=c(vars1ono,"Ch.9_circulatory"),
                      data=data.dm[!is.na(data.dm$ethnic3.onomap),], add.reflevel=TRUE)
temp1<-multivariate.dm.ono1[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp2<-multivariate.dm.ono2[c("Not diabetic","Type 1 diabetes","Type 2 diabetes","Other/unknown type"),]
temp3<-combine.tables2(table.dm.aug.ono,temp1)
temp4<-combine.tables2(temp3, temp2)
SevDEATHAge70<-temp4; colnames(SevDEATHAge70)<-c("Controls","Cases","u.ci","u.pvalue","m1.ci","m1.pvalue","m2.ci","m2.pvalue")

SevDEATHunimulti<-rbind(AllAgeSEVDEATH,SevDEATHAge50,SevDEATHAge5070,SevDEATHAge70)
rownames(SevDEATHunimulti)<-c("Not Diabetic","Type1", "Type 2","Other/unknown type",
                              "Not Diabetic <50","Type1 <50", "Type 2 <50","Other/unknown type <50",
                              "Not Diabetic 50-70","Type1 50-70", "Type 2 50-70","Other/unknown type 50-70",
                              "Not Diabetic >70","Type1 >70", "Type 2 >70","Other/unknown type >70")

colnames(SevDEATHunimulti)<-c("Controls","Cases","Uni OR (95% CI)","P-value","Multi1 OR (95% CI)","P-value","Multi2 OR (95% CI)","P-value")

SevDEATHunimulti[,4]<-ifelse(round(SevDEATHunimulti[,4],3)<0.001, paste("<0.001"), round(SevDEATHunimulti[,4],3))
SevDEATHunimulti[,6]<-ifelse(round(SevDEATHunimulti[,6],3)<0.001, paste("<0.001"), round(SevDEATHunimulti[,6],3))
SevDEATHunimulti[,8]<-ifelse(round(SevDEATHunimulti[,8],3)<0.001, paste("<0.001"), round(SevDEATHunimulti[,8],3))

write.csv(SevDEATHunimulti, "./Amanda files/DiabetesUniMultiAgeSEVDEATH_180620extract_260620.csv")

