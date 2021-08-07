library(data.table)
library(ggplot2)

source("helperfunctions.R")

year.to2020 <- function(x) {
    lubridate::year(x) <- 2020
    return(x)
}

year.to2021 <- function(x) {
    lubridate::year(x) <- 2021
    return(x)
}

linkdate <- "jul12"
linkdate <- "jul28"

if(linkdate == "jul12") {
        datadir <- "./data/2021-07-12/"
        ct.filename <-  paste0(datadir, "CC_CT_Test_Results_2021-07-12.rds")
        ecoss.tests.filename <- paste0(datadir, "CC_ecoss_tests_2021-07-12.rds")
} else if(linkdate=="jul28") {
    datadir <- "./data/2021-07-28/"
    ct.filename <-  paste0(datadir, "CC_CT_Test_Results_2021-07-28.rds")
    ecoss.tests.filename <- paste0(datadir, "CC_ecoss_tests_2021-07-28.rds")
}

################## read Ct table ####################################

ct <- RDStodt(ct.filename) # all rows are unique
## test_ch[1-4]_result_ct_value variables have been miscoded as character
ctvalue.cols <- grep("test_ch[1-4]_result_ct_value", names(ct), value=TRUE)
ct[, (ctvalue.cols) := lapply(.SD, as.numeric), .SDcols=ctvalue.cols]

## drop records that do not have a result on at least one PCR channel
ct <- ct[!is.na(ct$test_ch1_target_gene_result) |
         !is.na(ct$test_ch2_target_gene_result) |
         !is.na(ct$test_ch3_target_gene_result) |
         !is.na(ct$test_ch4_target_gene_result) ]
cat(nrow(ct), "records in ct with at least one non-empty target gene result\n")
setnames(ct, "test_result", "ct.result")

## jan 6 linkage has no variable date_appointment
setnames(ct, "date_ecoss_specimen", "date_appointment", skip_absent=TRUE)

ct[, date_appointment := as.Date(date_appointment)] # mostly missing
ct[, date_onset_of_symptoms := as.Date(date_onset_of_symptoms)] 
ct[, date_reporting := as.Date(date_reporting)] # only 3 missing values and these have all fields missing except ID

#print(ct[date_appointment > Sys.Date(), .(date_appointment, date_onset_of_symptoms, date_reporting)])
#print(ct[date_appointment < as.Date("2020-03-01"), .(date_appointment, date_onset_of_symptoms, date_reporting)])

ct[date_appointment > Sys.Date() & lubridate::month(date_appointment) > 9,
   date_appointment := year.to2020(date_appointment)]

ct[date_appointment < as.Date("2020-03-01") & lubridate::month(date_appointment) < 3,
   date_appointment := year.to2021(date_appointment)]

ct[date_appointment > Sys.Date() | date_appointment < as.Date("2020-03-01"),
   date_appointment := NA]

ct[date_onset_of_symptoms < as.Date("2020-02-01"), date_onset_of_symptoms := NA]

ct[TRUE, ct.id := .I] # set an ID for the record in the Ct table

############### read ECOSS tests table ####################

ecoss <- RDStodt(ecoss.tests.filename, keyname="anon_id")
ecoss <- unique(ecoss) ## about 1% of rows are duplicated
setnames(ecoss, "SPECIMENDATE", "SpecimenDate_ecosswrong", skip_absent=TRUE) # where did this field come from? 
setnames(ecoss, "date_ecoss_specimen", "SpecimenDate", skip_absent=TRUE)
setnames(ecoss, "Type", "ecoss.result", skip_absent=TRUE)
setnames(ecoss, "test_result", "ecoss.result", skip_absent=TRUE)
setnames(ecoss, "ecoss_submitting_laboratory", "SourceLab", skip_absent=TRUE)
ecoss[, SpecimenDate := as.Date(SpecimenDate)]

## SpecimenType field is tissue / compartment sampled
## Tests field is free text comments
## Category field contains occupational groups care home etc. 
ecoss[, ecoss.result := car::recode(ecoss.result,
                                    "'INSUFFICIENT'=NA;
                                         'VOID'=NA;
                                         'NEGATIVE'='Negative';
                                         'POSITIVE'='Positive'", 
                                    levels=c("Positive", "Negative"))]
## order levels Positive before Negative so that sorting in ascending order of level and dropping duplicated anon_ids will drop negatives
## many of those with only negative results are test-positive cases 
#ecoss[!(anon_id %in% cc.all$anon_id), .N, by=ecoss.result]
#ecoss <- ecoss[anon_id %in% cc.all$anon_id] # all unmatched IDs have negative or missing result
ecoss[SourceLab=="", SourceLab := "Not recorded in ECOSS"]

ecoss.pos <- ecoss[ecoss.result=="Positive"]
setkey(ecoss.pos, SpecimenDate)
ecoss.pos <- ecoss.pos[!duplicated(anon_id)]

setkey(ecoss, anon_id, SpecimenDate)

#############################################################

cat(nrow(ct), "records in Ct table\n")
cat(length(unique(ct$anon_id)), "unique anon_ids\n")
cat("Table of Ct anon_ids in ECOSS", table(ct$anon_id %in% ecoss$anon_id), "\n")
cat("Table of Ct anon_ids in ECOSS test-positives",
    table(ct$anon_id %in% ecoss[ecoss.result=="Positive", anon_id]), "\n")

#cat("Table of test-positive cases (defined in cc.all by exclusion of other ascertainment) in ECOSS test-positives",
#    table(cc.all[testpositive.case==TRUE, anon_id] %in%
#          ecoss[ecoss.result=="Positive", anon_id]), "\n")

## 4 channels: genes are ORF1ab, N, S, MS2
## all negative results have Ct value coded missing

genes <- c("ORF1ab", "N", "S", "MS2") # MS2 is a bacteriophage QC
ct.newnames <- paste0(genes, "_Ct")
test.newnames <- paste0(genes, "_result")

ct.names <- grep("value$", names(ct), value=TRUE)
test.names <- grep("^test_ch.+result$", names(ct), value=TRUE)
setnames(ct, old=ct.names, new=ct.newnames)
setnames(ct, old=test.names, new=test.newnames)

for(j in ct.newnames) set(x=ct, i=NULL, j = j, value = as.numeric(ct[[j]]))

ct[N_Ct > 50, N_Ct := NA]
ct[S_Ct > 50, S_Ct := NA]
ct[ORF1ab_Ct > 50, ORF1ab_Ct := NA]

ct[, av2channels := 0.5 * (ORF1ab_Ct + N_Ct)]
ct[, diff2channels := N_Ct - ORF1ab_Ct]

## could rewrite this more succinctly using set
ct[ORF1ab_result=="", ORF1ab_result := NA] 
ct[N_result=="", N_result := NA] 
ct[S_result=="" | S_result == "INCONCLUSIVE", S_result := NA]
ct[MS2_result=="", MS2_result := NA] 

## set negative to 40 so that medians can be calculated
ct[ORF1ab_result=="NEGATIVE", ORF1ab_Ct := 40] 
ct[N_result=="NEGATIVE", N_Ct := 40] 
ct[S_result=="NEGATIVE", S_Ct := 40] 
ct[MS2_result=="NEGATIVE", MS2_Ct := 40] 

## sgtf.orig defined as S gene negative & ORF1ab Ct < 31  & N Ct <31
ct[S_result=="POSITIVE", sgtf.orig := "No dropout"]
ct[S_result=="NEGATIVE" & ORF1ab_result=="POSITIVE" &  N_result=="POSITIVE" &
   ORF1ab_Ct < 31 & N_Ct < 31,  sgtf.orig := "Definite dropout"]
ct[S_result=="NEGATIVE" & ORF1ab_result=="POSITIVE" & N_result=="POSITIVE" &
   (ORF1ab_Ct >= 31 | N_Ct >= 31),
   sgtf.orig := "Undetermined"]
ct[, sgtf.orig := factor(sgtf.orig,
                         levels = c("No dropout", "Undetermined", "Definite dropout"))]

p.dropoutct <- ggplot(data=ct[ORF1ab_result=="POSITIVE" & N_result=="POSITIVE" & !is.na(S_result)],
       aes(x=ORF1ab_Ct, y=N_Ct, color=sgtf.orig)) +
    geom_point() +
    scale_color_manual(values=c("blue", "orange", "red")) +
    scale_x_reverse(breaks=c(40, 30, 20, 10),
                    labels=c("40", "30", "20", "10")) +
    scale_y_reverse(breaks=c(40, 30, 20, 10),
                    labels=c("40", "30", "20", "10")) +
    geom_segment(aes(x=31, y=31, xend=31, yend=10), linetype=3, color="black") +
    geom_segment(aes(x=31, y=31, xend=10, yend=31), linetype=3, color="black") +
    geom_abline(slope=1, intercept=-1, color="green") +
    ggtitle("Relation of S gene dropout (original definition with both other channels < 31) to ORF and N gene Ct values")
p.dropoutct

## sgtf.fix1 restricts definite dropout to those with N - ORF > -1.5
ct[, sgtf.fix1 := sgtf.orig]
ct[ORF1ab_Ct - N_Ct < -1.5 & sgtf.orig=="Definite dropout", sgtf.fix1 := "Undetermined"]
ct[, sgtf.fix1 := factor(sgtf.fix1,
                         levels = c("No dropout", "Undetermined", "Definite dropout"))]

with(ct, table(sgtf.orig, sgtf.fix1))
             
p.sepline <- ggplot(data=ct[ORF1ab_result=="POSITIVE" & N_result=="POSITIVE" & !is.na(S_result)],
       aes(x=ORF1ab_Ct, y=N_Ct, color=sgtf.fix1)) +
    geom_point() +
    scale_color_manual(values=c("blue", "orange", "red")) +
    scale_x_reverse(breaks=c(40, 30, 20, 10),
                    labels=c("40", "30", "20", "10")) +
    scale_y_reverse(breaks=c(40, 30, 20, 10),
                    labels=c("40", "30", "20", "10")) +
    geom_segment(aes(x=31, y=31, xend=29.5, yend=31), linetype=3, color="black") +
    geom_segment(aes(x=31, y=31, xend=31, yend=10), linetype=3, color="black") +
    geom_segment(aes(x=29.5, y=31, xend=8.5, yend=10), linetype=3, color="black") +
    ggtitle("Separation line at N - ORF = 1.5")
p.sepline

## use sgtf.fix1 to assign true neg / pos status using replicates where possible 
ct[, max.sgtf := max(as.integer(sgtf.fix1), na.rm=TRUE), by=anon_id]
ct[, min.sgtf := min(as.integer(sgtf.fix1), na.rm=TRUE), by=anon_id]
ct[is.infinite(max.sgtf), max.sgtf := NA]
ct[is.infinite(-min.sgtf), min.sgtf := NA]
ct[max.sgtf==3, truepos := "Other test definite dropout"]
ct[min.sgtf==1, truepos := "Other test no dropout"]

p.testvalid <-ggplot(data=ct[sgtf.fix1=="Undetermined" & !is.na(truepos)], 
       aes(x=ORF1ab_Ct, y=N_Ct, color=truepos)) +
    geom_point() +
    scale_color_manual(values=c("red", "blue")) +
    scale_x_reverse(breaks=c(40, 30, 20, 10),
                    labels=c("40", "30", "20", "10")) +
    scale_y_reverse(breaks=c(40, 30, 20, 10),
                    labels=c("40", "30", "20", "10")) +
    # horiz line from 33 to 28 at y=31
    geom_segment(aes(x=33, y=31, xend=28, yend=31), linetype=3, color="black") +
    # vertical line from 31 to 10 at x=31
    geom_segment(aes(x=31, y=31, xend=31, yend=10), linetype=3, color="black") +
    geom_segment(aes(x=29.5, y=31, xend=8.5, yend=10), linetype=3, color="black") +
    # vertical line from 31 to 10 at x=31
    geom_segment(aes(x=33, y=31, xend=33, yend=10), color="green") +
    geom_segment(aes(x=28, y=31, xend=7, yend=10), color="green") +
    ggtitle("Tests with S gene dropout classified as Undetermined with separation line at N - ORF = 1.5, by status at repeat test. New separation lines at N - ORF = 3 and ORF = 33")
p.testvalid

## assign S_gene dropout using difference between N and ORF signals
## refine this by restricting diff2channels to <=2 and N Ct to < 30
## these values were chosen after plotting test results by truepos status

ct[, Sgene.dropout := sgtf.fix1]
ct[S_result=="NEGATIVE" & ORF1ab_result=="POSITIVE" &  N_result=="POSITIVE" &
   ORF1ab_Ct < 33 & N_Ct < 31 & diff2channels <= 3, Sgene.dropout := "Definite dropout"]

ct[S_result=="NEGATIVE" & ORF1ab_result=="POSITIVE" & N_result=="POSITIVE" &
   (ORF1ab_Ct >= 33 | N_Ct >= 31 | diff2channels > 3),
   Sgene.dropout := "Undetermined"]

ct[, Sgene.dropout := factor(Sgene.dropout,
                             levels = c("No dropout", "Undetermined", "Definite dropout"))]
ct[is.na(S_result), Sgene.dropout := NA]

# with(ct, print(table(sgtf.orig, Sgene.dropout, exclude=NULL)))

## use sgene.dropout to generate variables that can be used to reassign undetermined values where possible 
ct[, max.Sgene.dropout := max(as.integer(Sgene.dropout), na.rm=TRUE), by=anon_id]
ct[, min.Sgene.dropout := min(as.integer(Sgene.dropout), na.rm=TRUE), by=anon_id]
ct[is.infinite(max.Sgene.dropout), max.Sgene.dropout := NA]
ct[is.infinite(-min.Sgene.dropout), min.Sgene.dropout := NA]
ct[max.Sgene.dropout==3, truepos := "Other test definite dropout"]
ct[min.Sgene.dropout==1, truepos := "Other test no dropout"]

########################  merge Ct with ECOSS ###########################################

## this renames the date_appointment field in Ct with SpecimenDate
# cat("Table of missingness status for original date_appointment in Ct table: \n")
# print(table(is.na(ct$date_appointment)))
## 177 appointment dates are missing

ct.nodate <- ct[is.na(date_appointment), .(anon_id, ct.id, date_appointment, date_reporting)]
ecoss.withdate <- ecoss[SourceLab == "NHS:COV", .(anon_id, SpecimenDate, SourceLab)]
setkey(ecoss.withdate, anon_id)
setkey(ct.nodate, anon_id)

## left join ct.nodate with ecoss.withdate
ct.nodate <- ecoss.withdate[ct.nodate] 

## restrict to no more than 14 days from SpecimenDate to date_reporting
ct.nodate <- ct.nodate[date_reporting >= SpecimenDate & date_reporting - SpecimenDate < 14]
ct.nodate <- ct.nodate[order(date_reporting - SpecimenDate), ]
ct.nodate <- ct.nodate[!duplicated(ct.id)]
ct.nodate <- ct.nodate[, .(ct.id, SpecimenDate)]
setnames(ct.nodate, "SpecimenDate", "matched.date_reporting")

## now merge these records back into ct
setkey(ct, ct.id)
setkey(ct.nodate, ct.id)
## left join ct with ct.unmatched on ct.id
ct <- ct.nodate[ct]
## records not matched on SpecimenDate in ECOSS have missing SourceLab
## for records with missing SpecimenDate and missing SourceLab, impute the SpecimenDate from ecossSpecimenDate 
ct[is.na(date_appointment), date_appointment := matched.date_reporting]

cat("Table of missingness status for date_appointment in Ct table after matching in ECOSS: \n")
print(table(is.na(ct$date_appointment)))
## now only 12 missing
setnames(ct, "date_appointment", "SpecimenDate") 
maxdate.ct <- max(ct$SpecimenDate, na.rm=TRUE)

setkey(ct, anon_id, SpecimenDate) ## ct contains only unique values of key
setkey(ecoss, anon_id, SpecimenDate) ## about 1% of key values are duplicated

## this step imports the flag_lighthouse_labs_testing variable into ecoss
ecoss <- ct[, .(anon_id, SpecimenDate, flag_lighthouse_labs_testing, ct.result)][ecoss]

ecoss[, flag_lighthouse_labs_testing := factor(car::recode(flag_lighthouse_labs_testing, 
                                                      "NA='NHS test';
                                                       0='Ct record, not Lighthouse';
                                                       1='Lighthouse test'"))]

cat("ECOSS results by lighthouse / NHS status and month\n") 
ecoss.table <-
    with(ecoss, 
         table(flag_lighthouse_labs_testing, lubridate::month(SpecimenDate, label=TRUE))
         )
print(paste.colpercent(ecoss.table)[, c(9:12, 1:6)])

save(ecoss, file=paste0(datadir, "ecoss.RData"))
save(ct, file=paste0(datadir, "ct.RData"))


rm(ecoss.withdate)





 
