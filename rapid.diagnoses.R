library(data.table)
source("helperfunctions.R")

datadir <- "./data/2021-09-22/"

rapid.filename <- paste0(datadir, "CC_RAPID_ANON_2021-09-22.rds")
diagnoses.filename <- paste0(datadir, "CC_SMR01_ICD10_2021-09-22.rds")

rapid <- RDStodt(rapid.filename, keyname="anon_id")
setnames(rapid, "admission_date", "AdmissionDate.rapid", skip_absent=TRUE)
setnames(rapid, "discharge_date", "DischargeDate.rapid", skip_absent=TRUE)
rapid <- unique(rapid, by=c("anon_id", "AdmissionDate.rapid", "DischargeDate.rapid"))

diagnoses <- RDStodt(diagnoses.filename, keyname="anon_id")
diagnoses <- unique(diagnoses)
diagnoses <- diagnoses[inpatient_daycase_identifier == "I"]
setnames(diagnoses, "admission_date", "admissiondate", skip_absent=TRUE)
setnames(diagnoses, "discharge_date", "dischargedate", skip_absent=TRUE)
diagnoses[, admissiondate := as.Date(admissiondate)]
diagnoses[, dischargedate := as.Date(dischargedate)]
diagnoses <- diagnoses[admissiondate >= as.Date("2021-03-01")]
## set missing discharge dates to last date in table
lastdate.diagnoses <-  max(c(diagnoses$admissiondate, diagnoses$dischargedate), na.rm=TRUE)
diagnoses[is.na(dischargedate), dischargedate := lastdate.diagnoses]

## each record with same anon_id and cis_marker is a diagnosis (?or procedure)
## records with same anon_id, admission date, and discharge date are episodes
## records with same cis_marker are stays
setorder(diagnoses, anon_id, admissiondate, dischargedate)

## number episodes within stays
diagnoses[, episode := as.integer(c(0, diff(dischargedate)) > 0), by=c("anon_id", "cis_marker")]
diagnoses[, episode := 1 + cumsum(episode), by=c("anon_id", "cis_marker")]

#diagnoses[, numepisodes := max(episode), by=c("anon_id", "cis_marker")]

diagnoses[, AdmissionDate.rapid := min(admissiondate), by=c("anon_id", "cis_marker")]
diagnoses[, DischargeDate.rapid := max(dischargedate), by=c("anon_id", "cis_marker")]


episode1 <- unique(diagnoses[episode==1, .(anon_id, AdmissionDate.rapid, DischargeDate.rapid,
                                           specialty)])
setorder(episode1, specialty)
episode1 <- unique(episode1, by=c("anon_id", "AdmissionDate.rapid", "DischargeDate.rapid"))
    
setkey(episode1, anon_id, AdmissionDate.rapid, DischargeDate.rapid)
setkey(rapid, anon_id, AdmissionDate.rapid, DischargeDate.rapid)

rapid <- episode1[rapid]

#covid.specialtycodes <- "^A[16BQ]"
#noncovid.specialtycodes <- "^A[289DGHMR]|^[CEFHJ]"
#diagnoses[, noncovid.specialty := as.integer(grepl(noncovid.specialtycodes, specialty))]
#diagnoses[, covid.specialty := as.integer(grepl(covid.specialtycodes, specialty))]
## 1840 unique SMR01 admissions where positive test on or day after date of admission.

save(rapid, file=paste0(datadir, "rapid.RData"))
