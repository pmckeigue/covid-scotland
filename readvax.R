library(data.table)

source("helperfunctions.R")

linkdate <- "sep22"

datadir <- "./data/2021-09-22/"
vacc.filename <- paste0(datadir, "CC_VACC_2021-09-22.rds")

############# vaccination database ############################################

vacc <- RDStodt(vacc.filename, key="anon_id") # ## long format, one record per dose
vacc <- vacc[vacc_status == "completed"]
vacc <- unique(vacc)
setnames(vacc, "vacc_occurence_time", "vaxdate")
vacc[, vaxdate := as.Date(vaxdate)]
vacc[vaxdate < as.Date("2020-12-01"), vaxdate := as.Date(vacc_event_created_at)]
vacc <- vacc[grep("^Covid-19", vacc_product_name)]
vacc[, vacc_product_name := factor(as.character(vacc_product_name))]
with(vacc, table(toupper(substr(vacc_batch_number, 1, 1)),
                 vacc_product_name)[, -1], exclude=NULL)

## Pfizer vaccines have batchnum beginning with E or F
## AZ vaccines have batchnum beginning with 4, A, K or P
## unknown batchnums coded as UNKNOWN or unknown

## clean up miscoded product using batch number
vacc[toupper(substr(vacc_batch_number, 1, 2)) == "UN",
         vacc_batch_number := NA]
vacc[toupper(substr(vacc_batch_number, 1, 1)) %in% c("E", "F"),
         productcorrect := "Covid-19 mRNA Vaccine Pfizer"]
vacc[toupper(substr(vacc_batch_number, 1, 1)) %in% c("4", "A", "K", "P"),
         productcorrect := "Covid-19 Vaccine AstraZeneca"]
vacc[!is.na(productcorrect), vacc_product_name := productcorrect]

table(vacc$vacc_product_name, exclude=NULL)

vacc <- vacc[, .(anon_id, vacc_dose_number, vaxdate, vacc_product_name,
                     vacc_batch_number)]
setkey(vacc, anon_id, vacc_dose_number) 
vacc <- unique(vacc, by=key(vacc))
vacc.wide <- dcast(vacc, anon_id ~ vacc_dose_number,
                   value.var=c("vaxdate", "vacc_product_name",
                               "vacc_batch_number"))
vacc.wide[is.na(vaxdate_1), vaxdate_1 := vaxdate_2]
vacc.wide[vaxdate_2 <= vaxdate_1, vaxdate_2 := NA]
vacc.wide[vaxdate_3 <= vaxdate_2, vaxdate_3 := NA]
vacc.wide <- unique(vacc.wide)
vacc.wide[, weekdose1 := floor(as.integer(vaxdate_1 - as.Date("2020-12-01")) / 7)]
vacc.wide[, weekdose1 := as.Date("2010-12-01") + 7 * weekdose1]
setkey(vacc.wide, anon_id)

save(vacc.wide, file=paste0(datadir, "vacc.wide.RData"))
