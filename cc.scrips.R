library(data.table)

linkdate <- "jul28"
datadir <- "./data/2021-07-28/"
 
scripsobject.filename <- paste0(datadir, "scrips.pruned.RData")
scrips.firsttime <- FALSE
if(scrips.firsttime) {
    source("readscrips.R") # saves pruned scrips table with anon_id, dispensed date, item code
}
cat("Loading saved scripsbnf file ...")
load(paste0(datadir, "scripsbnf.RData"))
cat("done\n")

load(paste0(datadir, "cc.specimendate.RData"))
## create table cc.scripsbnf with 240-day lookback from specimen date for each anon_id / specimendate
cat("creating table cc.scripsbnf ...")
cc.specimendate[, lookback240day := specimen_date - 240]
setkey(cc.specimendate, anon_id, lookback240day, specimen_date)
scripsbnf[, enddate := dispensed_date + 1]
setkey(scripsbnf, anon_id, dispensed_date, enddate)
cc.scripsbnf <- foverlaps(scripsbnf,
                          cc.specimendate[, .(anon_id, specimen_date, lookback240day)],
                          type="within", nomatch=NULL)
rm(scripsbnf)
gc()

cc.scripsbnf <- cc.scripsbnf[as.integer(specimen_date - dispensed_date) >= 15,
                             .(anon_id, specimen_date, dispensed_date, bnf_item_code)]
cc.scripsbnf <- unique(cc.scripsbnf, by=c("anon_id", "specimen_date", "bnf_item_code"))
setkey(cc.scripsbnf, anon_id, specimen_date)
cat("done\n")
cat("cc.scripsbnf object uses", object.size(cc.scripsbnf) * 1E-6, "MB\n")
gc()
save(cc.scripsbnf, file=paste0(datadir, "cc.scripsbnf.RData"))
