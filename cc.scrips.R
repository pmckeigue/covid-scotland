library(data.table)

linkdate <- "jul28"
datadir <- "./data/2021-07-28/"

linkdate <- "sep02"
datadir <- "./data/2021-09-02/"

linkdate <- "sep22"
datadir <- "./data/2021-09-22/"

#scripsobject.filename <- paste0(datadir, "scrips.pruned.RData")

scrips.firsttime <- FALSE
if(scrips.firsttime) {
    source("readscrips.R") # saves scripsbnf.RData with tables scripsonly, itemcodes, chemsubsts
}
cat("Loading saved scripsbnf file ...")
load(paste0(datadir, "scripsbnf.RData"))
cat("done\n")

load(paste0(datadir, "cc.specimendate.RData"))
## create table cc.scripsbnf with (240 + 14)-day lookback from specimen date for each anon_id / specimendate
cat("creating table cc.scripsbnf ...")
cc.specimendate[, lookback255day := specimen_date - 240 - 14]
setkey(cc.specimendate, anon_id, lookback255day, specimen_date)
scripsonly[, enddate := dispensed_date + 1]
setkey(scripsonly, anon_id, dispensed_date, enddate)
cc.scripsbnf <- foverlaps(scripsonly,
                          cc.specimendate[, .(anon_id, specimen_date, lookback255day)],
                          type="within", nomatch=NULL)
rm(scripsonly)
gc()

## drop scrips dispensed in last 14 days before specimendate
cc.scripsbnf <- cc.scripsbnf[as.integer(specimen_date - dispensed_date) >= 15,
                             .(anon_id, specimen_date, dispensed_date, itemcodekey)]
cc.scripsbnf <- unique(cc.scripsbnf, by=c("anon_id", "specimen_date", "itemcodekey"))

## import bnf_item_code
setkey(cc.scripsbnf, itemcodekey)
setkey(itemcodes, itemcodekey)
cc.scripsbnf <- itemcodes[, .(itemcodekey, bnf_item_code)][cc.scripsbnf]

setkey(cc.scripsbnf, anon_id, specimen_date)
cat("done\n")
cat("cc.scripsbnf object uses", object.size(cc.scripsbnf) * 1E-6, "MB\n")
gc()
save(cc.scripsbnf, file=paste0(datadir, "cc.scripsbnf.RData"))

## FIXME 62k missing values of itemcodekey in cc.scripsbnf
