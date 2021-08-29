
library(data.table)

linkdate <- "jul28"
datadir <- "./data/2021-07-28/"
 
cat("Loading saved scrips.pruned table ...")
load(paste0(datadir, "scrips.pruned.RData"))

load(paste0(datadir, "cc.specimendate.RData"))
cc.specimendate[, lookback240day := specimen_date - 240 - 14]
setkey(cc.specimendate, anon_id, lookback240day, specimen_date)

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

summary(scrips)

weeklydose <- function(chemsubstcode, drugname) {
    drug <- scrips[grep(paste0("^", chemsubstcode), substr(bnf_item_code, 1, 9))]
    drug[, enddate := dispensed_date + 1]
    setkey(drug, anon_id, dispensed_date, enddate)
    cc.drug <- foverlaps(drug,
                          cc.specimendate[, .(anon_id, specimen_date, lookback240day)],
                         type="within", nomatch=NULL)
    drug.dose.weekly <- cc.drug[specimen_date - dispensed_date <= 120 &
                                 specimen_date - dispensed_date >= 1, 
                                 list(weeklydose=sum(item_strength * quantity) * 7 / 120),
                              by=c("anon_id", "specimen_date")]
    setnames(drug.dose.weekly, "weeklydose", paste0(drugname, ".weeklydose"))
    setkey(drug.dose.weekly, anon_id, specimen_date)
    return(drug.dose.weekly)
}

methotrexate <- weeklydose(chemsubstcode="1001030U0", drugname="methotrexate")
prednisolone <- weeklydose(chemsubstcode="0603020T0", drugname="prednisolone")

drugdoses <- merge(methotrexate, prednisolone, all=TRUE)

save(drugdoses, file=paste0(datadir, "drugdoses.RData"))

gc()

