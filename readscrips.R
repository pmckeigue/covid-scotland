library(data.table)
source("helperfunctions.R")

datadir <- "./data/2021-07-12/"
scrips.filename <- paste0(datadir, "CC_PIS_x15_ANON_2021-07-12.rds")
scrips.last15.filename <- paste0(datadir, "CC_PIS_15_ANON_2021-07-12.rds")
  
scripsobject.filename <- paste0(datadir, "scrips.pruned.RData")
scrips <- RDStodt(scrips.filename)
bnf.lookup <- unique(scrips[, c("bnf_item_code",
                                "bnf_item_description",
                                "bnf_paragraph_code",
                                "bnf_paragraph_description")])
table(duplicated(bnf.lookup$bnf_item_code))
save(bnf.lookup, file="bnf.lookup.RData")

scrips[, paid_date := NULL]        
scrips[, bnf_item_description := NULL]        
scrips[, bnf_paragraph_code := NULL]        
scrips[, bnf_paragraph_description := NULL]
save(scrips, file=scripsobject.filename)
scripsbnf <- scrips[, .(anon_id, dispensed_date, bnf_item_code)]
save(scripsbnf,
     file=paste0(datadir, "scripsbnf.RData"))

rm(scrips)
gc()

