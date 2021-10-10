library(data.table)
source("helperfunctions.R")


linkdate <- "jul28"
datadir <- "./data/2021-07-28/"

linkdate <- "sep22"
datadir <- "./data/2021-09-22/"

if(!exists("scripsonly")) {
    load(paste0(datadir, "scripsbnf.RData"))
}

load(paste0(datadir, "cc.specimendate.RData"))
cc.specimendate[, lookback240day := specimen_date - 240 - 14]
setkey(cc.specimendate, anon_id, lookback240day, specimen_date)

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

weeklydose <- function(chemsubstcode.arg, drugname, lookback.days=120) {
    keys <- itemcodes[chemsubstcode==chemsubstcode.arg,
                      itemcodekey]
    scrips.matched <- scripsonly[itemcodekey %in% keys]
    setkey(scrips.matched, itemcodekey)
    setkey(itemcodes, itemcodekey)
    drug <- scrips.matched[itemcodes]

    drug <- drug[!is.na(anon_id)]
    drug[, enddate := dispensed_date + 1]
    setkey(drug, anon_id, dispensed_date, enddate)
    cc.drug <- foverlaps(drug,
                          cc.specimendate[, .(anon_id, specimen_date, lookback240day)],
                         type="within", nomatch=NULL)
    drug.dose.weekly <- cc.drug[specimen_date - dispensed_date <= lookback.days &
                                specimen_date - dispensed_date >= 1 &
                                formulation_code=="TABS", 
                                list(weeklymg=sum(item_strength * quantity *
                                                    ifelse(item_strength_uom=="MICROGRAMS" |
                                                           item_strength_uom=="MCG", 0.001, 1)) *
                                         7 / lookback.days),
                                by=c("anon_id", "specimen_date")]
    
    setnames(drug.dose.weekly, "weeklymg", paste0(drugname, ".weeklymg"))
    setkey(drug.dose.weekly, anon_id, specimen_date)
    return(drug.dose.weekly)
}

#bnf.lookup[, chemsubstcode := substr(bnf_item_code, 1, 9)]
#bnf.chemsubst <- unique(bnf.lookup, by="chemsubstcode")


#https://cks.nice.org.uk/topics/corticosteroids-oral/background-information/equivalent-anti-inflammatory-doses/

#Drug	Dose equivalent to 5 mg of prednisolone
#Betamethasone	750 micrograms
#Cortisone acetate	25 mg
#Deflazacort	6 mg
#Dexamethasone	750 micrograms
#Hydrocortisone	20 mg
#Methylprednisolone	4 mg
#Prednisone	5 mg
#Triamcinolone	4 mg

if(FALSE) {
unique(scrips[grep("betamethasone", approved_name, ignore.case=TRUE)]) # no tabs
unique(scrips[grep("cortisone acetate", approved_name, ignore.case=TRUE)]) # tabs micrograms
unique(scrips[grep("deflazacort", approved_name, ignore.case=TRUE)]) # tabs mg
unique(scrips[grep("dexamethasone", approved_name, ignore.case=TRUE)]) # tabs mg or micrograms, also drops
unique(scrips[grep("hydrocortisone", approved_name, ignore.case=TRUE)]) # no tabs
unique(scrips[grep("methylprednisolone", approved_name, ignore.case=TRUE)]) # tabs or inj, mg or mg/ml
unique(scrips[grep("prednisone", approved_name, ignore.case=TRUE)]) # tabs mg
unique(scrips[grep("prednisolone", approved_name, ignore.case=TRUE)]) # tabs mg
unique(scrips[grep("triamcinolone", approved_name, ignore.case=TRUE)]) # inj or spray
}

#Betamethasone	0603020B0 not in tablet form
#Cortisone acetate	0603020F0
#Deflazacort	0603020I0
#Dexamethasone	0603020G0 tabs or drops
#Hydrocortisone	0603020J0 not in tablet form
#Methylprednisolone	0603020S0 tabs or inj
#Prednisone	0603020X0
#triamcinolone not in tablet form
                                        
methotrexate <- weeklydose(chemsubstcode.arg="1001030U0", drugname="methotrexate")

prednisolone <- weeklydose(chemsubstcode.arg="0603020T0", drugname="prednisolone")
cortisone_acetate <- weeklydose(chemsubstcode.arg="0603020F0", drugname="cortisone_acetate")
deflazacort <- weeklydose(chemsubstcode.arg="0603020I0", drugname="deflazacort")
dexamethasone <- weeklydose(chemsubstcode.arg="0603020G0", drugname="dexamethasone")
methylprednisolone <- weeklydose(chemsubstcode.arg="0603020S0", drugname="methylprednisolone")
prednisone <- weeklydose(chemsubstcode.arg="0603020X0", drugname="prednisone")

glucocorticoids <- merge(cortisone_acetate, prednisolone, all=TRUE)
glucocorticoids <- merge(deflazacort, glucocorticoids, all=TRUE)
glucocorticoids <- merge(dexamethasone, glucocorticoids, all=TRUE)
glucocorticoids <- merge(methylprednisolone, glucocorticoids, all=TRUE)
glucocorticoids <- merge(prednisone, glucocorticoids, all=TRUE)
setnafill(glucocorticoids, cols=c("prednisolone.weeklymg", "cortisone_acetate.weeklymg",
                                  "deflazacort.weeklymg", "dexamethasone.weeklymg",
                                  "methylprednisolone.weeklymg", "prednisone.weeklymg"), fill=0)

glucocorticoids[, prednisolone.equiv := prednisolone.weeklymg +
                      cortisone_acetate.weeklymg * 5 / 25 +
                      deflazacort.weeklymg * 5 / 6 +
                      dexamethasone.weeklymg * 5 / 0.75 +
                      methylprednisolone.weeklymg  * 5 / 4 +
                      prednisone.weeklymg]
glucocorticoids <- glucocorticoids[, .(anon_id, specimen_date, prednisolone.equiv)] 
                      
drugdoses <- merge(methotrexate, glucocorticoids, all=TRUE)
setkey(cc.specimendate, anon_id, specimen_date)
drugdoses <- drugdoses[cc.specimendate[, .(anon_id, specimen_date)]]
setnafill(drugdoses, cols=c("methotrexate.weeklymg", "prednisolone.equiv"), fill=0)
save(drugdoses, file=paste0(datadir, "drugdoses.RData"))

gc()

