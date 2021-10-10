library(data.table)

RDStodt <- function(rds.filename, keyname=NULL) {
    ## read RDS file 
    dt <- setDT(readRDS(rds.filename)) # modifies data frame by reference

    ## recode <32 value character cols as factor 
    factor.cols <- names(which(unlist(sapply(dt,
                                             function(x) is.character(x) &
                                                         length(unique(x)) < 32))))
    if(length(factor.cols) > 0) {
        for (col in factor.cols) {
            set(dt, j=col, value=as.factor(dt[[col]]))
        }
    }
    gc()
    
    ## recode integer-valued numeric cols as integer
    numeric.cols <- names(which(unlist(sapply(dt, is.numeric))))
    if(length(numeric.cols) > 0) {
        integer.cols <- names(which(unlist(sapply(dt[, ..numeric.cols],
                                                  function(x) is.numeric(x) &
                                                              isTRUE(all.equal(x, floor(x)))))))
        if(length(integer.cols) > 0) {
            for (col in integer.cols) {
                set(dt, j=col, value=as.integer(dt[[col]]))
            }
        }
    }
    
    if(!is.null(keyname)) {
        setkeyv(dt, keyname)
    }
    return(dt)
}

#linkdate <- "jul12"
linkdate <- "jul28"
linkdate <- "sep02"
## all saved data files should be saved to datadir defined by linkdate


if(linkdate=="jul12") {
    datadir <- "./data/2021-07-12/"
    scrips.filename <- paste0(datadir, "CC_PIS_x15_ANON_2021-07-12.rds")
    scrips.last15.filename <- paste0(datadir, "CC_PIS_15_ANON_2021-07-12.rds")
} else if(linkdate=="jul28") {
      datadir <- "./data/2021-07-28/"
      scrips.filename <- paste0(datadir, "CC_PIS_ANON_2021-07-28.rds")
} else if(linkdate=="sep02") {
      datadir <- "./data/2021-09-02/"
datadir <- ""
      scrips.filename <- paste0(datadir, "CC_PIS_ANON_2021-09-02.rds")
}

scrips <- RDStodt(scrips.filename) # 8.95 GB in memory

## relational format
## bnf_item_codes should specify unique formulation code and item strength
## but there are a few duplicated item codes with different formulations or strengths
## so we define a primary key based on fields itemcode, formulation, strength

## table itemcodes: bnf_item_code, bnf_item_description, formulation_code, item_strength, item_strength_uom, approved_name, chemsubstcode

itemcodes <- unique(scrips[, .(bnf_item_code,  bnf_item_description,
                               formulation_code, item_strength,
                               item_strength_uom,
                               approved_name, bnf_paragraph_code,
                               bnf_paragraph_description)])
itemcodes <- itemcodes[!is.na(bnf_item_code)]
itemcodes[, chemsubstcode := substr(bnf_item_code, 1, 9)]
length(unique(itemcodes$bnf_item_code))
itemcodes[, num.withcode := .N, by=bnf_item_code]
itemcodes[num.withcode > 1]
setkey(itemcodes, bnf_item_code, formulation_code, item_strength) # primary key
itemcodes <- unique(itemcodes, by=key(itemcodes))
itemcodes[, itemcodekey := .I]
setkey(itemcodes, bnf_item_code, formulation_code, item_strength) 

## table scrips: anon_id, dispensed_date, item_code_no, no_items, quantity
scripsonly <- scrips[, .(anon_id, dispensed_date, bnf_item_code, formulation_code,
                         item_strength, no_items, quantity)]
setkey(scripsonly, bnf_item_code, formulation_code, item_strength) 

itemcodes.join <- itemcodes[, .(bnf_item_code, formulation_code,
                                item_strength, itemcodekey)]
setkey(itemcodes.join, bnf_item_code, formulation_code, item_strength) 
scripsonly <- itemcodes.join[scripsonly]

scripsonly <- scripsonly[, .(anon_id, dispensed_date, itemcodekey, no_items, quantity)]
setkey(scripsonly, itemcodekey) 

## table chemsubst: chemsubstcode (primary key), approved_name
chemsubsts <- unique(itemcodes[, .(chemsubstcode, approved_name,
                                bnf_paragraph_code,
                                bnf_paragraph_description)])
setkey(chemsubsts, chemsubstcode)

# drop from itemcodes columns duplicated in chemsubsts
itemcodes <- itemcodes[, .(itemcodekey, bnf_item_code, bnf_item_description,
                               formulation_code, item_strength,
                               item_strength_uom,
                               chemsubstcode)]

save(scripsonly, itemcodes, chemsubsts, file="scripsbnf.RData")

#scrips[, paid_date := NULL]        
#scrips[, bnf_item_description := NULL]        
#scrips[, bnf_paragraph_code := NULL]        
#scrips[, bnf_paragraph_description := NULL]
rm(scrips)
objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))
gc()

