library(data.table)
source("helperfunctions.R")

RDStodt <- function(rds.filename, keyname=NULL) {
    ## read RDS file 
    dt <- setDT(readRDS(rds.filename)) # modifies data frame by reference
    cat(rds.filename, "loaded\n")
    
    ## recode <32 value character cols as factor 
    factor.cols <- names(which(unlist(sapply(dt,
                                             function(x) is.character(x) &
                                                         length(unique(x)) < 32))))
    if(length(factor.cols) > 0) {
        for (col in factor.cols) {
            set(dt, j=col, value=as.factor(dt[[col]]))
        }
        #dt[, (factor.cols) := lapply(.SD, as.factor), .SDcols = factor.cols]
    }
    gc()
    cat("recoding cols as factor completed\n")
    
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
            #dt[, (integer.cols) := lapply(.SD, as.integer), .SDcols = integer.cols]
        }
    }
    
    if(!is.null(keyname)) {
        setkeyv(dt, keyname)
    }
    return(dt)
}

linkdate <- "jul12"
linkdate <- "jul28"
## all saved data files should be saved to datadir defined by linkdate

datadir <- ""

if(linkdate=="jul12") {
    datadir <- "./data/2021-07-12/"
    scrips.filename <- paste0(datadir, "CC_PIS_x15_ANON_2021-07-12.rds")
    scrips.last15.filename <- paste0(datadir, "CC_PIS_15_ANON_2021-07-12.rds")
} else if(linkdate=="jul28") {
      datadir <- "./data/2021-07-28/"
      scrips.filename <- paste0(datadir, "CC_PIS_ANON_2021-07-28.rds")
}

scrips <- RDStodt(scrips.filename) # 8.95 GB in memory

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

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

objmem <- 1E-6 * sort( sapply(ls(), function(x) {object.size(get(x))}))
print(tail(objmem))

## scrips now reduced to 6.4 GB
scripsobject.filename <- paste0(datadir, "scrips.pruned.RData")
save(scrips, file=scripsobject.filename)

scripsbnf <- scrips[, .(anon_id, dispensed_date, bnf_item_code)]
save(scripsbnf,
     file=paste0(datadir, "scripsbnf.RData"))
rm(scrips)
gc()

