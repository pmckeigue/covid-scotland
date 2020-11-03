
################### short names for ICD chapters ########################
## icd10_chapters is a list of chapters with start and end ICD codes

## transform this to a data frame named icdchapters
icdchapters <- data.frame(names(icd::icd10_chapters),
                             t(matrix(as.character(unlist(icd10_chapters)), nrow=2)))
colnames(icdchapters) <- c("name", "start", "end")
## shorten name field
icdchapters$shortname <- gsub("Diseases of the ", "", icdchapters$name)
icdchapters$shortname <- gsub("Certain ", "", icdchapters$shortname)
icdchapters$shortname <- gsub("conditions originating in the ", "",
                              icdchapters$shortname)
icdchapters$shortname <- gsub("Factors influencing ", "", icdchapters$shortname)
truncate.at <- StrPos(icdchapters$shortname, " |,|-") - 1
truncate.at[is.na(truncate.at)] <- nchar(icdchapters$shortname[is.na(truncate.at)])
icdchapters$shortname <- substr(icdchapters$shortname, 1, truncate.at) 
icdchapters$shortname <- gsub("^health$", "Health_factors", icdchapters$shortname)
icdchapters$start <- as.character(icdchapters$start)
icdchapters$end <- as.character(icdchapters$end)
icdchapters$start[icdchapters$start == "F01"] <- "F00" # mental 
icdchapters$end[icdchapters$end == "T88"] <- "T98" # injuries 

## icd10_subchapters is a list of subchapters with start and end ICD codes
## transform this to a data frame
icdsubchapters <- data.frame(names(icd::icd10_sub_chapters),
                             t(matrix(as.character(unlist(icd10_sub_chapters)), nrow=2)))
colnames(icdsubchapters) <- c("name", "start", "end")
icdsubchapters$start <- as.character(icdsubchapters$start)
icdsubchapters$end <- as.character(icdsubchapters$end)

## unexplained missing entries in R package icd table
icdsubchapters$end[icdsubchapters$end == "B20"] <- "B24" # codes B21 to B24 for complications of HIV are not in original ICD-10
icdsubchapters$end[icdsubchapters$end == "B97"] <- "B98" # code B98 for other diseases caused by H.pylori or Vibrio vulnificus
icdsubchapters$end[icdsubchapters$end == "H05"] <- "H06" # lacrimal system
icdsubchapters$end[icdsubchapters$end == "H11"] <- "H13" # conjunctiva
icdsubchapters$end[icdsubchapters$end == "H44"] <- "H45" # globe
icdsubchapters$end[icdsubchapters$end == "H57"] <- "H58" # other eye
icdsubchapters$start[icdsubchapters$start == "T07"] <- "T00" # injuries 
icdsubchapters$end[icdsubchapters$end == "T07"] <- "T14" # injuries 
icdsubchapters$start[icdsubchapters$start == "T30"] <- "T29" # injuries 
icdsubchapters$end[icdsubchapters$end == "T34"] <- "T35" # injuries 
icdsubchapters$end[icdsubchapters$end == "X08"] <- "X09" # injuries 
icdsubchapters$end[icdsubchapters$end == "X50"] <- "X51" # injuries 
icdsubchapters$end[icdsubchapters$end == "X58"] <- "X59" # injuries 
icdsubchapters$end[icdsubchapters$end == "X08"] <- "X09" # injuries 
icdsubchapters$start[icdsubchapters$start == "X71"] <- "X60" # self-harm 
icdsubchapters$end[icdsubchapters$end == "X83"] <- "X84" # self-harm 
icdsubchapters$start[icdsubchapters$start == "X92"] <- "X85" # complications of medical and surgical care 
icdsubchapters$start[icdsubchapters$start == "Y62"] <- "Y40" # complications of medical and surgical care 
icdsubchapters$end[icdsubchapters$end == "Z53"] <- "Z54" # self-harm 
icdsubchapters$start[icdsubchapters$start == "Y21"] <- "Y10" #  
icdsubchapters$end[icdsubchapters$end == "Y33"] <- "Y34" #  
icdsubchapters <- rbind(icdsubchapters,
                        data.frame(name=c("Insulin-dependent diabetes mellitus", 
                                          "Non-insulin-dependent diabetes mellitus", 
                                          "Malnutrition-related diabetes mellitus", 
                                          "Other specified diabetes mellitus",  
                                          "Unspecified diabetes mellitus",
                                          "Dementia in Alzheimer's disease",
                                          "Sequelae of external causes",
                                          "Provisional assignment or emergency use",
                                          "Resistance to antimicrobial and antineoplastic drugs",
                                          "Other accidental threats to breathing",
                                          "Contact with venomous animals and plants",
                                          "Accidental poisoning by and exposure to noxious substances",
                                          "Sequelae of external causes"), 
                                   start=c("E10", "E11", "E12", "E13", "E14", "F00", "T90", "U00", "U82", "W75", "X20", "X40", "Y85"),
                                     end=c("E10", "E11", "E12", "E13", "E14", "F00", "T98", "U80", "U89", "W84", "X29", "X49", "Y89"))) 
icdsubchapters <- icdsubchapters[order(icdsubchapters$start), ] # ensures that subchapters will be correctly numbered after overlap join

## now do an interval join so that we have the chapter and subchapter
icdsubchapters <- as.data.table(icdsubchapters)
icdchapters <- as.data.table(icdchapters)

icdsubchapters[, startnum := icdToInt(start)]
icdsubchapters[, endnum := icdToInt(end)]
icdsubchapters[, subchnum := .I]
icdchapters[, startnum := icdToInt(start)]
icdchapters[, endnum := icdToInt(end)]
icdchapters[, chnum := .I]
setkey(icdsubchapters, startnum, endnum)
setkey(icdchapters, startnum, endnum)

icd.joined <- foverlaps(icdsubchapters, icdchapters, type="any", which=TRUE)

setkey(icdsubchapters, subchnum)
setkey(icd.joined, xid)
icd.subchapters <- icdsubchapters[icd.joined]
setnames(icd.subchapters, "yid", "chnum")
## fill in missing values of chnum by carrying forward
icd.subchapters[, chnum := zoo::na.locf(chnum)]

