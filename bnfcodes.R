## BNF chapter codes: up to 2 digits, leading 0 may be stripped
## section codes: 2 digits

## 103 is proton pump inhibitors 01 03

## paragraph codes: 2 digits
## subpara codes: 1 digit

## 1001030 is rheumatic disease suppressants 10 01 03 0

## chemical substance code: 2 characters, may include a letter

## 1001030C0 is hydroxychloroquine
## 1001030U0 is methotrexate

bnfcodes <- read_excel("./BNF_Code_Information.xlsx", sheet=4)
colnames(bnfcodes) <- c("chaptername", "chapternum", "sectionname", "sectioncode")

bnfparacodes <- read_excel("./BNF_Code_Information.xlsx", sheet=3)
colnames(bnfparacodes) <- c("chaptername", "chapternum",
                            "sectionname", "sectioncode",
                            "paraname", "paracode")
bnfparacodes$paracode <- as.integer(bnfparacodes$paracode)

bnfsubparacodes <- read_excel("./BNF_Code_Information.xlsx", sheet=1)[, 1:8]
colnames(bnfsubparacodes) <- c("chaptername", "chapternum",
                               "sectionname", "sectioncode",
                               "paraname", "paracode",
                               "subparaname", "subparacode")
bnfsubparacodes <- bnfsubparacodes[!duplicated(bnfsubparacodes$subparacode), ]
bnfsubparacodes$subparacode <- as.integer(bnfsubparacodes$subparacode)

bnfchapters <- bnfcodes[, 1:2]
bnfchapters <- base::unique(bnfchapters)
truncate.at <- StrPos(bnfchapters$chaptername, " |,|-") -1
truncate.at[is.na(truncate.at)] <- nchar(bnfchapters$chaptername[is.na(truncate.at)])
bnfchapters$shortname <- substr(bnfchapters$chaptername, 1, truncate.at) 

## recode BNF chapternums > 14 as 14 and label this category "Other"
bnfchapters$chapternum[bnfchapters$chapternum > 14] <- 14
bnfchapters$shortname[bnfchapters$chapternum==14] <- "Other" 
bnfchapters$shortname[bnfchapters$chapternum==4] <- "Nervous" 
