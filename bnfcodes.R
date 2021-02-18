
bnfcodes <- as.data.table(read_excel("./BNF_Code_Information.xlsx", sheet=4))
colnames(bnfcodes) <- c("chaptername", "chapternum", "sectionname", "sectioncode")
bnfcodes[, `:=`(chapternum=as.integer(chapternum), sectioncode=as.integer(sectioncode))]

bnfparacodes <- as.data.table(read_excel("./BNF_Code_Information.xlsx", sheet=3))
colnames(bnfparacodes) <- c("chaptername", "chapternum",
                            "sectionname", "sectioncode",
                            "paraname", "paracode")
bnfparacodes[, `:=`(chapternum=as.integer(chapternum), sectioncode=as.integer(sectioncode),
                    paracode=as.integer(paracode))]

bnfsubparacodes <- as.data.table(read_excel("./BNF_Code_Information.xlsx", sheet=1)[, 1:8])
colnames(bnfsubparacodes) <- c("chaptername", "chapternum",
                               "sectionname", "sectioncode",
                               "paraname", "paracode",
                               "subparaname", "subparacode")
bnfsubparacodes <- bnfsubparacodes[!duplicated(subparacode)]
bnfsubparacodes[, `:=`(chapternum=as.integer(chapternum), sectioncode=as.integer(sectioncode),
                       paracode=as.integer(paracode), subparacode=as.integer(subparacode))]
bnfsubparacodes[substring(subparaname, 1, 5) == "DUMMY", subparaname := sectionname]

bnfchapters <- bnfcodes[, 1:2]
bnfchapters <- base::unique(bnfchapters)
truncate.at <- StrPos(bnfchapters$chaptername, " |,|-") -1
truncate.at[is.na(truncate.at)] <- nchar(bnfchapters$chaptername[is.na(truncate.at)])
bnfchapters$shortname <- substr(bnfchapters$chaptername, 1, truncate.at) 

## recode BNF chapternums > 14 as 14 and label this category "Other"
bnfchapters$chapternum[bnfchapters$chapternum > 14] <- 14
bnfchapters$shortname[bnfchapters$chapternum==14] <- "Other" 
bnfchapters$shortname[bnfchapters$chapternum==4] <- "Nervous" 

## chemical codes
bnfchemicalcodes <- read_excel("./BNF_Code_Information.xlsx", sheet=1)[, 1:10]
colnames(bnfchemicalcodes) <- c(colnames(bnfsubparacodes),
                                "chemicalname", "chemicalcode")
# bnfchemicalcodes <- subset(bnfchemicalcodes, subset=!duplicated(chemicalcode))

