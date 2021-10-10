library(data.table)
library(readxl)
source("helperfunctions.R")

## datadir should already be set when this script is sourced from casecontrol.R
datadir <- "data/2021-09-22/"
smr06.filename <- paste0(datadir, "CC_SMR06_ICD10_ANON_2021-09-22.rds")

# source("cc.scrips.R")
load(paste0(datadir, "cc.specimendate.RData"))
load(paste0(datadir, "cc.scripsbnf.RData"))
load(paste0(datadir, "cc.diagnoses.RData"))
load(paste0(datadir, "cc.procedures.RData"))

## cc.specimendate, cc.diagnoses and cc.procedures should be loaded in memory

## BNF chapter codes: up to 2 digits, leading 0 may be stripped
## section codes: 2 digits
## paragraph codes: 2 digits
## subpara codes: 1 digit
## chemical substance code: 2 characters, may include a letter

bnfcodes <- as.data.table(read_excel("./BNF_Code_Information.xlsx", sheet=4))
colnames(bnfcodes) <- c("chaptername", "chapternum", "sectionname", "sectioncode")
bnfcodes[, `:=`(chapternum=as.integer(chapternum), sectioncode=as.integer(sectioncode))]
bnfchapters <- bnfcodes[, 1:2]
bnfchapters <- base::unique(bnfchapters)
 
## number of distinct BNF subparas with at least one prescription, excluding chapter 02 cardiovascular
#setkey(cc.all, anon_id, specimen_date)
#cc.all <- cc.numdrugs.notcv[cc.all]
#setnafill(cc.all, cols="numdrugs.notcv", fill=0) 

#cc.numdrugs.cns <- cc.scripsbnf[substr(bnf_item_code, 1, 2) == "04",
#                                  .(numdrugs.cns = uniqueN(substr(bnf_item_code, 1, 7))),
#                                  by=c("anon_id", "specimen_date")]
#setorder(cc.numdrugs.cns, -numdrugs.cns)
#setkey(cc.numdrugs.cns, anon_id, specimen_date)
#cc.numdrugs.cns <- unique(cc.numdrugs.cns, by=key(cc.numdrugs.cns))
#setkey(cc.all, anon_id, specimen_date)
#cc.all <- cc.numdrugs.cns[cc.all]
#setnafill(cc.all, cols="numdrugs.cns", fill=0) 

## 103 is proton pump inhibitors 01 03
## 1001030 is rheumatic disease suppressants 10 01 03 0
## 1001030C0 is hydroxychloroquine
## 1001030U0 is methotrexate

## recode scrips$bnf.chapter values > 14 or NA to 14
#scrips[is.na(chapternum), chapternum := 14]
#scrips[chapternum > 14, chapternum := 14] 
## colchicine
## BNF Subparagraph Code  1001040
## BNF Chemical Substance 1001040G0
## always prescribed as 500 mcg tablets, quantity usually 28, 56 or 100
## dose 2-4 tablets/day so at 2 tabs/day that corresponds to 14, 28 or 50 days supply
## gabapentin 0408010G0 pregabalin 0408010AE - classified in subpara 0408010 Control of epilepsy

## number of distinct BNF subparas with at least one prescription in chapter 02 cardiovascular
cc.numdrugs.cv <- cc.scripsbnf[substr(bnf_item_code, 1, 2) == "02",
                                  .(numdrugs.cv = uniqueN(substr(bnf_item_code, 1, 7))),
                                  by=c("anon_id", "specimen_date")]
setorder(cc.numdrugs.cv, -numdrugs.cv)
setkey(cc.numdrugs.cv, anon_id, specimen_date)
cc.numdrugs.cv <- unique(cc.numdrugs.cv, by=key(cc.numdrugs.cv))

cc.numdrugs.notcv <- cc.scripsbnf[substr(bnf_item_code, 1, 2) != "02",
                                  .(numdrugs.notcv = uniqueN(substr(bnf_item_code, 1, 7))),
                                  by=c("anon_id", "specimen_date")]
setorder(cc.numdrugs.notcv, -numdrugs.notcv)
setkey(cc.numdrugs.notcv, anon_id, specimen_date)
cc.numdrugs.notcv <- unique(cc.numdrugs.notcv, by=key(cc.numdrugs.notcv))
cc.numdrugs <- merge(cc.numdrugs.cv, cc.numdrugs.notcv, all=TRUE)
rm(cc.numdrugs.cv)
rm(cc.numdrugs.notcv)

## extract comorbid conditions from diagnoses, bnf item codes and procedures 
## FIXME: dayslookback should apply to diagnotes and bnf codes argument also
dt.diagorscrip <- function(diag.regex=NULL, bnf.regex=NULL, procedure.regex=NULL, dayslookback=NULL,
                           diag.name) {
    ## cc.specimendate keyed on anon_id, specimendate
    ## cc.diagnoses keyed on anon_id, specimendate
    ## extract records from cc.diagnoses that are unique on key, where icd10 matches regex
    ## create new column for this diagnosis, and retain key columns plus new column
    ## left join cc.specimendate with these records, creating a new column 
    diag.matched <- NULL
    scrips.matched <- NULL
    procedures.matched <- NULL
    if(!is.null(diag.regex)) {
        diag.matched <- unique(cc.diagnoses[grep(diag.regex, icd10), .(anon_id, specimen_date)])
    }
    if(!is.null(bnf.regex)) {
        scrips.matched <- unique(cc.scripsbnf[grep(bnf.regex, substr(bnf_item_code, 1, 6)),
                                              .(anon_id, specimen_date)])
    }
    if(!is.null(procedure.regex)) {
        if(!is.null(dayslookback)) {
            cc.procedures.lookback <- cc.procedures[specimen_date - discharge_date <= dayslookback &
                                                    specimen_date - discharge_date >= 0]
        } else cc.procedures.lookback <- cc.procedures
        procedures.matched <- unique(cc.procedures.lookback[grep(procedure.regex, main_operation),
                                                            .(anon_id, specimen_date)])
        
    }
    ## rbind but retain only records unique on anon_id, specimen_date
    diagorscrip <- unique(rbind(diag.matched, scrips.matched, procedures.matched))
    setkey(diagorscrip, anon_id, specimen_date)
    diagorscrip[, eval(diag.name) := 1]
    cc.diagorscrip <- diagorscrip[cc.specimendate[, .(anon_id, specimen_date)]]
    setnafill(cc.diagorscrip, fill=0, cols=3)
    diagorscrip[, eval(diag.name) := as.integer(eval(diag.name))]
    setkey(cc.diagorscrip, anon_id, specimen_date)
    cc.diagorscrip <- unique(cc.diagorscrip, by=key(cc.diagorscrip))
    return(cc.diagorscrip)
}

smr06 <- RDStodt(smr06.filename)
setnames(smr06, "icd10s_cancer_site", "icd10")
smr06[, incidence_date := as.Date(incidence_date)]
smr06[, joindate := incidence_date]
setkey(cc.specimendate, anon_id, lookback5yr, specimen_date)
setkey(smr06, anon_id, incidence_date, joindate)
cc.smr06 <- foverlaps(smr06,
                          cc.specimendate[, .(anon_id, specimen_date, lookback5yr)],
                          type="within", nomatch=NULL)
cc.smr06 <- cc.smr06[, .(anon_id, specimen_date, incidence_date, icd10)]
rm(smr06)
setkey(cc.smr06, anon_id, specimen_date)

##  diabetes ################################################

cc.diabetes <- dt.diagorscrip(diag.regex="^E1[0-4]", bnf.regex="^0601", diag.name="diabetes")

################  IHD ###################################
cc.ihd <- dt.diagorscrip(diag.regex="^I2[0-5]", bnf.regex="^020601",
                         procedure.regex="^K4[012349]|^K50",  diag.name="ihd")

## nitrates are BNF code 020601
## procedure codes for CABG and PTCA

##### other heart disease ####################################
## heart disease is I05 to I52
## 02023 anti-arrythmics
cc.heart.other <- dt.diagorscrip(diag.regex="^I0[01256789]|^I1[0-5]|^I2[6-8]|^I3[0-9]|^I4[0-9]|^I5[0-2]",
                                 bnf.regex="^0203",  procedure.regex="^K57", 
                                 diag.name="heart.other")

#######################################################################

## other circulatory disease is I60 to I99
cc.circulatory.other <- dt.diagorscrip(diag.regex="^I[6-9]|^Z95", diag.name="circulatory.other")

###############################################

## codes for solid organ transplants https://datadictionary.nhs.uk/Covid19PRA
## not M011 which is autotransplantation of kidney

solidtransplant.regex <- "^M01[23459]?|^B17[18]?|^E53[12389]?|^G26[19]?|^G688?|^G788|^J01[123589]?|^J54[12489]?|^K01[12]?|^K02[123489]?|^Y01[45689]"

cc.solidtransplants <- cc.procedures[grep(solidtransplant.regex, main_operation)]
setorder(cc.solidtransplants, -discharge_date)
cc.lastsolidtransplant <- unique(cc.solidtransplants, by=c("anon_id", "specimen_date"))
rm(cc.solidtransplants)
cc.lastsolidtransplant <- cc.lastsolidtransplant[, .(anon_id, specimen_date, discharge_date, main_operation)]
setnames(cc.lastsolidtransplant, "discharge_date", "date.lasttransplant")
setnames(cc.lastsolidtransplant, "main_operation", "code.lasttransplant")
setkey(cc.lastsolidtransplant, anon_id, specimen_date)

cc.solidtransplant <- dt.diagorscrip(procedure.regex=solidtransplant.regex,
                                     diag.name="solidtransplant")

############# chronic kidney disease ##########################
## includes CKD stage 4
cc.ckd <- dt.diagorscrip(diag.regex="^N18[45]|^Z49[0-2]|^Z94[02]",
                         procedure.regex="^X40[12]", diag.name="ckd")

##### asthma and chronic lower respiratory disease #################

#cc.asthma <- unique(diagnoses$ANON_ID[grep("^J4[56]", diagnoses$ICD10)])
cc.chronresp <- dt.diagorscrip(diag.regex="^J4[01234567]|^J6[0-9]|^J70|^J8[0-6]|^J9[0-9]|^G47\\.?3",
                               bnf.regex="^030(1|3)",
                               diag.name="chronresp")

##### Neurological disorders #######################

## include all Nervous chapter except G40 "Episodic and Paroxysmal Disorders"
## omit G0 meningitis and encephalitis, and G5 local neuropathies

## also include F03 dementia NOS
cc.neuro <- dt.diagorscrip(diag.regex="^F03|^G[1236789]",
                           bnf.regex="^(409)|(411)",
                           diag.name="neuro")

## drugs listed by HPS pharmacist as used for multiple sclerosis
## interferon beta 080204M, Glatiramer acetate 0802040U0, Natalizumab 0802040W0
## Dimethyl fumar 0802040AK, Teriflunomide 0802040AL, Alemtuzumab 0802030

############## Liver disease #############################################

liver.grep.string <- "^C22\\.?0|^I85\\.?0|^I98\\.?3|^K70\\.?[234|^K71\\.?7|^K72\\.?[019]|^K72\\.?[019|^K73|^K74\\.?[023456]|^K76\\.?7|^R18"

cc.liver <- dt.diagorscrip(diag.regex=liver.grep.string,
                           diag.name="liver")
 
#### Immunodeficiency and immunosuppression #################################

cc.immune <- dt.diagorscrip(diag.regex="^B2[0-3|^D8[0-9]",
                            bnf.regex="^(0802)|(050301)",
                            diag.name="immune")

############# autoimmune rheumatology
	
## RA codes
##     M05 Seropositive rheumatoid arthritis  
##     M06 Other rheumatoid arthritis
##     M08 Juvenile arthritis  
##     M09 Juvenile arthritis in diseases classified elsewhere  
##     M12.3 Palindromic rheumatism
##     M13 Other arthritis

## AXSPA codes
##      M07 Psoriatic and enteropathic arthropathies
##	M45 Ankylosing spondylitis  
##	M46 Other inflammatory spondylopathies  

## CTD codes
##	M30 Polyarteritis nodosa and related conditions  
##      M31 Other necrotizing vasculopathies
##	M32 Systemic lupus erythematosus  
##	M33 Dermatopolymyositis  
##	M34 Systemic sclerosis   
##      M35 Other systemic involvement of connective tissue

cc.ra <- dt.diagorscrip(diag.regex="^M0[5689]|^M123|^M13",
                            diag.name="ra")
cc.axspa <- dt.diagorscrip(diag.regex="^M07|^M4[56]",
                           diag.name="axspa")
cc.ctd <- dt.diagorscrip(diag.regex="^M3[0-5]",
                           diag.name="ctd")

cc.rheumatol <- cc.axspa[cc.ra]
cc.rheumatol <- cc.ctd[cc.rheumatol]
rm(cc.ra)
rm(cc.axspa)
rm(cc.ctd)

cc.rheumatol[, rheumatol.diag := "No rheumatologic diagnosis"]
cc.rheumatol[ctd==1, rheumatol.diag := "CTD"]
cc.rheumatol[axspa==1, rheumatol.diag := "AXSPA"]
cc.rheumatol[ra==1, rheumatol.diag := "RA"]
cc.rheumatol <- cc.rheumatol[, .(anon_id, specimen_date, rheumatol.diag)]

##############################################

cc.methotrexate <- unique(cc.scripsbnf[grep("^1001030U0", substr(bnf_item_code, 1, 9)),
                                       .(anon_id, specimen_date)])
cc.methotrexate[, methotrexate := 1]
cc.methotrexate <- cc.methotrexate[cc.specimendate[, .(anon_id, specimen_date)]]
setnafill(cc.methotrexate, fill=0, cols=3)
setkey(cc.methotrexate, anon_id, specimen_date)

cc.hydroxychloroquine <- unique(cc.scripsbnf[grep("^1001030C0", substr(bnf_item_code, 1, 9)),
                                       .(anon_id, specimen_date)])
cc.hydroxychloroquine[, hydroxychloroquine := 1]
cc.hydroxychloroquine <- cc.hydroxychloroquine[cc.specimendate[, .(anon_id, specimen_date)]]
setnafill(cc.hydroxychloroquine, fill=0, cols=3)
setkey(cc.hydroxychloroquine, anon_id, specimen_date)

cc.leflunomide <- unique(cc.scripsbnf[grep("^1001030L0", substr(bnf_item_code, 1, 9)),
                                       .(anon_id, specimen_date)])
cc.leflunomide[, leflunomide := 1]
cc.leflunomide <- cc.leflunomide[cc.specimendate[, .(anon_id, specimen_date)]]
setnafill(cc.leflunomide, fill=0, cols=3)
setkey(cc.leflunomide, anon_id, specimen_date)

## prednisolone ch 6 0603020T0
cc.prednisolone <- unique(cc.scripsbnf[grep("^0603020T0", substr(bnf_item_code, 1, 9)),
                                       .(anon_id, specimen_date)])
cc.prednisolone[, prednisolone := 1]
cc.prednisolone <- cc.prednisolone[cc.specimendate[, .(anon_id, specimen_date)]]
setnafill(cc.prednisolone, fill=0, cols=3)
setkey(cc.prednisolone, anon_id, specimen_date)

## sulfasalazine 1001040E0, mesalazine 0105010B0, Olsalazine 105010C0
cc.sulfasalazine <- unique(cc.scripsbnf[grep("^0105010E0", substr(bnf_item_code, 1, 9)),
                                       .(anon_id, specimen_date)])
cc.sulfasalazine[, sulfasalazine := 1]
cc.sulfasalazine <- cc.sulfasalazine[cc.specimendate[, .(anon_id, specimen_date)]]
setnafill(cc.sulfasalazine, fill=0, cols=3)
setkey(cc.sulfasalazine, anon_id, specimen_date)

## any 5-ASA
cc.aminosalicylate <- unique(cc.scripsbnf[grep("^0105010B0|^0105010C0|^0105010E0", substr(bnf_item_code, 1, 9)),
                                       .(anon_id, specimen_date)])
cc.aminosalicylate[, aminosalicylate := 1]
cc.aminosalicylate <- cc.aminosalicylate[cc.specimendate[, .(anon_id, specimen_date)]]
setnafill(cc.aminosalicylate, fill=0, cols=3)
setkey(cc.aminosalicylate, anon_id, specimen_date)

############# IBD

cc.ibd <- dt.diagorscrip(diag.regex="^K5[01]",
                                  diag.name="ibd")

############# neoplasms ################

cc.neoplasm <- dt.diagorscrip(diag.regex="^C[0-9]|^D[0-4]",
                              bnf.regex="^0801",
                              diag.name="neoplasm")

## blood cancers identified from hospital diagnoses and cancer registrations
cc.bloodcancer <- dt.diagorscrip(diag.regex="^C8[1-8]|^C9[0-6]",
                                 diag.name="bloodcancer")
cc.bloodcancer.smr06 <- cc.smr06[grep("^C8[1-8]|^C9[0-6]", icd10)]
cc.bloodcancer <- cc.bloodcancer.smr06[cc.bloodcancer]
with(cc.bloodcancer, table(bloodcancer, !is.na(icd10)))
cc.bloodcancer[!is.na(icd10), bloodcancer := 1]
cc.bloodcancer <- cc.bloodcancer[, .(anon_id, specimen_date, bloodcancer)]
setorder(cc.bloodcancer, -bloodcancer)
cc.bloodcancer <- unique(cc.bloodcancer, by=c("anon_id", "specimen_date"))
setkey(cc.bloodcancer, anon_id, specimen_date)

###### disorders of esophagus, stomach and duodenum ############################

cc.esoph.stomach.duod <-  dt.diagorscrip(diag.regex="^K2[0-9]|^K3[01]",
                                         diag.name="esoph.stomach.duod")

###################################################################
cat(nrow(cc.ihd), nrow(cc.heart.other), "\n")
cc.comorbid <- cc.heart.other[cc.ihd]

cat(nrow(cc.comorbid), nrow(cc.circulatory.other), "\n") 
cc.comorbid <- cc.circulatory.other[cc.comorbid]

cat(nrow(cc.comorbid), nrow(cc.ckd), "\n") 
cc.comorbid <- cc.ckd[cc.comorbid]

cat(nrow(cc.comorbid), nrow(cc.chronresp), "\n") 
cc.comorbid <- cc.chronresp[cc.comorbid]
cat(nrow(cc.comorbid), "\n")

cc.comorbid <- cc.neuro[cc.comorbid]
cc.comorbid <- cc.liver[cc.comorbid]
cc.comorbid <- cc.immune[cc.comorbid]

cat(nrow(cc.comorbid), nrow(cc.neoplasm), "\n") 
cc.comorbid <- cc.neoplasm[cc.comorbid]

cc.comorbid <- cc.bloodcancer[cc.comorbid]

cat(nrow(cc.comorbid), nrow(cc.esoph.stomach.duod), "\n") 
cc.comorbid <- cc.esoph.stomach.duod[cc.comorbid]

cat(nrow(cc.comorbid), nrow(cc.rheumatol), "\n") 
cc.comorbid <- cc.rheumatol[cc.comorbid]

cat(nrow(cc.comorbid), nrow(cc.methotrexate), "\n") ## duplicate key rows in cc.methotrexate
cc.methotrexate <- unique(cc.methotrexate, by=key(cc.methotrexate))
cat(nrow(cc.comorbid), nrow(cc.methotrexate), "\n") 
cc.comorbid <- cc.methotrexate[cc.comorbid]

cat(nrow(cc.comorbid), nrow(cc.hydroxychloroquine), "\n") 
cc.comorbid <- cc.hydroxychloroquine[cc.comorbid]

cc.comorbid <- cc.sulfasalazine[cc.comorbid]
cc.comorbid <- cc.aminosalicylate[cc.comorbid]
cc.comorbid <- cc.leflunomide[cc.comorbid]
cc.comorbid <- cc.prednisolone[cc.comorbid]
cc.comorbid <- cc.ibd[cc.comorbid]
cc.comorbid <- cc.lastsolidtransplant[cc.comorbid]
cc.comorbid <- cc.solidtransplant[cc.comorbid]
                            
rm(cc.ihd)
rm(cc.heart.other)
rm(cc.circulatory.other)
rm(cc.ckd)
rm(cc.chronresp)
rm(cc.neuro)
rm(cc.liver)
rm(cc.immune)
rm(cc.neoplasm)
rm(cc.esoph.stomach.duod)
rm(cc.bloodcancer)
rm(cc.rheumatol)
rm(cc.ibd)
rm(cc.smr06)
rm(cc.methotrexate)
rm(cc.hydroxychloroquine)
rm(cc.sulfasalazine)
rm(cc.solidtransplant)
rm(cc.leflunomide)
gc()

save(cc.comorbid, file=paste0(datadir, "cc.comorbid.RData"))
save(cc.numdrugs, file=paste0(datadir, "cc.numdrugs.RData"))
