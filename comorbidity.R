## extract comorbid conditions from diagnoses and bnf item codes 

dt.diagorscrip <- function(diag.regex=NULL, bnf.regex=NULL, diag.name) {
    ## cc.specimendate keyed on anon_id, specimendate
    ## cc.diagnoses keyed on anon_id, specimendate
    ## extract records from cc.diagnoses that are unique on key, where icd10 matches regex
    ## create new column for this diagnosis, and retain key columns plus new column
    ## left join cc.specimendate with these records, creating a new column 
    diag.matched <- NULL
    scrips.matched <- NULL
    if(!is.null(diag.regex)) {
       diag.matched <- unique(cc.diagnoses[grep(diag.regex, icd10), .(anon_id, specimen_date)])
    }
    if(!is.null(bnf.regex)) {
        scrips.matched <- unique(cc.scripsbnf[grep(bnf.regex, substr(bnf_item_code, 1, 6)),
                                              .(anon_id, specimen_date)])
    }
    diagorscrip <- unique(rbind(diag.matched, scrips.matched))
    setkey(diagorscrip, anon_id, specimen_date)
    diagorscrip[, eval(diag.name) := 1]
    cc.diagorscrip <- diagorscrip[cc.specimendate[, .(anon_id, specimen_date)]]
    setnafill(cc.diagorscrip, fill=0, cols=3)
    setkey(cc.diagorscrip, anon_id, specimen_date)
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
cc.ihd <- dt.diagorscrip(diag.regex="^I2[0-5]", bnf.regex="^020601", diag.name="ihd")

## nitrates are BNF code 020601
## FIXME: procedure codes for CABG and PTCA
                                        #ids.procedures.IHD <- unique(procedures$ANON_ID[grep("^K4[012349]|^K50", procedures$MAIN_OPERATION)])

##### other heart disease ####################################
## heart disease is I05 to I52
## 02023 anti-arrythmics
cc.heart.other <- dt.diagorscrip(diag.regex="^I0[01256789]|^I1[0-5]|^I2[6-8]|^I3[0-9]|^I4[0-9]|^I5[0-2]",
                                 bnf.regex="^0203", 
                                 diag.name="heart.other")

## ids.procedures.heart.other <- unique(procedures$ANON_ID[grep("^K57", procedures$MAIN_OPERATION)])

#######################################################################

## other circulatory disease is I60 to I99
cc.circulatory.other <- dt.diagorscrip(diag.regex="^I[6-9]|^Z95", diag.name="circulatory.other")
    
############# chronic kidney disease ##########################
## includes CKD stage 4
cc.ckd <- dt.diagorscrip(diag.regex="^N18[45]|^Z49[0-2]|^Z94[02]", diag.name="ckd")

## FIXME
#ids.kidneytransplant <- unique(procedures$ANON_ID[grep("^M01[1234589]",
#                                                       procedures$MAIN_OPERATION)])

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

## these drugs listed by HPS pharmacist as used for multiple sclerosis
## interferon beta 080204M, Glatiramer acetate 0802040U0, Natalizumab 0802040W0
## Dimethyl fumar 0802040AK, Teriflunomide 0802040AL, Alemtuzumab 0802030
## no records in scrips for these drugs
## 526 records in scrips[substr(scrips$bnf_paragraph_code, 1, 5) == "08020", ]

############## Liver disease #############################################

liver.grep.string <- "^C22\\.?0|^I85\\.?0|^I98\\.?3|^K70\\.?[234|^K71\\.?7|^K72\\.?[019]|^K72\\.?[019|^K73|^K74\\.?[023456]|^K76\\.?7|^R18"

cc.liver <- dt.diagorscrip(diag.regex=liver.grep.string,
                           diag.name="liver")
 
#### Immunodeficiency and immunosuppression #################################

cc.immune <- dt.diagorscrip(diag.regex="^B2[0-3|^D8[0-9]",
                            bnf.regex="^(0802)|(050301)",
                            diag.name="immune")

############# autoimmune rheumatology

cc.rheumatol <- dt.diagorscrip(diag.regex="^M0[5-8]|^M3[0-5]",
                                  diag.name="rheumatol")

setkey(rheumatol.wide, anon_id)
setkey(cc.specimendate, anon_id)
cc.rheumatol.list <- rheumatol.wide[cc.specimendate]
setkey(cc.specimendate, anon_id, specimen_date)

############# IBD

cc.ibd <- dt.diagorscrip(diag.regex="^K5[01]",
                                  diag.name="ibd")

############# neoplasms ################

cc.neoplasm <- dt.diagorscrip(diag.regex="^C[0-9]|^D[0-4]",
                              bnf.regex="^0801",
                              diag.name="neoplasm")

cc.bloodcancer <- dt.diagorscrip(diag.regex="^C8[1-8]|^C9[0-6]",
                                 diag.name="bloodcancer")

cc.bloodcancer.smr06 <- cc.smr06[grep("^C8[1-8]|^C9[0-6]", icd10)]
cc.bloodcancer <- cc.bloodcancer.smr06[cc.bloodcancer]
with(cc.bloodcancer, table(bloodcancer, !is.na(icd10)))
cc.bloodcancer[!is.na(icd10), bloodcancer := 1]
cc.bloodcancer <- cc.bloodcancer[, .(anon_id, specimen_date, bloodcancer)]

###### disorders of esophagus, stomach and duodenum ############################

cc.esoph.stomach.duod <-  dt.diagorscrip(diag.regex="^K2[0-9]|^K3[01]",
                                         diag.name="esoph.stomach.duod")

cc.comorbid <- cc.heart.other[cc.ihd]
cc.comorbid <- cc.circulatory.other[cc.comorbid]
cc.comorbid <- cc.ckd[cc.comorbid]
cc.comorbid <- cc.chronresp[cc.comorbid]
cc.comorbid <- cc.neuro[cc.comorbid]
cc.comorbid <- cc.liver[cc.comorbid]
cc.comorbid <- cc.immune[cc.comorbid]
cc.comorbid <- cc.neoplasm[cc.comorbid]
cc.comorbid <- cc.esoph.stomach.duod[cc.comorbid]
cc.comorbid <- cc.neoplasm[cc.comorbid]
cc.comorbid <- cc.bloodcancer[cc.comorbid]
cc.comorbid <- cc.rheumatol[cc.comorbid]
cc.comorbid <- cc.ibd[cc.comorbid]
                            
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

gc()
