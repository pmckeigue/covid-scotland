dt.diagorscrip
function(diag.regex=NULL, bnf.regex=NULL, diag.name) {
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

################  IHD ###################################
## nitrates are BNF code 020601

cc.ihd <- dt.diagorscrip(diag.regex="^I2[0-5]", bnfregex="02060(1|3)", diag.name="ihd")

## procedure codes for CABG and PTCA
ids.procedures.IHD <- unique(procedures$ANON_ID[grep("^K4[012349]|^K50", procedures$MAIN_OPERATION)])

##### other heart disease ####################################
## heart disease is I05 to I52

ids.icd.heart.other <- unique(diagnoses$ANON_ID[grep("^I0[01256789]|^I1[0-5]|^I2[6-8]|^I3[0-9]|^I4[0-9]|^I5[0-2]",
                                  diagnoses$ICD10)])
ids.bnf.heart.other <-
    unique(scrips$ANON_ID[substr(as.character(scrips$bnf_paragraph_code), 1, 4) == "0203"])  # anti-arrhythmics
table(ids.bnf.heart.other %in% ids.icd.heart.other)

ids.procedures.heart.other <- unique(procedures$ANON_ID[grep("^K57", procedures$MAIN_OPERATION)])

ids.heart.other <- unique(c(ids.icd.heart.other, ids.bnf.heart.other, ids.procedures.heart.other))
cc.all$heart.other.any <- as.factor(as.integer(cc.all$ANON_ID %in% ids.heart.other))
cc.severe$heart.other.any <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.heart.other))

#######################################################################
## other circulatory disease is I60 to I99

ids.icd.circulatory.other <-
    unique(diagnoses$ANON_ID[grep("^I[6-9]|^Z95", diagnoses$ICD10)])

cc.all$circulatory.other <-
    as.factor(as.integer(cc.all$ANON_ID %in% ids.icd.circulatory.other))
cc.severe$circulatory.other <-
  as.factor(as.integer(cc.severe$ANON_ID %in% ids.icd.circulatory.other))

############# chronic kidney disease ##########################

## this includes CKD stage 4
ids.icd.ckd <- unique(diagnoses$ANON_ID[grep("^N18[45]|^Z49[0-2]|^Z94[02]",
                                             diagnoses$ICD10)])

ids.kidneytransplant <- unique(procedures$ANON_ID[grep("^M01[1234589]",
                                                       procedures$MAIN_OPERATION)])
table(ids.kidneytransplant %in% ids.icd.ckd)
ids.ckd.any <- unique(c(ids.icd.ckd, ids.kidneytransplant))
cc.all$ckd.any <-  as.factor(as.integer(cc.all$ANON_ID %in% ids.ckd.any))
cc.severe$ckd.any <-  as.factor(as.integer(cc.severe$ANON_ID %in% ids.ckd.any))

##### asthma and chronic lower respiratory disease #################

ids.icd.asthma <- unique(diagnoses$ANON_ID[grep("^J4[56]", diagnoses$ICD10)])
ids.icd.chronresp <- unique(diagnoses$ANON_ID[grep("^J4[012347]|^J6[0-9]|^J70|^J8[0-6]|^J9[0-9]|^G47\\.?3",
                                                   diagnoses$ICD10)])
ids.bnf.broncho <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) >= 301 &
                                         as.integer(scrips$sectioncode) <= 303])
table(ids.icd.asthma %in% ids.bnf.broncho)
table(ids.icd.chronresp %in% ids.bnf.broncho)
ids.oad.any <- unique(c(ids.icd.asthma, ids.icd.chronresp, ids.bnf.broncho))
cc.all$oad.any <- as.factor(as.integer(cc.all$ANON_ID %in% ids.oad.any))
cc.severe$oad.any <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.oad.any))

##### Neurological disorders #######################

## include all Nervous chapter except G40 "Episodic and Paroxysmal Disorders"
## Helen says leave out G0 meningitis and encephalitis, and G5 local neuropathies

## also include F03 dementia NOS
ids.icd.neuro <- unique(diagnoses$ANON_ID[grep("^F03|^G[1236789]", diagnoses$ICD10)])
ids.bnf.neuro <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 409 |
                                         as.integer(scrips$sectioncode) == 411])
table(ids.bnf.neuro %in% ids.icd.neuro)
## these drugs listed by HPS pharmacist as used for multiple sclerosis
## interferon beta 080204M, Glatiramer acetate 0802040U0, Natalizumab 0802040W0
## Dimethyl fumar 0802040AK, Teriflunomide 0802040AL, Alemtuzumab 0802030
## no records in scrips for these drugs
## 526 records in scrips[substr(scrips$bnf_paragraph_code, 1, 5) == "08020", ]
ids.neuro.any <- unique(c(ids.icd.neuro, ids.bnf.neuro))
cc.all$neuro.any <- as.factor(as.integer(cc.all$ANON_ID %in% ids.neuro.any))
cc.severe$neuro.any <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.neuro.any))

############## Liver disease #############################################

liver.grep.string <- "^C22\\.?0|^I85\\.?0|^I98\\.?3|^K70\\.?[234|^K71\\.?7|^K72\\.?[019]|^K72\\.?[019|^K73|^K74\\.?[023456]|^K76\\.?7|^R18"
table(grep(liver.grep.string, diagnoses$ICD10, value=TRUE))
ids.icd.liver <- unique(diagnoses$ANON_ID[grep(liver.grep.string, diagnoses$ICD10)])
cc.all$liver.any <- as.factor(as.integer(cc.all$ANON_ID %in% ids.icd.liver))
cc.severe$liver.any <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.icd.liver))

#### Immunodeficiency and immunosuppression #################################

## immune.any includes primary immunodeficiency and secondary immunosuppression
ids.icd.immune <- unique(diagnoses$ANON_ID[grep("^B2[0-3|^D8[0-9]", diagnoses$ICD10)])

## 802 other immunomodulating drugs
## Methotrexate and chloroquine appear in musculoskeletal chapter 
ids.bnf.immune <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 802 |
                                        as.integer(scrips$paracode) == 50301])

ids.immune.any <- unique(c(ids.icd.immune, ids.bnf.immune))
cc.all$immune.any <- as.factor(as.integer(cc.all$ANON_ID %in% ids.immune.any))
cc.severe$immune.any <- as.factor(as.integer(cc.severe$ANON_ID %in% ids.immune.any))

############# neoplasms ################

ids.icd.neoplasm <- unique(diagnoses$ANON_ID[grep("^C[0-9]|^D[0-4]", diagnoses$ICD10)])
ids.bnf.neoplasm <- unique(scrips$ANON_ID[as.integer(scrips$sectioncode) == 801])
ids.neoplasm.any <- unique(c(ids.icd.neoplasm, ids.bnf.neoplasm))

cc.all <- mutate(cc.all,
                    neoplasm.any = as.factor(as.integer(ANON_ID %in% ids.neoplasm.any |
                                             can.reg==1)))
cc.severe <- mutate(cc.severe,
                 neoplasm.any = as.factor(as.integer(ANON_ID %in% ids.neoplasm.any |
                                                       can.reg==1)))

###### disorders of esophagus, stomach and duodenum ############################

ids.esoph.stomach.duod <-  unique(diagnoses$ANON_ID[grep("^K2[0-9]|^K3[01]", diagnoses$ICD10)])
cc.all$esoph.stomach.duod <-
    as.factor(as.integer(cc.all$ANON_ID %in% ids.esoph.stomach.duod))
cc.severe$esoph.stomach.duod <-
  as.factor(as.integer(cc.severe$ANON_ID %in% ids.esoph.stomach.duod))


