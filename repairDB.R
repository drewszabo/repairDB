library(tidyverse)
library(webchem)
library(rcdk)
library(rJava)

repairDB <- function(filename, reference, susdat) {
  
  compoundDB <- read.csv(filename, check.names = FALSE)
  susDatDB <- read.csv(susdat, check.names = FALSE)
  compoundDB[compoundDB == ""] <- NA

  #-----Create Backup-----
  if (askYesNo("Create a backup?")) {
    write.csv(compoundDB, file = paste("backup_", as.character(format(Sys.time(), "%Y-%m-%d")),".csv", sep = ""))
  } else (NA)
  
  #-----Remove invalid entries-----
  compoundDB <- compoundDB %>%
    rowwise() %>%
    mutate(SMILES = ifelse(is.smiles(SMILES),
                           as.character(SMILES),
                           NA),
           CAS_RN = ifelse(is.cas(CAS_RN),
                           as.character(CAS_RN),
                           NA)) %>%
             ungroup()
  
  #-----Grab PubChem_CID from Reference (SusDat)-----
  if(reference == "SMILES"){
    compoundDB$PubChem_CID[is.na(compoundDB$PubChem_CID)] <- 
      susDatDB$PubChem_CID[match(compoundDB$SMILES, susDatDB$SMILES_Dashboard)[which(is.na(compoundDB$PubChem_CID))]]
    compoundDB[compoundDB == ""] <- NA
    compoundDB <- compoundDB %>%
      rowwise() %>%
      mutate(PubChem_CID = ifelse(is.na(PubChem_CID) & !is.na(SMILES),
                                  unlist(get_cid(SMILES, from = "smiles", domain = "compound", match = "na", verbose = TRUE)[[2]]),
                                  as.character(PubChem_CID))) %>%
      ungroup()
  } else if(reference == "InChI"){
    compoundDB$PubChem_CID[is.na(compoundDB$PubChem_CID)] <- 
      susDatDB$PubChem_CID[match(compoundDB$InChI, susDatDB$StdInChI)[which(is.na(compoundDB$PubChem_CID))]]
    compoundDB[compoundDB == ""] <- NA
    compoundDB <- compoundDB %>%
      rowwise() %>%
      mutate(PubChem_CID = ifelse(is.na(PubChem_CID) & !is.na(InChI),
                                  unlist(get_cid(InChI, from = "inchi", domain = "compound", match = "na", verbose = TRUE)[[2]]),
                                  as.character(PubChem_CID))) %>%
      ungroup()
  } else if(reference == "InChIKey"){
    compoundDB$PubChem_CID[is.na(compoundDB$PubChem_CID)] <- 
      susDatDB$PubChem_CID[match(compoundDB$InChIKey, susDatDB$StdInChIKey)[which(is.na(compoundDB$PubChem_CID))]]
    compoundDB[compoundDB == ""] <- NA
    compoundDB <- compoundDB %>%
      rowwise() %>%
      mutate(PubChem_CID = ifelse(is.na(PubChem_CID) & !is.na(InChIKey),
                                  unlist(get_cid(InChIKey, from = "inchikey", domain = "compound", match = "na", verbose = TRUE)[[2]]),
                                  as.character(PubChem_CID))) %>%
      ungroup()
  } else if(reference == "CAS_RN"){
    compoundDB$PubChem_CID[is.na(compoundDB$PubChem_CID)] <- 
      susDatDB$PubChem_CID[match(compoundDB$CAS_RN, susDatDB$CAS_RN_Dashboard)[which(is.na(compoundDB$PubChem_CID))]]
    compoundDB[compoundDB == ""] <- NA
    compoundDB <- compoundDB %>%
      rowwise() %>%
      mutate(PubChem_CID = ifelse(is.na(PubChem_CID) & !is.na(CAS_RN),
                                  unlist(get_cid(CAS_RN, from = "xref/rn", domain = "compound", match = "na", verbose = TRUE)[[2]]),
                                  as.character(PubChem_CID))) %>%
      ungroup()
  } else {NA}
  
  compoundDB <- compoundDB %>%
    rowwise() %>%
    mutate(PubChem_CID = ifelse(is.na(PubChem_CID), 
                                unlist(get_cid(Name, from = "name", match = "na", domain = "compound", verbose = TRUE)[[2]]),
                                as.character(PubChem_CID))) %>%
    ungroup()
  compoundDB[compoundDB == ""] <- NA
  
  #-----Query PubChem for information-----
  pubChemProp <- pc_prop(compoundDB$PubChem_CID)
  
  #-----Match information with compoundDB-----
  
  #Molecular Formula
  compoundDB$Molecular_Formula[is.na(compoundDB$Molecular_Formula)] <- 
    pubChemProp$MolecularFormula[match(compoundDB$PubChem_CID, pubChemProp$CID)[which(is.na(compoundDB$Molecular_Formula))]]
  
  #IUPAC Name
  compoundDB$Name_IUPAC[is.na(compoundDB$Name_IUPAC)] <- 
    pubChemProp$IUPACName[match(compoundDB$PubChem_CID, pubChemProp$CID)[which(is.na(compoundDB$Name_IUPAC))]]
  
  #Monoisotopic Mass
  compoundDB$Monoiso_Mass[is.na(compoundDB$Monoiso_Mass)] <- 
    pubChemProp$MonoisotopicMass[match(compoundDB$PubChem_CID, pubChemProp$CID)[which(is.na(compoundDB$Monoiso_Mass))]]
  
  #SMILES
  compoundDB$SMILES[is.na(compoundDB$SMILES)] <- 
    pubChemProp$CanonicalSMILES[match(compoundDB$PubChem_CID, pubChemProp$CID)[which(is.na(compoundDB$SMILES))]]
  
  #InChI
  compoundDB$InChI[is.na(compoundDB$InChI)] <- 
    pubChemProp$InChI[match(compoundDB$PubChem_CID, pubChemProp$CID)[which(is.na(compoundDB$InChI))]]
  
  #logP - PubChem
  compoundDB$XlogP3[is.na(compoundDB$XlogP3)] <- 
    pubChemProp$XLogP[match(compoundDB$PubChem_CID, pubChemProp$CID)[which(is.na(compoundDB$XlogP3))]]
  
  #InChIKey
  compoundDB$InChIKey[is.na(compoundDB$InChIKey)] <- 
    pubChemProp$InChIKey[match(compoundDB$PubChem_CID, pubChemProp$CID)[which(is.na(compoundDB$InChIKey))]]
  
  #CAS_RN
  compoundDB$CAS_RN[is.na(compoundDB$CAS_RN)] <- 
    susDatDB$CAS_RN_Dashboard[match(compoundDB$PubChem_CID, susDatDB$PubChem_CID)[which(is.na(compoundDB$CAS_RN))]]
  
  #DTXSID
  compoundDB$DTXSID[is.na(compoundDB$DTXSID)] <- 
    susDatDB$DTXSID[match(compoundDB$PubChem_CID, susDatDB$PubChem_CID)[which(is.na(compoundDB$DTXSID))]]
  
  #ChemSpiderID
  compoundDB$ChemSpiderID[is.na(compoundDB$ChemSpiderID)] <- 
    susDatDB$ChemSpiderID[match(compoundDB$PubChem_CID, susDatDB$PubChem_CID)[which(is.na(compoundDB$ChemSpiderID))]]
  
  #M+H+
  compoundDB$`M+H+`[is.na(compoundDB$`M+H+`)] <- 
    susDatDB$`M+H+`[match(compoundDB$PubChem_CID, susDatDB$PubChem_CID)[which(is.na(compoundDB$`M+H+`))]]
  
  #M-H-
  compoundDB$`M-H-`[is.na(compoundDB$`M-H-`)] <- 
    susDatDB$`M-H-`[match(compoundDB$PubChem_CID, susDatDB$PubChem_CID)[which(is.na(compoundDB$`M-H-`))]]
  
  #-----Removed duplicates based on Name-----
  compoundDB <- compoundDB %>%
    distinct(Name, .keep_all = TRUE)
  
  #-----Return compoundDB-----
  
  return(compoundDB)
  
}
