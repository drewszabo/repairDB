require(tidyverse)
require(webchem)
require(Rdisop)
require(plotly)
require(cowplot)

#From Scratch

kemiDB <- read.csv("KEMI_MarketList_29122017_fixedIKs.csv", check.names = FALSE)

kemiPriority <- kemiDB %>%
  mutate(PubChem_CID = NA,
         DTXSID = NA,
         logP_EPI = NA,
         InChI = NA,
         Molecular_Formula = NA,
         Name_IUPAC = NA,
         Monoiso_Mass = NA,
         logP_PubChem = NA,
         "Prob. +ESI" = NA,
         "Prob. -ESI" = NA,
         "M+H+" = NA,
         "M-H-" = NA)

kemiPriority[kemiPriority == ""] <- NA

kemiPriority$PubChem_CID <- as.numeric(kemiPriority$PubChem_CID)
kemiPriority$`ExposureScore-Wat (max27)` <- as.numeric(kemiPriority$`ExposureScore-Wat (max27)`)
kemiPriority$`HazScore_EcoChronic(max9)` <- as.numeric(kemiPriority$`HazScore_EcoChronic(max9)`)


#Load data previously used for toxicity prediction

kemiPriority <- read.csv("kemiPriority.csv", check.names = FALSE)

#---------------------------
#Remove duplicates based on CAS, InChI and name (Working on based on presence of duplicates)
#---------------------------

kemiPriority <- kemiPriority %>%
  distinct(CASno, .keep_all = TRUE)

kemiPriority <- kemiPriority %>%
  distinct(Name, .keep_all = TRUE)

kemiPriority <- kemiPriority %>%
  distinct(InChIKey, .keep_all = TRUE)

#---------------------------
#Remove invalid CAS, InChIKey, SMILES & CID
#---------------------------

kemiPriority <- kemiPriority %>%
  rowwise() %>%
  mutate(CASno = ifelse(is.cas(CASno),
                        as.character(CASno),
                        NA),
         InChIKey = ifelse(is.inchikey_format(InChIKey),
                           as.character(InChIKey),
                           NA),
         Smiles = ifelse(is.smiles(Smiles),
                         as.character(Smiles),
                         NA)) %>%
  ungroup()


kemiPriority <- kemiPriority %>%
  rowwise() %>%
  mutate(PubChem_CID = ifelse(PubChem_CID == 0,
                              NA,
                              as.character(PubChem_CID))) %>%
  ungroup()


#---------------------------
#Grab EMPTY variables from NormanDB based on PubChem_CID
#---------------------------

#PubChem_CID
kemiPriority$PubChem_CID[is.na(kemiPriority$PubChem_CID)] <- 
  susDatDB$PubChem_CID[match(kemiPriority$Smiles, susDatDB$SMILES_Dashboard)[which(is.na(kemiPriority$PubChem_CID))]]

#DTXSID (USEPA CompTox)
kemiPriority$DTXSID[is.na(kemiPriority$DTXSID)] <- 
  susDatDB$DTXSID[match(kemiPriority$Smiles, susDatDB$SMILES_Dashboard)[which(is.na(kemiPriority$DTXSID))]]

#CAS_RN
kemiPriority$CASno[is.na(kemiPriority$CASno)] <- 
  susDatDB$CAS_RN_Dashboard[match(kemiPriority$Smiles, susDatDB$SMILES_Dashboard)[which(is.na(kemiPriority$CASno))]]

#logP - EPI-Suite
kemiPriority$logP_EPI[is.na(kemiPriority$logP_EPI)] <- 
  susDatDB$logKow_EPISuite[match(kemiPriority$Smiles, susDatDB$SMILES_Dashboard)[which(is.na(kemiPriority$logP_EPI))]]

#Prob. +ESI
kemiPriority$"Prob. +ESI"[is.na(kemiPriority$"Prob. +ESI")] <- 
  susDatDB$"Prob. +ESI"[match(kemiPriority$PubChem_CID, susDatDB$PubChem_CID)[which(is.na(kemiPriority$"Prob. +ESI"))]]

#Prob. -ESI
kemiPriority$"Prob. -ESI"[is.na(kemiPriority$"Prob. -ESI")] <- 
  susDatDB$"Prob. -ESI"[match(kemiPriority$PubChem_CID, susDatDB$PubChem_CID)[which(is.na(kemiPriority$"Prob. -ESI"))]]

#[M+H]+
kemiPriority$"M+H+"[is.na(kemiPriority$"M+H+")] <- 
  susDatDB$"M+H+"[match(kemiPriority$Smiles, susDatDB$SMILES_Dashboard)[which(is.na(kemiPriority$"M+H+"))]]

#[M-H]-
kemiPriority$"M-H-"[is.na(kemiPriority$"M-H-")] <- 
  susDatDB$"M-H-"[match(kemiPriority$Smiles, susDatDB$SMILES_Dashboard)[which(is.na(kemiPriority$"M-H-"))]]


#---------------------------
#Query empty pubChem for CID from InChIKey, CAS, then Name
#---------------------------

kemiPriority <- kemiPriority %>%
  rowwise() %>%
  mutate(PubChem_CID = ifelse(is.na(PubChem_CID), 
                              unlist(get_cid(InChIKey, from = "inchikey", match = "na", domain = "compound", verbose = TRUE)[[2]]),
                              as.character(PubChem_CID))) %>%
  ungroup()

kemiPriority <- kemiPriority %>%
  rowwise() %>%
  mutate(PubChem_CID = ifelse(is.na(PubChem_CID), 
                              unlist(get_cid(Smiles, from = "smiles", match = "na", domain = "compound", verbose = TRUE)[[2]]),
                              as.character(PubChem_CID))) %>%
  ungroup()

kemiPriority <- kemiPriority %>%
  rowwise() %>%
  mutate(PubChem_CID = ifelse(is.na(PubChem_CID) & !is.na(CASno), 
                              unlist(get_cid(CASno, from = "xref/rn", match = "na", domain = "compound", verbose = TRUE)[[2]]),
                              as.character(PubChem_CID))) %>%
  ungroup()

kemiPriority <- kemiPriority %>%
  rowwise() %>%
  mutate(PubChem_CID = ifelse(is.na(PubChem_CID), 
                              unlist(get_cid(Name, from = "name", match = "na", domain = "compound", verbose = TRUE)[[2]]),
                              as.character(PubChem_CID))) %>%
  ungroup()


#---------------------------
#Update pubChemProp list. 
#Unknown limit >10000. Slice data into <=10000 rows for pc_prop query
#---------------------------




pubChemProp1 <- kemiPriority %>%
  distinct(PubChem_CID) %>%
  slice(1:10000)

pubChemProp2 <- kemiPriority %>%
  distinct(PubChem_CID) %>%
  slice(10001:20000)

pubChemProp3 <- kemiPriority %>%
  distinct(PubChem_CID) %>%
  slice(20001:30000)

pubChemProp4 <- kemiPriority %>%
  distinct(PubChem_CID) %>%
  slice(30001:40000)

pubChemProp5 <- kemiPriority %>%
  distinct(PubChem_CID) %>%
  slice(40001:50000)


pubChemProp <- pc_prop(pubChemProp1$PubChem_CID)

pubChemProp <- pubChemProp %>%
  add_row(pc_prop(pubChemProp2$PubChem_CID))

pubChemProp <- pubChemProp %>%
  add_row(pc_prop(pubChemProp3$PubChem_CID))

pubChemProp <- pubChemProp %>%
  add_row(pc_prop(pubChemProp4$PubChem_CID))

pubChemProp <- pubChemProp %>%
  add_row(pc_prop(pubChemProp5$PubChem_CID))

#---------------------------
#Match compoundDB with pubChemProp
#---------------------------

#Molecular Formula
kemiPriority$Molecular_Formula[is.na(kemiPriority$Molecular_Formula)] <- 
  pubChemProp$MolecularFormula[match(kemiPriority$PubChem_CID, pubChemProp$CID)[which(is.na(kemiPriority$Molecular_Formula))]]

#IUPAC Name
kemiPriority$Name_IUPAC[is.na(kemiPriority$Name_IUPAC)] <- 
  pubChemProp$IUPACName[match(kemiPriority$PubChem_CID, pubChemProp$CID)[which(is.na(kemiPriority$Name_IUPAC))]]

#SMILES
kemiPriority$Smiles[is.na(kemiPriority$Smiles)] <- 
  pubChemProp$CanonicalSMILES[match(kemiPriority$PubChem_CID, pubChemProp$CID)[which(is.na(kemiPriority$Smiles))]]

#Monoisotopic Mass
kemiPriority$Monoiso_Mass[is.na(kemiPriority$Monoiso_Mass)] <- 
  pubChemProp$MonoisotopicMass[match(kemiPriority$PubChem_CID, pubChemProp$CID)[which(is.na(kemiPriority$Monoiso_Mass))]]

#InChI
kemiPriority$InChI[is.na(kemiPriority$InChI)] <- 
  pubChemProp$InChI[match(kemiPriority$PubChem_CID, pubChemProp$CID)[which(is.na(kemiPriority$InChI))]]

#logP - PubChem
kemiPriority$logP_PubChem[is.na(kemiPriority$logP_PubChem)] <- 
  pubChemProp$XLogP[match(kemiPriority$PubChem_CID, pubChemProp$CID)[which(is.na(kemiPriority$logP_PubChem))]]

#InChIKey
kemiPriority$InChIKey[is.na(kemiPriority$InChIKey)] <- 
  pubChemProp$InChIKey[match(kemiPriority$PubChem_CID, pubChemProp$CID)[which(is.na(kemiPriority$InChIKey))]]

#Charge
kemiPriority$Charge <- 
  pubChemProp$Charge[match(kemiPriority$PubChem_CID, pubChemProp$CID)]

#---------------------------
#Write kemiPriority file
#---------------------------

write.csv(kemiPriority, "kemiPriority.csv")


#---------------------------
#Filter based on distinct PubChem CID
#This also removes compounds not matched with the PubChem library
#49264 results
#---------------------------

kemiShortlist <- kemiPriority %>%
  filter(PubChem_CID > 0) %>%
  distinct(PubChem_CID, .keep_all = TRUE)

kemiPriority$RiskScore <- (kemiPriority$ExposureScore/27*0.5) + (kemiPriority$ToxScore/10*0.5)

kemiPriority <- kemiPriority %>%
  mutate(ToxScore = ifelse(LC50_predicted < -5.011754,
                           10,
                           ifelse(LC50_predicted < -4.097302,
                                  9,
                                  ifelse(LC50_predicted < -3.182849,
                                         8,
                                         ifelse(LC50_predicted < -2.26837,
                                                7,
                                                ifelse(LC50_predicted < -1.353944,
                                                       6,
                                                       ifelse(LC50_predicted < -0.4394921,
                                                              5,
                                                              ifelse(LC50_predicted < 0.4749602,
                                                                     4,
                                                                     ifelse(LC50_predicted < 1.389413,
                                                                            3,
                                                                            ifelse(LC50_predicted < 2.303865,
                                                                                   2,
                                                                                   ifelse(LC50_predicted < 3.218317,
                                                                                          1,
                                                                                          0))))))))))) %>%
  mutate(`Prob. ESI` = ifelse(`Prob. +ESI`>`Prob. -ESI`,
                              "Positive",
                              "Negative"))


#---------------------------
# Filter based on m/z >100, SMILE contains "C", ESI Prob >0.3
# Limited by SusDat list now because of ESI prediction
# 30941 results
#---------------------------

kemiShortlist <- kemiPriority %>%
  distinct(PubChem_CID, .keep_all = TRUE) %>%
  filter(PubChem_CID > 0,
         Monoiso_Mass > 100,
         `Prob. +ESI` > 0.5 | `Prob. -ESI` > 0.5,
         !is.na(`M+H+`))

#---------------------------
#Filter based on weighted RiskScore
#Selected Top 700 for +ve/-ve (little bit more/less if Score is the same for cutoff)
#---------------------------

kemiShortlist <- kemiShortlist %>%
  group_by(`Prob. ESI`) %>%
  top_n(n =700, wt =RiskScore)


#---------------------------
# Get m/z for M-Li, -Na, -K, -NH4 adducts
# Not in use as filtered by M+H+
#---------------------------

kemiShortlist <- kemiShortlist %>%
  rowwise() %>%
  mutate(`M-H-` = ifelse(is.na(`M-H-`) & `Prob. ESI` == "Negative" & grepl("Li", Molecular_Formula),
                         getMolecule(gsub("Li", "", Molecular_Formula), z = -1)[[3]],
                         `M-H-`)) %>%
  ungroup()

kemiShortlist <- kemiShortlist %>%
  rowwise() %>%
  mutate(`M-H-` = ifelse(is.na(`M-H-`) & `Prob. ESI` == "Negative" & grepl("Na", Molecular_Formula),
                         getMolecule(gsub("Na", "", Molecular_Formula), z = -1)[[3]],
                         `M-H-`)) %>%
  ungroup()

kemiShortlist <- kemiShortlist %>%
  rowwise() %>%
  mutate(`M-H-` = ifelse(is.na(`M-H-`) & `Prob. ESI` == "Negative" & grepl("K", Molecular_Formula),
                         getMolecule(gsub("K", "", Molecular_Formula), z = -1)[[3]],
                         `M-H-`)) %>%
  ungroup()


#---------------------------
#Write kemiShortlist file
#---------------------------

write.csv(kemiShortlist, "kemiShortlist.csv")

#---------------------------
#Visualise Groups
#---------------------------

#contains?

PFAS.test <- kemiShortlist %>%
  filter(grepl("Br", Molecular_Formula))

#basic plots
ggplot() +
  geom_density(aes(x = as.numeric(kemiShortlist$LC50_predicted), y = ..density.., color = kemiShortlist$`Prob. ESI`))
ggplot(data = kemiShortlist) +
  geom_histogram(aes(x = as.numeric(LC50_predicted), color = `Prob. ESI`),
                 position = "identity")

plot.mass <- ggplot() +
  geom_histogram(data = kemiPriority %>% filter(`Prob. ESI` == "Positive"),
               aes(x = Monoiso_Mass, y = ..density.., fill = `Prob. ESI`)) +
  geom_histogram(data = kemiPriority %>% filter(`Prob. ESI` == "Negative"),
               aes(x = Monoiso_Mass, y = -..density.., fill = `Prob. ESI`)) +
  scale_x_continuous(name = "Monoisotopic Mass") +
  scale_y_continuous(name = "",
                     limits = c(-0.005,0.005)) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

plot.mass

plot.logP <- ggplot() +
  geom_histogram(data = kemiPriority %>% filter(`Prob. ESI` == "Positive"),
                 aes(x = logP_PubChem, y = ..density.., fill = `Prob. ESI`)) +
  geom_histogram(data = kemiPriority %>% filter(`Prob. ESI` == "Negative"),
                 aes(x = logP_PubChem, y = -..density.., fill = `Prob. ESI`)) +
  scale_x_continuous(name = "XlogP") +
  scale_y_continuous(name = "",
                     limits = c(-0.25,0.25)) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

plot.logP

plot.exp <- ggplot() +
  geom_histogram(data = kemiPriority %>% filter(`Prob. ESI` == "Positive"),
                 aes(x = ExposureScore, y = ..density.., fill = `Prob. ESI`)) +
  geom_histogram(data = kemiPriority %>% filter(`Prob. ESI` == "Negative"),
                 aes(x = ExposureScore, y = -..density.., fill = `Prob. ESI`)) +
  scale_x_continuous(name = "Exposure Score") +
  scale_y_continuous(name = "",
                     limits = c(-0.5,0.5)) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

plot.exp

plot.tox <- ggplot() +
  geom_histogram(data = kemiPriority %>% filter(`Prob. ESI` == "Positive"),
                 aes(x = LC50_predicted, y = ..density.., fill = `Prob. ESI`)) +
  geom_histogram(data = kemiPriority %>% filter(`Prob. ESI` == "Negative"),
                 aes(x = LC50_predicted, y = -..density.., fill = `Prob. ESI`)) +
  scale_x_continuous(name = "Predicted LC50") +
  scale_y_continuous(name = "",
                     limits = c(-1,1)) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

plot.tox

plot.risk <- ggplot() +
  geom_histogram(data = kemiPriority %>% filter(`Prob. ESI` == "Positive"),
                 aes(x = RiskScore, y = ..density.., fill = `Prob. ESI`)) +
  geom_histogram(data = kemiPriority %>% filter(`Prob. ESI` == "Negative"),
                 aes(x = RiskScore, y = -..density.., fill = `Prob. ESI`)) +
  scale_x_continuous(name = "Risk Score") +
  scale_y_continuous(name = "",
                     limits = c(-25,25)) +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

plot.risk

cowplot::plot_grid(plot.mass + theme(legend.position = "none"),
          plot.logP + theme(legend.position = "none"),
          plot.exp + theme(legend.position = "none"),
          plot.tox + theme(legend.position = "none"),
          plot.risk + theme(legend.position = "none"), labels = "AUTO", nrow = 1)


ggsave("market.png",
       width = 10,
       height = 3)

#3DPlot
ggplot(data = kemiShortlist) +
  geom_point(aes(x = Monoiso_Mass, y = logP_PubChem, color = ifelse(grepl("F", Molecular_Formula),
                                                                    "red",
                                                                    `Prob. ESI`)))

fig <- plot_ly(kemiShortlist, x = ~Monoiso_Mass, y = ~logP_PubChem, z = ~LC50_predicted, color = ~`Prob. ESI`,
               text = ~paste('Compound:', Name, '<br>CAS:', CASno))
fig

#PCA?
kemiShortlist$ExposureScore <- as.numeric(kemiShortlist$ExposureScore)
kemiShortlist$HazardScore <- as.numeric(kemiShortlist$HazardScore)
kemiShortlist$logP_EPI <- as.numeric(kemiShortlist$logP_EPI)
kemiShortlist$logP_PubChem <- as.numeric(kemiShortlist$logP_PubChem)
kemiShortlist$Monoiso_Mass <- as.numeric(kemiShortlist$Monoiso_Mass)
kemiShortlist$LC50_predicted <- as.numeric(kemiShortlist$LC50_predicted)

kemiShortlist %>% select(ExposureScore, HazardScore, logP_PubChem, Monoiso_Mass, LC50_predicted)

pca.kemi <- prcomp(kemiShortlist %>% select(ExposureScore, HazardScore, logP_PubChem, Monoiso_Mass, LC50_predicted) %>% drop_na())
summary(pca.kemi)
ggbiplot(pca.kemi, groups = "Prob. ESI")


#Tox distribution score

shapiro.test(kemiPriority$LC50_predicted)

ggplot(data = kemiPriority) +
  geom_histogram(aes(x = LC50_predicted, y = ..density..)) +
  geom_vline(xintercept = mean(kemiPriority$LC50_predicted, na.rm = TRUE)) +
  geom_vline(xintercept = mean(kemiPriority$LC50_predicted, na.rm = TRUE) + sd(kemiPriority$LC50_predicted, na.rm = TRUE)*c(1,2,3,4), linetype = 2) +
  geom_vline(xintercept = mean(kemiPriority$LC50_predicted, na.rm = TRUE) - sd(kemiPriority$LC50_predicted, na.rm = TRUE)*c(1,2,3,4), linetype = 2) +
  scale_x_continuous(name = "Predicted LC50") +
  scale_y_continuous(name = "") +
  theme_classic() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

ggsave("test.png", height = 2, width = 5)
