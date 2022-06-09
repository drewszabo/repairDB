# repairDB

## Description
A solution to complete missing information from a list of chemicals.<br>
Currently, the format required for the .csv is very strict. The naming convention for each column follows the SusDat(S0) List as closely as possible. Columns must be named as follows:<br>
##### Name, Name_IUPAC, CAS_RN, PubChem_CID, ChemSpiderID, DTXSID, SMILES, InChI, InChIKey, Molecular_Formula, Monoiso_Mass, XlogP3<br>
A template is available on the repairDB repo: [template.csv](https://github.com/drewszabo/repairDB/template.csv)

## Dependancies
```
library(tidyverse)
library(webchem)
```

## Usage
```
repairDB(filename = "database.csv", reference = "InChIKey", susdat = "susdat.csv")
```

## Arguments
#### filename
Path and filename of .csv file for repair
#### reference
A column with complete information to reference agianst. Eg. CAS_RN, SMILES, InChI, or InChIKey
#### susdat
Path and filename of .csv file containing the S0 SusDat List (available from: https://www.norman-network.com/nds/SLE/)

## Value
A data.frame object.

## Author
Drew Szabo, PhD <br>
GitHub: [drewszabo](https://github.com/drewszabo) <br>
ORCID: [0000-0002-0089-9218](https://orcid.org/0000-0002-0089-9218)
