# repairDB

## Description
A solution to complete missing information from a list of chemicals.<p>
The naming convention for each column follows the SusDat(S0) List as closely as possible. Each column will be completed from the SusDat List where possible. However, unmatched chemicals will then be queried to the PubChem database and this can take some time. If the PubChem query returns multiple results, it will return a `NA` value to avoid false positive results. At this time, some chemicals will not find results so manual entry of the `PubChem_CID` may be neccesary to ensure complete table is returned. <p> 
Currently, the format required for the .csv is very strict. Columns must be named as follows:<br>
`Name`, `Name_IUPAC`, `CAS_RN`, `PubChem_CID`, `ChemSpiderID`, `DTXSID`, `SMILES`, `InChI`, `InChIKey`, `Molecular_Formula`, `Monoiso_Mass`, `XlogP3`<p>
A template is available on the repairDB repo: [template.csv](https://github.com/drewszabo/repairDB/blob/main/template.csv)

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
Path and filename of `.csv` file for repair
#### reference
A column with complete information to reference agianst. Eg. `CAS_RN`, `SMILES`, `InChI`, or `InChIKey`
#### susdat
Path and filename of `.csv` file containing the S0 SusDat List (available from: https://www.norman-network.com/nds/SLE/)

## Value
A data.frame object.

## Example
Download [example.csv](https://github.com/drewszabo/repairDB/blob/main/template.csv) from repo. <br>
Download [susdat_2022-04-09-194500.csv](https://www.norman-network.com/sites/default/files/files/suspectListExchange/010322Update/susdat_2022-01-18-104316.csv) from Norman SLE.

```
compoundDB <- repairDB("example.csv", reference = "InChIKey", susdat = "susdat_2022-04-09-194500.csv"
compoundDB
```

```
              
```

## Author
Drew Szabo <br>
GitHub: [drewszabo](https://github.com/drewszabo) <br>
ORCID: [0000-0002-0089-9218](https://orcid.org/0000-0002-0089-9218)

## Copyright
©️ KruveLab 2022. Stockholm University, Sweden 114 18.
