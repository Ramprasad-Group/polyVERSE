# Data
## Overview
This repository contains two JSON files: (1) `polyimides.json` and (2) `romp_polymers.json`.

File 1 contains five keys, named: (a) "zinc_id1" (b) "zinc_id2" (c) "polymer" (d) "reactant1" and (e) "reactant2". Key A contains
the ZINC or ChEMBL ID corresponding to a dianhydride. Key B contains the ZINC or ChEMBL ID corresponding to a diamine. Key C contains the
SMILES string of hypothetical polyimides. For example, the SMILES string in key C, value 0 corresponds to the hypothetical polymimide that would be
made by reacting the molecule in key A, value 0 with the molecule in key B, value 0. Key D contains the SMILES string of the molecule in Key A.
Key E contains the SMILES string of the molecule in B.

File 2 contains three keys, named: (a) "zinc_ids" (b) "polymer" and (c) "reactants". Key A contains
the ZINC or ChEMBL ID corresponding to a cyclic olefin. Key B contains the
SMILES string of hypothetical polymers. For example, the SMILES string in key B, value 0 corresponds to the hypothetical polymer that would be
made by performing Ring-opening metathesis polymerization (ROMP) on the molecule in key A, value 0. Key C contains the SMILES string of the
molecule in Key A.

## CSV files
The corresponding CSV files are: (1) `polyimides.csv` and (2) `romp_polymers.csv`.
