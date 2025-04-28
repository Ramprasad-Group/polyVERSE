# VFS ROP Results Data

This repository holds data associated with the paper "An Informatics Framework for the Design of Sustainable, Chemically Recyclable, Synthetically-Accessible and Durable Polymers"

## Accessing Polymer Designs
The hypothetical ROP polymers generated using VFS are separated into 20 spreadsheets in the `polymers` folder, labeled `split_##.csv` where i ranges from 00 to 19. They we separated due to size constraints. To merge these into a single file, you can use the following code (found in merge.py)
```Python
from pathlib import Path
import pandas as pd

folder = Path('polymers')

dfs = []
for i in range(20):
    fileid = str(i)
    if i < 10:
        fileid = '0' + fileid
    dfs.append(pd.read_csv(folder / f'split_{fileid}.csv'))

df = pd.concat(dfs, ignore_index=True)
```

- Each spreadsheet has the following columns:
    - `polymer_smiles`: polymer SMILES
    - `ROP of`: Which ring opening reaction is run on the monomer smiles
    - `monomer_smiles`: SMILES of the monomer that creates the polymer smiles
    - The following predicted properties
        - `Enthalpy of Polymerization (kJ/mol)`
        - `Glass Transition Temperature (K)`
        - `Specific Heat Capacity (J {gK}^{-1})`
        - `Tensile Strength at break (MPa)`
        - `Thermal Decomposition Temperature (K)`
        - `Youngs Modulus (GPa)`
    - `fitness`: The fitness value found using the following algorithm (also found in merge.py)
        ```Python
        props = {
            'Youngs Modulus (GPa)': 2,
            'Glass Transition Temperature (K)': 373,
            'Thermal Decomposition Temperature (K)': 473,
            'Specific Heat Capacity (J {gK}^{-1})': 1.24,
            'Tensile Strength at break (MPa)': 39,
            'Enthalpy of Polymerization (kJ/mol)': 10
        }

        df['fitness'] = 1
        for prop, target in props.items():
            fitkey = f"{prop}_fitness"
            # note, we would then drop the new fitness columns 
            # from the dataframe
            df[fitkey] = df[prop]
            # transform all enthalpies less than -20 so the distance is the same as if 
            # they were greater than -10 e.g., transform -40 to 10 since 10 is 
            # 20 greater than -10 and -40 is 20 less than -20
            if 'Enthalpy' in fitkey:
                df.loc[df[fitkey] < -20, fitkey] = (-20 - df[fitkey] + -10)
                df.loc[df[fitkey] < -10, fitkey] = -10
                # multiply by negative one so we can find all less than 10 instead of
                # greater than -10
                df[fitkey] *= -1
            min_ = df[fitkey].min()
            df[fitkey] = (df[fitkey] - min_) / (target - min_)
            df.loc[df[fitkey] > 1, fitkey] = 1
            df['fitness'] *= df[fitkey]
        ```

## Reaction Procedures
The `reaction_procedure.json` file contains all reaction procedures used to generate these polymers. The json file contains a list of key, value pairs that can be loaded into a Python dictionary. Each item contains a "class" key indicating which ROP class is supposed to be run, as well as a "reaction_procedure" key which points to the reaction procedure used to generate the polymer.

`react.ipynb` has some example code that explores how to use these reaction procedures to create new polymers, and models.py has some pydantic `models.py` that explain further different keywords used in the json file.


## Database
`schema.sql` holds a SQL representation of the database we used for this project. It can be visualized at [dbdiagram.io](https://dbdiagram.io/) by copy-pasting the schema text into the user interface.

## Prediction Models
All predicted properties except for Enthalpy are from the model in the paper "C. Kuenneth and R. Ramprasad, “polyBERT: a chemical language model to enable fully machine-driven ultrafast polymer informatics,” Nat Commun, vol. 14, no. 1, Art. no. 1, Jul. 2023, doi: 10.1038/s41467-023-39868-6.

The Enthalpy predictions came from the model in the paper "A. Toland et al., “Accelerated Scheme to Predict Ring-Opening Polymerization Enthalpy: Simulation-Experimental Data Fusion and Multitask Machine Learning,” J. Phys. Chem. A, Dec. 2023, doi: 10.1021/acs.jpca.3c05870.

## Installation
The script used in here can be installed using [poetry v1.8.3](https://python-poetry.org/) and Python v3.9

After installing poetry and configuring env to use Python 3.9, run `poetry install`
