# Ring-Opening Polymerization (ROP) Dataset

This repository contains a dataset of 1,087,564 ring-opening polymerization (ROP) polymer structures generated using VFS, polyBART, and POLYT5.

Due to GitHub file size limitations, the complete dataset is distributed across three Parquet files which can be concatenated row-wise. To merge these into a single file, you can use the following code:
```Python
from pathlib import Path
import pandas as pd

df1 = pd.read_parquet("data_split1.parquet")
df2 = pd.read_parquet("data_split2.parquet")
df3 = pd.read_parquet("data_split3.parquet")
df = pd.concat([df1, df2, df3], ignore_index=True)
```

- Each file has the following columns:
    - `smiles`: the polymer SMILES representation
    - `canonical_smiles`: the canonicalized polymer SMILES
    - `source`: the method by which the structure was generated or curated
    - `monomer_smiles`: the corresponding cyclic ROP monomer SMILES
    - `rop_type`: the monomer class assigned by the chemical heuristics filter
    - `mechanism_profile`: the polymerization mechanism profile predicted by the chemical heuristics filter
    - `commercial_availability`: a Boolean indicator of whether the corresponding monomer is commercially available
    - `SAscore`: the RDKit synthetic accessibility score of the hydrogen-capped polymer SMILES
