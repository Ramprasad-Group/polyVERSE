# Solvent diffusivity in polymers

This repository contains simulation code and datasets used in the paper _"Polymer Design for Solvent Separations: Integrating Simulations, Experiments, and Known Physics via Machine Learning"_.

## Abstract

This study focuses on guiding the discovery of sustainable and high-performance polymer membranes for organic binary solvent separation. We specifically focus on understanding solvent diffusivity in polymers, a key factor in quantifying solvent transport. Traditional experimental and computational methods for determining diffusivity rates are time and resource-intensive, while current machine learning (ML) models often lack accuracy outside their training domains. To address these limitations, we fuse experimental and simulated diffusivity data to train multi-task ML models that are further guided by physics-based laws, resulting in more robust predictions in unseen chemical spaces. These models outperform traditional single-task models, especially in data-limited scenarios. We then address the critical challenge of identifying optimal membranes for a given binary solvent separation problem, namely, toluene-heptane separation. Amongst 13,000 known polymers, polyvinyl chloride (PVC) is identified as an optimal polymer for this separation, consistent with literature findings, thereby validating our methodology. Moreover, we expand the screening to a broader space of synthetically accessible, virtually produced polymers, utilizing a database comprising 1 million publicly available candidates and 7 million chemically recyclable halogen-free candidates. We identify environmentally-friendly alternatives to PVC for this technologically important binary solvent separation.  This new capability is expected to effectively guide and advance membrane design for solvent separation technologies.

## Usage
### Simulation Details

The **Polymer Structure Predictor (PSP)** package is used to generate initial polymer structures and is freely available at [PSP GitHub Repository](https://github.com/Ramprasad-Group/PSP).

#### Input Structure Details

##### `data.lmps`
This file contains the input structure for the system, including:
- A polymer with the SMILES `COC(=O)C(*)C*`
- 150 atoms per chain, 12 repeat units, 32 chains (capped with `*C`)
- 12 solvents (`CC(=O)OC`)
- A total of 4996 atoms

##### `lmp.in`
A LAMMPS script used for minimization, equilibration, and production runs. Solvent diffusivity is analyzed through the mean square displacement in the production run.

#### Simulation Parameters
- **Temperature**: 308 K
- **Pressure**: 1 atm
- **Equilibration**: 21 steps, followed by 10 ns NPT and 200 ns NVT production runs

### Data Files

##### `master_diffusivity_dataset.csv`
This spreadsheet contains experimental and simulated data for solvent diffusivity in polymers as a function of temperature and concentration. The dataset includes the following columns:

1. **Polymer Canonical SMILES**
2. **Solvent Canonical SMILES**
3. **Temperature** at which solvent diffusivity is recorded
4. **Weight Fraction of Solvent** (mass of solvent / (mass of solvent + mass of polymer))
5. **Experimental / Simulated Selector**: One-hot encoded vectors to denote whether the data point is experimental or simulated
6. **Diffusivity(log10)**: Solvent diffusivity in polymer in log10 scale

##### Link to the sorption uptake dataset: [Sorption Uptake Dataset (Excel)](https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-023-40257-2/MediaObjects/41467_2023_40257_MOESM3_ESM.xlsx)

## Citation
If you use this repository in your work, please consider citing the original paper:
- _"Polymer Design for Solvent Separations: Integrating Simulations, Experiments, and Known Physics via Machine Learning"_

