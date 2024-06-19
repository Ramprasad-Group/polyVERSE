# Gas transport simulations

This repository contains the molecular dynamics (MD) and Monte Carlo (MC) code and data for the paper "Gas permeability, diffusivity, and solubility in polymers: Simulation-experiment data fusion and multi-task machine learning". It contains LAMMPS scripts for simulations of gas diffusion and solubility in amorphous polymers. In addition, the repo also contains the gas transport data generated from our simulations as well as those collected from literatures.



Abstract:

Machine learning (ML) models for predicting gas permeability through polymers have traditionally relied on experimental data. While these models exhibit robustness within familiar chemical domains, reliability wanes when applied to new spaces. To address this challenge, we present a multi-tiered multi-task learning framework empowered with advanced machine-crafted polymer fingerprinting algorithms and data fusion techniques. This framework combines scarce "high-fidelity" experimental data with abundant diverse "low-fidelity" simulation or synthetic data, resulting in predictive models that display a high level of generalizability across novel chemical spaces. Additionally, this multi-task scheme capitalizes on known physics and interrelated properties, such as gas diffusivity and solubility, both of which are closely tied to permeability. By amalgamating high-throughput generated simulation data with available experimental data for gas permeability, diffusivity, and solubility for various gases, we construct multi-task deep learning models. These models can simultaneously predict all three properties for all gases under consideration. With markedly enhanced predictive accuracy, particularly compared to traditional models reliant solely on experimental data for a singular property. This strategy underscores the potential of coupling high-throughput classical simulations with data fusion methodologies to yield state-of-the-art property predictors, especially when experimental data for targeted properties is scarce.

# Usage 
The simulation files in this repository can be run with LAMMPS. For installation and usage visit https://www.lammps.org/. The Polymer Structure Predictor (PSP) package to create simulation polymer structures is available free of charge at https://github.com/Ramprasad-Group/PSP. The model was trained using polyGNN and the code and usage details can be found at https://github.com/Ramprasad-Group/polygnn. 

**Gas diffusivity simulation files**:  \
[amorphous_polymer_diffusivity.lmps](amorphous_polymer_diffusivity.lmps) - a structure file for polyethylene, comprised of 27 polymer chains and 27 methane molecules. There are approximately ~150 atoms per polymer chain. \
\
[gas_diffusivity.in](gas_diffusivity.in) - A Lammps script to run a minimization, temperature ramping, equilibration, and production run for the polymer system. Gas diffusivity can be found via analysis of the production run. \
\
**Gas solubility simulation files**: \
[amorphous_polymer_solubility.lmps](amorphous_polymer_solubility.lmps) - a structure file for polyethylene, comprised of 27 polymer chains. There are approximately ~150 atoms per polymer chain. \
\
[widom_insertions.in](widom_insertions.in) - A lammps script to run 25,000 widom insertions for methane. Gas solubility can be found via analysis of the CH4xcp.txt output file.\
\
[CH4_TraPPE_FF.txt](CH4_TraPPE_FF.txt) - TraPPE forcefield parameters for methane that the widom_insertions.in script calls upon. \


**Datasets**: The datasets here were used to train our gas transport predictor model. The datasets provided are provided in two formats: long and wide. The long format organizes data with one observation per row, while the wide format arranges variables in separate columns with each row representing a unique case or observation. 

[Gas_permeability_solubility_diffusivity_long.csv](Gas_permeability_solubility_diffusivity_long.csv) 
The dataset comprises 3 columns:
1. **property**: Format follows the form; property name, data type, and gas type, separated by "_". Where permeability = "p", diffusivity = "d", solubility = "s", experimental data = "exp", and simulation data = "sim".
2. **value**: Property value in log(10) form. Units for each property are: Permeability- "(cc\*cm)/(cm^2\*s\*cmHg)", solubility- "cc/(cc\*cmHg)", and diffusivity- "cm^2/s".
3. **smiles_string**: SMILES representation of polymers

[Gas_permeability_solubility_diffusivity_wide.csv](Gas_permeability_solubility_diffusivity_wide.csv)
The dataset comprises 31 columns:
1. **smiles_string**: SMILES representation of polymers.    
<!-- -->
2-31. **property**: Format follows the form; property name, data type, and gas type, separated by "_". Where permeability = "p", diffusivity = "d", solubility = "s", experimental data = "exp", and simulation data = "sim". A full list of the properties are shown below. 

Property list:
'd_exp_CH4', 'd_exp_CO2', 'd_exp_H2', 'd_exp_He', 'd_exp_N2', 'd_exp_O2', 'd_sim_CH4', 'd_sim_CO2', 'd_sim_N2', 'd_sim_O2', 'p_exp_CH4', 'p_exp_CO2', 'p_exp_H2', 'p_exp_He', 'p_exp_N2', 'p_exp_O2', 'p_sim_CH4', 'p_sim_CO2', 'p_sim_N2', 'p_sim_O2', 's_exp_CH4', 's_exp_CO2', 's_exp_H2', 's_exp_He', 's_exp_N2', 's_exp_O2', 's_sim_CH4', 's_sim_CO2', 's_sim_N2', 's_sim_O2'

