# Gas transport simulations

This repository contains the most up-to-date gas transport dataset, combining experimental and simulation data. It includes data from the following publications: “Gas permeability, diffusivity, and solubility in polymers: Simulation-experiment data fusion and multi-task machine learning” and “polyBERT: a chemical language model to enable fully machine-driven ultrafast polymer informatics”, as well as newly generated or collected data. The dataset features a source column indicating the origin of each entry, with values corresponding to the sources: "Gas npj (2024)", "polyBERT", and "polyVERSE_2025_06_23".

\
[master_transport_2025_08_13.csv](master_transport_2025_08_13.csv) 
\
The dataset comprises 9 columns:
1. **property**: Format follows the form; property name, data type, and gas type, separated by "_". Where permeability = "p", diffusivity = "d", solubility = "s", experimental data = "exp", and simulation data = "sim".
2. **smiles_string**: SMILES representation of polymers
3. **p_csmiles**: canonicalized version of the smiles string, using "canonicalize_psmiles" package
4. **value**: Property value in log(10) form. Units for each property are: Permeability- "(cc\*cm)/(cm^2\*s\*cmHg)", solubility- "cc/(cc\*cmHg)", and diffusivity- "cm^2/s".
5. **unscaled_value**: Unscaled value
6. **Activity**: Activity at which measurement was taken
7. **Temp**: Temperature at which measurement was taken
8. **Solvent_SMILES**: The type of gas 
9. **source**: Source of the data

\
Property list:
'd_exp_CH4', 'd_exp_CO2', 'd_exp_H2', 'd_exp_He', 'd_exp_N2', 'd_exp_O2', 'd_exp_H2O', 'd_sim_CH4', 'd_sim_CO2', 'd_sim_N2', 'd_sim_O2', 'p_exp_CH4', 'p_exp_CO2', 'p_exp_H2', 'p_exp_He', 'p_exp_N2', 'p_exp_O2', 'p_exp_H2O, 'p_sim_CH4', 'p_sim_CO2', 'p_sim_N2', 'p_sim_O2', 's_exp_CH4', 's_exp_CO2', 's_exp_H2', 's_exp_He', 's_exp_N2', 's_exp_O2', 's_exp_H2O, 's_sim_CH4', 's_sim_CO2', 's_sim_N2', 's_sim_O2',  's_sim_H2O'

