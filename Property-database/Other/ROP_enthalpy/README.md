The data has been made available along with the following paper:

Aubrey Toland, Huan Tran, Lihua Chen, Yinghao Li, Chao Zhang, Will Gutekunst, and Rampi Ramprasad, [Accelerated Scheme to Predict Ring-Opening Polymerization Enthalpy: Simulation-Experimental Data Fusion and Multitask Machine Learning](https://doi.org/10.1021/acs.jpca.3c05870), J. Phys. Chem. A 2023, 127, 50, 10709–10716 (2023).

The data set is available in two different formats, a long format where every data point relating to an enthalpy value is explicitly in a different row, and a wide format, where all enthalpy data (experimental and produced by DFT) for a polymer are grouped in one row.

The long format dataset comprises of 6 columns:

1. **ID**: Serial Number (Arbitrarily Assigned for Easier Reference)
2. **smiles_polymer**: SMILES Representation of Polymer
3. **smiles_monomer**: SMILES Representation of Monomer
4. **1/length_non_normalized**: Inverse of the Loop Size for a given Calculation (0 for Experimental Data)
5. **source**: Indicator for if the Data Comes from Experiment or DFT Calculation (expt for Experiment, or dft for DFT calculation)
6. **roe_kj/mol**: Enthalpy Value in kJ/mol


The wide format dataset comprises of 10 columns:

1. **ID**: Serial Number (Arbitrarily Assigned for Easier Reference)
2. **smiles_polymer**: SMILES Representation of Polymer
3. **smiles_monomer**: SMILES Representation of Monomer
4.  **monomer_image**: 2D Image of the Monomer Chemical Structure
5.  **delta_H_exp_KJ/mol**: Experimental Enthalpy Value in kJ/mol
6.  **delta_H_1/l-3_KJ/mol**: DFT Calculated Enthalpy Value in kJ/mol for a polymer system of loop size 3
7.  **delta_H_1/l-4_KJ/mol**: DFT Calculated Enthalpy Value in kJ/mol for a polymer system of loop size 4
8.  **delta_H_1/l-5_KJ/mol**: DFT Calculated Enthalpy Value in kJ/mol for a polymer system of loop size 5
9. **delta_H_1/l-6_KJ/mol**: DFT Calculated Enthalpy Value in kJ/mol for a polymer system of loop size 6
10. **Reference**: DOI of Data Point (If the Polymer Has Experimental Data)
