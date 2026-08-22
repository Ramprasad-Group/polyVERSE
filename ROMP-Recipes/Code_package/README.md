# Representative ROMP code for peer review

This package provides representative code for three modules described in the manuscript:

1. a descriptor-based Gaussian process classifier for ROMP feasibility;
2. a structure-based heuristic explainer for potential ROMP incompatibility; and
3. a ring-aware retriever for related ROMP reaction-condition precedents.

The package is intended for confidential editorial and peer-review evaluation. It provides executable representative code, a small example dataset, installation instructions, and demo commands. The example dataset is included only to demonstrate code execution and data format, not to reproduce the full manuscript performance metrics.

## Package contents

```text
Code_package/
├── README.md
├── requirements.txt
├── gpc_training/
│   ├── DFT_GPC_MODEL.py
│   └── dft_feature_columns.json
├── heuristic_explainer/
│   └── romp_heuristic_explainer.py
├── ring_aware_retriever/
│   └── ring_aware_reaction_condition_retriever.py
└── example_data/
    └── ROMP_DFT_GPC_representative_subset_for_review.xlsx
```

## System requirements

The scripts were tested with Python 3.10 on Windows 11 and Linux. Required Python packages are listed in `requirements.txt`.

No non-standard hardware is required for the demo scripts. The provided examples are expected to run on a normal desktop or laptop computer.

Typical installation time on a normal desktop computer is approximately 5–10 minutes, depending on internet speed and package availability.

## Installation

From the `Code_package/` directory:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
```

On Windows:

```bash
python -m venv .venv
.venv\Scripts\activate
python -m pip install -r requirements.txt
```

## 1. DFT-GPC training

The GPC training script fits a binary Gaussian process classifier using the 55 model-compatible feature columns listed in `gpc_training/dft_feature_columns.json`.

Example interface:

```bash
python gpc_training/DFT_GPC_MODEL.py \
  --input example_data/ROMP_DFT_GPC_representative_subset_for_review.xlsx \
  --output_dir gpc_model_output \
  --label_col target_property
```

The script saves the fitted model package, feature-column definition, training-set medians used for missing-value imputation, run configuration, and an apparent in-sample training summary.

A small representative labeled DFT descriptor subset is included under `example_data/` for testing the script interface. This subset is intended only to demonstrate code execution and data format, not to reproduce the full manuscript performance metrics.

Expected output files are saved in the specified output directory and include:

```text
model_dft_gpc.joblib
dft_feature_columns.json
training_medians.json
run_config.json
production_training_summary.json
```

Expected terminal output includes a short training summary, such as the number of rows used, label counts, feature count, missing-descriptor summary, and saved output paths.

The DFT-GPC demo with the representative subset is expected to finish within a few minutes on a normal desktop computer.

To run the script on user-provided data, replace the `--input` path with a CSV or Excel file containing the required DFT descriptor columns and a binary label column. If the label column has a different name, update the `--label_col` argument accordingly.

## 2. Heuristic explainer

The heuristic explainer identifies structure-centered cautionary features and returns qualitative rationales and possible structure-modification suggestions. It is an interpretation aid rather than an experimental verdict.

Example:

```bash
python heuristic_explainer/romp_heuristic_explainer.py \
  --smiles "C1=CCCCCCC1"
```

The script also supports user-provided CSV or Excel files containing a SMILES column.

To run the heuristic explainer on user data, provide a SMILES string directly using `--smiles`, or provide a CSV/Excel file containing a SMILES column according to the script interface.

The expected output is a structure-centered interpretation report containing detected structural features, triggered heuristic warnings where applicable, qualitative rationales, and possible structure-modification suggestions.

The demo is expected to finish within seconds on a normal desktop computer.

## 3. Ring-aware reaction-condition retriever

The retriever uses a weighted score combining whole-molecule, reactive-alkene-local, alkene-containing-ring, and ring-metadata similarities. The default weights are 0.15, 0.30, 0.35, and 0.20, respectively, with a minimum weighted similarity of 0.6.

Example interface:

```bash
python ring_aware_retriever/ring_aware_reaction_condition_retriever.py \
  --input PATH_TO_REACTION_DATABASE.csv \
  --query-smiles "C1=CCCCCCC1" \
  --topk 3 \
  --mechanism ROMP \
  --min-similarity 0.6 \
  --require-ring-alkene
```

The input reaction table must contain a monomer SMILES column, such as `monomer_smiles`, `Monomer_smiles`, `Monomer_SMILES`, or `SMILES`. Additional reaction-condition fields are returned when present in the input table.

To run the retriever on user data, replace `PATH_TO_REACTION_DATABASE.csv` with a CSV file containing monomer SMILES and reaction-condition fields, and replace `--query-smiles` with the target monomer SMILES.

The expected output is a ranked list of structurally related ROMP precedent reactions, including similarity scores and available reaction-condition information.

Runtime depends on the size of the reaction database. For small demonstration tables, the retriever is expected to finish within seconds to a few minutes on a normal desktop computer.

## Reproducibility note

This reviewer code package is intended to demonstrate the implementation, data format, and execution of representative modules described in the manuscript. It is not intended to reproduce all quantitative results reported in the manuscript. 