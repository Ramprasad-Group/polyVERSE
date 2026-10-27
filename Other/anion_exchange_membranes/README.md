The data has been made available along with the following paper:

William Schertzer, Sonakshi Gupta, Rampi Ramprasad, "SALSA: Semi-Autonomous Literature Summarization Assistant", arXiv:2609.22210 (2026), https://doi.org/10.48550/arXiv.2609.22210

The dataset is split into two files that are joined on **DOI** + **Sample Name**:
* **aem_data_all.csv**: anion exchange membrane (AEM) property measurements (9,535 rows, 586 samples from 121 papers)
* **aem_data_all_chem.csv**: chemical composition of each sample (600 rows, one per sample)

The **batch** column in both files identifies the provenance of each row:
* `original`: data from the previously released AEM aging dataset (see `Conductivity_anionic_aging`)
* `SALSA`: new data extracted from the literature with SALSA

Rows with a negative property value have been removed.

## aem_data_all.csv

The file comprises 13 columns:
* 1: **batch**: `original` or `SALSA` (see above)
* 2: **DOI**: DOI of the source paper
* 3: **Sample Name**: sample name as used in the source paper
* 4: **time(h)**: the amount of time (in hours) elapsed since the beginning of the alkaline stability test
* 5: **rh(%)**: measurement relative humidity
* 6: **solvent type**: base used for the alkaline stability test (NaOH or KOH)
* 7: **solvent molarity (M)**: molarity of the solvent used for the alkaline stability test
* 8: **deg test temp (C)**: temperature at which the stability test is conducted
* 9: **prop test temp (C)**: temperature at which the property is measured
* 10: **property**: property name, one of `OHCond(mS/cm)` (hydroxide conductivity), `IEC` (ion exchange capacity, meq/g), `WU(wt%)` (water uptake), `SR(%)` (swelling ratio), `TS(MPa)` (tensile strength), `EB(%)` (elongation at break), `YM(MPa)` (Young's modulus)
* 11: **value**: recorded property value
* 12: **source kind**: where in the paper the value was extracted from (table or figure); empty for `original` rows
* 13: **source name**: identifier of the specific table or figure; empty for `original` rows and for some `SALSA` rows

Values suffixed with ` Not Found` (e.g. `100 Not Found`, `None Not Found`) were not stated in the source paper; the value before the suffix is the default that was assumed. `N/A` means the condition does not apply (e.g. a measurement taken before any stability test). A small number of `value` entries are reported as in the source paper with uncertainties or units (e.g. `14.6 ± 0.92`).

## aem_data_all_chem.csv

The file comprises 24 columns:
* 1: **batch**: `original` or `SALSA` (see above)
* 2: **DOI**: DOI of the source paper
* 3: **Sample Name**: sample name, matching aem_data_all.csv
* 4-7: **smilesx**: SMILES of repeat unit x of the (co)polymer, with `*` marking the attachment points
* 8-11: **cx**: relative composition of smilesx in the copolymer
* 12-13: **additive_name_x**: name of additive x (filler, crosslinker, small molecule) present in the sample
* 14-16: **additive_smilesx**: SMILES of additive x
* 17-19: **additivecx**: amount of additive x present in the sample
* 20: **cont_add**: 1 if the sample contains an additive, 0 otherwise (not filled for all samples)
* 21: **cont_cl**: 1 if the sample is crosslinked, 0 otherwise (not filled for all samples)
* 22: **confidence**: confidence in the composition: 1 = structure and composition fully specified by the paper, 0.5 = partly inferred (e.g. crosslinker structure or ratio), 0 = cannot be fully represented (e.g. filler-dominated composites, unknown functionalization degree)
* 23: **llm notes**: SALSA's justification for including the sample and for the assigned structure and composition (`SALSA` rows only)
* 24: **my notes**: manual curation notes (currently empty)

Unused repeat-unit and additive slots are left empty or filled with the placeholder `*CC*` and a composition of 0.
