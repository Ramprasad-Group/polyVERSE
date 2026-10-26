The data has been made available along with the following paper (to be published):

Schertzer, William; Al Otmi, Mohammed; Sampath, Janani; Lively, Ryan; Ramprasad, Rampid, "AI-Assisted Physics-Informed Predictions of Degradation Behavior of Polymeric Anion Exchange Membranes", The Journal of Physical Chemistry, Part B: Biophysics, Biomaterials, Liquids, Soft Matter, Section: B: Soft Matter, Fluid Interfaces, Colloids, Polymers, and Glassy Materials 

The dataset is split into two files that are joined on **DOI** + **Sample Name**:
* **aem_aging.csv**: property measurements over the course of alkaline stability tests (4,232 rows, 272 samples from 54 papers)
* **aem_aging_chem.csv**: chemical composition of each sample (304 rows, one per sample)

## aem_aging.csv

The file comprises 12 columns:
* 1: **DOI**: DOI of the source paper
* 2: **Sample Name**: sample name as used in the source paper
* 3: **time(h)**: the amount of time (in hours) elapsed since the beginning of the alkaline stability test
* 4: **rh(%)**: measurement relative humidity
* 5: **solvent type**: base used for the alkaline stability test (NaOH or KOH)
* 6: **solvent molarity (M)**: molarity of the solvent used for the alkaline stability test
* 7: **deg test temp (C)**: temperature at which the stability test is conducted
* 8: **prop test temp (C)**: temperature at which the property is measured
* 9: **property**: property name, one of `OHCond(mS/cm)` (hydroxide conductivity), `IEC` (ion exchange capacity, meq/g), `WU(wt%)` (water uptake), `SR(%)` (swelling ratio), `TS(MPa)` (tensile strength), `EB(%)` (elongation at break), `YM(MPa)` (Young's modulus)
* 10: **value**: recorded property value
* 11: **source kind**: where in the paper the value was extracted from (table or figure); currently empty
* 12: **source name**: identifier of the specific table or figure; currently empty

Values suffixed with ` Not Found` (e.g. `100 Not Found`, `None Not Found`) were not stated in the source paper; the value before the suffix is the default that was assumed. `N/A` means the condition does not apply (e.g. a measurement taken before any stability test).

## aem_aging_chem.csv

The file comprises 21 columns:
* 1: **DOI**: DOI of the source paper
* 2: **Sample Name**: sample name, matching aem_aging.csv
* 3-6: **smilesx**: SMILES of repeat unit x of the (co)polymer, with `*` marking the attachment points
* 7-10: **cx**: relative composition (%) of smilesx in the copolymer; c1-c4 sum to 100
* 11-12: **additive_name_x**: name of additive x (filler, crosslinker, small molecule) present in the sample
* 13-15: **additive_smilesx**: SMILES of additive x
* 16-18: **additivecx**: amount of additive x present in the sample
* 19: **cont_add**: 1 if the sample contains an additive, 0 otherwise (filled for only some samples)
* 20: **cont_cl**: 1 if the sample is crosslinked, 0 otherwise (currently empty)
* 21: **confidence**: confidence in the composition: 1 = structure and composition fully specified by the paper, 0.5 = partly inferred (e.g. crosslinker structure or ratio, random copolymer placement), 0 = cannot be fully represented (e.g. filler-dominated composites, unknown functionalization)

Unused repeat-unit and additive slots are left empty or filled with the placeholder `*CC*` and a composition of 0.
