# polyverse
This repo contains the code used to generate the polymers studied in the companion paper, [AI-assisted discovery of high-temperature dielectrics for energy storage](https://www.nature.com/articles/s41467-024-50413-x). At a high level, the code uses known molecules and know reaction rules (for ROMP and imide polycondensation) to generate plausible polymer repeat units.
## Installation
1. Install [poetry](https://python-poetry.org/) on your machine.
2. If Python3.9 is installed on your machine skip to step 3, if not you will need to install it.
3. Clone this repo on your machine.
4. Open a terminal at the root directory of this repository.
5. Run `poetry env use /path/to/python3.9/executable`.
6. Run `poetry install`.
## Examples
Examples of how to use the codebase are provided in the file named `tests/test_polyverse.py`. Basically, the user needs to first create a data frame of molecules. The data frame should have one column named "id-smiles". Each entry in this column should be a `tuple`, with the first element being a unique identifier for the molecule (I use the ZINC ID) and the second element being the SMILES string for the molecule. Then, the user can input these molecules into an instance of the `ROMP` or `Polyimide` classes (defined in `polyverse/polymerize.py`) to generate polymers.
## Testing
There are two tests in this repo. One is to test the `ROMP` class, the other to test the `Polyimide` class. To execute the tests, run `poetry run pytest`.
## Data
- `data/romp_polymers.json` contains the ROMP polymers designed in the companion paper.
- `data/polyimides.json` contains the polyimides designed in the companion paper.
- The data base of monomers used in this work may be found at [doi.org/10.5281/zenodo.12535176](https://doi.org/10.5281/zenodo.12535176).
## Reproducibility
The version of this codebase that was used in the companion paper is v0.1.0.
