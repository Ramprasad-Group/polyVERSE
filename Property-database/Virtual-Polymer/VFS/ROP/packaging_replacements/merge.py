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

# method to calculate new fitness based on differing targets
props = {
    'Youngs Modulus (GPa)': 2,
    'Glass Transition Temperature (K)': 373,
    'Thermal Decomposition Temperature (K)': 473,
    'Specific Heat Capacity (J {gK}^{-1})': 1.24,
    'Tensile Strength at break (MPa)': 39,
    'Enthalpy of Polymerization (kJ/mol)': 10
}

df['new_fitness'] = 1
for prop, target in props.items():
    fitkey = f"{prop}_fitness"
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
    df['new_fitness'] *= df[fitkey]
