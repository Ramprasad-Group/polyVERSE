import pandas as pd

df = pd.read_csv('combined_dft_exp_with_features.csv')

# 1. Identify the compositions which are in source 'exp'
compositions_in_exp = df[df['source'] == 'exp']['composition']

# 2. Identify compositions which are only in source 'exp' and type 'overlapped'
compositions_in_exp_overlapped = df[(df['source'] == 'exp') & (df['type'] == 'overlapped')]['composition']

# 3. Identify compositions which are only in source 'exp' and type 'only_exp'
compositions_in_exp_only_exp = df[(df['source'] == 'exp') & (df['type'] == 'only_exp')]['composition']

# 4. Identify the compositions in source except 'exp'
compositions_except_exp = df[df['source'] != 'exp']['composition'].unique()

# 5. Identify the compositions in source except 'exp' and type 'only_dft'
compositions_except_exp_only_dft = df[(df['source'] != 'exp') & (df['type'] == 'only_dft')]['composition']

# 6. Identify the compositions in source except 'exp' and in type 'overlapped'
compositions_except_exp_overlapped = df[(df['source'] != 'exp') & (df['type'] == 'overlapped')]['composition']

print("Compositions in source 'exp':", len(compositions_in_exp))
print("Compositions in source 'exp' and type 'overlapped':", len(compositions_in_exp_overlapped))
print("Compositions in source 'exp' and type 'only_exp':", len(compositions_in_exp_only_exp))
print("Compositions in source except 'exp':", len(compositions_except_exp))
print("Compositions in source except 'exp' and type 'only_dft':", len(compositions_except_exp_only_dft))
print("Compositions in source except 'exp' and type 'overlapped':", len(compositions_except_exp_overlapped))
