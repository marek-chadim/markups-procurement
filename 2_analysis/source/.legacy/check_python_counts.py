"""Diagnostic: count observations in the same filter sequence as Stata."""
import pandas as pd
import numpy as np

df = pd.read_stata('../input/data_rebuilt.dta')
df = df.sort_values(['id', 'year']).reset_index(drop=True)
print(f'data_rebuilt.dta: {len(df)} obs')

df['Lcogs'] = df.groupby('id')['cogs'].shift(1)
df['Lk'] = df.groupby('id')['k'].shift(1)
df['L2cogs'] = df.groupby('id')['cogs'].shift(2)
df['Lpp'] = df.groupby('id')['pp_dummy'].shift(1)

m1 = df[['k','cogs','go','Lk','Lcogs','Lpp']].notna().all(axis=1)
print(f'  + first lags + Lpp: {m1.sum()}')
m2 = df[['k','cogs','go','Lk','Lcogs','Lpp','L2cogs']].notna().all(axis=1)
print(f'  + L2cogs: {m2.sum()}')

print('\nBy NACE2 (first-lag sample):')
for nace in sorted(df['nace2'].unique()):
    print(f'  NACE {int(nace)}: CD={m1[df["nace2"]==nace].sum()}  TL={m2[df["nace2"]==nace].sum()}')
