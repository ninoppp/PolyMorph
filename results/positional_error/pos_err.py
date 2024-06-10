import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# --- cv ---
pos_offset = 19
data = pd.read_csv(f'data/positional_error_cv.csv')
data['readout_pos'] = data['readout_pos'] + pos_offset
grad_cv = data['cv'].unique()
thresholds = data['threshold'].unique()
pos_err = data.groupby(['threshold', 'cv'])['readout_pos'].std().reset_index()
print(pos_err.head(10))

plt.figure(figsize=(10, 5))
for thresh in thresholds: 
    subset = pos_err[pos_err['threshold'] == thresh]
    plt.plot(subset['cv'], subset['readout_pos'], marker='o', label=f'Threshold = {thresh}') # TODO add error bars

plt.title('Positional Error vs Gradient Coefficient-Variation')
plt.xscale('log')
plt.xlabel('Gradient CV (D,k,p)')
plt.ylabel('Positional Error')
plt.grid(True)
plt.legend()
plt.savefig("positional_error_cv.pdf")

# readout position with min and max values as error estimates
plt.figure(figsize=(10, 5))
for thresh in [0.01]:
    subset = data[data['threshold'] == thresh]
    readout_pos = subset.groupby('cv')['readout_pos'].mean()
    min = subset.groupby('cv')['readout_pos'].min()
    max = subset.groupby('cv')['readout_pos'].max()
    plt.plot(grad_cv, readout_pos, marker='o', label=f'Threshold = {thresh}')
    plt.fill_between(grad_cv, min, max, alpha=0.2)

plt.title('Readout Position vs Gradient Coefficient-Variation')
plt.xscale('log')
plt.xlabel('Gradient CV (D,k,p)')
plt.ylabel('Readout Position (distance to source)')
plt.grid(True)
plt.legend()
plt.savefig("readout_position.pdf")

# --- width ---


