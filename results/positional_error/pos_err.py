import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
import pandas as pd
from math import sqrt

n = 100
def SE(stddev):
    return stddev / sqrt(2*(n-1)) # standard error of the standard deviation

# --- cv ---
pos_offset = 19
data = pd.read_csv(f'data/positional_error_cv.csv')
data['readout_pos'] = data['readout_pos'] + pos_offset
grad_cv = data['cv'].unique()
thresholds = data['threshold'].unique()
pos_err = data.groupby(['threshold', 'cv'])['readout_pos'].std().reset_index()

print(pos_err.head(10))
print(len(thresholds))

viridis = colormaps.get_cmap('viridis')  # Get the viridis colormap
colors = viridis(np.linspace(0, 1, len(thresholds)))  # Get colors from viridis
plt.rcParams['text.usetex'] = True

plt.figure(figsize=(10, 5))
for i, thresh in enumerate(thresholds): 
    subset = pos_err[pos_err['threshold'] == thresh]
    plt.plot(subset['cv'], subset['readout_pos'], marker='o', label=f'Threshold c0 = {thresh}', color=colors[i])
    plt.errorbar(subset['cv'], subset['readout_pos'], yerr=SE(subset['readout_pos']), color=colors[i])

plt.title('Positional error vs gradient variability CV')
plt.xscale('log')
plt.xlabel('Gradient CV (D,k,p)')
plt.ylabel('Positional error')
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
    plt.plot(grad_cv, readout_pos, marker='o', label=r'Threshold C_{\theta} = '+f'{thresh}')
    plt.fill_between(grad_cv, min, max, alpha=0.2)

plt.title('Readout Position vs Gradient Coefficient-Variation')
plt.xscale('log')
plt.xlabel('Gradient CV (D,k,p)')
plt.ylabel('Readout Position (distance to source)')
plt.grid(True)
plt.legend()
plt.savefig("readout_position.pdf")

# --- width ---
data = pd.read_csv(f'data/positional_error_width.csv')
data['readout_pos'] = data['readout_pos'] + pos_offset
widths = data['width'].unique()
thresholds = data['threshold'].unique()
pos_err = data.groupby(['threshold', 'width'])['readout_pos'].std().reset_index()

plt.figure(figsize=(10, 5))
for i, thresh in enumerate(thresholds): 
    subset = pos_err[pos_err['threshold'] == thresh]
    plt.plot(subset['width'], subset['readout_pos'], marker='o', label=f'Threshold c0 = {thresh}', color=colors[i])
    plt.errorbar(subset['width'], subset['readout_pos'], yerr=SE(subset['readout_pos']), color=colors[i])

plt.title('Positional error vs domain width')
plt.xlabel('Domain width (units ~ avg. cell radius)')
plt.ylabel('Positional error')
plt.grid(True)
plt.legend()
plt.savefig("positional_error_width.pdf")