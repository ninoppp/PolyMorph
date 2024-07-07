import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps
import pandas as pd
from math import sqrt
from scipy import stats

# constants (simulation specific)
pos_offset = 19 # source position (data is in absolute coordinates)
n = 100 # number of samples
cell_diameter = 1.31 # avg. cell diameter
mu_lambda = 4.31 # lambda / cell diameter = 5.65 / 1.31

def stderr(stddev):
    return stddev / sqrt(2*(n-1)) # standard error of the standard deviation

# colormap & latex
viridis = colormaps.get_cmap('viridis')  # Get the viridis colormap
plt.rcParams['text.usetex'] = True
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],  # This line may need to be adjusted
})

# Continue as previously demonstrated


# --- cv ---
# positional error
data = pd.read_csv('data/positional_error_cv.csv')
data['readout_pos'] = (data['readout_pos'] + pos_offset) / cell_diameter
grad_cv = data['cv'].unique()
thresholds = data['threshold'].unique()
thresh_labels = [r'$c_\theta=$' + f'{thresh:.0e}' for thresh in thresholds]
pos_err = data.groupby(['threshold', 'cv'])['readout_pos'].std().reset_index()
colors = viridis(np.linspace(0, 1, len(thresholds)+1))

plt.figure(figsize=(8, 5))
for i, thresh in enumerate(thresholds): 
    mu_x = data[data['threshold'] == thresh]['readout_pos'].mean()
    subset = pos_err[pos_err['threshold'] == thresh]
    plt.plot(subset['cv'], subset['readout_pos'], marker='o', label=thresh_labels[i], color=colors[i])
    plt.errorbar(subset['cv'], subset['readout_pos'], yerr=stderr(subset['readout_pos']), color=colors[i])

#plt.title('Positional error vs gradient variability')
plt.xscale('log')
plt.xlabel('Gradient variability CV (D,k,p)')
plt.ylabel('Positional error / avg. cell diameter')
plt.grid(True)
plt.legend()
plt.savefig("positional_error_cv.pdf")

# readout position
plt.figure(figsize=(8, 5))
for i, thresh in enumerate(thresholds):
    subset = data[data['threshold'] == thresh]
    mu_x = subset['readout_pos'].mean()
    readout_pos = subset.groupby('cv')['readout_pos'].mean()
    min = subset.groupby('cv')['readout_pos'].min()
    max = subset.groupby('cv')['readout_pos'].max()
    plt.plot(grad_cv, readout_pos, marker='o', label=thresh_labels[i], color=colors[i])
    plt.fill_between(grad_cv, min, max, alpha=0.2, color=colors[i]) # min and max values as error estimates

#plt.title('Readout position vs gradient variability')
plt.xscale('log')
plt.xlabel('Gradient variability CV (D,k,p)')
plt.ylabel('Readout position / avg. cell diameter')
plt.grid(True)
plt.legend()
plt.savefig("readout_position.pdf")

    
# --- width ---
# positional error
data = pd.read_csv('data/positional_error_width.csv')
data['readout_pos'] = (data['readout_pos'] + pos_offset) / cell_diameter
widths = data['width'].unique() / cell_diameter
thresholds = data['threshold'].unique()
thresh_labels = [r'$c_\theta=$' + f'{thresh:.0e}' for thresh in thresholds]
pos_err = data.groupby(['threshold', 'width'])['readout_pos'].std().reset_index()
colors = viridis(np.linspace(0, 1, len(thresholds+2)))

plt.figure(figsize=(8, 5))
for i, thresh in enumerate(thresholds): 
    if thresh == 0.0001:
        continue # skip this because too close to opposing boundary
    mu_x = data[data['threshold'] == thresh]['readout_pos'].mean()
    print(mu_x)
    subset = pos_err[pos_err['threshold'] == thresh]
    plt.plot(widths, subset['readout_pos'], marker='o', label=thresh_labels[i] +  r',$\hspace{0.3cm}$' + r'$\mu_{x_\theta}$' + f'={mu_x/mu_lambda:.1f}'  r'$\mu_\lambda$', color=colors[i])
    plt.errorbar(widths, subset['readout_pos'], yerr=stderr(subset['readout_pos']), color=colors[i])

# plot reference line 1/sqrt(width)
x = np.linspace(widths.min(), widths.max(), 100)
y = 1.0 / np.sqrt(x)
plt.plot(x, y, label=r'$1 / \sqrt{width}$', color='gray', linestyle='--')

#plt.title('Positional error vs domain width')
plt.xlabel('Domain width / avg. cell diameter')
plt.ylabel('Positional error / avg. cell diameter')
plt.grid(True)
plt.legend()
plt.savefig("positional_error_width.pdf")

# print average readout position for threshold 0.1
#subset = data[data['threshold'] == 0.0001]
#readout_pos = subset.groupby('width')['readout_pos'].mean()
# mean over all widths
#print(readout_pos) #6.1, 16, 
#print(min) 
#print(max)
#print(readout_pos.mean())