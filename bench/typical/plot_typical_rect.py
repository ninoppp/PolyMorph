import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np

# area in rect tissue = 1.348
# gridpoints per cell = 1.348 / dx^2
cell_area = 1.348 # dense tissue simulation
viridis = colormaps.get_cmap('viridis')
colors = viridis(np.linspace(0, 1, 4))

df = pd.read_csv('benchmark_typical_rect.csv')
df['ratio'] = df['time_all']/df['time_ensemble']
min = df.groupby('dx')['ratio'].min()
max = df.groupby('dx')['ratio'].max()
stddev = df.groupby('dx').std()
df = df.groupby('dx').mean()

plt.figure()
plt.plot(cell_area / (df.index**2), df['ratio'], marker='o', color=colors[0])
plt.errorbar(cell_area / df.index**2, df['ratio'], yerr=stddev['ratio'], color=colors[0])
#plt.fill_between(cell_area / df.index**2, min, max, color=colors[0], alpha=0.2)

plt.xlabel('Number of grid points per cell')
plt.ylabel('Time coupled / time mechanical ')
#plt.xscale('log', base=2)
plt.grid()
plt.savefig('plot_typical_rect.pdf')