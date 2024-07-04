import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np

viridis = colormaps.get_cmap('viridis')
colors = viridis(np.linspace(0, 1, 2))

df = pd.read_csv('benchmark_typical_test.csv')
gp_per_cell = [1, 2, 4, 8, 16, 24, 32] #64
relative_increase = df['time_all'] / df['time_ensemble']

plt.plot(gp_per_cell, relative_increase, marker='o', label=f'', color=colors[0])

#plt.plot(gp_per_cell, df['time_all'], marker='o', label=f'Full simulation', color=colors[1])
#plt.errorbar(subset[''], subset[''], yerr=0, color=colors[i])

#plt.title('Relative increase in runtime of coupled simulation compared to\n mechanical simulation for a typical test case')
plt.xlabel('Number of gridpoints per cell')
plt.ylabel('Time coupled / time mechanical')
#plt.xscale('log', base=2)
plt.grid(True)
#plt.legend()
plt.savefig("benchmark_typical.pdf")