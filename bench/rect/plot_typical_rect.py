import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np

# area in rect tissue = 1.348
# gridpoints per cell = 1.348 / dx^2
cell_area = 1.348

viridis = colormaps.get_cmap('viridis')
colors = viridis(np.linspace(0, 1, 4))
data = pd.read_csv('benchmark_rect_sorted.csv')

df = data[data['width'] == 30]
df = df.groupby('dx').mean()
#drop dx=0.125
#df = df.drop(0.125)
df['ratio'] = df['all']/df['ensemble']

print(df)

plt.figure()
plt.plot(cell_area / (df.index**2), df['ratio'], marker='o', color=colors[0])

plt.xlabel('Number of grid points per cell')
plt.ylabel('Time coupled / time mechanical ')
#plt.xscale('log', base=2)
plt.grid()
plt.savefig('plot_typical_rect.pdf')