import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np

viridis = colormaps.get_cmap('viridis')
colors = viridis(np.linspace(0, 1, 4))
data = pd.read_csv('benchmark_rect_sorted.csv')

# over number of polygons
fig, ax = plt.subplots()
# group by width
df = data
df['solver'] /= df['solver'][0]
df['scatter'] /= df['scatter'][0]
df['gather'] /= df['gather'][0]
df['ensemble'] /= df['ensemble'][0]

df = df.groupby('polygons').mean().reset_index()
# normalize solver by number of gridpoitns
# smooth the data by averaging
#df = df.rolling(window=3, center=True).mean()
df.plot(x='polygons', y=['ensemble', 'solver', 'scatter', 'gather'], ax=ax, marker='o', color=colors)
plt.xlabel('Number of Polygons')
plt.ylabel('Runtime (s)')
plt.title('Runtime of Ensemble, Solver, Scatter, Gather')
plt.savefig('benchmark_rect_polygons.png')


# over number of grid points
fig, ax = plt.subplots()
df = data[data['polygons'] == 871] # 1743, 1322, 871
df = df.groupby('gridpoints').mean().reset_index()
# normalize
df['solver'] /= df['solver'][0]
df['scatter'] /= df['scatter'][0]
df['gather'] /= df['gather'][0]
df['ensemble'] /= df['ensemble'][0]

df.plot(x='gridpoints', y=['ensemble', 'solver', 'scatter', 'gather'], ax=ax, marker='o', color=colors)

plt.xlabel('Number of grid points')
plt.ylabel('Relative increase in runtime (T / T0)')
plt.xscale('log')
plt.yscale('log')
plt.title('Runtime vs. number of grid points')
plt.savefig('benchmark_rect_gridpoints.png')