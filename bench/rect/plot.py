import pandas as pd
import matplotlib.pyplot as plt

data = pd.read_csv('benchmark_rect_sorted.csv')

# plot runtime of ensemble, solver, scatter, gather over number of polygons
fig, ax = plt.subplots()
# group by width
df = data.groupby('polygons').mean().reset_index()
# normalize solver by number of gridpoitns
df['solver'] /= df['gridpoints']
df['scatter'] /= df['gridpoints']
df['gather'] /= df['gridpoints']
# smooth the data by averaging
df = df.rolling(window=3, center=True).mean()
df.plot(x='polygons', y=['ensemble', 'solver', 'scatter', 'gather'], ax=ax)
plt.xlabel('Number of Polygons')
plt.ylabel('Runtime (s)')
plt.title('Runtime of Ensemble, Solver, Scatter, Gather')
plt.savefig('benchmark_rect_polygons.png')

# plot runtime of ensemble, solver, scatter, gather over number of grid points
fig, ax = plt.subplots()
# group by gridpoints
df = data.groupby('gridpoints').mean().reset_index()
# normalize everything by number of polygons
df['ensemble'] /= df['polygons']
# smooth the data by averaging
df = df.rolling(window=3, center=True).mean()
df.plot(x='gridpoints', y=['ensemble', 'solver', 'scatter', 'gather'], ax=ax)
plt.xlabel('Number of Grid Points')
plt.ylabel('Runtime (s)')
plt.title('Runtime of Ensemble, Solver, Scatter, Gather')
plt.savefig('benchmark_rect_gridpoints.png')