import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np

viridis = colormaps.get_cmap('viridis')
colors = viridis(np.linspace(0, 1, 3))
length = 300

df = pd.read_csv('convergence.csv')

plt.figure(figsize=(10, 5))

plt.plot(df['dx'], df['inf_norm'], label=f'Infinity norm', marker='o', color=colors[1])
plt.plot(df['dx'], df['rmse'], label=f'RMSE', marker='o', color=colors[0])

# reference (order 2)
dx = df['dx'].values
rmse = df['rmse'].values
inf_norm = df['inf_norm'].values
#plt.plot(dx, rmse[0]*(dx/dx[0])**2, label='Reference 2nd order', linestyle='--', color='gray')
plt.plot(dx, 1/(length/df['dx'])**2, label='Reference: 1/N^2', linestyle=':', color='gray')

plt.gca().invert_xaxis()
plt.title('')
plt.xlabel('dx')
plt.ylabel('numerical error')
plt.xscale('log')
plt.yscale('log')
plt.xticks(df['dx'], [f'{dx:.{len(str(dx))-2}f}' for dx in df['dx']])
#plt.grid(True)
plt.legend()
plt.savefig("convergence.pdf")