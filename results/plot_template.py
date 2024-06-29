import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np

viridis = colormaps.get_cmap('viridis')
colors = viridis(np.linspace(0, 1, NUMBER))

df = pd.read_csv('data.csv')
df = df.groupby('').mean().reset_index()

plt.figure(figsize=(10, 5))
for i, var in enumerate(VARIABLES):
    subset = df[df['threshold'] == var]
    plt.plot(subset[''], subset[''], marker='o', label=f'Label = {1}', color=colors[i])
    plt.errorbar(subset[''], subset[''], yerr=0, color=colors[i])

plt.title('Title')
plt.xlabel('xlabel')
plt.ylabel('ylabel')
plt.grid(True)
plt.legend()
plt.savefig("plot.pdf")