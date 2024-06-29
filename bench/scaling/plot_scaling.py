import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np

df = pd.read_csv('scaling.csv')
resolutions = df['dx'].unique()

#df.sort_values(['node_id', 'dx', 'gridpoints', 'threads'], ascending=[True, False, True, True], inplace=True)
#df.to_csv('scaling.csv', index=False)
#exit()
viridis = colormaps.get_cmap('viridis')  # Get the viridis colormap
colors = viridis(np.linspace(0, 1, resolutions.size))

plt.figure(figsize=(10, 5))
for i, dx in enumerate(resolutions):
    # average over dx and threads
    df_dx = df[df['dx'] == dx].groupby('threads').mean()
    plt.plot(df_dx.index, df_dx['time_all'], label=f'dx={dx}', marker='o', color=colors[i])
    df_dx['time_all_std'] = df[df['dx'] == dx].groupby('threads')['time_all'].std()
    plt.errorbar(df_dx.index, df_dx['time_all'], yerr=df_dx['time_all_std'], fmt='o', capsize=5, color=colors[i])

plt.xlabel('Threads')
plt.ylabel('Time (s)')
plt.xscale('log', base=2)
#plt.yscale('log', base=10)
plt.title('Scaling')
plt.legend()
plt.grid()
plt.savefig('scaling.png')

# TODO: make scaling and not runtime plot. do scaling measurement multiple times for better statistics


