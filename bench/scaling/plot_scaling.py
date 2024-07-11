import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np

df = pd.read_csv('scaling.csv')
resolutions = df['dx'].unique()
threads = df['threads'].unique()

def amdahl(P, t):
    return 1.0 / (1.0 - P + P / t)

#df.sort_values(['node_id', 'dx', 'gridpoints', 'threads'], ascending=[True, False, True, True], inplace=True)
#df.to_csv('scaling.csv', index=False)
#exit()
viridis = colormaps.get_cmap('viridis')  # Get the viridis colormap
colors = viridis(np.linspace(0, 1, resolutions.size+1))

plt.figure(figsize=(7, 5))
for i, dx in enumerate(resolutions):
    df_dx = df[df['dx'] == dx].copy()
    df_dx.reset_index(drop=True, inplace=True)
    df_dx.loc[:, 'speedup'] = df_dx['time_all'][0] / df_dx['time_all']
    speedup = df_dx.groupby('threads')['speedup'].mean()
    plt.plot(threads, speedup, label=f'Δx={dx}', marker='o', color=colors[i])
    #plt.errorbar(df_dx['threads'].unique(), speedup, yerr=speedup.std(), fmt='o', capsize=5, color=colors[i])
    #T0 = df_dx['time_all'][1]
    #speedup = T0 / df_dx['time_all']
    #plt.plot(df_dx.index, speedup, label=f'dx={dx}', marker='o', color=colors[i])
    #plt.plot(df_dx.index, df_dx['time_all'], label=f'dx={dx}', marker='o', color=colors[i])
    #df_dx['time_all_std'] = df[df['dx'] == dx].groupby('threads')['time_all'].std()
    #plt.errorbar(df_dx.index, df_dx['time_all'], yerr=df_dx['time_all_std'], fmt='o', capsize=5, color=colors[i])

speedup_reference = amdahl(0.90, threads)
plt.plot(threads, speedup_reference, label='Amdahl\'s Law, P=90%', color='gray', linestyle='--')

plt.xlabel('Number of threads p')
plt.ylabel('Speedup S = T(1) / T(p)')
plt.xscale('log', base=2)
plt.yscale('log', base=10)
# add more 10^1 ticks
plt.gca().set_yticks([1, 10])
plt.title('')
plt.legend()
#plt.grid()
plt.savefig('scaling.pdf')

