import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import colormaps

data = pd.read_csv('gradient.csv')
length = 60 # interesting part of domain
viridis = colormaps.get_cmap('coolwarm')
colors = viridis(np.linspace(0, 1, 4))
data['i'] = data['i'] * 0.3 - 10
data['u'] = data['u'] + 1e-6 # avoid log(0)

fig, ax = plt.subplots()
line1, = ax.plot([], [], color=colors[0], label='y = width / 4')  # Red line for analytic solution
line2, = ax.plot([], [], color=colors[2], label='y = width / 2')  # Blue line for simulated solution
line3, = ax.plot([], [], color=colors[3], label='y = 3 * width / 2')  # Blue line for simulated solution
ax.set_xlim(data['i'].min(), data['i'].max())  # Set x-axis limits
ax.set_ylim(min(data['u'].min(), data['u'].min()), max(data['u'].max(), data['u'].max()))  # Set y-axis limits
ax.legend()
# use log scale for y
#ax.set_yscale('log')
ax.set_xlabel('x')
ax.set_ylabel('Concentration c(x)')
#ax.set_title('Convergence of numerical solution (dx=0.3, dt=1e-4)')

def update(frame):
    df = data[data['frame'] == frame]
    df1 = df[df['j'] == 25]
    line1.set_data(df1['i'], df1['u'])
    df2 = df[df['j'] == 50]
    line2.set_data(df2['i'], df2['u'])
    df3 = df[df['j'] == 75]
    line3.set_data(df3['i'], df3['u'])
    return line1, line2, line3

ani = FuncAnimation(fig, update, frames=np.unique(data['frame']), blit=True, interval=83, repeat=True)
ani.save('gradient_animation.mp4', writer='ffmpeg')

