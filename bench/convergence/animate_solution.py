import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib import colormaps

data = pd.read_csv('solution.csv')
length = 50 # interesting part of domain
viridis = colormaps.get_cmap('viridis')
colors = viridis(np.linspace(0, 1, 2))

fig, ax = plt.subplots()
line1, = ax.plot([], [], color=colors[0], label='Analytic solution')  # Red line for analytic solution
line2, = ax.plot([], [], color=colors[1], label='Simulated solution')  # Blue line for simulated solution
ax.set_xlim(data['x'].min(), length)  # Set x-axis limits
ax.set_ylim(min(data['analytic'].min(), data['simulated'].min()), max(data['analytic'].max(), data['simulated'].max()))  # Set y-axis limits
ax.legend()
ax.set_xlabel('Position x')
ax.set_ylabel('Concentration c(x)')
ax.set_title('Convergence of numerical solution (dx=0.3, dt=1e-4)')

def update(frame):
    df = data[data['t'] == frame]
    line1.set_data(df['x'], df['analytic'])
    line2.set_data(df['x'], df['simulated'])
    return line1, line2

ani = FuncAnimation(fig, update, frames=np.unique(data['t']), blit=True, interval=100, repeat=True)
ani.save('animation.mp4', writer='ffmpeg')

