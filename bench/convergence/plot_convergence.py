import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

data = pd.read_csv('solution.csv')

fig, ax = plt.subplots()
line1, = ax.plot([], [], 'r-', label='Analytic')  # Red line for analytic solution
line2, = ax.plot([], [], 'b-', label='Simulated')  # Blue line for simulated solution
ax.set_xlim(data['x'].min(), data['x'].max())  # Set x-axis limits
ax.set_ylim(min(data['analytic'].min(), data['simulated'].min()), max(data['analytic'].max(), data['simulated'].max()))  # Set y-axis limits
ax.legend()
ax.set_xlabel('Position (x)')
ax.set_ylabel('Solution Value')
ax.set_title('Convergence of Analytic and Simulated Solutions')

def update(frame):
    df = data[data['t'] == frame]
    line1.set_data(df['x'], df['analytic'])
    line2.set_data(df['x'], df['simulated'])
    return line1, line2

ani = FuncAnimation(fig, update, frames=np.unique(data['t']), blit=True, interval=100, repeat=True)
ani.save('animation.mp4', writer='ffmpeg')
