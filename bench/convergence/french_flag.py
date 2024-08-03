import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colormaps
import numpy as np
import math
import matplotlib.patches as patches

# colormap & latex
coolwarm = colormaps.get_cmap('coolwarm')  # Get the viridis colormap
plt.rcParams['text.usetex'] = True
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],  # This line may need to be adjusted
})
colors = coolwarm(np.linspace(0, 1, 3))

c0 = 1
x_theta_1 = 1/3
x_theta_2 = 2/3
x_ticks = np.array([0, x_theta_1, x_theta_2])
alpha = 0.8

def f(x):
    return c0 * np.exp(-3*x)

x = np.linspace(0, 1, 1000)
c = f(x)

plt.figure(figsize=(7, 5))
plt.plot(x, c, color='black')

plt.fill_between(x, c, where=(x <= x_theta_1), color=colors[0], alpha=alpha)
plt.fill_between(x, c, where=(x > x_theta_1) & (x < x_theta_2), color=colors[1], alpha=alpha)
plt.fill_between(x, c, where=(x >= x_theta_2), color=colors[2], alpha=alpha)

plt.xticks([])
plt.yticks([])
plt.xticks(x_ticks, [r'$0$', r'$x_{\theta_1}$', r'$x_{\theta_2}$'])  # Place the tick at x=5, labeled 'x_theta'
plt.yticks(f(x_ticks), [r'$c_0$', r'$c_{\theta_1}$', r'$c_{\theta_2}$'])  # Place the tick at y=0, labeled 'c_theta'

# Add cell boxes below the plot
ax = plt.gca()
for i in range(24):
    x_box = i * (1/24)
    y_box = -0.2
    color = color=colors[0] if i < 8 else colors[2] if i >= 16 else colors[1]
    rect = patches.Rectangle((x_box, y_box), 1/24, 0.1, linewidth=1, edgecolor='black', facecolor=color, alpha=alpha)
    ax.add_patch(rect)

# Add label for cells
plt.text(0.5, -0.25, r'cells', ha='center')

equation = r"$c(x) = c_0 e^{-x/\lambda}$"
plt.text(0.5, 0.6, equation, fontsize=15, ha='center')

# Adjust plot limits
plt.ylim(-0.3, 1)
plt.xlim(0, 1)

plt.grid(True)
plt.xlabel(r'Distance $x$ from source')
plt.ylabel(r'Concentration $c(x)$')
plt.savefig('french_flag.pdf')
