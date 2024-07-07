import matplotlib.pyplot as plt
import numpy as np

# Function to create a convex polygon
def generate_convex_polygon(num_vertices):
    # Random angles
    angles = np.sort(np.random.rand(num_vertices) * 2 * np.pi)
    # Random radii
    radii = np.random.rand(num_vertices) * 5 + 1  # radius between 1 and 6
    # Cartesian coordinates
    points = np.vstack((radii * np.cos(angles), radii * np.sin(angles))).T
    # Shift center to middle of the grid
    points[:, 0] += 5
    points[:, 1] += 5
    return points

# Generate convex polygon vertices
np.random.seed(0)
vertices = generate_convex_polygon(10)

# Set up the figure and axis
fig, ax = plt.subplots(figsize=(6,6))
ax.set_xlim([0, 10])
ax.set_ylim([0, 10])
ax.set_xticks(np.arange(0, 11, 1))
ax.set_yticks(np.arange(0, 11, 1))
ax.grid(True)

# Plot the convex polygon
polygon = plt.Polygon(vertices, closed=True, fill=None, edgecolor='black')
ax.add_patch(polygon)

# Hide axes' labels
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.xaxis.set_tick_params(size=0)
ax.yaxis.set_tick_params(size=0)

# Show the plot
plt.savefig("illu.pdf")
