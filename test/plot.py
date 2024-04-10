import matplotlib.pyplot as plt
import math

c0 = 1.0
D = 0.03
dx = 0.01
k = 0.01
L = 100.0
lam = 1.0
if k != 0:
    lam = math.sqrt(D / k)

# analytical solution
def u(x):
    return c0 * math.exp(-x/lam)

def u_linear(x):
    return c0 * (L-x)/L

# Read the values from final_frame.txt
with open('test/final_frame.txt', 'r') as file:
    values = [float(line.strip()) for line in file]
    values = values[1:int(len(values)/5)]

# Generate x values for the analytical solution
x_values = [i * dx for i in range(len(values))]

# Calculate the analytical solution for each x value
analytical_values = [u(x) for x in x_values]

difference = [v - a for v, a in zip(values, analytical_values)]

# Plot the analytical solution
plt.plot(x_values, analytical_values, label='Analytical Solution')
# Plot the polymorph solution
plt.plot(x_values, values, label='Numerical Solution')
# Plot the difference
plt.plot(x_values, difference, label='Difference')
plt.ylim(min(difference), max(difference))

plt.xlabel('x')
plt.ylabel('concentration u')
# Set the legend
plt.legend()
# Save the plot as a PDF file
plt.savefig('/home/uranus/PolyMorph/test/solution_plot.pdf')
