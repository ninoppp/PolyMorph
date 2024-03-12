import matplotlib.pyplot as plt
import math

c0 = 1.0
D = 0.03
dx = 0.01
k = 0.01
lam = math.sqrt(D / k)

# analytical solution
def u(x):
    return c0 * math.exp(-x/lam)

# Read the values from final_frame.txt
with open('final_frame.txt', 'r') as file:
    values = [float(line.strip()) for line in file]

# Generate x values for the analytical solution
x_values = [i * dx for i in range(len(values))]

# Calculate the analytical solution for each x value
analytical_values = [u(x) for x in x_values]

# Plot the analytical solution
plt.plot(x_values, analytical_values, label='Analytical Solution')
# Plot the polymorph solution
plt.plot(x_values, values, label='Numerical Solution')

# Set the legend
plt.legend()
# Save the plot as a PDF file
plt.savefig('/home/uranus/PolyMorph/test/analytical_solution.pdf')
