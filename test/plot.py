import matplotlib.pyplot as plt
import math

c0 = 1
D = 0.03
j = 0.01 # influx
dx = 0.01
k = 0.01
L = 100.0
lam = 1.0
if k != 0:
    lam = math.sqrt(D / k)

# analytical solution of standard model
def u_standard(x):
    return c0 * math.exp(-x/lam)

# solution of double dirichlet
def u_linear(x):
    return c0 * (L-x)/L

# solution of double neumann
def u_2xneumann(x):
    mu = 1/lam
    return j/mu * (math.exp(mu * (-2 * L + x)) + math.exp(-mu * x)) / (1 - math.exp(2 * mu))

# - select the analytical solution - #
u = u_standard

# Read the values from final_frame.txt
with open('test/final_frame.txt', 'r') as file:
    values = [float(line.strip()) for line in file]
    values = values[0:int(len(values)/5)]

# Generate x values for the analytical solution
x_values = [i * dx for i in range(len(values))]

# Calculate the analytical solution for each x value
analytical_values = [u(x) for x in x_values]

difference = [abs(v - a) for v, a in zip(values, analytical_values)]
print(values[0], analytical_values[0], difference[0])

# Plot the analytical solution
#plt.plot(x_values, analytical_values, ':', label='Analytical Solution')
# Plot the polymorph solution
#plt.plot(x_values, values, ':', label='Numerical Solution')
# Plot the difference
plt.plot(x_values[::10], difference[::10], 'o', ms=3, label='absolute numerical error |u - u_analytical|')

#plt.ylabel('concentration u')
#plt.ylim(min(difference), max(difference))
plt.yscale('log')
plt.xlabel('x')

# Set the legend
plt.legend()
# Save the plot as a PDF file
plt.savefig('/home/uranus/PolyMorph/test/solution_plot.pdf')
