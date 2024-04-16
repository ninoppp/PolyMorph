import matplotlib.pyplot as plt
import math

c0 = 1
D = 0.03
j = 0.03 # influx
dx = 0.01
k = 0.05
L = 1.0
lam = 1.0
if k != 0:
    lam = math.sqrt(D / k)

# analytical solution of standard model u(x->inf) = 0
def u_standard(x):
    return c0 * math.exp(-x/lam)

# solution of double dirichlet
def u_linear(x):
    return c0 * (L-x)/L

# solution of double neumann
def u_2xneumann(x):
    mu = 1/lam
    return j/mu * (math.exp(mu * (2 - x)) + math.exp(mu * x)) / (math.exp(2 * mu) - 1)

# - select the analytical solution - #
u = u_standard

# Read the values from final_frame.txt
with open('test/final_frame.txt', 'r') as file:
    values = [float(line.strip()) for line in file]
    values = values[0:int(len(values)/5)]

x_values = [i * dx for i in range(len(values))]
analytical_values = [u(x) for x in x_values]
difference = [abs(v - a) for v, a in zip(values, analytical_values)]

#plt.plot(x_values[0::4], analytical_values[0::4], ':', label='Analytical Solution')
#plt.plot(x_values[2::4], values[2::4], ':', label='Numerical Solution')
plt.plot(x_values[::10], difference[::10], 'o', ms=2, label='absolute numerical error |u - u_analytical|')

plt.xlabel('x')
#plt.ylabel('concentration u')
plt.yscale('log')

plt.legend()
plt.savefig('/home/uranus/PolyMorph/test/solution_plot.pdf')
