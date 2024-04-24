import numpy as np
import matplotlib.pyplot as plt

CV = [0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 1.0]
measurements = {
    0.1: [3.39008, 6.12622],
    0.2: [7.77738, 10.1409],
    0.3: [9.78458, 15.3179],
    0.4: [14.8297],
    0.5: [14.8039, 16.6108],
    0.7: [20.2002],
    1.0: [24.8346]
}

avg_sharpness = [np.mean(measurements[cv]) for cv in CV]

plt.plot(CV, avg_sharpness, 'o-')
plt.xlabel('CV value for D, k, p, readout')
plt.ylabel('border width [units]')
plt.savefig('sharpness.pdf')