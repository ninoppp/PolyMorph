import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# by tcv
data = pd.read_csv('sharpness.csv', sep='\s+')
threshCV = data['ThreshCV'].unique()
avg_sharpness_by_threshCV = data.groupby('ThreshCV')['sharpness'].mean().reset_index()
plt.figure(figsize=(10, 5))
plt.plot(threshCV, avg_sharpness_by_threshCV['sharpness'], marker='o', label='Average')

#for gcv in data['GradCV'].unique():
#    data_by_gcv = data[data['GradCV'] == gcv]
#    avg_sharpness_by_tcv = data_by_gcv.groupby('ThreshCV')['sharpness'].mean().reset_index()
#    plt.plot(avg_sharpness_by_tcv['ThreshCV'], avg_sharpness_by_tcv['sharpness'], marker='o', label=f'Gradient CV = {gcv}')

plt.title('Sharpness vs Threshold Coefficient-Variation')
plt.xlabel('Threshold CV')
plt.ylabel('Border width [units]')
plt.grid(True)
plt.legend()
plt.savefig("sharpness_threshold.pdf")

# by gcv
plt.figure(figsize=(10, 5))
for tcv in threshCV:
    data_by_tcv = data[data['ThreshCV'] == tcv]
    avg_sharpness_by_gcv = data_by_tcv.groupby('GradCV')['sharpness'].mean().reset_index()
    plt.plot(avg_sharpness_by_gcv['GradCV'], avg_sharpness_by_gcv['sharpness'], marker='o', label=f'Threshold CV = {tcv}')

plt.title('Sharpness vs Gradient Coefficient-Variation')
plt.xlabel('Gradient CV (D,k,p)')
plt.ylabel('Border width [units]')
plt.grid(True)
plt.legend()
plt.savefig("sharpness_gradient.pdf")