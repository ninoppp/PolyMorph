import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

data = pd.read_csv('positional_error.csv')
grad_cv = data['grad_cv'].unique()
thresh_cv = data['thresh_cv'].unique()

grouped = data.groupby(['thresh_cv', 'grad_cv'])['readout_pos'].std().reset_index()
print(grouped.head(10))

# by tcv
plt.figure(figsize=(10, 5))
for tcv in thresh_cv: 
    subset = grouped[grouped['thresh_cv'] == tcv]
    plt.plot(subset['grad_cv'], subset['readout_pos'], marker='o', label=f'Threshold CV = {tcv}')
plt.title('Positional Error vs Gradient Coefficient-Variation')
plt.xlabel('Gradient CV (D,k,p)')
plt.ylabel('Positional Error')
plt.grid(True)
plt.legend()
plt.savefig("positional_error.pdf")

# by gcv
plt.figure(figsize=(10, 5))
for gcv in grad_cv: 
    subset = grouped[grouped['grad_cv'] == gcv]
    plt.plot(subset['thresh_cv'], subset['readout_pos'], marker='o', label=f'Gradient CV = {gcv}')
plt.title('Positional Error vs Threshold Coefficient-Variation')
plt.xlabel('Threshold CV')
plt.ylabel('Positional Error')
plt.grid(True)
plt.legend()
plt.savefig("positional_error_2.pdf")

# readout position with min and max values as error estimates
plt.figure(figsize=(10, 5))
for tcv in [0.01, 0.3]:
    subset = data[data['thresh_cv'] == tcv]
    readout_pos = subset.groupby('grad_cv')['readout_pos'].mean()
    min = subset.groupby('grad_cv')['readout_pos'].min()
    max = subset.groupby('grad_cv')['readout_pos'].max()
    plt.plot(grad_cv, readout_pos, marker='o', label=f'Threshold CV = {tcv}')
    plt.fill_between(grad_cv, min, max, alpha=0.2)
plt.title('Readout Position vs Gradient Coefficient-Variation')
plt.xlabel('Gradient CV (D,k,p)')
plt.ylabel('Readout Position')
plt.grid(True)
plt.legend()
plt.savefig("positional_error_3.pdf")


