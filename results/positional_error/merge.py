import os
import pandas as pd

# merge all csv files in the folder. for both experiments separately
source = 'data'
folders = ['thresh_01'] #, 'thresh_001', 'thresh_0001']
nodes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
experiments = ['cv', 'width']

for experiment in experiments:
    data = pd.DataFrame()
    for folder in folders:
        for node in nodes:
            file = f'{source}/{folder}/positional_error_{experiment}_{node}.csv'
            if os.path.exists(file):
                df = pd.read_csv(file)
                data = pd.concat([data, df])
            else:
                print(f'Error: {file} does not exist')
    data.to_csv(f'{source}/positional_error_{experiment}.csv', index=False)
    # sort
    sortby = 'grad_cv' if experiment == 'cv' else 'width'
    df.sort_values(['threshold', sortby, 'seed'], ascending=[False, False])
    print(f'{folder}/{experiment}.csv saved')