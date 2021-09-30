import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

df = pd.read_csv('SP500.csv')

returns = [1]
for i in range(len(df)-1):
    values = df.loc[i:i+1, 'SP500']
    if values[i+1] == '.' or values[i] == '.':
        rate = np.NaN
    else:
        rate = float(values[i+1]) / float(values[i])
    returns.append(rate)
    
df['returns'] = returns

df['r-1'] = df['returns'] - 1

df['log_r'] = np.log(df['returns'])

df.to_pickle("transformed.pkl")

print(df['log_r'])

