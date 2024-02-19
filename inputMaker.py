import pandas as pd
import numpy as np
from itertools import product

def createModCopy(filename1, filename2, var, val):
    with open(filename1, 'r') as file:
        lines = file.readlines()

    with open(filename2, 'w') as file:
        for line in lines:
            if line.startswith(var):
                file.write(f"{var} = {val}\n")
            else:
                file.write(line)

filename_base = "Input/inputData.py"

Hw      = [1, 2, 4]
H_CB    = [0.25, 0.4]
LDR_CB  = [3.0, 5.0]
bf      = [0.2]
RhoW    = [0.05]
fpc     = [40, 80]
Fy      = [360, 420]
lsr     = [24]
H_typical = [3.]
n_story = [8, 15, 20]
LoadG   = [0]

parameters = [Hw, H_CB, LDR_CB, fpc, Fy, n_story, bf, RhoW, lsr, H_typical, LoadG]

size = np.prod([len(param) for param in parameters])
Array = np.zeros((size, 12))
df = pd.DataFrame(Array, columns=['FileName', 'Hw', 'H_CB', 'LDR_CB', 'fpc', 'Fy', 'n_story', 'bf', 'RhoW', 'lsr', 'H_typical', 'LoadG'])

for i, combination in enumerate(product(*parameters), start=1):
    filename2 = f"Input/inputData{i}.py"
    createModCopy(filename_base, filename2, 'Hw',           combination[0])
    createModCopy(filename2, filename2,     'H_CB',         combination[1])
    createModCopy(filename2, filename2,     'LDR_CB',       combination[2])
    createModCopy(filename2, filename2,     'fpc',          combination[3])
    createModCopy(filename2, filename2,     'Fy',           combination[4])
    createModCopy(filename2, filename2,     'n_story',      combination[5])
    createModCopy(filename2, filename2,     'bf',           combination[6])
    createModCopy(filename2, filename2,     'RhoW',         combination[7])
    createModCopy(filename2, filename2,     'lsr',          combination[8])
    createModCopy(filename2, filename2,     'H_typical',    combination[9])
    createModCopy(filename2, filename2,     'LoadG',        combination[10])

    df.at[i-1, 'FileName'] = filename2[6:-3]
    df.iloc[i-1, 1:] = combination

df.to_excel("Input/inputDataTable.xlsx", sheet_name='table')















