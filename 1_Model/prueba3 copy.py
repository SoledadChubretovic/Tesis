#%% otro
import numpy as np
import tabulate

N_VARIABLES = 19
resF = [[4.19593172, 1.75351547], [8.35110426, 1.21019621]]
resX = [[ 11,  13,  25,  24,  22,  32,  11,  11,  40,  14,  13,  13,  19,  23,  14,  46,  11,  78, 61], [26,  17,  21,  17,  34,  29,  13,  26,  12,  14,  17 , 11,  44 , 15 , 24 , 21 , 12 ,145, 61]]
data1 = resF
data2 = resX

table_data = []
for i in range(0,len(data1)):
    linedata = []
    for j in data1[i]:
        linedata.append(j)
    for k in data2[i]:
        linedata.append(k)
    linedata[-1] = linedata[-1]/100
    linedata[-2] = linedata[-2]*2
    table_data.append(linedata)

# column tags
column_names = ["Weight [kg]","U_value [W/m2K]"]
for i in range(0,N_VARIABLES - 3):
    column_names = column_names + (["l" + str(i + 1) + " [mm]"])
column_names = column_names + (["w [mm]"])
column_names = column_names + (["W [mm]"])
column_names = column_names + (["λ_clay [W/mK]"])

#%%
from tabulate import tabulate
result_table = tabulate(table_data, column_names)
print(result_table)
# %%
results = open(r"C:\Users\nchubretovic\OneDrive - Entel\Escritorio\Sole\Tesis\OptPython\OptPython_Log.txt", "w")
results.write(result_table)
results.close()

# %%
